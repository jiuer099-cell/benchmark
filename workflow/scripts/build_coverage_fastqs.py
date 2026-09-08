#!/usr/bin/env python3
"""Build deterministic paired-FASTQ coverage replicates shared by every tool."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import sys
from itertools import zip_longest
from pathlib import Path
from typing import BinaryIO, Iterator, Sequence


class CoverageError(ValueError):
    """Raised when paired input or coverage parameters are invalid."""


def _open(path: Path, mode: str) -> BinaryIO:
    if path.name.endswith(".gz"):
        return gzip.open(path, mode)
    return path.open(mode)


def _records(path: Path) -> Iterator[tuple[bytes, bytes, bytes, bytes]]:
    with _open(path, "rb") as handle:
        while True:
            record = tuple(handle.readline() for _ in range(4))
            if record[0] == b"":
                return
            if any(part == b"" for part in record) or not record[0].startswith(b"@"):
                raise CoverageError(f"malformed FASTQ: {path}")
            yield record  # type: ignore[misc]


def _name(record: tuple[bytes, bytes, bytes, bytes]) -> bytes:
    value = record[0][1:].split(maxsplit=1)[0].rstrip(b"\r\n")
    return value[:-2] if value.endswith((b"/1", b"/2")) else value


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def build(
    *, r1: Path, r2: Path, output_dir: Path, source_coverage: float,
    coverages: Sequence[int], seeds: Sequence[int],
) -> dict:
    if source_coverage <= 0 or any(value <= 0 or value > source_coverage for value in coverages):
        raise CoverageError("requested coverage must be within the source coverage")
    output_dir.mkdir(parents=True, exist_ok=True)
    outputs: dict[
        tuple[int, int], tuple[BinaryIO, BinaryIO, BinaryIO, BinaryIO, Path, Path]
    ] = {}
    for coverage in coverages:
        for seed in seeds:
            directory = output_dir / f"{coverage}x" / f"seed-{seed}"
            directory.mkdir(parents=True, exist_ok=True)
            out1, out2 = directory / "reads.R1.fastq.gz", directory / "reads.R2.fastq.gz"
            raw1, raw2 = out1.open("wb"), out2.open("wb")
            outputs[(coverage, seed)] = (
                gzip.GzipFile(filename="", mode="wb", fileobj=raw1, mtime=0),
                gzip.GzipFile(filename="", mode="wb", fileobj=raw2, mtime=0),
                raw1, raw2, out1, out2,
            )
    counts = {key: 0 for key in outputs}
    bases = {key: 0 for key in outputs}
    total = 0
    total_bases = 0
    try:
        for pair_index, pair in enumerate(zip_longest(_records(r1), _records(r2)), start=1):
            left, right = pair
            if left is None or right is None or _name(left) != _name(right):
                raise CoverageError(f"paired FASTQ mismatch at pair {pair_index}")
            total += 1
            pair_bases = len(left[1].rstrip(b"\r\n")) + len(right[1].rstrip(b"\r\n"))
            total_bases += pair_bases
            name = _name(left)
            for (coverage, seed), (writer1, writer2, _, _, _, _) in outputs.items():
                value = int.from_bytes(hashlib.sha256(str(seed).encode() + b"\0" + name).digest()[:8], "big")
                if value / 2**64 < coverage / source_coverage:
                    writer1.write(b"".join(left))
                    writer2.write(b"".join(right))
                    counts[(coverage, seed)] += 1
                    bases[(coverage, seed)] += pair_bases
    finally:
        for writer1, writer2, raw1, raw2, _, _ in outputs.values():
            writer1.close()
            writer2.close()
            raw1.close()
            raw2.close()
    if total == 0:
        raise CoverageError("paired FASTQ is empty")
    replicates = []
    for (coverage, seed), (_, _, _, _, out1, out2) in sorted(outputs.items()):
        replicates.append({
            "coverage_x": coverage, "seed": seed, "read_pairs": counts[(coverage, seed)],
            "read_count": counts[(coverage, seed)] * 2,
            "read_bases": bases[(coverage, seed)],
            "fastq_r1": str(out1), "fastq_r2": str(out2),
            "fastq_r1_sha256": _sha256(out1), "fastq_r2_sha256": _sha256(out2),
        })
    return {
        "schema_version": 1, "contract": "pgbench_coverage_replicates_v1",
        "source_coverage_x": source_coverage,
        "source_fastq_r1_sha256": _sha256(r1), "source_fastq_r2_sha256": _sha256(r2),
        "selection": "sha256(seed,read_name)_threshold_v1", "replicates": replicates,
        "full_depth": {
            "coverage_x": source_coverage,
            "seed": None,
            "read_pairs": total,
            "read_count": total * 2,
            "read_bases": total_bases,
            "fastq_r1": str(r1),
            "fastq_r2": str(r2),
            "fastq_r1_sha256": _sha256(r1),
            "fastq_r2_sha256": _sha256(r2),
        },
    }


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fastq-r1", required=True, type=Path)
    parser.add_argument("--fastq-r2", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--source-coverage", required=True, type=float)
    parser.add_argument("--coverage", action="append", type=int, default=[])
    parser.add_argument("--seed", action="append", type=int, default=[])
    parser.add_argument("--manifest", required=True, type=Path)
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        result = build(
            r1=args.fastq_r1, r2=args.fastq_r2, output_dir=args.output_dir,
            source_coverage=args.source_coverage,
            coverages=args.coverage or [10, 20, 30], seeds=args.seed or [1701, 1702, 1703],
        )
        args.manifest.parent.mkdir(parents=True, exist_ok=True)
        args.manifest.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    except (OSError, CoverageError) as exc:
        print(f"build_coverage_fastqs: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
