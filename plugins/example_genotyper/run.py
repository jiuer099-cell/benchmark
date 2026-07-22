#!/usr/bin/env python3
"""Deterministic mock genotyper for the PGBench external-tool example."""

from __future__ import annotations

import gzip
import hashlib
import io
import os
import tempfile
from collections.abc import Iterator
from pathlib import Path
from typing import TextIO


def _open_text(path: Path) -> TextIO:
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("r", encoding="utf-8")


def _deterministic_gt(candidate_id: str, sample_id: str) -> str:
    value = hashlib.sha256(f"{candidate_id}\0{sample_id}".encode()).digest()[0]
    return ("0/0", "0/1", "0/1", "1/1")[value % 4]


def _output_lines(candidate_vcf: Path, sample_id: str) -> Iterator[str]:
    saw_format = False
    saw_header = False
    with _open_text(candidate_vcf) as source:
        for line in source:
            if line.startswith("##FORMAT=<ID=GT,"):
                saw_format = True
                yield line if line.endswith("\n") else f"{line}\n"
                continue
            if line.startswith("##"):
                yield line if line.endswith("\n") else f"{line}\n"
                continue
            if line.startswith("#CHROM"):
                if not saw_format:
                    yield (
                        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
                    )
                yield "##source=PGBench-example-genotyper-0.1.0\n"
                header = line.rstrip("\n").split("\t")[:8]
                yield "\t".join([*header, "FORMAT", sample_id]) + "\n"
                saw_header = True
                continue
            if not line.strip():
                continue
            if not saw_header:
                raise ValueError("candidate VCF is missing a #CHROM header")
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                raise ValueError("candidate VCF record has fewer than 8 columns")
            candidate_id = fields[2]
            if not candidate_id or candidate_id == ".":
                raise ValueError("candidate VCF record has no blinded candidate ID")
            gt = _deterministic_gt(candidate_id, sample_id)
            yield "\t".join([*fields[:8], "GT", gt]) + "\n"
    if not saw_header:
        raise ValueError("candidate VCF is missing a #CHROM header")


def write_deterministic_vcf(
    candidate_vcf: Path,
    output_vcf: Path,
    sample_id: str,
) -> None:
    output_vcf.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        dir=output_vcf.parent,
        prefix=f".{output_vcf.name}.",
        suffix=".tmp",
    )
    temporary = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "wb") as raw:
            if output_vcf.suffix == ".gz":
                with gzip.GzipFile(
                    filename="",
                    mode="wb",
                    fileobj=raw,
                    mtime=0,
                ) as compressed:
                    with io.TextIOWrapper(
                        compressed, encoding="utf-8", newline=""
                    ) as text:
                        text.writelines(_output_lines(candidate_vcf, sample_id))
            else:
                text = io.TextIOWrapper(raw, encoding="utf-8", newline="")
                text.writelines(_output_lines(candidate_vcf, sample_id))
                text.flush()
                text.detach()
            raw.flush()
            os.fsync(raw.fileno())
        os.replace(temporary, output_vcf)
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise


def _required_environment(name: str) -> str:
    value = os.environ.get(name)
    if not value:
        raise RuntimeError(f"required environment variable is absent: {name}")
    return value


def main() -> int:
    candidate_vcf = Path(_required_environment("PGBENCH_CANDIDATE_VCF"))
    output_vcf = Path(_required_environment("PGBENCH_OUTPUT_VCF"))
    output_dir = Path(_required_environment("PGBENCH_OUTPUT_DIR")).resolve()
    sample_id = _required_environment("PGBENCH_SAMPLE_ID")
    try:
        output_vcf.resolve().relative_to(output_dir)
    except ValueError as exc:
        raise RuntimeError("PGBENCH_OUTPUT_VCF escapes PGBENCH_OUTPUT_DIR") from exc
    write_deterministic_vcf(candidate_vcf, output_vcf, sample_id)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
