#!/usr/bin/env python3
"""Fail-closed lock for the one BAM shared by alignment-consuming adapters."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import tempfile
from pathlib import Path


class SharedAlignmentError(ValueError):
    pass


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def paired_fastq_sha256(r1: Path, r2: Path) -> str:
    """Hash the ordered pair, never an ambiguous concatenation of FASTQs."""

    payload = {"r1_sha256": sha256(r1), "r2_sha256": sha256(r2)}
    return hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


def validate(args: argparse.Namespace) -> dict[str, object]:
    declared = {
        "source_fastq_sha256": args.source_fastq_sha256,
        "reference_sha256": args.reference_sha256,
        "bam_sha256": args.bam_sha256,
        "bai_sha256": args.bai_sha256,
    }
    actual = {
        "source_fastq_sha256": paired_fastq_sha256(args.fastq_r1, args.fastq_r2),
        "reference_sha256": sha256(args.reference),
        "bam_sha256": sha256(args.bam),
        "bai_sha256": sha256(args.bai),
    }
    mismatches = [name for name, expected in declared.items() if actual[name] != expected]
    if mismatches:
        raise SharedAlignmentError("shared alignment hash mismatch: " + ", ".join(mismatches))
    if len(args.command_sha256) != 64 or any(char not in "0123456789abcdef" for char in args.command_sha256.casefold()):
        raise SharedAlignmentError("command_sha256 must be a SHA-256 hex digest")
    if not args.aligner or not args.aligner_version:
        raise SharedAlignmentError("aligner and aligner_version must be declared")
    return {
        "contract": "pgbench_shared_shortread_alignment_v1",
        "status": "verified",
        "source_fastq_sha256": actual["source_fastq_sha256"],
        "reference_sha256": actual["reference_sha256"],
        "aligner": args.aligner,
        "aligner_version": args.aligner_version,
        "command_sha256": args.command_sha256,
        "bam_sha256": actual["bam_sha256"],
        "bai_sha256": actual["bai_sha256"],
        "bam": str(args.bam.resolve()),
        "bai": str(args.bai.resolve()),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    for name in ("fastq-r1", "fastq-r2", "reference", "bam", "bai"):
        parser.add_argument(f"--{name}", required=True, type=Path)
    for name in ("source-fastq-sha256", "reference-sha256", "command-sha256", "bam-sha256", "bai-sha256"):
        parser.add_argument(f"--{name}", required=True)
    parser.add_argument("--aligner", required=True)
    parser.add_argument("--aligner-version", required=True)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    try:
        for path in (args.fastq_r1, args.fastq_r2, args.reference, args.bam, args.bai):
            if not path.is_file() or path.stat().st_size == 0:
                raise SharedAlignmentError(f"required shared-alignment input is missing or empty: {path}")
        payload = validate(args)
        args.output.parent.mkdir(parents=True, exist_ok=True)
        descriptor, temporary_name = tempfile.mkstemp(dir=args.output.parent, prefix=".shared-alignment-", text=True)
        with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
            json.dump(payload, handle, indent=2, sort_keys=True)
            handle.write("\n")
        Path(temporary_name).replace(args.output)
    except SharedAlignmentError as exc:
        parser.error(str(exc))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
