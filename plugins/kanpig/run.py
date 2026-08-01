#!/usr/bin/env python3
"""PGBench adapter for KanPIG germline SV genotyping."""

from __future__ import annotations

import os
import subprocess
from pathlib import Path


def required(name: str) -> str:
    value = os.environ.get(name)
    if not value:
        raise RuntimeError(f"required environment variable is absent: {name}")
    return value


def main() -> int:
    candidate = Path(required("PGBENCH_CANDIDATE_VCF"))
    alignment = Path(required("PGBENCH_SHARED_ALIGNMENT"))
    reference = Path(required("PGBENCH_REFERENCE_FASTA"))
    output = Path(required("PGBENCH_OUTPUT_VCF"))
    threads = required("PGBENCH_THREADS")
    output.parent.mkdir(parents=True, exist_ok=True)

    command = [
        "kanpig",
        "gt",
        "--input",
        str(candidate),
        "--reads",
        str(alignment),
        "--reference",
        str(reference),
        "--out",
        str(output),
        "--threads",
        threads,
    ]
    completed = subprocess.run(command, check=False)
    return completed.returncode


if __name__ == "__main__":
    raise SystemExit(main())
