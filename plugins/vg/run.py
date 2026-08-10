#!/usr/bin/env python3
"""PGBench paired-short-read adapter for vg Giraffe + pack + call."""

from __future__ import annotations

import os
import subprocess
from pathlib import Path


def required(name: str) -> str:
    value = os.environ.get(name)
    if not value:
        raise RuntimeError(f"required environment variable is absent: {name}")
    return value


def run(command: list[str], *, stdout: Path | None = None) -> None:
    if stdout is None:
        subprocess.run(command, check=True)
        return
    with stdout.open("wb") as handle:
        subprocess.run(command, check=True, stdout=handle)


def main() -> int:
    reads_r1 = Path(required("PGBENCH_INPUT_FASTQ_R1"))
    reads_r2 = Path(required("PGBENCH_INPUT_FASTQ_R2"))
    graph_dir = Path(required("PGBENCH_GRAPH_DIR"))
    output_dir = Path(required("PGBENCH_OUTPUT_DIR"))
    output_vcf = Path(required("PGBENCH_OUTPUT_VCF"))
    sample = required("PGBENCH_SAMPLE_ID")
    threads = required("PGBENCH_THREADS")
    reference_path = required("PGBENCH_GRAPH_REFERENCE_PATH")

    gbz = graph_dir / "graph.gbz"
    minimizer = graph_dir / "graph.shortread.withzip.min"
    zipcodes = graph_dir / "graph.shortread.zipcodes"
    distance = graph_dir / "graph.dist"
    for path in (gbz, minimizer, zipcodes, distance):
        if not path.is_file():
            raise RuntimeError(f"required vg graph asset is absent: {path}")

    work = output_dir / "work"
    work.mkdir(parents=True, exist_ok=True)
    output_vcf.parent.mkdir(parents=True, exist_ok=True)
    gam = work / "HG002.gam"
    pack = work / "HG002.pack"
    snarls = work / "graph.snarls"

    # Giraffe is vg's haplotype-aware production mapper for paired short reads.
    # Explicit index paths prevent implicit, mutable index construction.
    run(
        [
            "vg",
            "giraffe",
            "-Z",
            str(gbz),
            "-d",
            str(distance),
            "-m",
            str(minimizer),
            "-z",
            str(zipcodes),
            "-f",
            str(reads_r1),
            "-f",
            str(reads_r2),
            "-t",
            threads,
            "-p",
        ],
        stdout=gam,
    )
    run(
        [
            "vg",
            "pack",
            "-x",
            str(gbz),
            "-g",
            str(gam),
            "-o",
            str(pack),
            "-t",
            threads,
        ]
    )
    run(["vg", "snarls", str(gbz), "-t", threads], stdout=snarls)
    run(
        [
            "vg",
            "call",
            str(gbz),
            "-r",
            str(snarls),
            "-k",
            str(pack),
            "-s",
            sample,
            "-z",
            "--ref-sample",
            reference_path,
        ],
        stdout=output_vcf,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
