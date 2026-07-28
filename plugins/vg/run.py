#!/usr/bin/env python3
"""PGBench CLR graph adapter for GraphAligner + vg pack + vg call."""

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
    reads = Path(required("PGBENCH_INPUT_FASTQ"))
    graph_dir = Path(required("PGBENCH_GRAPH_DIR"))
    output_dir = Path(required("PGBENCH_OUTPUT_DIR"))
    output_vcf = Path(required("PGBENCH_OUTPUT_VCF"))
    sample = required("PGBENCH_SAMPLE_ID")
    threads = required("PGBENCH_THREADS")

    gbz = graph_dir / "graph.gbz"
    for path in (
        gbz,
        graph_dir / "graph.xg",
        graph_dir / "graph.min",
        graph_dir / "graph.dist",
    ):
        if not path.is_file():
            raise RuntimeError(f"required vg graph asset is absent: {path}")

    work = output_dir / "work"
    work.mkdir(parents=True, exist_ok=True)
    output_vcf.parent.mkdir(parents=True, exist_ok=True)
    gfa = work / "graph.gfa"
    gam = work / "HG002.gam"
    pack = work / "HG002.pack"
    snarls = work / "graph.snarls"

    # GraphAligner is appropriate for noisy PacBio CLR reads and can emit GAM
    # with node identities preserved from the exact frozen graph.
    run(
        [
            "vg",
            "convert",
            "-f",
            "--no-translation",
            str(gbz),
        ],
        stdout=gfa,
    )
    run(
        [
            "GraphAligner",
            "-g",
            str(gfa),
            "-f",
            str(reads),
            "-a",
            str(gam),
            "-x",
            "vg",
            "-t",
            threads,
        ]
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
        ],
        stdout=output_vcf,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
