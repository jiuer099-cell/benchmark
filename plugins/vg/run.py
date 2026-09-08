#!/usr/bin/env python3
"""PGBench paired-short-read adapter for vg Giraffe + pack + call."""

from __future__ import annotations

import gzip
import os
import subprocess
from pathlib import Path
from typing import TextIO


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


def open_text(path: Path, mode: str) -> TextIO:
    return gzip.open(path, mode, encoding="utf-8") if path.name.endswith(".gz") else path.open(mode, encoding="utf-8")


def _allele_id(info: str) -> str | None:
    for item in info.split(";"):
        if item.startswith("PANGENOME_ALLELE_ID="):
            return item.split("=", 1)[1]
    return None


def _genotype(fields: list[str]) -> str:
    if len(fields) < 10:
        return "./."
    names = fields[8].split(":")
    values = fields[9].split(":")
    if "GT" not in names or names.index("GT") >= len(values):
        return "./."
    value = values[names.index("GT")]
    return value if value not in {"", "."} else "./."


def project_all_sites(native: Path, candidates: Path, destination: Path, sample: str) -> tuple[int, int, int]:
    """Project vg calls onto the complete frozen panel with explicit no-calls."""

    order: list[str] = []
    records: dict[str, list[str]] = {}
    by_allele: dict[str, str] = {}
    by_key: dict[tuple[str, str, str, str], str] = {}
    meta: list[str] = []
    with open_text(candidates, "rt") as handle:
        for line in handle:
            if line.startswith("##"):
                meta.append(line)
                continue
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8 or not fields[2]:
                raise RuntimeError("candidate panel contains a malformed record")
            candidate_id = fields[2]
            if candidate_id in records:
                raise RuntimeError(f"candidate panel duplicates {candidate_id}")
            records[candidate_id] = fields
            order.append(candidate_id)
            by_key[(fields[0], fields[1], fields[3], fields[4])] = candidate_id
            allele_id = _allele_id(fields[7])
            if allele_id:
                by_allele[allele_id] = candidate_id
    if not order:
        raise RuntimeError("candidate panel contains no records")
    calls: dict[str, str] = {}
    outside = 0
    with open_text(native, "rt") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                raise RuntimeError("vg output contains a malformed record")
            candidate_id = fields[2] if fields[2] in records else None
            if candidate_id is None:
                allele_id = _allele_id(fields[7])
                candidate_id = by_allele.get(allele_id) if allele_id else None
            if candidate_id is None:
                candidate_id = by_key.get((fields[0], fields[1], fields[3], fields[4]))
            if candidate_id is None:
                outside += 1
                continue
            if candidate_id in calls:
                raise RuntimeError(f"vg output duplicates {candidate_id}")
            calls[candidate_id] = _genotype(fields)
    destination.parent.mkdir(parents=True, exist_ok=True)
    with open_text(destination, "wt") as output:
        output.write("##fileformat=VCFv4.2\n")
        for line in meta:
            if not line.startswith("##fileformat") and not line.startswith("##FORMAT=<ID=GT,"):
                output.write(line)
        output.write("##source=PGBench-vg-all-sites-adapter\n")
        output.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        output.write(f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}\n")
        for candidate_id in order:
            output.write("\t".join(records[candidate_id][:8] + ["GT", calls.get(candidate_id, "./.")]) + "\n")
    return len(calls), len(order) - len(calls), outside


def main() -> int:
    reads_r1 = Path(required("PGBENCH_INPUT_FASTQ_R1"))
    reads_r2 = Path(required("PGBENCH_INPUT_FASTQ_R2"))
    graph_dir = Path(required("PGBENCH_GRAPH_DIR"))
    output_dir = Path(required("PGBENCH_OUTPUT_DIR"))
    output_vcf = Path(required("PGBENCH_OUTPUT_VCF"))
    sample = required("PGBENCH_SAMPLE_ID")
    candidates = Path(required("PGBENCH_CANDIDATE_VCF"))
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
    native_vcf = work / "native.calls.vcf"

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
        stdout=native_vcf,
    )
    matched, no_call, outside = project_all_sites(native_vcf, candidates, output_vcf, sample)
    print(f"vg candidate projection: matched={matched} no_call={no_call} outside={outside}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
