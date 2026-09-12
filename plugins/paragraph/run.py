#!/usr/bin/env python3
"""PGBench adapter for Paragraph from a frozen shared BAM through all-sites VCF."""

from __future__ import annotations

import os
import subprocess
from pathlib import Path
from typing import TextIO


def required(name: str) -> str:
    value = os.environ.get(name)
    if not value:
        raise RuntimeError(f"required environment variable is absent: {name}")
    return value


def open_text(path: Path, mode: str) -> TextIO:
    import gzip
    return gzip.open(path, mode, encoding="utf-8") if path.name.endswith(".gz") else path.open(mode, encoding="utf-8")


def read_length_from_bam(path: Path) -> int:
    """Read a representative sequenced-read length without reopening FASTQ."""

    result = subprocess.run(
        ["samtools", "view", str(path)], capture_output=True, text=True, check=True
    )
    for line in result.stdout.splitlines():
        fields = line.split("\t")
        if len(fields) >= 10 and fields[9] not in {"", "*"}:
            return len(fields[9])
    raise RuntimeError("frozen shared BAM contains no sequenced reads")


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


def project_all_sites(native: Path, candidates: Path, destination: Path, sample: str) -> tuple[int, int]:
    """Project native Paragraph genotypes onto every frozen candidate."""

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
            order.append(candidate_id)
            records[candidate_id] = fields
            by_key[(fields[0], fields[1], fields[3], fields[4])] = candidate_id
            allele_id = _allele_id(fields[7])
            if allele_id:
                by_allele[allele_id] = candidate_id
    if not order:
        raise RuntimeError("candidate panel contains no records")
    calls: dict[str, str] = {}
    with open_text(native, "rt") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                raise RuntimeError("Paragraph output contains a malformed record")
            candidate_id = fields[2] if fields[2] in records else None
            if candidate_id is None:
                allele_id = _allele_id(fields[7])
                candidate_id = by_allele.get(allele_id) if allele_id else None
            if candidate_id is None:
                candidate_id = by_key.get((fields[0], fields[1], fields[3], fields[4]))
            if candidate_id is not None:
                if candidate_id in calls:
                    raise RuntimeError(f"Paragraph output duplicates {candidate_id}")
                calls[candidate_id] = _genotype(fields)
    destination.parent.mkdir(parents=True, exist_ok=True)
    with open_text(destination, "wt") as output:
        output.write("##fileformat=VCFv4.2\n")
        for line in meta:
            if not line.startswith("##fileformat") and not line.startswith("##FORMAT=<ID=GT,"):
                output.write(line)
        output.write("##source=PGBench-Paragraph-all-sites-adapter\n")
        output.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        output.write(f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}\n")
        for candidate_id in order:
            output.write("\t".join(records[candidate_id][:8] + ["GT", calls.get(candidate_id, "./.")]) + "\n")
    return len(calls), len(order) - len(calls)


def main() -> int:
    bam = Path(required("PGBENCH_SHARED_ALIGNMENT_BAM"))
    bai = Path(required("PGBENCH_SHARED_ALIGNMENT_BAI"))
    reference = Path(required("PGBENCH_REFERENCE_FASTA"))
    candidates = Path(required("PGBENCH_CANDIDATE_VCF"))
    sample = required("PGBENCH_SAMPLE_ID")
    threads = required("PGBENCH_THREADS")
    output_dir = Path(required("PGBENCH_OUTPUT_DIR")).resolve()
    output_vcf = Path(required("PGBENCH_OUTPUT_VCF")).resolve()
    output_vcf.relative_to(output_dir)
    work = output_dir / "paragraph-work"
    work.mkdir(parents=True, exist_ok=True)
    if not bai.is_file() or bai.stat().st_size == 0:
        raise RuntimeError("Paragraph requires the frozen shared BAM index")
    subprocess.run(["samtools", "quickcheck", "-v", str(bam)], check=True)
    coverage = subprocess.check_output(
        ["samtools", "coverage", str(bam)], text=True
    ).splitlines()
    depths = [float(line.split("\t")[6]) for line in coverage if line and not line.startswith("#")]
    depth = sum(depths) / len(depths) if depths else 1.0
    manifest = work / "sample.tsv"
    manifest.write_text(
        "id\tpath\tread length\tdepth\n"
        f"{sample}\t{bam}\t{read_length_from_bam(bam)}\t{depth:.6f}\n",
        encoding="utf-8",
    )
    native = work / "native"
    subprocess.run(
        [
            "multigrmpy.py", "-i", str(candidates), "-m", str(manifest),
            "-r", str(reference), "-o", str(native), "-t", threads,
            "-M", str(max(20, round(depth * 20))),
        ],
        check=True,
    )
    source = native / "genotypes.vcf.gz"
    if not source.is_file():
        raise RuntimeError("Paragraph did not create genotypes.vcf.gz")
    matched, no_call = project_all_sites(source, candidates, output_vcf, sample)
    print(f"Paragraph candidate projection: matched={matched} no_call={no_call}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
