#!/usr/bin/env python3
"""PGBench external adapter for Varigraph pangenome SV genotyping."""

from __future__ import annotations

import gzip
import os
import subprocess
from pathlib import Path
from typing import TextIO


BENCHMARK_ALIASES = {"hg002", "na24385"}


def required(name: str) -> str:
    value = os.environ.get(name)
    if not value:
        raise RuntimeError(f"required environment variable is absent: {name}")
    return value


def open_text(path: Path, mode: str) -> TextIO:
    if path.suffix.casefold() == ".gz":
        return gzip.open(path, mode, encoding="utf-8")
    return path.open(mode, encoding="utf-8")


def run(command: list[str], *, cwd: Path) -> None:
    print("+ " + " ".join(command), flush=True)
    subprocess.run(command, cwd=cwd, check=True)


def validate_panel_header(path: Path) -> list[str]:
    """Require a multi-sample population panel without benchmark leakage."""

    with open_text(path, "rt") as handle:
        for raw_line in handle:
            if not raw_line.startswith("#CHROM"):
                continue
            fields = raw_line.rstrip("\n").split("\t")
            samples = fields[9:]
            if not samples:
                raise RuntimeError("Varigraph population panel has no samples")
            leaked = sorted(
                sample for sample in samples if sample.casefold() in BENCHMARK_ALIASES
            )
            if leaked:
                raise RuntimeError(
                    "Varigraph population panel contains benchmark aliases: "
                    + ", ".join(leaked)
                )
            return samples
    raise RuntimeError("Varigraph population panel has no #CHROM header")


def _allele_id(info: str) -> str | None:
    for item in info.split(";"):
        if item.startswith("PANGENOME_ALLELE_ID="):
            return item.split("=", 1)[1]
    return None


def _genotype(fields: list[str]) -> str:
    if len(fields) < 10:
        return "./."
    format_keys = fields[8].split(":")
    if "GT" not in format_keys:
        return "./."
    index = format_keys.index("GT")
    values = fields[9].split(":")
    if index >= len(values) or values[index] in {"", "."}:
        return "./."
    return values[index]


def project_to_candidates(
    generated: Path,
    candidates: Path,
    destination: Path,
    sample: str,
) -> tuple[int, int, int]:
    """Emit the blinded all-sites universe with Varigraph GT or explicit no-call."""

    order: list[str] = []
    records: dict[str, list[str]] = {}
    by_allele: dict[str, str] = {}
    by_key: dict[tuple[str, str, str, str], str] = {}
    candidate_headers: list[str] = []
    with open_text(candidates, "rt") as handle:
        for line_number, line in enumerate(handle, start=1):
            if line.startswith("##"):
                candidate_headers.append(line)
                continue
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8 or not fields[2].startswith("CAND_"):
                raise RuntimeError(
                    f"candidate VCF line {line_number} has an invalid record"
                )
            candidate_id = fields[2]
            key = (fields[0], fields[1], fields[3], fields[4])
            if candidate_id in records or key in by_key:
                raise RuntimeError("candidate VCF contains a duplicate record")
            allele_id = _allele_id(fields[7])
            if allele_id:
                if allele_id in by_allele:
                    raise RuntimeError(
                        f"candidate VCF has duplicate allele ID {allele_id}"
                    )
                by_allele[allele_id] = candidate_id
            by_key[key] = candidate_id
            records[candidate_id] = fields
            order.append(candidate_id)
    if not order:
        raise RuntimeError("candidate VCF contains no records")

    genotypes: dict[str, str] = {}
    outside = 0
    with open_text(generated, "rt") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                raise RuntimeError(
                    f"Varigraph output line {line_number} is malformed"
                )
            allele_id = _allele_id(fields[7])
            candidate_id = by_allele.get(allele_id) if allele_id else None
            if candidate_id is None:
                candidate_id = by_key.get((fields[0], fields[1], fields[3], fields[4]))
            if candidate_id is None:
                outside += 1
                continue
            if candidate_id in genotypes:
                raise RuntimeError(f"Varigraph output duplicates {candidate_id}")
            genotypes[candidate_id] = _genotype(fields)

    destination.parent.mkdir(parents=True, exist_ok=True)
    matched = no_calls = 0
    with destination.open("w", encoding="utf-8") as output:
        output.write("##fileformat=VCFv4.2\n")
        for header in candidate_headers:
            if not header.startswith("##fileformat") and not header.startswith(
                "##FORMAT=<ID=GT,"
            ):
                output.write(header)
        output.write("##source=PGBench-Varigraph-1.0.8-adapter\n")
        output.write(
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        )
        output.write(
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
            + sample
            + "\n"
        )
        for candidate_id in order:
            genotype = genotypes.get(candidate_id, "./.")
            if genotype in {".", "./.", ".|."}:
                no_calls += 1
            else:
                matched += 1
            output.write("\t".join(records[candidate_id][:8] + ["GT", genotype]) + "\n")
    return matched, outside, no_calls


def main() -> int:
    read1 = Path(required("PGBENCH_INPUT_FASTQ_R1"))
    read2 = Path(required("PGBENCH_INPUT_FASTQ_R2"))
    reference = Path(required("PGBENCH_REFERENCE_FASTA"))
    reference_index = Path(required("PGBENCH_REFERENCE_INDEX"))
    panel = Path(required("PGBENCH_PANEL_VCF"))
    candidates = Path(required("PGBENCH_CANDIDATE_VCF"))
    output = Path(required("PGBENCH_OUTPUT_VCF"))
    output_dir = Path(required("PGBENCH_OUTPUT_DIR"))
    sample = required("PGBENCH_SAMPLE_ID")
    threads = required("PGBENCH_THREADS")

    expected_index = Path(f"{reference}.fai")
    if reference_index.resolve() != expected_index.resolve():
        raise RuntimeError(
            "Varigraph requires the .fai adjacent to the registered reference"
        )
    validate_panel_header(panel)

    work = output_dir / "work"
    work.mkdir(parents=True, exist_ok=True)
    graph = work / "graph.bin"
    marker = work / "graph.bin.complete"
    if not graph.is_file() or not marker.is_file():
        graph.unlink(missing_ok=True)
        marker.unlink(missing_ok=True)
        run(
            [
                "varigraph",
                "construct",
                "-r",
                str(reference),
                "-v",
                str(panel),
                "--save-graph",
                str(graph),
                "-t",
                threads,
            ],
            cwd=work,
        )
        if not graph.is_file() or graph.stat().st_size == 0:
            raise RuntimeError("Varigraph construct produced no graph index")
        marker.write_text("varigraph=1.0.8\n", encoding="utf-8")

    samples = work / "samples.cfg"
    samples.write_text(
        f"{sample} {read1.resolve()} {read2.resolve()}\n", encoding="utf-8"
    )
    generated = work / f"{sample}.varigraph.vcf.gz"
    generated.unlink(missing_ok=True)
    run(
        [
            "varigraph",
            "genotype",
            "--load-graph",
            str(graph),
            "-s",
            str(samples),
            "--sv",
            "--use-depth",
            "-t",
            threads,
        ],
        cwd=work,
    )
    if not generated.is_file() or generated.stat().st_size == 0:
        raise RuntimeError(f"Varigraph produced no output at {generated}")
    matched, outside, no_calls = project_to_candidates(
        generated, candidates, output, sample
    )
    print(
        "Varigraph candidate projection: "
        f"matched={matched} outside_blinded_universe={outside} no_call={no_calls}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
