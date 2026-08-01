#!/usr/bin/env python3
"""PGBench adapter for PanGenie short-read pangenome genotyping."""

from __future__ import annotations

import gzip
import os
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import TextIO


def required(name: str) -> str:
    value = os.environ.get(name)
    if not value:
        raise RuntimeError(f"required environment variable is absent: {name}")
    return value


def run(command: list[str]) -> None:
    subprocess.run(command, check=True)


def stage_uncompressed(source: Path, destination: Path) -> Path:
    """Materialize an uncompressed input because PanGenie rejects gzip input."""

    if source.suffix.casefold() != ".gz":
        return source
    with gzip.open(source, "rb") as input_handle, destination.open("wb") as output_handle:
        shutil.copyfileobj(input_handle, output_handle, length=16 * 1024 * 1024)
    return destination


def concatenate_paired_reads(read1: Path, read2: Path, destination: Path) -> Path:
    """Create the deterministic combined k-mer stream required by PanGenie."""

    with destination.open("wb") as output_handle:
        for source in (read1, read2):
            opener = gzip.open if source.suffix.casefold() == ".gz" else Path.open
            if opener is gzip.open:
                input_handle = opener(source, "rb")
            else:
                input_handle = opener(source, "rb")
            with input_handle:
                shutil.copyfileobj(
                    input_handle,
                    output_handle,
                    length=16 * 1024 * 1024,
                )
    return destination


def validate_pangenie_panel(path: Path) -> None:
    """Fail early when the panel violates PanGenie's documented VCF contract."""

    samples: list[str] | None = None
    last_end: dict[str, int] = {}
    record_count = 0
    with path.open("r", encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            if raw_line.startswith("##") or not raw_line.strip():
                continue
            fields = raw_line.rstrip("\n").split("\t")
            if raw_line.startswith("#CHROM"):
                samples = fields[9:]
                if not samples:
                    raise RuntimeError("PanGenie panel must contain phased panel samples")
                leaked = sorted(
                    sample
                    for sample in samples
                    if sample.casefold() in {"hg002", "na24385"}
                )
                if leaked:
                    raise RuntimeError(
                        "PanGenie panel contains benchmark sample aliases: "
                        + ", ".join(leaked)
                    )
                continue
            if samples is None:
                raise RuntimeError("PanGenie panel has no #CHROM header")
            if len(fields) != 9 + len(samples):
                raise RuntimeError(
                    f"PanGenie panel line {line_number} has inconsistent sample columns"
                )
            chrom, pos_raw, _, ref, alt = fields[:5]
            if any(token.startswith("<") and token.endswith(">") for token in alt.split(",")):
                raise RuntimeError(
                    f"PanGenie panel line {line_number} contains a symbolic ALT"
                )
            if "[" in alt or "]" in alt or alt in {"", "."}:
                raise RuntimeError(
                    f"PanGenie panel line {line_number} is not sequence-resolved"
                )
            try:
                pos = int(pos_raw)
            except ValueError as exc:
                raise RuntimeError(
                    f"PanGenie panel line {line_number} has invalid POS"
                ) from exc
            end = pos + max(len(ref), 1) - 1
            for item in fields[7].split(";"):
                if item.startswith("END="):
                    try:
                        end = int(item.split("=", 1)[1].split(",", 1)[0])
                    except ValueError as exc:
                        raise RuntimeError(
                            f"PanGenie panel line {line_number} has invalid END"
                        ) from exc
            if pos <= last_end.get(chrom, 0):
                raise RuntimeError(
                    f"PanGenie panel line {line_number} overlaps a previous record"
                )
            last_end[chrom] = end

            format_fields = fields[8].split(":")
            if "GT" not in format_fields:
                raise RuntimeError(
                    f"PanGenie panel line {line_number} has no GT field"
                )
            gt_index = format_fields.index("GT")
            for sample, sample_value in zip(samples, fields[9:], strict=True):
                values = sample_value.split(":")
                genotype = values[gt_index] if gt_index < len(values) else ""
                if "|" not in genotype or "." in genotype:
                    raise RuntimeError(
                        f"PanGenie panel line {line_number} has unphased/missing "
                        f"genotype for {sample}"
                    )
            record_count += 1
    if samples is None:
        raise RuntimeError("PanGenie panel has no #CHROM header")
    if record_count == 0:
        raise RuntimeError("PanGenie panel contains no variant records")


def open_text(path: Path, mode: str) -> TextIO:
    if path.suffix.casefold() == ".gz":
        return gzip.open(path, mode, encoding="utf-8")
    return path.open(mode, encoding="utf-8")


def _allele_id(info: str) -> str | None:
    for item in info.split(";"):
        if item.startswith("PANGENOME_ALLELE_ID="):
            return item.split("=", 1)[1]
    return None


def remap_to_candidate_ids(
    generated: Path,
    candidate_vcf: Path,
    destination: Path,
) -> None:
    """Restore the blinded CAND_* namespace required by the benchmark."""

    by_allele: dict[str, str] = {}
    by_record: dict[tuple[str, str, str, str], str] = {}
    with open_text(candidate_vcf, "rt") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8 or not fields[2].startswith("CAND_"):
                raise RuntimeError(
                    f"candidate VCF line {line_number} has an invalid record"
                )
            key = (fields[0], fields[1], fields[3], fields[4])
            allele_id = _allele_id(fields[7])
            if allele_id:
                if allele_id in by_allele:
                    raise RuntimeError(
                        f"candidate VCF has duplicate allele ID {allele_id}"
                    )
                by_allele[allele_id] = fields[2]
            if key in by_record:
                raise RuntimeError(f"candidate VCF has duplicate record key {key}")
            by_record[key] = fields[2]
    if not by_record:
        raise RuntimeError("candidate VCF contains no records")

    observed: set[str] = set()
    with generated.open("r", encoding="utf-8") as source, destination.open(
        "w", encoding="utf-8"
    ) as output:
        for line_number, line in enumerate(source, start=1):
            if not line.strip() or line.startswith("#"):
                output.write(line)
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                raise RuntimeError(
                    f"PanGenie output line {line_number} is malformed"
                )
            allele_id = _allele_id(fields[7])
            candidate_id = by_allele.get(allele_id) if allele_id else None
            if candidate_id is None:
                candidate_id = by_record.get(
                    (fields[0], fields[1], fields[3], fields[4])
                )
            if candidate_id is None:
                raise RuntimeError(
                    "PanGenie output contains a record outside the blinded "
                    f"candidate universe at line {line_number}"
                )
            if candidate_id in observed:
                raise RuntimeError(
                    f"PanGenie output duplicates candidate {candidate_id}"
                )
            observed.add(candidate_id)
            fields[2] = candidate_id
            output.write("\t".join(fields) + "\n")
    missing = sorted(set(by_record.values()) - observed)
    if missing:
        raise RuntimeError(
            f"PanGenie output omitted {len(missing)} blinded candidate records"
        )


def main() -> int:
    read1 = Path(required("PGBENCH_INPUT_FASTQ_R1"))
    read2 = Path(required("PGBENCH_INPUT_FASTQ_R2"))
    reference = Path(required("PGBENCH_REFERENCE_FASTA"))
    panel = Path(required("PGBENCH_PANEL_VCF"))
    candidate = Path(required("PGBENCH_CANDIDATE_VCF"))
    output = Path(required("PGBENCH_OUTPUT_VCF"))
    sample = required("PGBENCH_SAMPLE_ID")
    threads = required("PGBENCH_THREADS")

    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(
        prefix="pangenie-",
        dir=required("TMPDIR"),
    ) as temporary_name:
        temporary = Path(temporary_name)
        staged_reads = concatenate_paired_reads(
            read1,
            read2,
            temporary / "reads.fastq",
        )
        staged_panel = stage_uncompressed(panel, temporary / "panel.vcf")
        staged_reference = stage_uncompressed(reference, temporary / "reference.fa")
        validate_pangenie_panel(staged_panel)
        index_prefix = temporary / "panel-index"
        result_prefix = temporary / "HG002"

        run(
            [
                "PanGenie-index",
                "-v",
                str(staged_panel),
                "-r",
                str(staged_reference),
                "-t",
                threads,
                "-o",
                str(index_prefix),
            ]
        )
        run(
            [
                "PanGenie",
                "-f",
                str(index_prefix),
                "-i",
                str(staged_reads),
                "-s",
                sample,
                "-j",
                threads,
                "-t",
                threads,
                "-u",
                "-o",
                str(result_prefix),
            ]
        )
        generated = Path(f"{result_prefix}_genotyping.vcf")
        if not generated.is_file():
            raise RuntimeError(f"PanGenie did not create expected VCF: {generated}")
        remapped = temporary / "candidate-genotypes.vcf"
        remap_to_candidate_ids(generated, candidate, remapped)
        os.replace(remapped, output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
