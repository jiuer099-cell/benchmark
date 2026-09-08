#!/usr/bin/env python3
"""PGBench adapter for PanGenie short-read pangenome genotyping."""

from __future__ import annotations

import gzip
import os
import shutil
import subprocess
import sys
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


def filter_pangenie_panel(source: Path, destination: Path) -> tuple[int, int]:
    """Drop records that cannot be indexed without inventing panel genotypes.

    PanGenie requires every population-panel genotype to be complete and phased.
    Population releases can contain a small number of missing or unphased calls.
    Imputing those calls would leak an unsupported assumption into the benchmark,
    so the formal adapter excludes the affected record from the private index and
    later emits an explicit no-call for its blinded candidate.
    """

    samples: list[str] | None = None
    kept = 0
    dropped = 0
    with source.open("r", encoding="utf-8") as input_handle, destination.open(
        "w", encoding="utf-8"
    ) as output_handle:
        for line_number, raw_line in enumerate(input_handle, start=1):
            if raw_line.startswith("##") or not raw_line.strip():
                output_handle.write(raw_line)
                continue
            fields = raw_line.rstrip("\n").split("\t")
            if raw_line.startswith("#CHROM"):
                samples = fields[9:]
                output_handle.write(raw_line)
                continue
            if samples is None:
                raise RuntimeError("PanGenie panel has no #CHROM header")
            if len(fields) != 9 + len(samples):
                raise RuntimeError(
                    f"PanGenie panel line {line_number} has inconsistent sample columns"
                )
            format_fields = fields[8].split(":")
            if "GT" not in format_fields:
                raise RuntimeError(
                    f"PanGenie panel line {line_number} has no GT field"
                )
            gt_index = format_fields.index("GT")
            usable = True
            for sample_value in fields[9:]:
                values = sample_value.split(":")
                genotype = values[gt_index] if gt_index < len(values) else ""
                if "|" not in genotype or "." in genotype:
                    usable = False
                    break
            if usable:
                output_handle.write(raw_line)
                kept += 1
            else:
                dropped += 1
    if kept == 0:
        raise RuntimeError(
            "PanGenie panel has no records with complete phased genotypes"
        )
    return kept, dropped


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
) -> tuple[int, int, int]:
    """Restore the blinded CAND_* namespace required by the benchmark.

    PanGenie genotypes the complete population panel, whereas the blinded
    challenge is deliberately restricted to the benchmark universe.  Records
    outside that universe are filtered at this boundary and never enter the
    scored VCF.  Candidate records not emitted by PanGenie remain explicit
    no-calls rather than being imputed as reference.
    """

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

    generated_headers: list[str] = []
    generated_by_candidate: dict[str, list[str]] = {}
    outside_candidate_count = 0
    with generated.open("r", encoding="utf-8") as source:
        for line_number, line in enumerate(source, start=1):
            if not line.strip() or line.startswith("#"):
                if line.strip():
                    generated_headers.append(line)
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
                outside_candidate_count += 1
                continue
            if candidate_id in generated_by_candidate:
                raise RuntimeError(
                    f"PanGenie output duplicates candidate {candidate_id}"
                )
            fields[2] = candidate_id
            generated_by_candidate[candidate_id] = fields

    matched_count = len(generated_by_candidate)
    no_call_count = 0
    with destination.open("w", encoding="utf-8") as output:
        output.writelines(generated_headers)
        with open_text(candidate_vcf, "rt") as candidates:
            for line_number, line in enumerate(candidates, start=1):
                if not line.strip() or line.startswith("#"):
                    continue
                candidate_fields = line.rstrip("\n").split("\t")
                candidate_id = candidate_fields[2]
                generated_fields = generated_by_candidate.pop(candidate_id, None)
                if generated_fields is not None:
                    output.write("\t".join(generated_fields) + "\n")
                    continue
                if len(candidate_fields) < 8:
                    raise RuntimeError(
                        f"candidate VCF line {line_number} is malformed"
                    )
                # Preserve the blinded candidate and explicitly represent the
                # genotype as unavailable instead of guessing a reference call.
                output.write("\t".join(candidate_fields[:8] + ["GT", "./."]) + "\n")
                no_call_count += 1
    if generated_by_candidate:
        raise RuntimeError(
            "PanGenie output contains candidates absent from the blinded panel: "
            + ", ".join(sorted(generated_by_candidate)[:5])
        )
    return matched_count, outside_candidate_count, no_call_count


def main() -> int:
    read1 = Path(required("PGBENCH_INPUT_FASTQ_R1"))
    read2 = Path(required("PGBENCH_INPUT_FASTQ_R2"))
    reference = Path(required("PGBENCH_REFERENCE_FASTA"))
    formal = os.environ.get("PGBENCH_EXECUTION_PURPOSE") == "formal"
    private_phased_raw = os.environ.get("PGBENCH_PANGENIE_PRIVATE_PHASED_PANEL")
    private_biallelic_raw = os.environ.get("PGBENCH_PANGENIE_PRIVATE_BIALLELIC_PANEL")
    converter_raw = os.environ.get("PGBENCH_PANGENIE_BIALLELIC_CONVERTER")
    projection_raw = os.environ.get("PGBENCH_CANONICAL_ALLELE_PROJECTION")
    if formal and not all((private_phased_raw, private_biallelic_raw, converter_raw, projection_raw)):
        raise RuntimeError(
            "formal PanGenie requires a private phased panel, biallelic "
            "projection resources, and a canonical allele projection"
        )
    private_phased_panel = Path(private_phased_raw or required("PGBENCH_PANEL_VCF"))
    private_biallelic_panel = Path(private_biallelic_raw) if private_biallelic_raw else None
    biallelic_converter = Path(converter_raw) if converter_raw else None
    canonical_projection = Path(projection_raw) if projection_raw else None
    candidate = Path(required("PGBENCH_CANDIDATE_VCF"))
    output = Path(required("PGBENCH_OUTPUT_VCF"))
    sample = required("PGBENCH_SAMPLE_ID")
    threads = required("PGBENCH_THREADS")

    output.parent.mkdir(parents=True, exist_ok=True)
    if canonical_projection is not None and not canonical_projection.is_file():
        raise RuntimeError("PanGenie canonical allele projection is missing")
    if biallelic_converter is not None and not biallelic_converter.is_file():
        raise RuntimeError("PanGenie biallelic converter is missing")
    with tempfile.TemporaryDirectory(
        prefix="pangenie-",
        dir=required("TMPDIR"),
    ) as temporary_name:
        temporary = Path(temporary_name)
        # Validate the benchmark-owned population haplotypes before expanding
        # the (potentially hundreds-of-GB) paired-read input.  A source panel
        # with no usable fully phased genotypes is an input-contract failure,
        # not a condition that can be repaired by consuming target reads.
        staged_panel_source = stage_uncompressed(
            private_phased_panel, temporary / "private.phased.source.vcf"
        )
        staged_panel = temporary / "panel.vcf"
        kept_records, dropped_records = filter_pangenie_panel(
            staged_panel_source,
            staged_panel,
        )
        print(
            "PanGenie panel genotype filter: "
            f"kept={kept_records} dropped_to_no_call={dropped_records}",
            flush=True,
        )
        validate_pangenie_panel(staged_panel)
        staged_reads = concatenate_paired_reads(
            read1,
            read2,
            temporary / "reads.fastq",
        )
        staged_reference = stage_uncompressed(reference, temporary / "reference.fa")
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
        converted = generated
        if biallelic_converter is not None and private_biallelic_panel is not None:
            converted = temporary / "native-biallelic-genotypes.vcf"
            with generated.open("rb") as source, converted.open("wb") as converted_output:
                subprocess.run(
                    [sys.executable, str(biallelic_converter), str(private_biallelic_panel)],
                    stdin=source,
                    stdout=converted_output,
                    check=True,
                )
        remapped = temporary / "candidate-genotypes.vcf"
        matched, outside, no_call = remap_to_candidate_ids(
            converted, candidate, remapped
        )
        print(
            "PanGenie candidate projection: "
            f"matched={matched} filtered_outside={outside} no_call={no_call}",
            flush=True,
        )
        os.replace(remapped, output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
