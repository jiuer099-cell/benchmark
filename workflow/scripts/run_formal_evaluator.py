#!/usr/bin/env python3
"""Run one frozen formal evaluator and emit an auditable semantic ledger."""

from __future__ import annotations

import argparse
import bisect
import csv
import gzip
import hashlib
import json
import os
import shutil
import subprocess
import sys
import tempfile
from collections import defaultdict
from pathlib import Path
from typing import Any, Callable, TextIO

import yaml  # type: ignore[import-untyped]

from sv_matching import (
    SvMatchError,
    SvRecord,
    genotypes_equal,
    gt_state,
    load_evaluator_profile,
    load_vcf,
    one_to_one_match,
    query_record_in_universe,
    truth_record_in_universe,
)


class FormalEvaluatorError(RuntimeError):
    """Raised when an evaluator cannot produce a complete semantic ledger."""


LEDGER_FIELDS = [
    "result_id",
    "evaluator",
    "scope_eligible",
    "event_eligible",
    "truth_event_id",
    "evaluator_accept",
    "evaluator_resolved",
    "mapping_status",
    "detection_correct",
    "genotype_scorable",
    "genotype_correct",
    "no_call",
    "gt_state",
    "query_svtype",
    "query_svlen",
    "match_score",
    "start_distance",
    "end_distance",
    "size_similarity",
    "sequence_similarity",
    "correct",
]


CANONICAL_VCF_HEADER_LINES = {
    "INFO:SVTYPE": (
        '##INFO=<ID=SVTYPE,Number=1,Type=String,'
        'Description="Type of structural variant">\n'
    ),
    "INFO:END": (
        '##INFO=<ID=END,Number=1,Type=Integer,'
        'Description="End position of the structural variant">\n'
    ),
    "INFO:SVLEN": (
        '##INFO=<ID=SVLEN,Number=.,Type=Integer,'
        'Description="Difference in length between REF and ALT alleles">\n'
    ),
    "FORMAT:FT": (
        '##FORMAT=<ID=FT,Number=1,Type=String,'
        'Description="Sample-level genotype filter">\n'
    ),
}


def open_text(path: Path) -> TextIO:
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("r", encoding="utf-8")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def vcf_records(path: Path) -> list[list[str]]:
    records: list[list[str]] = []
    with open_text(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                raise FormalEvaluatorError(f"malformed VCF record in {path}")
            records.append(fields)
    return records


def record_key(fields: list[str]) -> tuple[str, int, str, str]:
    return fields[0], int(fields[1]), fields[3], fields[4]


def aardvark_record_key(fields: list[str]) -> tuple[str, int, str, str]:
    """Return the representation emitted by Aardvark's default VCF writer.

    Aardvark rebuilds output records from internal variants, drops input IDs,
    and removes matching suffix bases while both alleles remain anchored.
    Mirroring that transformation keeps every evaluator decision attached to
    its canonical query event without silently treating an unmatched row as FP.
    """
    chrom, pos, ref, alt = record_key(fields)
    if "," in alt:
        raise FormalEvaluatorError(
            "Aardvark formal query must be biallelic, but encountered "
            f"{chrom}:{pos} {ref}>{alt}"
        )
    while len(ref) > 1 and len(alt) > 1 and ref[-1] == alt[-1]:
        ref = ref[:-1]
        alt = alt[:-1]
    return chrom, pos, ref, alt


def vcfdist_record_keys(fields: list[str]) -> list[tuple[str, int, str, str]]:
    """Mirror vcfdist's internal variant representation used in query.tsv.

    vcfdist writes zero-based positions.  For length-changing alleles it trims
    their common prefix and suffix and advances the position by the prefix
    length; either emitted allele may consequently be empty.  Remaining
    complex alleles are emitted as separate insertion and deletion components.
    """
    chrom, vcf_pos, ref, alt = record_key(fields)
    if "," in alt:
        raise FormalEvaluatorError(
            "vcfdist formal query must be biallelic, but encountered "
            f"{chrom}:{vcf_pos} {ref}>{alt}"
        )
    position = vcf_pos - 1
    ref_length = len(ref)
    alt_length = len(alt)
    complex_variant = False
    if alt_length != ref_length:
        prefix_limit = min(ref_length, alt_length)
        prefix = 0
        while prefix < prefix_limit and ref[prefix] == alt[prefix]:
            prefix += 1

        ref_end = ref_length
        alt_end = alt_length
        while (
            ref_end > prefix
            and alt_end > prefix
            and ref[ref_end - 1] == alt[alt_end - 1]
        ):
            ref_end -= 1
            alt_end -= 1

        position += prefix
        ref = ref[prefix:ref_end]
        alt = alt[prefix:alt_end]
        complex_variant = bool(ref and alt)
    elif ref_length > 1 and ref[1:] == alt[1:]:
        ref = ref[0]
        alt = alt[0]
    elif ref_length > 1:
        complex_variant = True

    if complex_variant:
        return [
            (chrom, position, "", alt),
            (chrom, position, ref, ""),
        ]
    return [(chrom, position, ref, alt)]


def query_universe(
    path: Path,
    *,
    key_fn: Callable[[list[str]], tuple[str, int, str, str]] = record_key,
) -> tuple[list[str], dict[str, tuple[str, int, str, str]]]:
    order: list[str] = []
    keys: dict[str, tuple[str, int, str, str]] = {}
    for fields in vcf_records(path):
        result_id = fields[2]
        if not result_id or result_id == ".":
            raise FormalEvaluatorError("canonical query VCF contains an empty ID")
        if result_id in keys:
            raise FormalEvaluatorError(f"duplicate canonical query ID: {result_id}")
        order.append(result_id)
        keys[result_id] = key_fn(fields)
    return order, keys


def _resolve_votes(
    query: Path,
    decisions_by_id: dict[str, bool],
    decisions_by_key: dict[tuple[str, int, str, str], list[bool]],
    *,
    key_fn: Callable[[list[str]], tuple[str, int, str, str]] = record_key,
) -> tuple[list[str], dict[str, bool]]:
    order, keys = query_universe(query, key_fn=key_fn)
    votes: dict[str, bool] = {}
    for result_id in order:
        if result_id in decisions_by_id:
            votes[result_id] = decisions_by_id[result_id]
            continue
        key_votes = decisions_by_key.get(keys[result_id], [])
        if not key_votes:
            # Native evaluators may intentionally omit records they cannot
            # represent (for example calls on auxiliary contigs or very large
            # complex alleles).  Preserve that distinction in the formal
            # ledger as evaluator_resolved=0 instead of aborting the complete
            # benchmark.  Missing decisions are deliberately not converted to
            # false votes here.
            continue
        if len(set(key_votes)) != 1:
            raise FormalEvaluatorError(
                f"evaluator produced conflicting classifications for {result_id}"
            )
        votes[result_id] = key_votes[0]
    return order, votes


def parse_truvari(query: Path, artifacts: Path) -> tuple[list[str], dict[str, bool]]:
    by_id: dict[str, bool] = {}
    by_key: dict[tuple[str, int, str, str], list[bool]] = defaultdict(list)
    artifact_candidates = (
        (("tp-comp.vcf.gz", "tp-comp.vcf", "tp-call.vcf.gz", "tp-call.vcf"), True),
        (("fp.vcf.gz", "fp.vcf"), False),
    )
    for names, correct in artifact_candidates:
        path = next(
            (artifacts / name for name in names if (artifacts / name).is_file()),
            None,
        )
        if path is None:
            searched = ", ".join(str(artifacts / name) for name in names)
            raise FormalEvaluatorError(
                f"Truvari output is absent; searched: {searched}"
            )
        for fields in vcf_records(path):
            if fields[2] and fields[2] != ".":
                by_id[fields[2]] = correct
            by_key[record_key(fields)].append(correct)
    return _resolve_votes(query, by_id, by_key)


def _format_value(fields: list[str], name: str) -> str:
    if len(fields) < 10:
        raise FormalEvaluatorError("evaluator VCF has no sample column")
    keys = fields[8].split(":")
    values = fields[9].split(":")
    if name not in keys:
        raise FormalEvaluatorError(f"evaluator VCF is missing FORMAT/{name}")
    index = keys.index(name)
    if index >= len(values):
        raise FormalEvaluatorError(f"evaluator VCF has no FORMAT/{name} value")
    return values[index]


def parse_aardvark(query: Path, artifacts: Path) -> tuple[list[str], dict[str, bool]]:
    path = artifacts / "query.vcf.gz"
    if not path.is_file():
        raise FormalEvaluatorError(f"Aardvark output is absent: {path}")
    by_id: dict[str, bool] = {}
    by_key: dict[tuple[str, int, str, str], list[bool]] = defaultdict(list)
    for fields in vcf_records(path):
        decision = _format_value(fields, "BD")
        if decision not in {"TP", "FP"}:
            raise FormalEvaluatorError(f"unexpected Aardvark BD={decision!r}")
        correct = decision == "TP"
        if fields[2] and fields[2] != ".":
            by_id[fields[2]] = correct
        by_key[aardvark_record_key(fields)].append(correct)
    return _resolve_votes(query, by_id, by_key, key_fn=aardvark_record_key)


def parse_vcfdist(
    query: Path,
    artifacts: Path,
    *,
    credit_threshold: float,
    diagnostics: dict[str, Any] | None = None,
) -> tuple[list[str], dict[str, bool]]:
    path = artifacts / "query.tsv"
    if not path.is_file():
        raise FormalEvaluatorError(f"vcfdist output is absent: {path}")
    credits: dict[tuple[str, int, str, str], list[float]] = defaultdict(list)
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"CONTIG", "POS", "REF", "ALT", "CREDIT"}
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise FormalEvaluatorError("vcfdist query.tsv has an unsupported header")
        for row in reader:
            key = (
                row["CONTIG"],
                int(row["POS"]),
                row["REF"],
                row["ALT"],
            )
            credits[key].append(float(row["CREDIT"]))
    order: list[str] = []
    votes: dict[str, bool] = {}
    status_counts: dict[str, int] = defaultdict(int)
    key_owners: dict[tuple[str, int, str, str], str] = {}
    for fields in vcf_records(query):
        result_id = fields[2]
        if not result_id or result_id == ".":
            raise FormalEvaluatorError("canonical query VCF contains an empty ID")
        if result_id in votes:
            raise FormalEvaluatorError(f"duplicate canonical query ID: {result_id}")

        component_values: list[float] = []
        expected_components = vcfdist_record_keys(fields)
        mapped_components = 0
        for key in expected_components:
            owner = key_owners.setdefault(key, result_id)
            if owner != result_id:
                raise FormalEvaluatorError(
                    "vcfdist normalization ambiguously maps query results "
                    f"{owner} and {result_id} to {key}"
                )
            values = credits.get(key, [])
            if not values:
                continue
            else:
                mapped_components += 1
                component_values.extend(values)

        order.append(result_id)
        if mapped_components == len(expected_components):
            status_counts["exact"] += 1
        elif mapped_components > 0:
            # vcfdist can absorb one normalized component of a complex allele
            # into its supercluster representation. The remaining component
            # is still an explicit evaluator decision for the source event.
            status_counts["partial_complex"] += 1
        else:
            status_counts["unresolved"] += 1
            continue
        votes[result_id] = sum(component_values) / len(component_values) >= credit_threshold
    if diagnostics is not None:
        resolved = len(votes)
        diagnostics.update(
            {
                "algorithm": "vcfdist_component_event_mapping_v2",
                "query_events": len(order),
                "resolved_events": resolved,
                "unresolved_events": len(order) - resolved,
                "mapping_coverage": resolved / len(order) if order else 1.0,
                "status_counts": dict(sorted(status_counts.items())),
            }
        )
    return order, votes


def parse_vcfdist_native_summary(artifacts: Path) -> dict[str, Any]:
    path = artifacts / "precision-recall-summary.tsv"
    if not path.is_file():
        raise FormalEvaluatorError(f"vcfdist native summary is absent: {path}")
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if row.get("VAR_TYPE") == "ALL" and row.get("THRESHOLD") == "NONE":
                return {
                    "unit": "vcfdist_native_haplotype_components",
                    "truth_true_positive": int(row["TRUTH_TP"]),
                    "query_true_positive": int(row["QUERY_TP"]),
                    "truth_false_negative": int(row["TRUTH_FN"]),
                    "query_false_positive": int(row["QUERY_FP"]),
                    "precision": float(row["PREC"]),
                    "recall": float(row["RECALL"]),
                    "f1": float(row["F1_SCORE"]),
                }
    raise FormalEvaluatorError("vcfdist native summary has no ALL/NONE row")


def load_regions(path: Path) -> dict[str, tuple[list[int], list[int]]]:
    raw: dict[str, list[tuple[int, int]]] = defaultdict(list)
    with path.open("r", encoding="utf-8") as handle:
        for number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                raise FormalEvaluatorError(f"malformed BED record at {path}:{number}")
            start, end = int(fields[1]), int(fields[2])
            if start < 0 or end <= start:
                raise FormalEvaluatorError(f"invalid BED interval at {path}:{number}")
            raw[fields[0]].append((start, end))
    merged: dict[str, tuple[list[int], list[int]]] = {}
    for chrom, intervals in raw.items():
        compact: list[list[int]] = []
        for start, end in sorted(intervals):
            if compact and start <= compact[-1][1]:
                compact[-1][1] = max(compact[-1][1], end)
            else:
                compact.append([start, end])
        merged[chrom] = (
            [interval[0] for interval in compact],
            [interval[1] for interval in compact],
        )
    if not merged:
        raise FormalEvaluatorError(f"benchmark BED contains no intervals: {path}")
    return merged


def fully_contained(
    regions: dict[str, tuple[list[int], list[int]]], record: SvRecord
) -> bool:
    indexed = regions.get(record.chrom)
    if indexed is None:
        return False
    starts, ends = indexed
    start = record.pos - 1
    end = max(record.pos, record.end)
    index = bisect.bisect_right(starts, start) - 1
    return index >= 0 and ends[index] >= end


def reference_contig_headers(reference: Path) -> list[str]:
    """Return vcfdist-compatible contig declarations from a FASTA index."""

    fai = Path(f"{reference}.fai")
    if not fai.is_file():
        raise FormalEvaluatorError(
            f"reference FASTA index is required for evaluator input: {fai}"
        )
    headers: list[str] = []
    seen: set[str] = set()
    with fai.open("r", encoding="utf-8") as handle:
        for index, line in enumerate(handle):
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 2 or not fields[0]:
                raise FormalEvaluatorError(f"malformed reference FASTA index: {fai}")
            contig = fields[0]
            if contig in seen:
                raise FormalEvaluatorError(
                    f"duplicate contig {contig!r} in reference FASTA index: {fai}"
                )
            seen.add(contig)
            try:
                length = int(fields[1])
            except ValueError as error:
                raise FormalEvaluatorError(
                    f"invalid length for contig {contig!r} in {fai}"
                ) from error
            if length <= 0:
                raise FormalEvaluatorError(
                    f"non-positive length for contig {contig!r} in {fai}"
                )
            headers.append(
                f"##contig=<ID={contig},length={length},IDX={index}>\n"
            )
    if not headers:
        raise FormalEvaluatorError(f"reference FASTA index contains no contigs: {fai}")
    return headers


def write_variant_query(
    source: Path,
    destination: Path,
    allowed_ids: set[str],
    *,
    reference: Path | None = None,
    phase_unphased_genotypes: bool = False,
    project_detection_genotypes: bool = False,
) -> None:
    """Write evaluator input with canonical definitions for fields we emit.

    Upstream callers sometimes omit INFO declarations or declare FORMAT/FT as
    an integer flag.  htslib may tolerate those records while formal
    evaluators reject them, so the benchmark owns and freezes these four
    definitions in its evaluator-facing VCF.  For ``variant_sites`` tools,
    the evaluator-facing copy may additionally represent a no-call as a
    heterozygous *detection proxy*.  This is solely for site matching: the
    canonical VCF is never changed and its no-call remains a separate
    genotype diagnostic in the score ledger.
    """

    replaced = set(CANONICAL_VCF_HEADER_LINES)
    contig_headers = reference_contig_headers(reference) if reference else []
    with open_text(source) as reader, destination.open("w", encoding="utf-8") as writer:
        for line in reader:
            if contig_headers and line.startswith("##contig=<"):
                # Rebuild all contig declarations from the exact reference used
                # by the evaluator. vcfdist requires both length and IDX.
                continue
            if line.startswith("##INFO=<ID=SVTYPE,"):
                replaced.discard("INFO:SVTYPE")
                writer.write(CANONICAL_VCF_HEADER_LINES["INFO:SVTYPE"])
                continue
            if line.startswith("##INFO=<ID=END,"):
                replaced.discard("INFO:END")
                writer.write(CANONICAL_VCF_HEADER_LINES["INFO:END"])
                continue
            if line.startswith("##INFO=<ID=SVLEN,"):
                replaced.discard("INFO:SVLEN")
                writer.write(CANONICAL_VCF_HEADER_LINES["INFO:SVLEN"])
                continue
            if line.startswith("##FORMAT=<ID=FT,"):
                replaced.discard("FORMAT:FT")
                writer.write(CANONICAL_VCF_HEADER_LINES["FORMAT:FT"])
                continue
            if line.startswith("#CHROM"):
                for key in sorted(replaced):
                    writer.write(CANONICAL_VCF_HEADER_LINES[key])
                replaced.clear()
                for header in contig_headers:
                    writer.write(header)
                if phase_unphased_genotypes:
                    writer.write(
                        '##pgbench_vcfdist_detection_phase="deterministic 0|1; '
                        'phase is ignored for detection voting"\n'
                    )
                if project_detection_genotypes:
                    writer.write(
                        '##pgbench_detection_genotype_projection="no-call '
                        'projected to heterozygous only in evaluator input; '
                        'canonical genotype semantics are unchanged"\n'
                    )
                writer.write(line)
                continue
            if line.startswith("#"):
                writer.write(line)
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 3 and fields[2] in allowed_ids:
                if (
                    phase_unphased_genotypes or project_detection_genotypes
                ) and len(fields) >= 10:
                    format_keys = fields[8].split(":")
                    sample_values = fields[9].split(":")
                    if "GT" in format_keys:
                        gt_index = format_keys.index("GT")
                        if gt_index < len(sample_values):
                            gt = sample_values[gt_index]
                            if (
                                project_detection_genotypes
                                and gt in {".", "./.", ".|."}
                            ):
                                sample_values[gt_index] = (
                                    "0|1"
                                    if phase_unphased_genotypes
                                    else "0/1"
                                )
                                gt = sample_values[gt_index]
                            if phase_unphased_genotypes and gt in {"0/1", "1/0"}:
                                sample_values[gt_index] = "0|1"
                            fields[9] = ":".join(sample_values)
                writer.write("\t".join(fields) + "\n")


def command_prefix(config: dict[str, Any], evaluator: str) -> list[str]:
    value = config["evaluation"]["commands"][evaluator]
    if not isinstance(value, list) or not all(
        isinstance(item, str) and item for item in value
    ):
        raise FormalEvaluatorError(f"invalid command prefix for {evaluator}")
    command = list(value)
    # `mamba run -n ...` may contend on shared wrapper/lock state when the
    # three formal evaluators launch concurrently. Resolve the conventional
    # named-environment form to its immutable executable whenever available.
    if len(command) >= 5 and command[:3] == ["mamba", "run", "-n"]:
        environment_name = command[3]
        executable = command[4]
        roots = []
        configured_root = os.environ.get("MAMBA_ROOT_PREFIX")
        if configured_root:
            roots.append(Path(configured_root))
        roots.append(Path.home() / ".local" / "share" / "mamba")
        for root in roots:
            candidate = root / "envs" / environment_name / "bin" / executable
            if candidate.is_file():
                return [str(candidate.resolve()), *command[5:]]
    return command


def build_command(
    evaluator: str,
    prefix: list[str],
    extra_args: list[str],
    *,
    query: Path,
    truth: Path,
    reference: Path,
    regions: Path,
    artifacts: Path,
    threads: int,
) -> list[str]:
    if evaluator == "truvari":
        return [
            *prefix,
            "bench",
            "-b",
            str(truth),
            "-c",
            str(query),
            "-f",
            str(reference),
            "--includebed",
            str(regions),
            "-o",
            str(artifacts),
            *extra_args,
        ]
    if evaluator == "aardvark":
        return [
            *prefix,
            "compare",
            "--reference",
            str(reference),
            "--truth-vcf",
            str(truth),
            "--query-vcf",
            str(query),
            "--regions",
            str(regions),
            "--output-dir",
            str(artifacts),
            "--threads",
            str(threads),
            *extra_args,
        ]
    if evaluator == "vcfdist":
        return [
            *prefix,
            str(query),
            str(truth),
            str(reference),
            "-b",
            str(regions),
            "-p",
            str(artifacts) + os.sep,
            "--max-threads",
            str(threads),
            *extra_args,
        ]
    raise FormalEvaluatorError(f"unsupported evaluator: {evaluator}")


def evaluator_version(
    prefix: list[str], version_args: list[str]
) -> dict[str, Any]:
    command = [*prefix, *version_args]
    completed = subprocess.run(
        command,
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
    )
    text = completed.stdout.strip()
    if not text:
        raise FormalEvaluatorError("evaluator version command returned no text")
    return {
        "command": command,
        "output": text,
        "sha256": hashlib.sha256(text.encode("utf-8")).hexdigest(),
    }


def native_evaluator_cache_key(
    *,
    evaluator: str,
    prepared_query: Path,
    truth: Path,
    reference: Path,
    regions: Path,
    version_sha256: str,
    prefix: list[str],
    extra_args: list[str],
    threads: int,
) -> tuple[str, dict[str, Any]]:
    """Key expensive native execution independently from parser code.

    Parser/report changes must not rerun a multi-hour evaluator.  The cache is
    content-addressed by every biological input, the frozen executable
    fingerprint and native parameters; no parser implementation hash enters
    this key.
    """

    payload = {
        "schema_version": 1,
        "evaluator": evaluator,
        "prepared_query_sha256": sha256_file(prepared_query),
        "truth_sha256": sha256_file(truth),
        "reference_sha256": sha256_file(reference),
        "regions_sha256": sha256_file(regions),
        "version_sha256": version_sha256,
        "command_prefix": prefix,
        "extra_args": extra_args,
        "threads": threads,
    }
    canonical = json.dumps(payload, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(canonical.encode("utf-8")).hexdigest(), payload


def store_native_cache(source_artifacts: Path, source_log: Path, cache: Path, manifest: dict[str, Any]) -> None:
    cache.parent.mkdir(parents=True, exist_ok=True)
    temporary = cache.with_name(f".{cache.name}.tmp-{os.getpid()}")
    shutil.rmtree(temporary, ignore_errors=True)
    temporary.mkdir()
    shutil.copytree(source_artifacts, temporary / "artifacts")
    shutil.copy2(source_log, temporary / "evaluator.log")
    (temporary / "native-run.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    try:
        temporary.replace(cache)
    except FileExistsError:
        shutil.rmtree(temporary, ignore_errors=True)


def event_query_indices(
    queries: list[SvRecord],
    *,
    scope_eligible_ids: set[str],
    candidate_output_contract: str,
) -> list[int]:
    """Select submitted detection events without conflating GT and detection."""

    if candidate_output_contract not in {"all_sites", "variant_sites"}:
        raise FormalEvaluatorError(
            f"invalid candidate_output_contract: {candidate_output_contract!r}"
        )
    selected: list[int] = []
    for index, record in enumerate(queries):
        if record.record_id not in scope_eligible_ids:
            continue
        state = gt_state(record.gt)
        if candidate_output_contract == "variant_sites":
            if state == "hom_ref":
                raise FormalEvaluatorError(
                    "variant_sites output contains an in-scope hom-ref record: "
                    f"{record.record_id}"
                )
            # Discovery VCFs commonly omit GT. Presence in a variant-sites VCF
            # is the detection assertion; no-call remains a separate GT
            # diagnostic and must not erase the detection.
            selected.append(index)
        elif state == "variant":
            selected.append(index)
    return selected


def write_ledger(
    path: Path,
    evaluator: str,
    queries: list[SvRecord],
    truths: list[SvRecord],
    assignments: dict[int, Any],
    evaluator_votes: dict[str, bool],
    *,
    scope_eligible_ids: set[str],
    event_eligible_ids: set[str],
    require_phase: bool,
) -> dict[str, int]:
    counts = {
        "rows": 0,
        "scope_eligible": 0,
        "scope_excluded": 0,
        "events": 0,
        "detection_correct": 0,
        "evaluator_resolved": 0,
        "evaluator_unresolved": 0,
        "genotype_scorable": 0,
        "genotype_correct": 0,
        "no_call": 0,
    }
    seen_truth: set[str] = set()
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=LEDGER_FIELDS,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for index, query in enumerate(queries):
            state = gt_state(query.gt)
            scope_eligible = query.record_id in scope_eligible_ids
            event_eligible = query.record_id in event_eligible_ids
            if event_eligible and not scope_eligible:
                raise FormalEvaluatorError(
                    f"out-of-scope query was marked as a detection: {query.record_id}"
                )
            match = assignments.get(index) if event_eligible else None
            truth = truths[match.truth_index] if match is not None else None
            truth_id = truth.record_id if truth is not None else ""
            if truth_id:
                if truth_id in seen_truth:
                    raise FormalEvaluatorError(
                        f"truth event was assigned more than once: {truth_id}"
                    )
                seen_truth.add(truth_id)
            evaluator_resolved = (
                query.record_id in evaluator_votes if event_eligible else True
            )
            accepted = evaluator_votes.get(query.record_id, False) if event_eligible else False
            detection_correct = bool(accepted and truth is not None)
            genotype_scorable = truth is not None and state != "no_call"
            genotype_correct = bool(
                genotype_scorable
                and genotypes_equal(
                    query.gt,
                    truth.gt,
                    require_phase=require_phase,
                )
            )
            no_call = scope_eligible and state == "no_call"
            writer.writerow(
                {
                    "result_id": query.record_id,
                    "evaluator": evaluator,
                    "scope_eligible": int(scope_eligible),
                    "event_eligible": int(event_eligible),
                    "truth_event_id": truth_id,
                    "evaluator_accept": int(accepted),
                    "evaluator_resolved": int(evaluator_resolved),
                    "mapping_status": (
                        "resolved" if evaluator_resolved else "unresolved"
                    ),
                    "detection_correct": int(detection_correct),
                    "genotype_scorable": int(genotype_scorable),
                    "genotype_correct": int(genotype_correct),
                    "no_call": int(no_call),
                    "gt_state": state,
                    "query_svtype": query.svtype,
                    "query_svlen": query.svlen,
                    "match_score": "" if match is None else f"{match.score:.12g}",
                    "start_distance": "" if match is None else match.start_distance,
                    "end_distance": "" if match is None else match.end_distance,
                    "size_similarity": (
                        "" if match is None else f"{match.size_similarity:.12g}"
                    ),
                    "sequence_similarity": (
                        ""
                        if match is None or match.sequence_similarity is None
                        else f"{match.sequence_similarity:.12g}"
                    ),
                    # Compatibility alias: the primary vote remains detection only.
                    "correct": int(detection_correct),
                }
            )
            counts["rows"] += 1
            counts["scope_eligible"] += int(scope_eligible)
            counts["scope_excluded"] += int(not scope_eligible)
            counts["events"] += int(event_eligible)
            counts["detection_correct"] += int(detection_correct)
            counts["evaluator_resolved"] += int(event_eligible and evaluator_resolved)
            counts["evaluator_unresolved"] += int(event_eligible and not evaluator_resolved)
            counts["genotype_scorable"] += int(genotype_scorable)
            counts["genotype_correct"] += int(genotype_correct)
            counts["no_call"] += int(no_call)
    return counts


def run_evaluator(args: argparse.Namespace) -> None:
    config = yaml.safe_load(args.config.read_text(encoding="utf-8"))
    if not isinstance(config, dict):
        raise FormalEvaluatorError("configuration must be a YAML mapping")
    profile = load_evaluator_profile(args.evaluator_profile)
    tool_manifest = yaml.safe_load(args.tool_manifest.read_text(encoding="utf-8"))
    if not isinstance(tool_manifest, dict) or not isinstance(
        tool_manifest.get("outputs"), dict
    ):
        raise FormalEvaluatorError("tool manifest has no outputs contract")
    candidate_output_contract = tool_manifest["outputs"].get(
        "candidate_output_contract"
    )
    evaluator_profile = profile["evaluators"][args.evaluator]
    prefix = command_prefix(config, args.evaluator)
    version = evaluator_version(prefix, list(evaluator_profile["version_args"]))
    if version["sha256"] != evaluator_profile["expected_version_sha256"]:
        raise FormalEvaluatorError(
            f"{args.evaluator} executable fingerprint differs from the frozen "
            "evaluator profile"
        )
    regions = load_regions(args.regions)
    all_queries = load_vcf(args.query, prefix="query", allow_empty=True)
    queries = [
        record for record in all_queries if fully_contained(regions, record)
    ]
    all_truths = load_vcf(args.truth, prefix="truth")
    outside_truth = [
        record.record_id
        for record in all_truths
        if not fully_contained(regions, record)
    ]
    if outside_truth:
        raise FormalEvaluatorError(
            "materialized truth contains events not fully contained in the "
            f"benchmark BED: {', '.join(outside_truth[:5])}"
        )
    out_of_scope_truth = [
        record.record_id
        for record in all_truths
        if not truth_record_in_universe(record, profile)
    ]
    if out_of_scope_truth:
        raise FormalEvaluatorError(
            "materialized truth violates the frozen SV universe: "
            f"{', '.join(out_of_scope_truth[:5])}"
        )
    truths = all_truths
    if not truths:
        raise FormalEvaluatorError("truth VCF has no events in the benchmark regions")
    if len({record.record_id for record in queries}) != len(queries):
        raise FormalEvaluatorError("query result IDs must be unique")
    scope_eligible_ids = {
        record.record_id
        for record in queries
        if query_record_in_universe(record, profile)
    }
    detection_indices = event_query_indices(
        queries,
        scope_eligible_ids=scope_eligible_ids,
        candidate_output_contract=str(candidate_output_contract),
    )
    event_queries = [queries[index] for index in detection_indices]
    event_assignments = one_to_one_match(event_queries, truths, profile)
    assignments = {
        detection_indices[event_index]: match
        for event_index, match in event_assignments.items()
    }
    event_ids = {queries[index].record_id for index in detection_indices}

    args.output_dir.parent.mkdir(parents=True, exist_ok=True)
    work = Path(
        tempfile.mkdtemp(prefix=f".{args.output_dir.name}.", dir=args.output_dir.parent)
    )
    command: list[str] = []
    try:
        evaluator_votes: dict[str, bool] = {}
        native_cache_hit = False
        native_cache_key_value: str | None = None
        mapping_diagnostics: dict[str, Any] = {
            "algorithm": "native_record_identity",
            "query_events": len(event_ids),
            "resolved_events": len(event_ids),
            "unresolved_events": 0,
            "mapping_coverage": 1.0,
            "status_counts": {"exact": len(event_ids)},
        }
        if event_ids:
            plain = work / "input" / "query.events.vcf"
            plain.parent.mkdir(parents=True)
            write_variant_query(
                args.query,
                plain,
                event_ids,
                reference=args.reference,
                phase_unphased_genotypes=args.evaluator == "vcfdist",
                project_detection_genotypes=(
                    str(candidate_output_contract) == "variant_sites"
                ),
            )
            prepared = work / "input" / "query.events.vcf.gz"
            bcftools = command_prefix(config, "bcftools")
            tabix = command_prefix(config, "tabix")
            subprocess.run(
                [*bcftools, "view", "-Oz", "-o", str(prepared), str(plain)],
                check=True,
            )
            subprocess.run([*tabix, "-f", "-p", "vcf", str(prepared)], check=True)
            artifacts = work / "artifacts"
            command = build_command(
                args.evaluator,
                prefix,
                list(evaluator_profile.get("extra_args", [])),
                query=prepared,
                truth=args.truth,
                reference=args.reference,
                regions=args.regions,
                artifacts=artifacts,
                threads=args.threads,
            )
            native_cache_key_value, native_manifest = native_evaluator_cache_key(
                evaluator=args.evaluator,
                prepared_query=prepared,
                truth=args.truth,
                reference=args.reference,
                regions=args.regions,
                version_sha256=version["sha256"],
                prefix=prefix,
                extra_args=list(evaluator_profile.get("extra_args", [])),
                threads=args.threads,
            )
            native_cache = (
                args.output_dir.parent
                / f".{args.evaluator}.native-cache"
                / native_cache_key_value
            )
            if (
                (native_cache / "native-run.json").is_file()
                and (native_cache / "artifacts").is_dir()
                and (native_cache / "evaluator.log").is_file()
            ):
                shutil.copytree(native_cache / "artifacts", artifacts)
                shutil.copy2(native_cache / "evaluator.log", work / "evaluator.log")
                native_cache_hit = True
            else:
                with (work / "evaluator.log").open("w", encoding="utf-8") as log:
                    try:
                        subprocess.run(
                            command,
                            check=True,
                            stdout=log,
                            stderr=subprocess.STDOUT,
                            text=True,
                        )
                    except subprocess.CalledProcessError as error:
                        log.flush()
                        evaluator_log = (work / "evaluator.log").read_text(
                            encoding="utf-8", errors="replace"
                        )
                        tail = "\n".join(evaluator_log.splitlines()[-80:])
                        raise FormalEvaluatorError(
                            f"{args.evaluator} exited with status {error.returncode}. "
                            f"Evaluator log tail:\n{tail or '<empty>'}"
                        ) from error
                store_native_cache(
                    artifacts,
                    work / "evaluator.log",
                    native_cache,
                    native_manifest,
                )
            if args.evaluator == "truvari":
                _, evaluator_votes = parse_truvari(prepared, artifacts)
            elif args.evaluator == "aardvark":
                _, evaluator_votes = parse_aardvark(prepared, artifacts)
            else:
                _, evaluator_votes = parse_vcfdist(
                    prepared,
                    artifacts,
                    credit_threshold=float(evaluator_profile["credit_threshold"]),
                    diagnostics=mapping_diagnostics,
                )
        else:
            (work / "evaluator.log").write_text(
                "No non-reference query events in benchmark regions; evaluator skipped.\n",
                encoding="utf-8",
            )

        ledger = work / "votes.tsv"
        counts = write_ledger(
            ledger,
            args.evaluator,
            queries,
            truths,
            assignments,
            evaluator_votes,
            scope_eligible_ids=scope_eligible_ids,
            event_eligible_ids=event_ids,
            require_phase=bool(
                profile["semantics"]["genotype"].get("require_phase", False)
            ),
        )
        complete = {
            "schema_version": 3,
            "evaluator": args.evaluator,
            "command": command,
            "version": version,
            "evaluator_profile_id": profile["profile"]["id"],
            "evaluator_profile_sha256": profile["_sha256"],
            "evaluator_profile_path": str(args.evaluator_profile),
            "query_sha256": sha256_file(args.query),
            "truth_sha256": sha256_file(args.truth),
            "benchmark_bed_sha256": sha256_file(args.regions),
            "reference_sha256": sha256_file(args.reference),
            "vote_ledger_sha256": sha256_file(ledger),
            "matching": {
                "algorithm": (
                    "deterministic_max_cardinality_max_weight_one_to_one_v1"
                ),
                "truth_assignments": len(assignments),
                "truth_ids_unique": True,
            },
            "event_mapping": mapping_diagnostics,
            "native_metrics": (
                parse_vcfdist_native_summary(artifacts)
                if args.evaluator == "vcfdist" and event_ids
                else None
            ),
            "native_execution_cache": {
                "key": native_cache_key_value,
                "hit": native_cache_hit,
                "contract": "content_addressed_native_evaluator_v1",
            },
            "semantics": {
                "primary_vote": "detection",
                "candidate_output_contract": candidate_output_contract,
                "genotype_separate": True,
                "no_call_separate": True,
            },
            "universe": {
                **profile["universe"],
                "query_records_outside_regions": len(all_queries) - len(queries),
            },
            **counts,
        }
        (work / ".complete.json").write_text(
            json.dumps(complete, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        backup = args.output_dir.with_name(f".{args.output_dir.name}.previous")
        if backup.exists():
            shutil.rmtree(backup)
        if args.output_dir.exists():
            args.output_dir.replace(backup)
        try:
            work.replace(args.output_dir)
        except BaseException:
            if backup.exists() and not args.output_dir.exists():
                backup.replace(args.output_dir)
            raise
        shutil.rmtree(backup, ignore_errors=True)
    except BaseException:
        # Preserve the newest failed evaluator workspace for diagnosis.  Formal
        # evaluators can run for hours, and deleting query.tsv/evaluator.log on
        # failure makes a reproducible parser or tool incompatibility needlessly
        # expensive to investigate.  The hidden sibling is replaced atomically
        # by the next failure and is never treated as a completed result.
        failed = args.output_dir.with_name(f".{args.output_dir.name}.failed")
        shutil.rmtree(failed, ignore_errors=True)
        if work.exists():
            try:
                work.replace(failed)
            except OSError:
                shutil.rmtree(work, ignore_errors=True)
        raise


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--evaluator",
        required=True,
        choices=["truvari", "aardvark", "vcfdist"],
    )
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--evaluator-profile", required=True, type=Path)
    parser.add_argument("--tool-manifest", required=True, type=Path)
    parser.add_argument("--query", required=True, type=Path)
    parser.add_argument("--truth", required=True, type=Path)
    parser.add_argument("--reference", required=True, type=Path)
    parser.add_argument("--regions", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--threads", type=int, default=4)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    try:
        run_evaluator(parse_args(argv))
    except (
        OSError,
        ValueError,
        KeyError,
        TypeError,
        subprocess.CalledProcessError,
        FormalEvaluatorError,
        SvMatchError,
    ) as exc:
        print(f"run_formal_evaluator: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
