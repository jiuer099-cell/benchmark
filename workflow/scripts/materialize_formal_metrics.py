#!/usr/bin/env python3
"""Fuse frozen formal evaluator ledgers without double-counting truth events."""

from __future__ import annotations

import argparse
import bisect
import csv
import gzip
import hashlib
import json
import math
import random
import statistics
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any, TextIO

import yaml  # type: ignore[import-untyped]

from build_pangenome_manifest import parse_info
from sv_matching import (
    SvRecord,
    genotypes_equal,
    gt_state,
    load_evaluator_profile,
    load_vcf,
    sha256_file,
    truth_record_in_universe,
)


class FormalMetricError(ValueError):
    """Raised when evaluator ledgers do not describe one identical universe."""


def load_json_object(path: Path, label: str) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise FormalMetricError(f"cannot load {label}: {path}") from exc
    if not isinstance(value, dict):
        raise FormalMetricError(f"{label} must contain a JSON object")
    return value


EVALUATORS = ("truvari", "aardvark", "vcfdist")
COMPARISON_TASKS = {"panel_genotyping"}
REQUIRED_CONTEXT_STRATA = frozenset(
    {
        "non_repeat",
        "tandem_repeat",
        "segmental_duplication",
        "low_complexity",
        "low_mappability",
        "other_difficult",
    }
)
PROJECT_ROOT = Path(__file__).resolve().parents[2]
EXTENDED_FIELDS = {
    "scope_eligible",
    "event_eligible",
    "truth_event_id",
    "detection_correct",
    "genotype_scorable",
    "genotype_correct",
    "no_call",
    "gt_state",
    "query_svtype",
    "query_svlen",
}


def open_text(path: Path) -> TextIO:
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("r", encoding="utf-8")


def semantic_asset_hash(
    path: Path,
    *,
    excluded_root_keys: frozenset[str] = frozenset(),
) -> str:
    """Hash a YAML/JSON resource contract without volatile filesystem metadata.

    Pangenome and graph lock manifests include absolute or run-scoped paths, and
    the pangenome manifest also records its generation time.  Those fields do
    not change the biological resource.  Removing only those presentation
    fields gives cross-run comparisons a stable content identity while all
    embedded checksums, profiles, exclusions, and sample counts remain locked.
    """

    try:
        payload = yaml.safe_load(path.read_text(encoding="utf-8"))
    except yaml.YAMLError as exc:
        raise FormalMetricError(
            f"cannot parse semantic asset manifest {path}: {exc}"
        ) from exc
    if not isinstance(payload, dict):
        raise FormalMetricError(
            f"semantic asset manifest must be a mapping: {path}"
        )

    volatile_keys = {"asset_root", "generated_at", "path"}

    def normalize(value: Any, *, root: bool = False) -> Any:
        if isinstance(value, dict):
            return {
                str(key): normalize(item)
                for key, item in sorted(value.items(), key=lambda pair: str(pair[0]))
                if key not in volatile_keys
                and (not root or key not in excluded_root_keys)
            }
        if isinstance(value, list):
            return [normalize(item) for item in value]
        return value

    canonical = json.dumps(
        normalize(payload, root=True),
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=False,
    ).encode("utf-8")
    return hashlib.sha256(canonical).hexdigest()


def comparison_track_contract(
    *,
    tool_manifest: Path,
    official_score_mode: str,
    primary_truth_profile: str,
    score_profile: Path,
    evaluator_profile_sha256: str,
    asset_hashes: dict[str, str | None],
    evidence: dict[str, Any] | None,
) -> tuple[dict[str, Any], str]:
    """Freeze the complete boundary that makes formal results comparable."""

    try:
        manifest = yaml.safe_load(tool_manifest.read_text(encoding="utf-8"))
    except (OSError, yaml.YAMLError) as exc:
        raise FormalMetricError(
            f"cannot load comparison contract from {tool_manifest}"
        ) from exc
    if not isinstance(manifest, dict):
        raise FormalMetricError("tool manifest must be a mapping")
    task = manifest.get("comparison_task")
    if task not in COMPARISON_TASKS:
        raise FormalMetricError(
            "tool manifest comparison_task must declare one supported task"
        )
    outputs = manifest.get("outputs")
    if not isinstance(outputs, dict):
        raise FormalMetricError("tool manifest has no outputs contract")
    candidate_contract = outputs.get("candidate_output_contract")
    if candidate_contract != "all_sites":
        raise FormalMetricError("invalid candidate_output_contract")
    if not isinstance(evidence, dict):
        raise FormalMetricError(
            "formal comparison track requires a frozen evidence profile"
        )
    if evidence.get("official_score_mode") != official_score_mode:
        raise FormalMetricError(
            "evidence profile official mode differs from the score mode"
        )
    actual_technology = evidence.get("actual_technology")
    input_evidence_sha256 = evidence.get("resolved_inputs_sha256")
    if not isinstance(actual_technology, str) or not actual_technology:
        raise FormalMetricError(
            "comparison track requires actual sequencing technology"
        )
    if (
        not isinstance(input_evidence_sha256, str)
        or len(input_evidence_sha256) != 64
    ):
        raise FormalMetricError(
            "comparison track requires resolved input evidence SHA-256"
        )

    evaluation_universe_sha256 = (
        asset_hashes["challenge_hidden_ledger"]
        if task == "panel_genotyping"
        else asset_hashes["truth_vcf"]
    )
    boundaries = {
        "contract": "pgbench_comparison_track_v1",
        "task": task,
        "official_score_mode": official_score_mode,
        "candidate_output_contract": candidate_contract,
        "actual_technology": actual_technology,
        "primary_truth_profile": primary_truth_profile,
        "score_profile_sha256": sha256_file(score_profile),
        "evaluator_profile_sha256": evaluator_profile_sha256,
        "reference_sha256": asset_hashes["reference"],
        "truth_vcf_sha256": asset_hashes["truth_vcf"],
        "benchmark_bed_sha256": asset_hashes["benchmark_bed"],
        "pangenome_manifest_sha256": asset_hashes["pangenome_manifest"],
        "tool_manifest_sha256": asset_hashes["tool_manifest"],
        "evaluation_universe_sha256": evaluation_universe_sha256,
        "input_evidence_sha256": input_evidence_sha256,
    }
    canonical = json.dumps(
        boundaries,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=False,
    ).encode("utf-8")
    digest = hashlib.sha256(canonical).hexdigest()
    track = {
        "id": (
            f"{task}.{official_score_mode}.{actual_technology}."
            f"{candidate_contract}.{digest[:12]}"
        ),
        **boundaries,
    }
    return track, digest


def load_benchmark_regions(path: Path) -> dict[str, tuple[list[int], list[int]]]:
    raw: dict[str, list[tuple[int, int]]] = defaultdict(list)
    with open_text(path) as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                raise FormalMetricError(
                    f"malformed BED record at {path}:{line_number}"
                )
            start, end = int(fields[1]), int(fields[2])
            if start < 0 or end <= start:
                raise FormalMetricError(
                    f"invalid BED interval at {path}:{line_number}"
                )
            raw[fields[0]].append((start, end))
    if not raw:
        raise FormalMetricError(f"benchmark BED contains no regions: {path}")
    merged: dict[str, tuple[list[int], list[int]]] = {}
    for contig, intervals in raw.items():
        compact: list[list[int]] = []
        for start, end in sorted(intervals):
            if compact and start <= compact[-1][1]:
                compact[-1][1] = max(compact[-1][1], end)
            else:
                compact.append([start, end])
        merged[contig] = (
            [interval[0] for interval in compact],
            [interval[1] for interval in compact],
        )
    return merged


def fully_contained_in_regions(
    regions: dict[str, tuple[list[int], list[int]]],
    contig: str,
    start: int,
    end: int,
) -> bool:
    interval_index = regions.get(contig)
    if interval_index is None:
        return False
    starts, ends = interval_index
    index = bisect.bisect_right(starts, start) - 1
    return index >= 0 and ends[index] >= end


def eligible_truth_records(
    truth_vcf: Path,
    benchmark_bed: Path,
    *,
    evaluator_profile: dict[str, Any] | None = None,
) -> list[SvRecord]:
    regions = load_benchmark_regions(benchmark_bed)
    all_records = load_vcf(truth_vcf, prefix="truth")
    records = [
        record
        for record in all_records
        if fully_contained_in_regions(
            regions,
            record.chrom,
            record.pos - 1,
            max(record.pos, record.end),
        )
    ]
    if evaluator_profile is not None:
        outside = [
            record.record_id
            for record in all_records
            if not fully_contained_in_regions(
                regions,
                record.chrom,
                record.pos - 1,
                max(record.pos, record.end),
            )
        ]
        if outside:
            raise FormalMetricError(
                "materialized truth contains events not fully contained in the "
                f"benchmark BED: {', '.join(outside[:5])}"
            )
        invalid = [
            record.record_id
            for record in records
            if not truth_record_in_universe(record, evaluator_profile)
        ]
        if invalid:
            raise FormalMetricError(
                "materialized truth violates the frozen evaluator universe: "
                f"{', '.join(invalid[:5])}"
            )
    if not records:
        raise FormalMetricError(
            "primary truth VCF has no records fully contained in the benchmark BED"
        )
    if len({record.record_id for record in records}) != len(records):
        raise FormalMetricError("truth event IDs are not unique")
    return records


def _binary(row: dict[str, str], field: str, path: Path) -> bool:
    value = row.get(field)
    if value not in {"0", "1"}:
        raise FormalMetricError(f"invalid binary {field} in {path}")
    return value == "1"


def load_ledger(
    path: Path, expected_evaluator: str
) -> tuple[dict[str, dict[str, Any]], bool]:
    rows: dict[str, dict[str, Any]] = {}
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"result_id", "evaluator", "correct", *EXTENDED_FIELDS}
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise FormalMetricError(f"unsupported vote ledger header: {path}")
        extended = True
        for raw in reader:
            if raw["evaluator"] != expected_evaluator:
                raise FormalMetricError(
                    f"{path} contains evaluator {raw['evaluator']!r}"
                )
            result_id = raw["result_id"]
            if not result_id or result_id in rows:
                raise FormalMetricError(
                    f"empty or duplicate result ID {result_id!r} in {path}"
                )
            rows[result_id] = {
                "scope_eligible": _binary(raw, "scope_eligible", path),
                "event_eligible": _binary(raw, "event_eligible", path),
                "truth_event_id": raw["truth_event_id"] or None,
                "evaluator_resolved": (
                    _binary(raw, "evaluator_resolved", path)
                    if raw.get("evaluator_resolved") not in {None, ""}
                    else True
                ),
                "mapping_status": raw.get("mapping_status", "resolved"),
                "detection_correct": _binary(raw, "detection_correct", path),
                "genotype_scorable": _binary(raw, "genotype_scorable", path),
                "genotype_correct": _binary(raw, "genotype_correct", path),
                "no_call": _binary(raw, "no_call", path),
                "gt_state": raw["gt_state"],
                "query_svtype": raw["query_svtype"],
                "query_svlen": int(raw["query_svlen"]),
                "correct": _binary(raw, "correct", path),
            }
            if rows[result_id]["correct"] != rows[result_id]["detection_correct"]:
                raise FormalMetricError(
                    f"primary vote is not detection semantics in {path}"
                )
    if not rows:
        raise FormalMetricError(f"empty vote ledger: {path}")
    truth_ids = [
        row["truth_event_id"]
        for row in rows.values()
        if row["truth_event_id"] is not None
    ]
    if len(truth_ids) != len(set(truth_ids)):
        raise FormalMetricError(
            f"{path} assigns one truth event to multiple query results"
        )
    return rows, extended


def provenance_bundle(paths: list[Path]) -> str:
    digest = hashlib.sha256()
    for path in sorted(paths, key=lambda item: str(item)):
        digest.update(str(path).encode("utf-8"))
        digest.update(path.read_bytes())
    return digest.hexdigest()


def content_seed_bundle(paths: list[tuple[str, Path]]) -> str:
    """Return a path-independent seed derived from the scored input contents."""

    digest = hashlib.sha256()
    for label, path in sorted(paths, key=lambda item: item[0]):
        digest.update(label.encode("utf-8"))
        digest.update(b"\0")
        digest.update(hashlib.sha256(path.read_bytes()).digest())
    return digest.hexdigest()


def load_completions(
    paths: list[Path], evaluator_profile_sha256: str
) -> dict[str, dict[str, Any]]:
    completions: dict[str, dict[str, Any]] = {}
    for path in paths:
        try:
            payload = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as exc:
            raise FormalMetricError(f"cannot read evaluator completion {path}") from exc
        evaluator = payload.get("evaluator")
        if evaluator not in EVALUATORS or evaluator in completions:
            raise FormalMetricError(f"invalid/duplicate evaluator completion: {path}")
        if payload.get("evaluator_profile_sha256") != evaluator_profile_sha256:
            raise FormalMetricError(
                f"{evaluator} used a different evaluator profile"
            )
        if payload.get("schema_version") not in {3, 4}:
            raise FormalMetricError(
                f"{evaluator} completion schema is not a supported frozen version"
            )
        version = payload.get("version")
        if not isinstance(version, dict) or not version.get("sha256"):
            raise FormalMetricError(f"{evaluator} version was not recorded")
        completions[evaluator] = payload
    if set(completions) != set(EVALUATORS):
        raise FormalMetricError("all three evaluator completion records are required")
    return completions


def _length_bin(length: int) -> str:
    value = abs(length)
    if value < 50:
        return "lt50"
    if value < 100:
        return "50_99"
    if value < 500:
        return "100_499"
    if value < 1000:
        return "500_999"
    if value < 10000:
        return "1000_9999"
    return "ge10000"


def load_hidden_truth_labels(path: Path) -> dict[str, dict[str, str]]:
    """Load candidate-level truth labels without exposing them to adapters."""

    rows: dict[str, dict[str, str]] = {}
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"candidate_id", "truth_label", "truth_event_id"}
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise FormalMetricError("hidden truth ledger has an invalid header")
        for row in reader:
            candidate_id = row["candidate_id"]
            if not candidate_id or candidate_id in rows:
                raise FormalMetricError(
                    f"empty or duplicate candidate in hidden truth ledger: {candidate_id!r}"
                )
            label = row["truth_label"]
            if label not in {"positive", "negative", "unscorable"}:
                raise FormalMetricError(
                    f"invalid truth label for {candidate_id}: {label!r}"
                )
            rows[candidate_id] = row
    if not rows:
        raise FormalMetricError("hidden truth ledger is empty")
    return rows


def load_addressability_rows(path: Path) -> dict[str, dict[str, str]]:
    rows: dict[str, dict[str, str]] = {}
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            "candidate_id",
            "tool_native_status",
            "output_status",
            "failure_reason",
        }
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise FormalMetricError("addressability TSV has an invalid header")
        for row in reader:
            candidate_id = row["candidate_id"]
            if not candidate_id or candidate_id in rows:
                raise FormalMetricError(
                    f"empty or duplicate candidate in addressability TSV: {candidate_id!r}"
                )
            rows[candidate_id] = row
    return rows


def load_candidate_af(path: Path) -> dict[str, float | None]:
    """Read frozen population AF from the canonical all-sites candidate VCF."""

    values: dict[str, float | None] = {}
    with open_text(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                raise FormalMetricError("canonical query VCF has an invalid row")
            candidate_id = fields[2]
            if not candidate_id or candidate_id == "." or candidate_id in values:
                raise FormalMetricError(
                    f"invalid/duplicate canonical candidate ID: {candidate_id!r}"
                )
            info = parse_info(fields[7])
            raw = info.get("PANEL_AF", info.get("AF"))
            if raw is None or raw is True or str(raw) in {"", "."}:
                values[candidate_id] = None
                continue
            try:
                value = float(str(raw).split(",", 1)[0])
            except ValueError as exc:
                raise FormalMetricError(
                    f"invalid population AF for {candidate_id}: {raw!r}"
                ) from exc
            if not math.isfinite(value) or not 0.0 <= value <= 1.0:
                raise FormalMetricError(
                    f"population AF for {candidate_id} is outside [0, 1]"
                )
            values[candidate_id] = value
    return values


def _af_bin(value: float | None) -> str:
    if value is None or value <= 0.0:
        return "unknown"
    if value < 0.001:
        return "ultrarare"
    if value < 0.01:
        return "rare"
    if value < 0.05:
        return "low_frequency"
    return "common"


def _overlaps_regions(
    regions: dict[str, tuple[list[int], list[int]]], record: SvRecord
) -> bool:
    interval_index = regions.get(record.chrom)
    if interval_index is None:
        return False
    starts, ends = interval_index
    event_start = record.pos - 1
    event_end = max(record.pos, record.end)
    index = bisect.bisect_right(starts, event_end - 1) - 1
    return index >= 0 and ends[index] > event_start


def _stratum_metrics(
    *,
    candidate_ids: set[str],
    positive_ids: set[str],
    event_ids: set[str],
    ledgers: dict[str, dict[str, dict[str, Any]]],
    addressability: dict[str, dict[str, str]],
) -> dict[str, Any]:
    truth_ids = candidate_ids & positive_ids
    query_ids = candidate_ids & event_ids
    output_counts: dict[str, int] = defaultdict(int)
    native_counts: dict[str, int] = defaultdict(int)
    for candidate_id in candidate_ids:
        row = addressability[candidate_id]
        output_counts[row["output_status"]] += 1
        native_counts[row["tool_native_status"]] += 1
    result: dict[str, Any] = {
        "status": "defined" if truth_ids else "excluded_zero_truth",
        "canonical_candidate_count": len(candidate_ids),
        "truth_positive_count": len(truth_ids),
        "query_event_count": len(query_ids),
        "called_count": output_counts["called"],
        "explicit_no_call_count": output_counts["explicit_no_call"],
        "missing_output_count": output_counts["missing_output"],
        "unsupported_representation_count": native_counts[
            "unsupported_representation"
        ],
        "linking_failure_count": native_counts["linking_failure"],
        "addressable_count": native_counts["addressable"],
        "addressability_rate": (
            native_counts["addressable"] / len(candidate_ids)
            if candidate_ids
            else None
        ),
        "evaluator_metrics": {},
        "me_f1": None,
    }
    if not truth_ids:
        return result
    for evaluator in EVALUATORS:
        tp = sum(
            int(
                ledgers[evaluator][candidate_id]["detection_correct"]
                and ledgers[evaluator][candidate_id]["genotype_scorable"]
                and ledgers[evaluator][candidate_id]["genotype_correct"]
            )
            for candidate_id in query_ids
        )
        fp = len(query_ids) - tp
        fn = len(truth_ids) - tp
        precision = tp / (tp + fp) if tp + fp else 0.0
        recall = tp / (tp + fn)
        f1 = (
            2.0 * precision * recall / (precision + recall)
            if precision + recall
            else 0.0
        )
        result["evaluator_metrics"][evaluator] = {
            "tp": tp,
            "fp": fp,
            "fn": fn,
            "precision": precision,
            "recall": recall,
            "f1": f1,
        }
    result["me_f1"] = (
        statistics.fmean(
            float(metrics["f1"])
            for metrics in result["evaluator_metrics"].values()
        )
        * 100.0
    )
    return result


def genotype_stratified_summary(
    *,
    query_records: dict[str, SvRecord],
    hidden_truth: dict[str, dict[str, str]],
    ledgers: dict[str, dict[str, dict[str, Any]]],
    event_ids: set[str],
    candidate_af: dict[str, float | None],
    context_beds: dict[str, Path],
    addressability: dict[str, dict[str, str]],
) -> dict[str, Any]:
    candidate_ids = set(hidden_truth)
    if set(query_records) != candidate_ids:
        raise FormalMetricError(
            "canonical query VCF and hidden truth ledger candidate sets differ"
        )
    if set(addressability) != candidate_ids:
        raise FormalMetricError(
            "addressability TSV and canonical candidate sets differ"
        )
    if not event_ids.issubset(candidate_ids):
        raise FormalMetricError("evaluator ledger contains off-panel result IDs")
    if set(candidate_af) != candidate_ids:
        raise FormalMetricError("population AF and canonical candidate sets differ")
    positive_ids = {
        candidate_id
        for candidate_id, row in hidden_truth.items()
        if row["truth_label"] == "positive"
    }
    dimensions: dict[str, dict[str, set[str]]] = {
        "svtype": defaultdict(set),
        "length_bin": defaultdict(set),
        "population_af": defaultdict(set),
        "genome_context": defaultdict(set),
    }
    contexts = {
        name: load_benchmark_regions(path)
        for name, path in sorted(context_beds.items())
    }
    for name in contexts:
        dimensions["genome_context"][name] = set()
    for candidate_id, record in query_records.items():
        dimensions["svtype"][record.svtype].add(candidate_id)
        dimensions["length_bin"][_length_bin(record.svlen)].add(candidate_id)
        dimensions["population_af"][_af_bin(candidate_af[candidate_id])].add(
            candidate_id
        )
        for context, regions in contexts.items():
            if _overlaps_regions(regions, record):
                dimensions["genome_context"][context].add(candidate_id)
    return {
        dimension: {
            stratum: _stratum_metrics(
                candidate_ids=members,
                positive_ids=positive_ids,
                event_ids=event_ids,
                ledgers=ledgers,
                addressability=addressability,
            )
            for stratum, members in sorted(strata.items())
        }
        for dimension, strata in dimensions.items()
    }


def _score(soft_tp: float, query_total: int, truth_total: int) -> float:
    if truth_total <= 0 or query_total + truth_total <= 0:
        return 0.0
    if soft_tp < 0 or soft_tp > float(truth_total) + 1e-9:
        raise FormalMetricError(
            "soft true-positive credit exceeds the one-to-one truth universe"
        )
    return 2.0 * soft_tp / (query_total + truth_total) * 100.0


def evaluator_genotype_f1(
    ledgers: dict[str, dict[str, dict[str, Any]]],
    event_ids: set[str],
    panel_truth_positive: int,
) -> dict[str, dict[str, float | int]]:
    """Recompute a common GT-aware F1 from every evaluator ledger.

    A query is a TP only when the evaluator accepts the biological match and
    the diploid genotype is exact.  A genotype mismatch therefore contributes
    one FP and one FN.  Every evaluator uses the same submitted-event count and
    the same frozen panel-addressable truth denominator.
    """

    if panel_truth_positive <= 0:
        raise FormalMetricError(
            "panel-addressable truth contains no positive candidate"
        )
    metrics: dict[str, dict[str, float | int]] = {}
    for evaluator in EVALUATORS:
        rows = ledgers[evaluator]
        true_positive = sum(
            int(
                rows[result_id]["detection_correct"]
                and rows[result_id]["genotype_scorable"]
                and rows[result_id]["genotype_correct"]
            )
            for result_id in event_ids
        )
        if true_positive > panel_truth_positive:
            raise FormalMetricError(
                f"{evaluator} GT true positives exceed panel-addressable truth"
            )
        false_positive = len(event_ids) - true_positive
        false_negative = panel_truth_positive - true_positive
        precision = (
            true_positive / (true_positive + false_positive)
            if true_positive + false_positive
            else 0.0
        )
        recall = true_positive / panel_truth_positive
        f1 = (
            2.0 * precision * recall / (precision + recall)
            if precision + recall
            else 0.0
        )
        metrics[evaluator] = {
            "tp": true_positive,
            "fp": false_positive,
            "fn": false_negative,
            "precision": precision,
            "recall": recall,
            "f1": f1,
        }
    return metrics


def percentile(values: list[float], probability: float) -> float:
    if not values:
        raise FormalMetricError("cannot take percentile of empty values")
    ordered = sorted(values)
    position = (len(ordered) - 1) * probability
    lower = math.floor(position)
    upper = math.ceil(position)
    if lower == upper:
        return ordered[lower]
    fraction = position - lower
    return ordered[lower] * (1.0 - fraction) + ordered[upper] * fraction


def benchmark_block_keys(path: Path, block_size_bp: int) -> list[tuple[str, int]]:
    """Enumerate a tool-independent set of fixed genomic blocks from a BED."""

    if block_size_bp < 1:
        raise FormalMetricError("bootstrap block_size_bp must be positive")
    keys: set[tuple[str, int]] = set()
    with path.open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                raise FormalMetricError(
                    f"malformed BED record at {path}:{line_number}"
                )
            start, end = int(fields[1]), int(fields[2])
            if start < 0 or end <= start:
                raise FormalMetricError(
                    f"invalid BED interval at {path}:{line_number}"
                )
            for block in range(start // block_size_bp, (end - 1) // block_size_bp + 1):
                keys.add((fields[0], block))
    if not keys:
        raise FormalMetricError("benchmark BED yields no bootstrap blocks")
    return sorted(keys)


def _record_block(record: SvRecord, block_size_bp: int) -> tuple[str, int]:
    return record.chrom, (record.pos - 1) // block_size_bp


def bootstrap_score_interval(
    truth_records: list[SvRecord],
    query_records: dict[str, SvRecord],
    representative_rows: dict[str, dict[str, Any]],
    vote_counts: dict[str, int],
    *,
    benchmark_bed: Path,
    block_size_bp: int,
    replicates: int,
    seed_material: str,
) -> dict[str, Any]:
    """Jointly resample truth, matched queries, and unmatched FPs by block."""

    seed = int(hashlib.sha256(seed_material.encode("utf-8")).hexdigest()[:16], 16)
    generator = random.Random(seed)
    samples: list[float] = []
    block_keys = benchmark_block_keys(benchmark_bed, block_size_bp)
    allowed_blocks = set(block_keys)
    block_stats: dict[tuple[str, int], list[float]] = {
        key: [0.0, 0.0, 0.0] for key in block_keys
    }
    truth_by_id = {record.record_id: record for record in truth_records}
    for truth in truth_records:
        key = _record_block(truth, block_size_bp)
        if key not in allowed_blocks:
            raise FormalMetricError(
                f"truth event falls outside bootstrap BED blocks: {truth.record_id}"
            )
        block_stats[key][0] += 1.0
    for result_id, votes in vote_counts.items():
        query = query_records.get(result_id)
        if query is None:
            raise FormalMetricError(
                f"event query is absent from canonical VCF: {result_id}"
            )
        truth_id = representative_rows[result_id]["truth_event_id"]
        matched_truth = truth_by_id.get(truth_id) if truth_id else None
        key = _record_block(matched_truth or query, block_size_bp)
        if key not in allowed_blocks:
            raise FormalMetricError(
                f"query event falls outside bootstrap BED blocks: {result_id}"
            )
        block_stats[key][1] += 1.0
        block_stats[key][2] += votes / 3.0
    for _ in range(replicates):
        sampled_truth = 0
        sampled_query = 0
        sampled_credit = 0.0
        for _block in block_keys:
            truth_count, query_count, credit = block_stats[
                generator.choice(block_keys)
            ]
            sampled_truth += int(truth_count)
            sampled_query += int(query_count)
            sampled_credit += credit
        samples.append(_score(sampled_credit, sampled_query, sampled_truth))
    return {
        "lower": percentile(samples, 0.025),
        "upper": percentile(samples, 0.975),
        "level": 0.95,
        "method": "paired_genomic_block_bootstrap",
        "bootstrap_unit": "fixed_genomic_block",
        "block_size_bp": block_size_bp,
        "block_count": len(block_keys),
        "replicates": replicates,
        "seed_sha256": hashlib.sha256(seed_material.encode("utf-8")).hexdigest(),
    }


def bootstrap_me_f1_interval(
    *,
    truth_records: list[SvRecord],
    query_records: dict[str, SvRecord],
    hidden_truth: dict[str, dict[str, str]],
    ledgers: dict[str, dict[str, dict[str, Any]]],
    event_ids: set[str],
    benchmark_bed: Path,
    block_size_bp: int,
    replicates: int,
    seed_material: str,
) -> dict[str, Any]:
    """Bootstrap the actual ME-F1 with tool-shared genomic block draws.

    The seed depends only on frozen benchmark assets. Therefore every tool in
    the same comparison track receives identical draws and the stored
    replicate vector can be subtracted pairwise without pretending the tools
    were evaluated on independent samples.
    """

    block_keys = benchmark_block_keys(benchmark_bed, block_size_bp)
    allowed_blocks = set(block_keys)
    # truth, query, then one exact-GT TP count per evaluator
    block_stats: dict[tuple[str, int], list[int]] = {
        key: [0, 0, 0, 0, 0] for key in block_keys
    }
    truth_by_id = {record.record_id: record for record in truth_records}
    positive_candidates = {
        candidate_id: row["truth_event_id"]
        for candidate_id, row in hidden_truth.items()
        if row["truth_label"] == "positive"
    }
    for candidate_id, truth_id in positive_candidates.items():
        truth = truth_by_id.get(truth_id)
        if truth is None:
            raise FormalMetricError(
                f"panel-positive candidate {candidate_id} has no eligible truth event"
            )
        key = _record_block(truth, block_size_bp)
        if key not in allowed_blocks:
            raise FormalMetricError(
                f"panel-positive truth falls outside bootstrap BED: {truth_id}"
            )
        block_stats[key][0] += 1
    for candidate_id in event_ids:
        query = query_records.get(candidate_id)
        if query is None:
            raise FormalMetricError(
                f"event query is absent from canonical VCF: {candidate_id}"
            )
        truth_id = hidden_truth[candidate_id]["truth_event_id"]
        matched_truth = truth_by_id.get(truth_id) if truth_id else None
        key = _record_block(matched_truth or query, block_size_bp)
        if key not in allowed_blocks:
            raise FormalMetricError(
                f"query event falls outside bootstrap BED: {candidate_id}"
            )
        block_stats[key][1] += 1
        for index, evaluator in enumerate(EVALUATORS, start=2):
            row = ledgers[evaluator][candidate_id]
            block_stats[key][index] += int(
                row["detection_correct"]
                and row["genotype_scorable"]
                and row["genotype_correct"]
            )

    seed_hash = hashlib.sha256(seed_material.encode("utf-8")).hexdigest()
    generator = random.Random(int(seed_hash[:16], 16))
    evaluator_samples = {name: [] for name in EVALUATORS}
    me_samples: list[float] = []
    for _ in range(replicates):
        selected = [generator.choice(block_keys) for _key in block_keys]
        truth_count = sum(block_stats[key][0] for key in selected)
        query_count = sum(block_stats[key][1] for key in selected)
        values = []
        for index, evaluator in enumerate(EVALUATORS, start=2):
            tp = sum(block_stats[key][index] for key in selected)
            value = _score(float(tp), query_count, truth_count)
            evaluator_samples[evaluator].append(value)
            values.append(value)
        me_samples.append(statistics.fmean(values))
    return {
        "lower": percentile(me_samples, 0.025),
        "upper": percentile(me_samples, 0.975),
        "level": 0.95,
        "method": "paired_genomic_block_bootstrap",
        "bootstrap_unit": "fixed_genomic_block",
        "block_size_bp": block_size_bp,
        "block_count": len(block_keys),
        "replicates": replicates,
        "seed_sha256": seed_hash,
        "replicate_me_f1": me_samples,
        "replicate_evaluator_f1": evaluator_samples,
    }


def stratified_summary(
    truth_records: list[SvRecord],
    representative_rows: dict[str, dict[str, Any]],
    vote_counts: dict[str, int],
) -> dict[str, list[dict[str, Any]]]:
    result: dict[str, list[dict[str, Any]]] = {}
    truth_by_id = {record.record_id: record for record in truth_records}
    for dimension, truth_key, query_key in (
        ("svtype", lambda record: record.svtype, lambda row: row["query_svtype"]),
        (
            "length_bin",
            lambda record: _length_bin(record.svlen),
            lambda row: _length_bin(int(row["query_svlen"])),
        ),
    ):
        truth_counts: dict[str, int] = defaultdict(int)
        query_counts: dict[str, int] = defaultdict(int)
        soft_tp: dict[str, float] = defaultdict(float)
        agreement_counts: dict[str, list[int]] = defaultdict(lambda: [0, 0, 0, 0])
        for truth in truth_records:
            truth_counts[str(truth_key(truth))] += 1
        for result_id, row in representative_rows.items():
            if not row["event_eligible"]:
                continue
            truth_id = row["truth_event_id"]
            matched_truth = truth_by_id.get(truth_id) if truth_id is not None else None
            # A matched pair belongs to the truth stratum. This keeps the TP
            # numerator and truth denominator aligned when equivalent query
            # representations cross an SVTYPE or length-bin boundary. Only
            # unmatched query events (false positives) use the query stratum.
            stratum = str(
                truth_key(matched_truth)
                if matched_truth is not None
                else query_key(row)
            )
            query_counts[stratum] += 1
            votes = vote_counts[result_id]
            soft_tp[stratum] += votes / 3.0
            agreement_counts[stratum][votes] += 1
        entries = []
        for stratum in sorted(set(truth_counts) | set(query_counts)):
            entries.append(
                {
                    "stratum": stratum,
                    "truth_events": truth_counts[stratum],
                    "query_events": query_counts[stratum],
                    "soft_true_positive_count": soft_tp[stratum],
                    "whole_truth_recovery": _score(
                        soft_tp[stratum],
                        query_counts[stratum],
                        truth_counts[stratum],
                    ),
                    "all_three_correct": agreement_counts[stratum][3],
                    "exactly_two_correct": agreement_counts[stratum][2],
                    "exactly_one_correct": agreement_counts[stratum][1],
                    "none_correct": agreement_counts[stratum][0],
                }
            )
        result[dimension] = entries
    return result


def _genotype_class(genotype: str) -> str:
    """Return an auditable ploidy-aware class for a VCF GT string."""

    if gt_state(genotype) == "no_call":
        return "no_call"
    alleles = genotype.replace("|", "/").split("/")
    if len(alleles) == 1:
        return "haploid_ref" if alleles[0] == "0" else "haploid_alt"
    if all(allele == "0" for allele in alleles):
        return "hom_ref"
    if len(alleles) == 2 and "0" in alleles and any(
        allele != "0" for allele in alleles
    ):
        return "heterozygous"
    if len(set(alleles)) == 1 and alleles[0] != "0":
        return "hom_alt"
    return "other_nonref"


def candidate_genotype_summary(
    *,
    query_vcf: Path,
    hidden_truth_ledger: Path,
    tool_manifest: Path,
    require_phase: bool,
    no_call_is_incorrect: bool = True,
) -> dict[str, Any]:
    """Score the blinded candidate universe without affecting detection credit."""

    manifest = yaml.safe_load(tool_manifest.read_text(encoding="utf-8"))
    if not isinstance(manifest, dict) or not isinstance(
        manifest.get("outputs"), dict
    ):
        raise FormalMetricError("tool manifest has no outputs contract")
    output_contract = manifest["outputs"]
    candidate_output_contract = output_contract.get("candidate_output_contract")
    absence_semantics = output_contract.get("absence_semantics")
    if candidate_output_contract != "all_sites":
        raise FormalMetricError("invalid candidate_output_contract")
    if absence_semantics != "no_call":
        raise FormalMetricError("invalid absence_semantics")

    hidden: dict[str, dict[str, str]] = {}
    by_allele_id: dict[str, str] = {}
    with hidden_truth_ledger.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            "candidate_id",
            "pangenome_allele_id",
            "truth_gt",
            "truth_label",
            "truth_scorable",
            "truth_event_id",
        }
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise FormalMetricError("hidden candidate ledger has invalid header")
        for row in reader:
            candidate_id = row["candidate_id"]
            if not candidate_id or candidate_id in hidden:
                raise FormalMetricError(
                    f"empty or duplicate hidden candidate ID: {candidate_id!r}"
                )
            allele_id = row["pangenome_allele_id"]
            if not allele_id or allele_id in by_allele_id:
                raise FormalMetricError(
                    "empty or duplicate hidden pangenome allele ID: "
                    f"{allele_id!r}"
                )
            if row["truth_scorable"] not in {"0", "1"}:
                raise FormalMetricError(
                    f"invalid truth_scorable for candidate {candidate_id}"
                )
            if row["truth_scorable"] == "1":
                if row["truth_label"] not in {"positive", "negative"}:
                    raise FormalMetricError(
                        f"scorable candidate {candidate_id} has invalid truth label"
                    )
            elif row["truth_label"] != "unscorable":
                raise FormalMetricError(
                    f"unscorable candidate {candidate_id} has invalid truth label"
                )
            hidden[candidate_id] = row
            by_allele_id[allele_id] = candidate_id
    if not hidden:
        raise FormalMetricError("hidden candidate ledger is empty")

    observed: dict[str, str] = {}
    query_record_count = 0
    unmapped_query_records = 0
    direct_candidate_records = 0
    linked_candidate_records = 0
    candidate_link_conflicts = 0
    rejected_link_records = 0
    duplicate_candidate_records = 0
    conflicting_duplicate_candidates: set[str] = set()
    observed_priority: dict[str, int] = {}
    accepted_link_statuses = {"in_panel_exact", "in_panel_equivalent"}
    with open_text(query_vcf) as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                raise FormalMetricError("linked query VCF has a malformed record")
            query_record_count += 1
            info = parse_info(fields[7])
            original_id = info.get("ORIG_ID")
            direct_candidate_id = (
                str(original_id)
                if isinstance(original_id, str) and original_id in hidden
                else (fields[2] if fields[2] in hidden else None)
            )
            linked_allele_id = info.get("PANGENOME_LINKED_ID")
            link_status = info.get("PANGENOME_LINK_STATUS")
            linked_candidate_id = None
            if isinstance(linked_allele_id, str):
                if (
                    isinstance(link_status, str)
                    and link_status in accepted_link_statuses
                ):
                    linked_candidate_id = by_allele_id.get(linked_allele_id)
                else:
                    rejected_link_records += 1
            if (
                direct_candidate_id is not None
                and linked_candidate_id is not None
                and direct_candidate_id != linked_candidate_id
            ):
                candidate_link_conflicts += 1
            candidate_id = direct_candidate_id or linked_candidate_id
            if candidate_id is None:
                unmapped_query_records += 1
                continue
            if direct_candidate_id is not None:
                direct_candidate_records += 1
            else:
                linked_candidate_records += 1
            if len(fields) < 10:
                query_gt = "./."
            else:
                format_keys = fields[8].split(":")
                sample_values = fields[9].split(":")
                query_gt = (
                    sample_values[format_keys.index("GT")]
                    if "GT" in format_keys
                    and format_keys.index("GT") < len(sample_values)
                    else "./."
                )
            # Discovery callers can emit several equivalent records for one
            # panel allele (for example, two nearby graph traversals for the
            # same repeat insertion).  A direct candidate ID is stronger than
            # a positional/sequence link.  Equal-priority duplicate calls are
            # accepted when their genotypes agree; conflicting calls become a
            # no-call so the benchmark cannot cherry-pick the favourable one.
            priority = 1 if direct_candidate_id is not None else 0
            if candidate_id in observed:
                duplicate_candidate_records += 1
                previous_priority = observed_priority[candidate_id]
                if priority > previous_priority:
                    observed[candidate_id] = query_gt
                    observed_priority[candidate_id] = priority
                    conflicting_duplicate_candidates.discard(candidate_id)
                elif priority == previous_priority and not genotypes_equal(
                    observed[candidate_id],
                    query_gt,
                    require_phase=require_phase,
                ):
                    observed[candidate_id] = "./."
                    conflicting_duplicate_candidates.add(candidate_id)
                continue
            observed[candidate_id] = query_gt
            observed_priority[candidate_id] = priority

    missing = set(hidden) - set(observed)
    if candidate_output_contract == "all_sites" and missing:
        raise FormalMetricError(
            "all-sites output is missing hidden candidate IDs after normalization"
        )
    inferred_gt = "0/0" if absence_semantics == "hom_ref" else "./."
    counts: dict[str, int] = {
        "candidate_count": len(hidden),
        "truth_scorable": 0,
        "truth_positive": 0,
        "truth_negative": 0,
        "truth_unscorable": 0,
        "observed_candidates": len(observed),
        "observed_truth_scorable_candidates": 0,
        "inferred_absent_candidates": len(missing),
        "inferred_absent_truth_scorable_candidates": 0,
        "genotype_scorable": 0,
        "called_candidates": 0,
        "genotype_correct": 0,
        "no_call": 0,
        "hom_ref": 0,
        "variant": 0,
    }
    confusion: dict[str, dict[str, int]] = defaultdict(
        lambda: defaultdict(int)
    )
    no_call_by_truth_class: dict[str, int] = defaultdict(int)
    binary_confusion = {
        "true_positive": 0,
        "false_positive": 0,
        "true_negative": 0,
        "false_negative": 0,
        "no_call_truth_positive": 0,
        "no_call_truth_negative": 0,
    }
    for candidate_id, truth_row in hidden.items():
        truth_is_scorable = truth_row["truth_scorable"] == "1"
        if not truth_is_scorable:
            counts["truth_unscorable"] += 1
            continue
        counts["truth_scorable"] += 1
        counts["observed_truth_scorable_candidates"] += int(
            candidate_id in observed
        )
        counts["inferred_absent_truth_scorable_candidates"] += int(
            candidate_id in missing
        )
        truth_gt = truth_row["truth_gt"]
        truth_state = gt_state(truth_gt)
        truth_nonref = truth_row["truth_label"] == "positive"
        if truth_nonref:
            counts["truth_positive"] += 1
        else:
            counts["truth_negative"] += 1
        query_gt = observed.get(candidate_id, inferred_gt)
        query_state = gt_state(query_gt)
        counts[query_state] += 1
        query_nonref = query_state == "variant"
        if truth_nonref and query_nonref:
            binary_confusion["true_positive"] += 1
        elif truth_nonref:
            binary_confusion["false_negative"] += 1
            binary_confusion["no_call_truth_positive"] += int(
                query_state == "no_call"
            )
        elif query_nonref:
            binary_confusion["false_positive"] += 1
        elif query_state == "hom_ref":
            binary_confusion["true_negative"] += 1
        else:
            binary_confusion["no_call_truth_negative"] += 1

        # Detection truth and genotype truth are separate contracts.  A
        # partially missing truth GT (for example .|1 on chrX) still provides
        # a valid positive/negative candidate label, but it cannot support an
        # exact genotype comparison or enter genotype denominators.
        if truth_state == "no_call":
            continue
        counts["genotype_scorable"] += 1
        truth_class = _genotype_class(truth_gt)
        query_class = _genotype_class(query_gt)
        confusion[truth_class][query_class] += 1
        no_call_by_truth_class[truth_class] += int(query_class == "no_call")
        called = query_state != "no_call"
        counts["called_candidates"] += int(called)
        counts["genotype_correct"] += int(
            called
            and genotypes_equal(
                query_gt,
                truth_gt,
                require_phase=require_phase,
            )
        )
    denominator = (
        counts["genotype_scorable"]
        if no_call_is_incorrect
        else counts["called_candidates"]
    )
    called_denominator = counts["called_candidates"]
    genotype_classes = sorted(confusion)
    all_prediction_classes = sorted(
        {
            predicted
            for predictions in confusion.values()
            for predicted in predictions
        }
        | set(genotype_classes)
        | {"no_call"}
    )
    confusion_matrix = {
        truth_class: {
            predicted_class: confusion[truth_class].get(predicted_class, 0)
            for predicted_class in all_prediction_classes
        }
        for truth_class in genotype_classes
    }
    class_metrics: dict[str, dict[str, float | int | None]] = {}
    for genotype_class in genotype_classes:
        true_positive = confusion[genotype_class].get(genotype_class, 0)
        truth_count = sum(confusion[genotype_class].values())
        predicted_count = sum(
            predictions.get(genotype_class, 0)
            for predictions in confusion.values()
        )
        precision = (
            true_positive / predicted_count if predicted_count > 0 else None
        )
        recall = true_positive / truth_count if truth_count > 0 else None
        f1 = (
            2.0 * precision * recall / (precision + recall)
            if precision is not None
            and recall is not None
            and precision + recall > 0
            else 0.0
        )
        class_metrics[genotype_class] = {
            "truth_count": truth_count,
            "predicted_count": predicted_count,
            "precision": precision,
            "recall": recall,
            "f1": f1,
        }
    macro_f1 = (
        sum(float(metrics["f1"]) for metrics in class_metrics.values())
        / len(class_metrics)
        if class_metrics
        else None
    )
    nonref_precision_denominator = (
        binary_confusion["true_positive"] + binary_confusion["false_positive"]
    )
    nonref_precision = (
        binary_confusion["true_positive"] / nonref_precision_denominator
        if nonref_precision_denominator > 0
        else None
    )
    nonref_sensitivity = (
        binary_confusion["true_positive"] / counts["truth_positive"]
        if counts["truth_positive"] > 0
        else None
    )
    specificity = (
        binary_confusion["true_negative"] / counts["truth_negative"]
        if counts["truth_negative"] > 0
        else None
    )
    balanced_accuracy = (
        (nonref_sensitivity + specificity) / 2.0
        if nonref_sensitivity is not None and specificity is not None
        else None
    )
    nonref_f1 = (
        2.0
        * nonref_precision
        * nonref_sensitivity
        / (nonref_precision + nonref_sensitivity)
        if nonref_precision is not None
        and nonref_sensitivity is not None
        and nonref_precision + nonref_sensitivity > 0
        else 0.0
    )
    return {
        "contract": "hidden_candidate_genotype_v3",
        "candidate_output_contract": candidate_output_contract,
        "absence_semantics": absence_semantics,
        "no_call_is_incorrect": no_call_is_incorrect,
        "query_record_count": query_record_count,
        "unmapped_query_records": unmapped_query_records,
        "direct_candidate_records": direct_candidate_records,
        "linked_candidate_records": linked_candidate_records,
        "candidate_link_conflicts": candidate_link_conflicts,
        "duplicate_candidate_records": duplicate_candidate_records,
        "conflicting_duplicate_candidates": len(conflicting_duplicate_candidates),
        "rejected_link_records": rejected_link_records,
        **counts,
        "genotype_accuracy": (
            counts["genotype_correct"] / denominator
            if denominator > 0
            else None
        ),
        "called_only_genotype_accuracy": (
            counts["genotype_correct"] / called_denominator
            if called_denominator > 0
            else None
        ),
        "genotype_confusion_matrix": confusion_matrix,
        "genotype_class_metrics": class_metrics,
        "genotype_macro_f1": macro_f1,
        "binary_nonref_confusion": binary_confusion,
        "nonref_precision": nonref_precision,
        "nonref_sensitivity": nonref_sensitivity,
        "specificity": specificity,
        "balanced_accuracy": balanced_accuracy,
        "nonref_f1": nonref_f1,
        "no_call_by_truth_class": {
            genotype_class: no_call_by_truth_class.get(genotype_class, 0)
            for genotype_class in genotype_classes
        },
    }


def resource_summary(path: Path | None, cache_policy: str) -> dict[str, Any]:
    summary: dict[str, Any] = {
        "cache_policy": cache_policy,
        "measurement_file": str(path) if path is not None else None,
        "repeat_count": 0,
        "median": {},
        "iqr": {},
    }
    if path is None or not path.is_file() or path.stat().st_size == 0:
        return summary
    rows: list[dict[str, Any]] = []
    text = path.read_text(encoding="utf-8").strip()
    if not text:
        return summary
    if text.lstrip().startswith("{"):
        for line in text.splitlines():
            loaded = json.loads(line)
            if isinstance(loaded, dict):
                rows.append(loaded)
    else:
        dialect = csv.Sniffer().sniff(text[:4096], delimiters="\t,")
        rows.extend(csv.DictReader(text.splitlines(), dialect=dialect))
    aliases = {
        "wall_seconds": ("s", "seconds", "walltime", "wall_time"),
        "cpu_seconds": ("cpu_time", "cpu_seconds"),
        "max_rss_mb": ("max_rss", "max_rss_mb"),
        "max_vms_mb": ("max_vms", "max_vms_mb"),
    }
    values: dict[str, list[float]] = defaultdict(list)
    for row in rows:
        for output_name, candidates in aliases.items():
            for candidate in candidates:
                value = row.get(candidate)
                if value not in {None, ""}:
                    values[output_name].append(float(value))
                    break
    summary["repeat_count"] = len(rows)
    for name, samples in values.items():
        summary["median"][name] = statistics.median(samples)
        if len(samples) >= 2:
            quartiles = statistics.quantiles(samples, n=4, method="inclusive")
            summary["iqr"][name] = quartiles[2] - quartiles[0]
        else:
            summary["iqr"][name] = 0.0
    return summary


def evidence_profile(
    path: Path,
    *,
    technology: str,
    library_id: str,
    source_evidence_id: str,
    coverage_x: float | None,
    read_count: int | None,
    read_bases: int | None,
    downsampling_seed: int | None,
) -> dict[str, Any]:
    """Freeze the actual read/alignment evidence used by one tool run."""

    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise FormalMetricError(f"cannot read resolved inputs: {path}") from exc
    if not isinstance(payload, dict) or payload.get("schema_version") != 1:
        raise FormalMetricError("resolved-input contract is not schema v1")
    inputs = payload.get("inputs")
    if not isinstance(inputs, list):
        raise FormalMetricError("resolved-input contract has no inputs list")
    evidence_names = {"short_fastq_r1", "short_fastq_r2"}
    selected: dict[str, dict[str, Any]] = {}
    for item in inputs:
        if not isinstance(item, dict):
            raise FormalMetricError("resolved-input entry is not an object")
        name = item.get("name")
        if name not in evidence_names:
            continue
        if name in selected:
            raise FormalMetricError(f"duplicate evidence input: {name}")
        digest = item.get("sha256")
        size_bytes = item.get("size_bytes")
        if not isinstance(digest, str) or len(digest) != 64:
            raise FormalMetricError(f"evidence input {name} has invalid SHA-256")
        if not isinstance(size_bytes, int) or size_bytes < 1:
            raise FormalMetricError(f"evidence input {name} has invalid size")
        selected[str(name)] = {
            "sha256": digest,
            "size_bytes": size_bytes,
            "path_type": item.get("path_type"),
        }
    primary = set(selected)
    paired = {"short_fastq_r1", "short_fastq_r2"}
    if primary != paired:
        raise FormalMetricError("resolved short-read evidence is not a complete pair")
    evidence_kind = "paired_fastq"
    read_validation = payload.get("read_validation", {})
    measured_read_count: int | None = None
    measured_read_bases: int | None = None
    if isinstance(read_validation, dict):
        validation = read_validation.get("paired_fastq")
        if isinstance(validation, dict):
            candidate_count = validation.get("read_count")
            candidate_bases = validation.get("read_bases")
            if isinstance(candidate_count, int) and candidate_count > 0:
                measured_read_count = candidate_count
            if isinstance(candidate_bases, int) and candidate_bases > 0:
                measured_read_bases = candidate_bases
    if measured_read_count is not None:
        if read_count is not None and read_count != measured_read_count:
            raise FormalMetricError(
                "configured read_count disagrees with complete FASTQ validation"
            )
        read_count = measured_read_count
    if measured_read_bases is not None:
        if read_bases is not None and read_bases != measured_read_bases:
            raise FormalMetricError(
                "configured read_bases disagrees with complete FASTQ validation"
            )
        read_bases = measured_read_bases
    return {
        "contract": "pgbench_evidence_v1",
        "actual_technology": technology,
        "library_id": library_id,
        "source_evidence_id": source_evidence_id,
        "evidence_kind": evidence_kind,
        "official_score_mode": payload.get("mode"),
        "alignment_kind": payload.get("alignment_kind"),
        "input_assets": dict(sorted(selected.items())),
        "resolved_inputs_sha256": sha256_file(path),
        "coverage_x": coverage_x,
        "read_count": read_count,
        "read_bases": read_bases,
        "downsampling_seed": downsampling_seed,
        "quantitative_metadata_complete": all(
            value is not None for value in (coverage_x, read_count, read_bases)
        ),
    }


def _metric_record(
    metric_id: str,
    value: int,
    *,
    denominator: int,
    truth_metric: bool,
    provenance_id: str,
) -> dict[str, Any]:
    return {
        "metric_id": metric_id,
        "value_type": "count",
        "value": value,
        "status": "defined",
        "evaluator": "fusion",
        "numerator": value,
        "denominator": denominator,
        "eligible_count": denominator,
        "universe_id": (
            "primary_truth_benchmark_universe_v2"
            if truth_metric
            else "formal_detection_event_universe_v2"
        ),
        "undefined_reason": None,
        "parser_id": "three_evaluator_detection_vote_fusion_v2",
        "parser_source_field": (
            "truth_event_id" if truth_metric else "detection_correct"
        ),
        "provenance_manifest_id": provenance_id,
        "strata": {},
    }


def context_bed_arguments(values: list[str]) -> dict[str, Path]:
    parsed: dict[str, Path] = {}
    for value in values:
        name, separator, path = value.partition("=")
        if not separator or not name or not path or name in parsed:
            raise FormalMetricError(
                "--context-bed must be supplied once per NAME=PATH"
            )
        parsed[name] = Path(path)
    if set(parsed) != set(REQUIRED_CONTEXT_STRATA):
        missing = sorted(REQUIRED_CONTEXT_STRATA - set(parsed))
        extra = sorted(set(parsed) - REQUIRED_CONTEXT_STRATA)
        details = []
        if missing:
            details.append("missing " + ", ".join(missing))
        if extra:
            details.append("unexpected " + ", ".join(extra))
        raise FormalMetricError("context BED contract differs: " + "; ".join(details))
    for name, path in parsed.items():
        if not path.is_file():
            raise FormalMetricError(f"missing context BED {name}: {path}")
    return parsed


def materialize(args: argparse.Namespace) -> dict[str, Any]:
    loaded = {
        name: load_ledger(getattr(args, f"{name}_ledger"), name)
        for name in EVALUATORS
    }
    ledgers = {name: value[0] for name, value in loaded.items()}
    extended_flags = {value[1] for value in loaded.values()}
    if len(extended_flags) != 1:
        raise FormalMetricError("evaluator ledgers do not share one semantic format")
    extended = extended_flags.pop()
    first = set(ledgers["truvari"])
    for name, rows in ledgers.items():
        if set(rows) != first:
            raise FormalMetricError(
                f"{name} result universe differs from Truvari"
            )
    for result_id in first:
        contracts = {
            (
                rows[result_id]["scope_eligible"],
                rows[result_id]["event_eligible"],
                rows[result_id]["truth_event_id"],
                rows[result_id]["gt_state"],
                rows[result_id]["query_svtype"],
                rows[result_id]["query_svlen"],
            )
            for rows in ledgers.values()
        }
        if len(contracts) != 1:
            raise FormalMetricError(
                f"semantic/truth assignment differs for {result_id}"
            )

    evaluator_profile_path = getattr(
        args,
        "evaluator_profile",
        PROJECT_ROOT / "config" / "evaluator_profile.yaml",
    )
    evaluator_profile = load_evaluator_profile(evaluator_profile_path)
    truth_records = eligible_truth_records(
        args.truth_vcf,
        args.benchmark_bed,
        evaluator_profile=evaluator_profile if extended else None,
    )
    truth_total = len(truth_records)
    truth_ids = [record.record_id for record in truth_records]
    truth_id_set = set(truth_ids)
    representative = ledgers["truvari"]
    event_ids = sorted(
        result_id
        for result_id, row in representative.items()
        if row["event_eligible"]
    )
    counts = [0, 0, 0, 0]
    vote_counts: dict[str, int] = {}
    credits: dict[str, float] = {}
    for result_id in event_ids:
        vote_count = sum(
            int(ledgers[name][result_id]["detection_correct"])
            for name in EVALUATORS
        )
        vote_counts[result_id] = vote_count
        counts[vote_count] += 1
        truth_id = representative[result_id]["truth_event_id"]
        if extended and truth_id is not None:
            if truth_id not in truth_id_set:
                raise FormalMetricError(
                    f"ledger references an ineligible truth event: {truth_id}"
                )
            if truth_id in credits:
                raise FormalMetricError(
                    f"truth event credited more than once: {truth_id}"
                )
            credits[truth_id] = vote_count / 3.0

    profile = yaml.safe_load(args.score_profile.read_text(encoding="utf-8"))
    profile_id = profile["profile"]["id"]
    completions: dict[str, dict[str, Any]] = {}
    evaluator_completions = list(getattr(args, "evaluator_completions", []))
    if evaluator_completions:
        completions = load_completions(
            evaluator_completions,
            evaluator_profile["_sha256"],
        )
    provenance_paths = [
        *args.evaluator_manifests,
        *evaluator_completions,
        evaluator_profile_path,
    ]
    provenance_id = provenance_bundle(provenance_paths)
    categories = (
        ("benchmark.truth.eligible.count", truth_total, True),
        ("diagnostic.evaluator_agreement.all_three_correct.count", counts[3], False),
        ("diagnostic.evaluator_agreement.exactly_two_correct.count", counts[2], False),
        ("diagnostic.evaluator_agreement.exactly_one_correct.count", counts[1], False),
        ("diagnostic.evaluator_agreement.none_correct.count", counts[0], False),
    )
    records = [
        _metric_record(
            metric_id,
            value,
            denominator=truth_total if truth_metric else len(event_ids),
            truth_metric=truth_metric,
            provenance_id=provenance_id,
        )
        for metric_id, value, truth_metric in categories
    ]

    soft_tp = sum(vote_counts.values()) / 3.0
    mapping_summary = {
        name: {
            "event_count": len(event_ids),
            "resolved_events": sum(
                int(rows[result_id]["evaluator_resolved"])
                for result_id in event_ids
            ),
            "unresolved_events": sum(
                int(not rows[result_id]["evaluator_resolved"])
                for result_id in event_ids
            ),
        }
        for name, rows in ledgers.items()
    }
    for values in mapping_summary.values():
        values["mapping_coverage"] = (
            values["resolved_events"] / values["event_count"]
            if values["event_count"]
            else 1.0
        )
    max_unresolved_fraction = float(
        evaluator_profile["semantics"].get(
            "maximum_unresolved_mapping_fraction", 0.01
        )
    )
    global_recovery_valid = all(
        1.0 - float(values["mapping_coverage"]) <= max_unresolved_fraction
        for values in mapping_summary.values()
    )
    unresolved_failures = [
        (
            f"{name}={int(values['unresolved_events'])}/"
            f"{int(values['event_count'])} "
            f"({1.0 - float(values['mapping_coverage']):.2%})"
        )
        for name, values in sorted(mapping_summary.items())
        if 1.0 - float(values["mapping_coverage"]) > max_unresolved_fraction
    ]
    if unresolved_failures:
        formal_score_status = "invalid_evaluator_mapping"
        formal_score_reason = (
            "evaluator unresolved mapping exceeds "
            f"{max_unresolved_fraction:.2%}: " + ", ".join(unresolved_failures)
        )
    else:
        formal_score_status = "valid"
        formal_score_reason = None
    query_vcf = getattr(args, "query_vcf", None)
    if extended and query_vcf is None:
        raise FormalMetricError(
            "semantic formal scoring requires the canonical query VCF"
        )
    loaded_query_records = (
        load_vcf(query_vcf, prefix="query", allow_empty=True)
        if query_vcf is not None
        else []
    )
    query_records = {
        record.record_id: record for record in loaded_query_records
    }
    if len(query_records) != len(loaded_query_records):
        raise FormalMetricError("canonical query VCF has duplicate result IDs")
    context_beds = context_bed_arguments(args.context_bed)
    hidden_truth = load_hidden_truth_labels(args.hidden_truth_ledger)
    addressability_rows = load_addressability_rows(args.addressability_tsv)
    candidate_af = load_candidate_af(args.query_vcf)
    genotype_strata = genotype_stratified_summary(
        query_records=query_records,
        hidden_truth=hidden_truth,
        ledgers=ledgers,
        event_ids=set(event_ids),
        candidate_af=candidate_af,
        context_beds=context_beds,
        addressability=addressability_rows,
    )
    ci = (
        bootstrap_score_interval(
            truth_records,
            query_records,
            representative,
            vote_counts,
            benchmark_bed=args.benchmark_bed,
            block_size_bp=int(
                evaluator_profile["bootstrap"]["block_size_bp"]
            ),
            replicates=int(evaluator_profile["bootstrap"]["replicates"]),
            seed_material=content_seed_bundle(
                [
                    ("truth_vcf", args.truth_vcf),
                    ("benchmark_bed", args.benchmark_bed),
                    ("score_profile", args.score_profile),
                    ("evaluator_profile", evaluator_profile_path),
                    *(
                        [("reference", args.reference)]
                        if getattr(args, "reference", None) is not None
                        else []
                    ),
                ]
            ),
        )
        if extended
        else None
    )
    me_f1_ci = bootstrap_me_f1_interval(
        truth_records=truth_records,
        query_records=query_records,
        hidden_truth=hidden_truth,
        ledgers=ledgers,
        event_ids=set(event_ids),
        benchmark_bed=args.benchmark_bed,
        block_size_bp=int(evaluator_profile["bootstrap"]["block_size_bp"]),
        replicates=int(evaluator_profile["bootstrap"]["replicates"]),
        seed_material=content_seed_bundle(
            [
                ("truth_vcf", args.truth_vcf),
                ("benchmark_bed", args.benchmark_bed),
                ("score_profile", args.score_profile),
                ("evaluator_profile", evaluator_profile_path),
                ("reference", args.reference),
                ("hidden_truth_ledger", args.hidden_truth_ledger),
            ]
        ),
    )
    asset_hashes = {
        "truth_vcf": sha256_file(args.truth_vcf),
        "benchmark_bed": sha256_file(args.benchmark_bed),
        "reference": sha256_file(args.reference),
        "pangenome_manifest": semantic_asset_hash(
            args.pangenome_manifest,
            excluded_root_keys=frozenset(
                {"graph_assets", "graph_build_recipe_sha256"}
            ),
        ),
        "challenge_hidden_ledger": sha256_file(args.hidden_truth_ledger),
        "tool_manifest": (
            sha256_file(args.tool_manifest)
            if getattr(args, "tool_manifest", None) is not None
            else None
        ),
        "graph_asset_lock": (
            semantic_asset_hash(args.graph_asset_lock)
            if getattr(args, "graph_asset_lock", None) is not None
            else None
        ),
        "stratification_catalogue": sha256_file(args.stratification_catalogue),
        "context_beds": {
            name: sha256_file(path) for name, path in sorted(context_beds.items())
        },
    }
    candidate_summary = (
        candidate_genotype_summary(
            query_vcf=args.query_vcf,
            hidden_truth_ledger=args.hidden_truth_ledger,
            tool_manifest=args.tool_manifest,
            require_phase=bool(
                evaluator_profile["semantics"]["genotype"].get(
                    "require_phase", False
                )
            ),
            no_call_is_incorrect=bool(
                evaluator_profile["semantics"].get(
                    "no_call_is_incorrect", True
                )
            ),
        )
        if all(
            getattr(args, name, None) is not None
            for name in ("query_vcf", "hidden_truth_ledger", "tool_manifest")
        )
        else None
    )
    genotype_macro_f1 = (
        candidate_summary.get("genotype_macro_f1")
        if isinstance(candidate_summary, dict)
        else None
    )
    nonref_f1 = (
        candidate_summary.get("nonref_f1")
        if isinstance(candidate_summary, dict)
        else None
    )
    panel_truth_positive = (
        int(candidate_summary.get("truth_positive", 0))
        if isinstance(candidate_summary, dict)
        else 0
    )
    panel_coverage = (
        panel_truth_positive / truth_total if truth_total > 0 else None
    )
    all_sites_contract = (
        isinstance(candidate_summary, dict)
        and candidate_summary.get("candidate_output_contract") == "all_sites"
    )
    evaluator_gt_metrics = evaluator_genotype_f1(
        ledgers,
        event_ids,
        panel_truth_positive,
    )
    tool_contract = yaml.safe_load(args.tool_manifest.read_text(encoding="utf-8"))
    information_contract = (
        tool_contract.get("information_contract", {})
        if isinstance(tool_contract, dict)
        else {}
    )
    semantic_validation = load_json_object(
        args.semantic_validation, "evaluator semantic validation"
    )
    addressability_validation = load_json_object(
        args.addressability_audit, "addressability audit"
    )
    allowed_information = load_json_object(
        args.allowed_information, "allowed-information manifest"
    )
    parameter_manifest = load_json_object(
        args.parameter_manifest, "parameter manifest"
    )
    external_resources = load_json_object(
        args.external_resource_hashes, "external-resource manifest"
    )
    tuning_status = load_json_object(args.tuning_status, "tuning-status manifest")
    quality_gates = {
        "evaluator_semantic_contract_valid": (
            isinstance(semantic_validation, dict)
            and semantic_validation.get("contract")
            == "pgbench_genotype_evaluator_contract_v1"
            and semantic_validation.get("status") == "valid"
        ),
        "addressability_complete": bool(
            all_sites_contract
            and isinstance(candidate_summary, dict)
            and candidate_summary.get("observed_candidates")
            == candidate_summary.get("candidate_count")
            and addressability_validation.get("contract")
            == "pgbench_addressability_v1"
            and addressability_validation.get("canonical_candidate_count")
            == candidate_summary.get("candidate_count")
            and addressability_validation.get("silent_candidate_deletion_count") == 0
        ),
        "allowed_information_complete": bool(
            isinstance(information_contract, dict)
            and information_contract.get("allowed_inputs_manifested") is True
            and information_contract.get("target_truth_used_for_calling") is False
            and information_contract.get("target_assembly_used") is False
            and information_contract.get("target_family_genotypes_used") is False
            and information_contract.get("target_specific_external_callset") is False
            and allowed_information.get("status") == "valid"
            and parameter_manifest.get("status") == "frozen"
            and external_resources.get("status") == "frozen"
        ),
        "tuning_frozen": bool(
            isinstance(information_contract, dict)
            and information_contract.get("tuning_status")
            in {"official_configuration", "independent_validation_frozen"}
            and tuning_status.get("status") == "frozen_before_test"
            and tuning_status.get("target_truth_inspected") is False
        ),
        "context_stratification_complete": bool(
            set(context_beds) == set(REQUIRED_CONTEXT_STRATA)
            and set(genotype_strata["genome_context"])
            == set(REQUIRED_CONTEXT_STRATA)
        ),
        "population_af_stratification_complete": bool(
            genotype_strata["population_af"].get("unknown", {}).get(
                "canonical_candidate_count", 0
            )
            == 0
        ),
    }
    semantic_summary = {
        name: {
            "query_rows": len(rows),
            "scope_eligible": sum(
                int(row["scope_eligible"]) for row in rows.values()
            ),
            "scope_excluded": sum(
                int(not row["scope_eligible"]) for row in rows.values()
            ),
            "detection_events": len(event_ids),
            "detection_correct": sum(
                int(row["detection_correct"]) for row in rows.values()
            ),
            "genotype_scorable": sum(
                int(row["genotype_scorable"]) for row in rows.values()
            ),
            "genotype_correct": sum(
                int(row["genotype_correct"]) for row in rows.values()
            ),
            "no_call": sum(int(row["no_call"]) for row in rows.values()),
        }
        for name, rows in ledgers.items()
    }
    version_summary = {
        name: {
            "output": completion["version"]["output"],
            "sha256": completion["version"]["sha256"],
        }
        for name, completion in completions.items()
    }
    evaluator_native_metrics = {
        name: completion.get("native_metrics")
        for name, completion in completions.items()
        if completion.get("native_metrics") is not None
    }
    frozen_evidence = (
        evidence_profile(
            args.resolved_inputs,
            technology=args.sample_technology,
            library_id=args.library_id,
            source_evidence_id=args.source_evidence_id,
            coverage_x=args.coverage_x,
            read_count=args.read_count,
            read_bases=args.read_bases,
            downsampling_seed=args.downsampling_seed,
        )
        if getattr(args, "resolved_inputs", None) is not None
        else None
    )
    if extended:
        comparison_track, comparison_track_sha256 = comparison_track_contract(
            tool_manifest=args.tool_manifest,
            official_score_mode=args.official_score_mode,
            primary_truth_profile=args.primary_truth_profile,
            score_profile=args.score_profile,
            evaluator_profile_sha256=evaluator_profile["_sha256"],
            asset_hashes=asset_hashes,
            evidence=frozen_evidence,
        )
    else:
        comparison_track = None
        comparison_track_sha256 = None
    return {
        "schema_version": 1,
        "dictionary_id": "pgbench_metrics_v1",
        "tuple": {
            "run_id": args.run_id,
            "sample_id": args.sample_id,
            "tool_id": args.tool_id,
            "official_score_mode": args.official_score_mode,
            "primary_truth_profile": args.primary_truth_profile,
            "score_profile": profile_id,
        },
        "records": records,
        "analysis": {
            "contract_version": "formal_genotype_me_f1_v1",
            "evaluator_profile_id": evaluator_profile["profile"]["id"],
            "evaluator_profile_sha256": evaluator_profile["_sha256"],
            "evaluator_bundle_sha256": provenance_id,
            "evaluator_versions": version_summary,
            "evaluator_native_metrics": evaluator_native_metrics,
            "evaluator_gt_metrics": evaluator_gt_metrics,
            "quality_gates": quality_gates,
            "asset_hashes": asset_hashes,
            "evidence_profile": frozen_evidence,
            "comparison_track": comparison_track,
            "comparison_track_sha256": comparison_track_sha256,
            "universe": evaluator_profile["universe"],
            "matching": {
                "algorithm": (
                    "deterministic_max_cardinality_max_weight_one_to_one_v1"
                ),
                "unique_truth_accounting": extended,
                "credited_truth_events": len(credits),
                # These are diagnostic denominators for whole-truth recovery. They
                # remain available even when the formal score is fail-closed
                # (and its top-level counts are deliberately unset/zeroed).
                "eligible_truth_event_count": truth_total,
                "query_event_count": len(event_ids),
                "soft_true_positive_count": soft_tp,
            },
            "evaluator_event_mapping": mapping_summary,
            "maximum_unresolved_mapping_fraction": max_unresolved_fraction,
            "whole_truth_recovery_status": (
                "valid" if global_recovery_valid else "invalid_unresolved_mapping"
            ),
            "formal_score_status": formal_score_status,
            "formal_score_reason": formal_score_reason,
            "whole_truth_recovery_score": _score(soft_tp, len(event_ids), truth_total),
            # Genotype and non-reference F1 are meaningful only for an
            # all-sites genotyping contract and is diagnostic only.
            "pangenome_genotyping_score": (
                float(genotype_macro_f1) * 100.0
                if all_sites_contract and genotype_macro_f1 is not None
                else None
            ),
            "non_reference_f1_score": (
                float(nonref_f1) * 100.0
                if all_sites_contract and nonref_f1 is not None
                else None
            ),
            "panel_coverage": panel_coverage,
            "panel_truth_positive_events": panel_truth_positive,
            "panel_addressable_truth_count": panel_truth_positive,
            "me_f1": (
                sum(float(values["f1"]) for values in evaluator_gt_metrics.values())
                / len(evaluator_gt_metrics)
                * 100.0
            ),
            "evaluator_range": (
                max(float(values["f1"]) for values in evaluator_gt_metrics.values())
                - min(float(values["f1"]) for values in evaluator_gt_metrics.values())
            )
            * 100.0,
            "evaluator_sd": statistics.pstdev(
                float(values["f1"]) for values in evaluator_gt_metrics.values()
            ) * 100.0,
            "leave_one_evaluator_out": {
                excluded: sum(
                    float(values["f1"])
                    for name, values in evaluator_gt_metrics.items()
                    if name != excluded
                ) / 2.0 * 100.0
                for excluded in EVALUATORS
            },
            "evaluator_agreement_counts": {
                "all_three_correct": counts[3],
                "exactly_two_correct": counts[2],
                "exactly_one_correct": counts[1],
                "none_correct": counts[0],
            },
            "whole_truth_recovery_confidence_interval": ci,
            "me_f1_confidence_interval": me_f1_ci,
            "semantic_summary": semantic_summary,
            "candidate_genotype_summary": candidate_summary,
            "addressability_audit": addressability_validation,
            "stratified_summary": (
                stratified_summary(truth_records, representative, vote_counts)
                if extended
                else {}
            ),
            "genotype_stratified_summary": genotype_strata,
            "resource_summary": resource_summary(
                getattr(args, "resource_benchmark", None),
                getattr(args, "cache_policy", "isolated_empty_tool_cache"),
            ),
        },
    }


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--truvari-ledger", required=True, type=Path)
    parser.add_argument("--aardvark-ledger", required=True, type=Path)
    parser.add_argument("--vcfdist-ledger", required=True, type=Path)
    parser.add_argument(
        "--evaluator-manifest",
        action="append",
        dest="evaluator_manifests",
        required=True,
        type=Path,
    )
    parser.add_argument(
        "--evaluator-completion",
        action="append",
        dest="evaluator_completions",
        default=[],
        type=Path,
    )
    parser.add_argument("--evaluator-profile", required=True, type=Path)
    parser.add_argument("--semantic-validation", required=True, type=Path)
    parser.add_argument("--addressability-audit", required=True, type=Path)
    parser.add_argument("--addressability-tsv", required=True, type=Path)
    parser.add_argument("--stratification-catalogue", required=True, type=Path)
    parser.add_argument(
        "--context-bed", action="append", required=True, default=[]
    )
    parser.add_argument("--allowed-information", required=True, type=Path)
    parser.add_argument("--parameter-manifest", required=True, type=Path)
    parser.add_argument("--external-resource-hashes", required=True, type=Path)
    parser.add_argument("--tuning-status", required=True, type=Path)
    parser.add_argument("--score-profile", required=True, type=Path)
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--tool-id", required=True)
    parser.add_argument("--official-score-mode", required=True)
    parser.add_argument("--primary-truth-profile", required=True)
    parser.add_argument("--truth-vcf", required=True, type=Path)
    parser.add_argument("--benchmark-bed", required=True, type=Path)
    parser.add_argument("--reference", required=True, type=Path)
    parser.add_argument("--pangenome-manifest", required=True, type=Path)
    parser.add_argument("--graph-asset-lock", type=Path)
    parser.add_argument("--query-vcf", type=Path)
    parser.add_argument("--hidden-truth-ledger", required=True, type=Path)
    parser.add_argument("--tool-manifest", type=Path)
    parser.add_argument("--resolved-inputs", type=Path)
    parser.add_argument("--sample-technology")
    parser.add_argument("--library-id")
    parser.add_argument("--source-evidence-id")
    parser.add_argument("--coverage-x", type=float)
    parser.add_argument("--read-count", type=int)
    parser.add_argument("--read-bases", type=int)
    parser.add_argument("--downsampling-seed", type=int)
    parser.add_argument("--resource-benchmark", type=Path)
    parser.add_argument("--cache-policy", default="isolated_empty_tool_cache")
    parser.add_argument("--output", required=True, type=Path)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        payload = materialize(args)
        args.output.parent.mkdir(parents=True, exist_ok=True)
        temporary = args.output.with_name(f".{args.output.name}.tmp")
        temporary.write_text(
            json.dumps(payload, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        temporary.replace(args.output)
    except (
        OSError,
        KeyError,
        TypeError,
        ValueError,
        json.JSONDecodeError,
        FormalMetricError,
    ) as exc:
        print(f"materialize_formal_metrics: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
