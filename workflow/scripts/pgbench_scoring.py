"""Deterministic multi-evaluator genotype-F1 scoring.

The primary score is the arithmetic mean of three independently materialized,
genotype-aware evaluator F1 values. Evaluator agreement is diagnostic only and
can never affect the primary score.
"""

from __future__ import annotations

import hashlib
import math
import re
from collections.abc import Mapping
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

import yaml  # type: ignore[import-untyped]


class ScoreInputError(ValueError):
    """Raised when ME-F1 input violates the frozen contract."""


DEFAULT_SCORE_PROFILE_PATH = (
    Path(__file__).resolve().parents[2] / "config" / "me_f1_scoring.yaml"
)
SCORE_TUPLE_FIELDS = frozenset(
    {"run_id", "sample", "tool", "official_score_mode", "primary_truth_profile"}
)
EVALUATORS = ("truvari", "aardvark", "vcfdist")
EVALUATOR_FIELDS = frozenset({"tp", "fp", "fn", "precision", "recall", "f1"})
QUALITY_GATE_FIELDS = frozenset(
    {
        "evaluator_semantic_contract_valid",
        "addressability_complete",
        "allowed_information_complete",
        "tuning_frozen",
        "context_stratification_complete",
        "population_af_stratification_complete",
    }
)
SCORE_PAYLOAD_FIELDS = frozenset(
    {
        "tuple",
        "score_profile",
        "eligibility_status",
        "infrastructure_valid",
        "truth_eligible_count",
        "evaluator_f1",
        "quality_gates",
    }
)
NON_SCORABLE_STATUSES = {
    "not_applicable", "runtime_failed", "invalid_output", "incompatible_output",
    "provenance_failed", "mode_contract_violation", "invalid_evaluator_mapping",
    "invalid_empty_submission",
}
INVALID_SCORE_STATUSES = {"invalid_evaluator_mapping", "invalid_empty_submission"}
_SHA256_RE = re.compile(r"^[0-9a-f]{64}$")


def _profile_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_score_profile(path: Path = DEFAULT_SCORE_PROFILE_PATH) -> dict[str, Any]:
    try:
        profile = yaml.safe_load(path.read_text(encoding="utf-8"))
    except (OSError, yaml.YAMLError) as exc:
        raise ScoreInputError(f"cannot load ME-F1 profile {path}: {exc}") from exc
    if not isinstance(profile, dict) or profile.get("schema_version") != 3:
        raise ScoreInputError("ME-F1 profile schema_version must be 3")
    meta = profile.get("profile")
    if not isinstance(meta, dict) or meta.get("id") != "pgbench_me_f1_v1":
        raise ScoreInputError("unsupported ME-F1 profile id")
    primary = profile.get("primary_score")
    if not isinstance(primary, dict):
        raise ScoreInputError("ME-F1 profile is missing primary_score")
    if primary.get("formula") != (
        "(F1_truvari + F1_aardvark_gt + F1_vcfdist) / 3"
    ):
        raise ScoreInputError("unsupported ME-F1 formula")
    if primary.get("denominator") != "frozen_panel_addressable_truth":
        raise ScoreInputError("ME-F1 must use frozen panel-addressable truth")
    weights = primary.get("evaluator_weights")
    if weights != {name: 1 for name in EVALUATORS}:
        raise ScoreInputError("ME-F1 requires equal evaluator weights")
    profile["_source_path"] = str(path)
    profile["_sha256"] = _profile_sha256(path)
    return profile


SCORE_PROFILE = load_score_profile()


@dataclass(frozen=True)
class ScoreResult:
    tuple_key: Mapping[str, str]
    score_profile: str
    score_profile_sha256: str
    evaluation_mode: str
    score_status: str
    benchmark_score: float | None
    benchmark_score_raw: float | None
    me_f1: float | None
    evaluator_range: float | None
    evaluator_sd: float | None
    pangenome_genotyping_score: float | None
    non_reference_f1_score: float | None
    panel_coverage: float | None
    evaluator_agreement_counts: Mapping[str, int]
    total_evaluated: int
    truth_eligible_count: int
    point_breakdown: Mapping[str, float]
    evaluator_scores: Mapping[str, float] | None = None
    required_f1_metrics: Mapping[str, Mapping[str, Any]] | None = None
    formal_analysis: Mapping[str, Any] | None = None
    reason: str | None = None

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def _reject_ranking_fields(value: Any, path: str = "payload") -> None:
    if isinstance(value, Mapping):
        for key, item in value.items():
            normalized = str(key).casefold()
            if normalized in {"rank", "ranking", "leaderboard"} or normalized.startswith("rank_"):
                raise ScoreInputError(f"forbidden ranking field: {path}.{key}")
            _reject_ranking_fields(item, f"{path}.{key}")
    elif isinstance(value, (list, tuple)):
        for index, item in enumerate(value):
            _reject_ranking_fields(item, f"{path}[{index}]")


def _exact_fields(value: Mapping[str, Any], expected: set[str] | frozenset[str], context: str) -> None:
    if set(value) != set(expected):
        missing = sorted(set(expected) - set(value))
        extra = sorted(set(value) - set(expected))
        details = (["missing " + ", ".join(missing)] if missing else []) + (["unexpected " + ", ".join(extra)] if extra else [])
        raise ScoreInputError(f"{context} has invalid fields: " + "; ".join(details))


def _validated_tuple(value: Any, context: str) -> dict[str, str]:
    if not isinstance(value, Mapping):
        raise ScoreInputError(f"{context} must be a mapping")
    _exact_fields(value, SCORE_TUPLE_FIELDS, context)
    result = dict(value)
    if not all(isinstance(v, str) and v for v in result.values()):
        raise ScoreInputError(f"{context} fields must be non-empty strings")
    return result  # type: ignore[return-value]


def _count(value: Any, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        raise ScoreInputError(f"evaluator_agreement.{name} must be a non-negative integer")
    return value


def _unit_interval(value: Any, name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ScoreInputError(f"{name} must be numeric")
    numeric = float(value)
    if not math.isfinite(numeric) or not 0.0 <= numeric <= 1.0:
        raise ScoreInputError(f"{name} must be finite in [0, 1]")
    return numeric


def _validated_evaluator_metrics(value: Any) -> dict[str, dict[str, float | int]]:
    if not isinstance(value, Mapping) or set(value) != set(EVALUATORS):
        raise ScoreInputError(
            "evaluator_f1 must contain exactly truvari, aardvark, and vcfdist"
        )
    validated: dict[str, dict[str, float | int]] = {}
    for evaluator in EVALUATORS:
        raw = value[evaluator]
        if not isinstance(raw, Mapping):
            raise ScoreInputError(f"evaluator_f1.{evaluator} must be a mapping")
        _exact_fields(raw, EVALUATOR_FIELDS, f"evaluator_f1.{evaluator}")
        tp = _count(raw["tp"], f"{evaluator}.tp")
        fp = _count(raw["fp"], f"{evaluator}.fp")
        fn = _count(raw["fn"], f"{evaluator}.fn")
        precision = _unit_interval(raw["precision"], f"{evaluator}.precision")
        recall = _unit_interval(raw["recall"], f"{evaluator}.recall")
        f1 = _unit_interval(raw["f1"], f"{evaluator}.f1")
        expected_precision = tp / (tp + fp) if tp + fp else 0.0
        expected_recall = tp / (tp + fn) if tp + fn else 0.0
        expected_f1 = (
            2.0 * expected_precision * expected_recall
            / (expected_precision + expected_recall)
            if expected_precision + expected_recall
            else 0.0
        )
        for name, observed, expected in (
            ("precision", precision, expected_precision),
            ("recall", recall, expected_recall),
            ("f1", f1, expected_f1),
        ):
            if not math.isclose(observed, expected, rel_tol=0.0, abs_tol=1e-9):
                raise ScoreInputError(
                    f"evaluator_f1.{evaluator}.{name} does not match TP/FP/FN"
                )
        validated[evaluator] = {
            "tp": tp,
            "fp": fp,
            "fn": fn,
            "precision": precision,
            "recall": recall,
            "f1": f1,
        }
    return validated


def _validated_quality_gates(value: Any) -> dict[str, bool]:
    if not isinstance(value, Mapping):
        raise ScoreInputError("quality_gates must be a mapping")
    _exact_fields(value, QUALITY_GATE_FIELDS, "quality_gates")
    if not all(isinstance(value[field], bool) for field in QUALITY_GATE_FIELDS):
        raise ScoreInputError("all quality_gates values must be boolean")
    return {field: bool(value[field]) for field in QUALITY_GATE_FIELDS}


def calculate_me_f1(
    payload: Mapping[str, Any], *, expected_tuple: Mapping[str, str],
    evaluation_mode: str, score_profile_path: Path | None = None,
) -> ScoreResult:
    _reject_ranking_fields(payload)
    allowed = SCORE_PAYLOAD_FIELDS | {"reason", "traceability", "analysis"}
    missing = SCORE_PAYLOAD_FIELDS - set(payload)
    extra = set(payload) - allowed
    if missing or extra:
        details = (["missing " + ", ".join(sorted(missing))] if missing else []) + (["unexpected " + ", ".join(sorted(extra))] if extra else [])
        raise ScoreInputError("payload has invalid fields: " + "; ".join(details))
    trusted = _validated_tuple(expected_tuple, "expected_tuple")
    if _validated_tuple(payload.get("tuple"), "payload.tuple") != trusted:
        raise ScoreInputError("payload tuple does not match trusted expected_tuple")
    profile = load_score_profile(score_profile_path or DEFAULT_SCORE_PROFILE_PATH)
    profile_id = profile["profile"]["id"]
    if payload.get("score_profile") != profile_id:
        raise ScoreInputError("payload score_profile does not match selected profile")
    if evaluation_mode not in {"formal", "synthetic_smoke"}:
        raise ScoreInputError("unsupported evaluation mode")
    status = payload.get("eligibility_status")
    analysis = payload.get("analysis")
    if status in NON_SCORABLE_STATUSES:
        return ScoreResult(
            tuple_key=trusted,
            score_profile=profile_id,
            score_profile_sha256=profile["_sha256"],
            evaluation_mode=evaluation_mode,
            score_status="invalid" if status in INVALID_SCORE_STATUSES else str(status),
            benchmark_score=None,
            benchmark_score_raw=None,
            me_f1=None,
            evaluator_range=None,
            evaluator_sd=None,
            pangenome_genotyping_score=None,
            non_reference_f1_score=None,
            panel_coverage=None,
            evaluator_agreement_counts={},
            total_evaluated=0,
            truth_eligible_count=0,
            point_breakdown={},
            formal_analysis=(dict(analysis) if isinstance(analysis, Mapping) else None),
            reason=str(payload.get("reason") or status),
        )
    if status != "eligible" or payload.get("infrastructure_valid") is not True:
        raise ScoreInputError("eligible ME-F1 input requires valid infrastructure")
    evaluator_metrics = _validated_evaluator_metrics(payload.get("evaluator_f1"))
    quality_gates = _validated_quality_gates(payload.get("quality_gates"))
    truth_total = payload.get("truth_eligible_count")
    if (
        isinstance(truth_total, bool)
        or not isinstance(truth_total, int)
        or truth_total <= 0
    ):
        raise ScoreInputError("truth_eligible_count must be a positive integer")
    for evaluator, metrics in evaluator_metrics.items():
        if int(metrics["tp"]) + int(metrics["fn"]) != truth_total:
            raise ScoreInputError(
                f"evaluator_f1.{evaluator} does not use the frozen truth denominator"
            )
    decimals = int(profile["profile"].get("display_decimals", 2))
    f1_values = [float(evaluator_metrics[name]["f1"]) for name in EVALUATORS]
    me_f1_raw = sum(f1_values) / len(f1_values) * 100.0
    me_f1 = round(me_f1_raw, decimals)
    evaluator_range = round((max(f1_values) - min(f1_values)) * 100.0, decimals)
    evaluator_mean = sum(f1_values) / len(f1_values)
    evaluator_sd = round(
        math.sqrt(sum((value - evaluator_mean) ** 2 for value in f1_values) / 3)
        * 100.0,
        decimals,
    )
    total = max(
        int(metrics["tp"]) + int(metrics["fp"])
        for metrics in evaluator_metrics.values()
    )
    counts = (
        dict(analysis.get("evaluator_agreement_counts", {}))
        if isinstance(analysis, Mapping)
        and isinstance(analysis.get("evaluator_agreement_counts"), Mapping)
        else {}
    )
    panel_score: float | None = None
    nonref_score: float | None = None
    panel_coverage: float | None = None
    if isinstance(analysis, Mapping):
        candidate_summary = analysis.get("candidate_genotype_summary")
        if (
            isinstance(candidate_summary, Mapping)
            and candidate_summary.get("candidate_output_contract") != "all_sites"
        ):
            raise ScoreInputError("formal ME-F1 requires an all-sites output contract")
        for field, target in (
            ("pangenome_genotyping_score", "panel"),
            ("non_reference_f1_score", "nonref"),
            ("panel_coverage", "coverage"),
        ):
            value = analysis.get(field)
            if value is not None and (
                isinstance(value, bool) or not isinstance(value, (int, float))
            ):
                raise ScoreInputError(f"formal analysis {field} must be numeric or null")
            if value is not None:
                numeric = float(value)
                upper = 1.0 if target == "coverage" else 100.0
                if not math.isfinite(numeric) or not 0.0 <= numeric <= upper:
                    raise ScoreInputError(f"formal analysis {field} is out of range")
                if target == "panel":
                    panel_score = numeric
                elif target == "nonref":
                    nonref_score = numeric
                else:
                    panel_coverage = numeric
    score_status = "valid" if evaluation_mode == "formal" else "provisional"
    gate_failures = sorted(field for field, passed in quality_gates.items() if not passed)
    traceability = payload.get("traceability")
    provenance_failed = (
        isinstance(traceability, Mapping)
        and traceability.get("core_provenance_valid") is False
    )
    if gate_failures or provenance_failed:
            return ScoreResult(
                tuple_key=trusted,
                score_profile=profile_id,
                score_profile_sha256=profile["_sha256"],
                evaluation_mode=evaluation_mode,
                score_status="invalid",
                benchmark_score=None,
                benchmark_score_raw=None,
                me_f1=None,
                evaluator_range=None,
                evaluator_sd=None,
                pangenome_genotyping_score=None,
                non_reference_f1_score=None,
                panel_coverage=None,
                evaluator_agreement_counts=counts,
                total_evaluated=total,
                truth_eligible_count=truth_total,
                point_breakdown={},
                evaluator_scores={
                    name: float(metrics["f1"]) * 100.0
                    for name, metrics in evaluator_metrics.items()
                },
                required_f1_metrics=evaluator_metrics,
                formal_analysis=(dict(analysis) if isinstance(analysis, Mapping) else None),
                reason=(
                    "quality gates failed: " + ", ".join(gate_failures)
                    if gate_failures
                    else "core provenance validation failed"
                ),
            )
    if isinstance(traceability, Mapping):
        if not all(traceability.get(k) is True for k in ("hash_lineage_complete", "environment_complete", "run_context_complete")) or traceability.get("manifest_completeness") != 1.0:
            score_status = "provisional"
    breakdown = {
        f"{name}_f1": float(metrics["f1"]) * 100.0
        for name, metrics in evaluator_metrics.items()
    }
    return ScoreResult(
        tuple_key=trusted,
        score_profile=profile_id,
        score_profile_sha256=profile["_sha256"],
        evaluation_mode=evaluation_mode,
        score_status=score_status,
        benchmark_score=me_f1,
        benchmark_score_raw=me_f1_raw,
        me_f1=me_f1,
        evaluator_range=evaluator_range,
        evaluator_sd=evaluator_sd,
        pangenome_genotyping_score=(
            round(panel_score, decimals) if panel_score is not None else None
        ),
        non_reference_f1_score=(
            round(nonref_score, decimals) if nonref_score is not None else None
        ),
        panel_coverage=panel_coverage,
        evaluator_agreement_counts=counts,
        total_evaluated=total,
        truth_eligible_count=truth_total,
        point_breakdown=breakdown,
        evaluator_scores={
            name: round(float(metrics["f1"]) * 100.0, decimals)
            for name, metrics in evaluator_metrics.items()
        },
        required_f1_metrics=evaluator_metrics,
        formal_analysis=dict(analysis) if isinstance(analysis, Mapping) else None,
    )
