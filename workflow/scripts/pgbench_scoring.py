"""Deterministic, unweighted three-evaluator consensus scoring."""

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
    """Raised when consensus input violates the frozen contract."""


DEFAULT_SCORE_PROFILE_PATH = (
    Path(__file__).resolve().parents[2] / "config" / "consensus_scoring.yaml"
)
SCORE_TUPLE_FIELDS = frozenset(
    {"run_id", "sample", "tool", "official_score_mode", "primary_truth_profile"}
)
CONSENSUS_FIELDS = frozenset(
    {"all_three_correct", "exactly_two_correct", "exactly_one_correct", "none_correct"}
)
SCORE_PAYLOAD_FIELDS = frozenset(
    {
        "tuple",
        "score_profile",
        "eligibility_status",
        "infrastructure_valid",
        "truth_eligible_count",
        "consensus",
    }
)
NON_SCORABLE_STATUSES = {
    "not_applicable", "runtime_failed", "invalid_output", "incompatible_output",
    "provenance_failed", "mode_contract_violation",
}
_SHA256_RE = re.compile(r"^[0-9a-f]{64}$")


def _profile_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_score_profile(path: Path = DEFAULT_SCORE_PROFILE_PATH) -> dict[str, Any]:
    try:
        profile = yaml.safe_load(path.read_text(encoding="utf-8"))
    except (OSError, yaml.YAMLError) as exc:
        raise ScoreInputError(f"cannot load consensus profile {path}: {exc}") from exc
    if not isinstance(profile, dict) or profile.get("schema_version") != 2:
        raise ScoreInputError("consensus profile schema_version must be 2")
    meta = profile.get("profile")
    consensus = profile.get("consensus")
    if not isinstance(meta, dict) or meta.get("id") != "pgbench_consensus_v2":
        raise ScoreInputError("unsupported consensus profile id")
    if meta.get("produce_ranking") is not False:
        raise ScoreInputError("consensus profile must disable ranking")
    if not isinstance(consensus, dict):
        raise ScoreInputError("consensus profile is missing consensus settings")
    if consensus.get("evaluators") != ["truvari", "aardvark", "vcfdist"]:
        raise ScoreInputError("consensus evaluators must be Truvari, Aardvark, vcfdist")
    if consensus.get("formula") != "(3*n3 + 2*n2 + n1) / (3*N) * 100":
        raise ScoreInputError("unsupported consensus formula")
    comparability = profile.get("comparability")
    if (
        not isinstance(comparability, dict)
        or comparability.get("truth_metric") != "benchmark.truth.eligible.count"
        or comparability.get("formula") != "2*min(softTP,T)/(Q+T)*100"
    ):
        raise ScoreInputError("unsupported unified comparability contract")
    profile["_source_path"] = str(path)
    profile["_sha256"] = _profile_sha256(path)
    return profile


SCORE_PROFILE = load_score_profile()
# Kept as empty compatibility exports; weighted components no longer exist.
EVALUATOR_COMPONENTS: dict[str, dict[str, float]] = {}
EVALUATOR_WEIGHTS: dict[str, float] = {}
PANGENOME_COMPONENTS: dict[str, float] = {}
RESOURCE_COMPONENTS: dict[str, float] = {}
TRACEABILITY_COMPONENTS: dict[str, float] = {}
RESOURCE_PROFILES: dict[str, dict[str, tuple[float, float]]] = {}


@dataclass(frozen=True)
class ScoreResult:
    tuple_key: Mapping[str, str]
    score_profile: str
    score_profile_sha256: str
    evaluation_mode: str
    score_status: str
    pgbench_score: float | None
    pgbench_score_raw: float | None
    consensus_score: float | None
    comparable_score: float | None
    comparable_score_raw: float | None
    pangenome_genotyping_score: float | None
    non_reference_f1_score: float | None
    panel_coverage: float | None
    global_end_to_end_sv_recovery_score: float | None
    consensus_counts: Mapping[str, int]
    total_evaluated: int
    truth_eligible_count: int
    soft_true_positive_count: float | None
    comparable_precision: float | None
    comparable_recall: float | None
    unanimous_correct_rate: float | None
    majority_correct_rate: float | None
    point_breakdown: Mapping[str, float]
    evaluator_points: None = None
    pangenome_points: None = None
    resource_points: None = None
    traceability_points: None = None
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
        raise ScoreInputError(f"consensus.{name} must be a non-negative integer")
    return value


def resource_fraction(value: float, target: float, limit: float) -> float:
    """Removed weighted-resource helper retained only to fail clearly."""
    raise ScoreInputError("resource weighting was removed in pgbench_consensus_v2")


def calculate_pgbench_score(
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
    if status in NON_SCORABLE_STATUSES:
        return ScoreResult(
            tuple_key=trusted,
            score_profile=profile_id,
            score_profile_sha256=profile["_sha256"],
            evaluation_mode=evaluation_mode,
            score_status=str(status),
            pgbench_score=None,
            pgbench_score_raw=None,
            consensus_score=None,
            comparable_score=None,
            comparable_score_raw=None,
            pangenome_genotyping_score=None,
            non_reference_f1_score=None,
            panel_coverage=None,
            global_end_to_end_sv_recovery_score=None,
            consensus_counts={},
            total_evaluated=0,
            truth_eligible_count=0,
            soft_true_positive_count=None,
            comparable_precision=None,
            comparable_recall=None,
            unanimous_correct_rate=None,
            majority_correct_rate=None,
            point_breakdown={},
            reason=str(payload.get("reason") or status),
        )
    if status != "eligible" or payload.get("infrastructure_valid") is not True:
        raise ScoreInputError("eligible consensus input requires valid infrastructure")
    raw_counts = payload.get("consensus")
    if not isinstance(raw_counts, Mapping):
        raise ScoreInputError("consensus must be a mapping")
    _exact_fields(raw_counts, CONSENSUS_FIELDS, "consensus")
    counts = {name: _count(raw_counts[name], name) for name in CONSENSUS_FIELDS}
    total = sum(counts.values())
    n3, n2, n1 = counts["all_three_correct"], counts["exactly_two_correct"], counts["exactly_one_correct"]
    truth_total = payload.get("truth_eligible_count")
    if (
        isinstance(truth_total, bool)
        or not isinstance(truth_total, int)
        or truth_total <= 0
    ):
        raise ScoreInputError("truth_eligible_count must be a positive integer")
    vote_points = 3 * n3 + 2 * n2 + n1
    consensus_raw = vote_points / (3 * total) * 100.0 if total else 0.0
    soft_tp = vote_points / 3.0
    if soft_tp > float(truth_total) + 1e-9:
        raise ScoreInputError(
            "soft true-positive credit exceeds the one-to-one truth universe"
        )
    effective_tp = soft_tp
    comparable_raw = 2.0 * soft_tp / (total + truth_total) * 100.0
    decimals = int(profile["profile"].get("display_decimals", 2))
    consensus_score = round(consensus_raw, decimals)
    comparable_score = round(comparable_raw, decimals)
    analysis = payload.get("analysis")
    if isinstance(analysis, Mapping):
        independently_materialized = analysis.get("comparable_score")
        if (
            isinstance(independently_materialized, bool)
            or not isinstance(independently_materialized, (int, float))
            or not math.isclose(
                float(independently_materialized),
                comparable_raw,
                rel_tol=0.0,
                abs_tol=1e-9,
            )
        ):
            raise ScoreInputError(
                "formal analysis ComparableScore does not match consensus counts"
            )
    panel_score: float | None = None
    nonref_score: float | None = None
    panel_coverage: float | None = None
    detection_only_contract = False
    if isinstance(analysis, Mapping):
        candidate_summary = analysis.get("candidate_genotype_summary")
        if isinstance(candidate_summary, Mapping):
            # A discovery adapter may intentionally emit variant sites with
            # no GT assertion.  Its site-recovery score is comparable, but a
            # panel-genotyping macro-F1 is not applicable and must never be
            # promoted to the primary score as a misleading zero.
            detection_only_contract = (
                candidate_summary.get("candidate_output_contract")
                == "variant_sites"
            )
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
    if detection_only_contract:
        panel_score = None
    score_status = "valid" if evaluation_mode == "formal" else "provisional"
    traceability = payload.get("traceability")
    if isinstance(traceability, Mapping):
        if traceability.get("core_provenance_valid") is False:
            return ScoreResult(
                tuple_key=trusted,
                score_profile=profile_id,
                score_profile_sha256=profile["_sha256"],
                evaluation_mode=evaluation_mode,
                score_status="invalid",
                pgbench_score=None,
                pgbench_score_raw=None,
                consensus_score=None,
                comparable_score=None,
                comparable_score_raw=None,
                pangenome_genotyping_score=panel_score,
                non_reference_f1_score=nonref_score,
                panel_coverage=panel_coverage,
                global_end_to_end_sv_recovery_score=comparable_score,
                consensus_counts=counts,
                total_evaluated=total,
                truth_eligible_count=truth_total,
                soft_true_positive_count=soft_tp,
                comparable_precision=effective_tp / total if total else 0.0,
                comparable_recall=effective_tp / truth_total,
                unanimous_correct_rate=n3 / total if total else 0.0,
                majority_correct_rate=(n3 + n2) / total if total else 0.0,
                point_breakdown={},
                reason="core provenance validation failed",
            )
        if not all(traceability.get(k) is True for k in ("hash_lineage_complete", "environment_complete", "run_context_complete")) or traceability.get("manifest_completeness") != 1.0:
            score_status = "provisional"
    breakdown = {"all_three_correct": float(n3), "exactly_two_correct": float(n2), "exactly_one_correct": float(n1), "none_correct": float(counts["none_correct"])}
    return ScoreResult(
        tuple_key=trusted,
        score_profile=profile_id,
        score_profile_sha256=profile["_sha256"],
        evaluation_mode=evaluation_mode,
        score_status=score_status,
        # One primary definition is used for every tool contract: the mean of
        # the three frozen evaluators' binary detection votes.  Panel GT and
        # fixed-universe recovery remain explicit secondary diagnostics; they
        # are not silently substituted into the primary score.
        pgbench_score=consensus_score,
        pgbench_score_raw=consensus_raw,
        consensus_score=consensus_score,
        comparable_score=comparable_score,
        comparable_score_raw=comparable_raw,
        pangenome_genotyping_score=(
            round(panel_score, decimals) if panel_score is not None else None
        ),
        non_reference_f1_score=(
            round(nonref_score, decimals) if nonref_score is not None else None
        ),
        panel_coverage=panel_coverage,
        global_end_to_end_sv_recovery_score=comparable_score,
        consensus_counts=counts,
        total_evaluated=total,
        truth_eligible_count=truth_total,
        soft_true_positive_count=soft_tp,
        comparable_precision=effective_tp / total if total else 0.0,
        comparable_recall=effective_tp / truth_total,
        unanimous_correct_rate=n3 / total if total else 0.0,
        majority_correct_rate=(n3 + n2) / total if total else 0.0,
        point_breakdown=breakdown,
        evaluator_scores={},
        required_f1_metrics={},
        formal_analysis=dict(analysis) if isinstance(analysis, Mapping) else None,
    )
