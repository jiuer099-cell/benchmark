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
    {"tuple", "score_profile", "eligibility_status", "infrastructure_valid", "consensus"}
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
    consensus_counts: Mapping[str, int]
    total_evaluated: int
    unanimous_correct_rate: float | None
    majority_correct_rate: float | None
    point_breakdown: Mapping[str, float]
    evaluator_points: None = None
    pangenome_points: None = None
    resource_points: None = None
    traceability_points: None = None
    evaluator_scores: Mapping[str, float] | None = None
    required_f1_metrics: Mapping[str, Mapping[str, Any]] | None = None
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
    allowed = SCORE_PAYLOAD_FIELDS | {"reason", "traceability"}
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
        return ScoreResult(trusted, profile_id, profile["_sha256"], evaluation_mode, str(status), None, None, None, {}, 0, None, None, {}, reason=str(payload.get("reason") or status))
    if status != "eligible" or payload.get("infrastructure_valid") is not True:
        raise ScoreInputError("eligible consensus input requires valid infrastructure")
    raw_counts = payload.get("consensus")
    if not isinstance(raw_counts, Mapping):
        raise ScoreInputError("consensus must be a mapping")
    _exact_fields(raw_counts, CONSENSUS_FIELDS, "consensus")
    counts = {name: _count(raw_counts[name], name) for name in CONSENSUS_FIELDS}
    total = sum(counts.values())
    if total == 0:
        raise ScoreInputError("consensus must contain at least one evaluated result")
    n3, n2, n1 = counts["all_three_correct"], counts["exactly_two_correct"], counts["exactly_one_correct"]
    raw = (3 * n3 + 2 * n2 + n1) / (3 * total) * 100.0
    score = round(raw, int(profile["profile"].get("display_decimals", 2)))
    score_status = "valid" if evaluation_mode == "formal" else "provisional"
    traceability = payload.get("traceability")
    if isinstance(traceability, Mapping):
        if traceability.get("core_provenance_valid") is False:
            return ScoreResult(trusted, profile_id, profile["_sha256"], evaluation_mode, "invalid", None, None, None, counts, total, n3 / total, (n3 + n2) / total, {}, reason="core provenance validation failed")
        if not all(traceability.get(k) is True for k in ("hash_lineage_complete", "environment_complete", "run_context_complete")) or traceability.get("manifest_completeness") != 1.0:
            score_status = "provisional"
    breakdown = {"all_three_correct": float(n3), "exactly_two_correct": float(n2), "exactly_one_correct": float(n1), "none_correct": float(counts["none_correct"])}
    return ScoreResult(trusted, profile_id, profile["_sha256"], evaluation_mode, score_status, score, raw, score, counts, total, n3 / total, (n3 + n2) / total, breakdown, evaluator_scores={}, required_f1_metrics={})
