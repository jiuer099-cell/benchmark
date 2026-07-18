"""Deterministic PGBenchScore calculation.

The functions in this module deliberately do not compare tools with one
another. Every score is calculated from one fixed tool/run tuple and the
versioned ``pgbench_v1`` point profile.
"""

from __future__ import annotations

import hashlib
import math
import re
from collections.abc import Mapping
from copy import deepcopy
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

import yaml  # type: ignore[import-untyped]


class ScoreInputError(ValueError):
    """Raised when a score input violates the fixed metric contract."""


DEFAULT_SCORE_PROFILE_PATH = (
    Path(__file__).resolve().parents[2] / "config" / "score_weights.yaml"
)
SCORE_TUPLE_FIELDS = frozenset(
    {
        "run_id",
        "sample",
        "tool",
        "official_score_mode",
        "primary_truth_profile",
    }
)
SCORE_PAYLOAD_FIELDS = frozenset(
    {
        "tuple",
        "score_profile",
        "eligibility_status",
        "infrastructure_valid",
        "evaluators",
        "pangenome",
        "resources",
        "traceability",
        "required_f1_metrics",
    }
)
_SHA256_RE = re.compile(r"^[0-9a-f]{64}$")


def _load_yaml_mapping(path: Path) -> dict[str, Any]:
    try:
        with path.open("r", encoding="utf-8") as handle:
            loaded = yaml.safe_load(handle)
    except (OSError, yaml.YAMLError) as exc:
        raise ScoreInputError(f"cannot load score profile {path}: {exc}") from exc
    if not isinstance(loaded, dict):
        raise ScoreInputError("score profile must be a YAML mapping")
    return loaded


def _profile_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    try:
        with path.open("rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(chunk)
    except OSError as exc:
        raise ScoreInputError(f"cannot hash score profile {path}: {exc}") from exc
    return digest.hexdigest()


def load_score_profile(path: Path = DEFAULT_SCORE_PROFILE_PATH) -> dict[str, Any]:
    """Load and structurally validate the versioned point profile."""

    profile = _load_yaml_mapping(path)
    if profile.get("schema_version") != 1:
        raise ScoreInputError("score profile schema_version must be 1")
    profile_meta = profile.get("profile")
    layers = profile.get("layers")
    if not isinstance(profile_meta, dict) or not isinstance(layers, dict):
        raise ScoreInputError("score profile requires profile and layers mappings")
    if profile_meta.get("id") != "pgbench_v1":
        raise ScoreInputError("unsupported score profile id")
    if profile_meta.get("produce_ranking") is not False:
        raise ScoreInputError("score profile must set produce_ranking to false")
    profile["_source_path"] = str(path)
    profile["_sha256"] = _profile_sha256(path)
    return profile


def _component_maxima(layer: Mapping[str, Any], name: str) -> dict[str, float]:
    raw_components = layer.get("components")
    if not isinstance(raw_components, Mapping):
        raise ScoreInputError(f"{name}.components must be a mapping")
    result: dict[str, float] = {}
    for component_name, component in raw_components.items():
        if not isinstance(component_name, str) or not isinstance(component, Mapping):
            raise ScoreInputError(f"invalid component in {name}")
        maximum = component.get("max_points")
        if isinstance(maximum, bool) or not isinstance(maximum, int | float):
            raise ScoreInputError(f"{name}.{component_name}.max_points must be numeric")
        numeric = float(maximum)
        if not math.isfinite(numeric) or numeric < 0:
            raise ScoreInputError(
                f"{name}.{component_name}.max_points must be finite and non-negative"
            )
        result[component_name] = numeric
    return result


def _derive_profile_constants(
    profile: Mapping[str, Any],
) -> tuple[
    dict[str, dict[str, float]],
    dict[str, float],
    dict[str, float],
    dict[str, float],
    dict[str, float],
    dict[str, dict[str, tuple[float, float]]],
]:
    layers = profile["layers"]
    evaluator_layer = layers["evaluator_accuracy"]
    evaluators = evaluator_layer["evaluators"]
    evaluator_components: dict[str, dict[str, float]] = {}
    evaluator_weights: dict[str, float] = {}
    for evaluator, evaluator_profile in evaluators.items():
        evaluator_components[evaluator] = _component_maxima(
            evaluator_profile, f"evaluator_accuracy.{evaluator}"
        )
        evaluator_weights[evaluator] = float(evaluator_profile["fusion_weight"])

    pangenome_components = _component_maxima(
        layers["pangenome_robustness"], "pangenome_robustness"
    )
    resource_layer = layers["resource_efficiency"]
    resource_components = _component_maxima(resource_layer, "resource_efficiency")
    traceability_components = _component_maxima(layers["traceability"], "traceability")
    resource_profiles: dict[str, dict[str, tuple[float, float]]] = {}
    for mode, budget in resource_layer["budgets"].items():
        resource_profiles[mode] = {
            metric: (float(values["target"]), float(values["limit"]))
            for metric, values in budget.items()
        }
    return (
        evaluator_components,
        evaluator_weights,
        pangenome_components,
        resource_components,
        traceability_components,
        resource_profiles,
    )


SCORE_PROFILE = load_score_profile()
(
    EVALUATOR_COMPONENTS,
    EVALUATOR_WEIGHTS,
    PANGENOME_COMPONENTS,
    RESOURCE_COMPONENTS,
    TRACEABILITY_COMPONENTS,
    RESOURCE_PROFILES,
) = _derive_profile_constants(SCORE_PROFILE)

NON_SCORABLE_STATUSES = {
    "not_applicable",
    "runtime_failed",
    "invalid_output",
    "incompatible_output",
    "provenance_failed",
    "mode_contract_violation",
}


@dataclass(frozen=True)
class ScoreResult:
    """One score result for one formal benchmark tuple."""

    tuple_key: Mapping[str, str]
    score_profile: str
    score_profile_sha256: str
    evaluation_mode: str
    score_status: str
    pgbench_score: float | None
    pgbench_score_raw: float | None
    evaluator_points: float | None
    pangenome_points: float | None
    resource_points: float | None
    traceability_points: float | None
    point_breakdown: Mapping[str, float]
    evaluator_scores: Mapping[str, float]
    required_f1_metrics: Mapping[str, Mapping[str, Any]]
    reason: str | None = None

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-serializable representation without ranking fields."""

        return asdict(self)


def _unit_interval(value: Any, name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, int | float):
        raise ScoreInputError(f"{name} must be a numeric value in [0, 1]")
    numeric = float(value)
    epsilon = 1e-12
    if not math.isfinite(numeric):
        raise ScoreInputError(f"{name} must be finite")
    if numeric < -epsilon or numeric > 1.0 + epsilon:
        raise ScoreInputError(f"{name}={numeric} is outside [0, 1]")
    return min(1.0, max(0.0, numeric))


def _positive_number(value: Any, name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, int | float):
        raise ScoreInputError(f"{name} must be a positive number")
    numeric = float(value)
    if not math.isfinite(numeric) or numeric <= 0:
        raise ScoreInputError(f"{name} must be finite and greater than zero")
    return numeric


def _required_mapping(
    parent: Mapping[str, Any], key: str, *, context: str
) -> Mapping[str, Any]:
    value = parent.get(key)
    if not isinstance(value, Mapping):
        raise ScoreInputError(f"{context}.{key} must be a mapping")
    return value


def _ranking_field_name(value: Any) -> bool:
    if not isinstance(value, str):
        return False
    normalized = value.casefold()
    return (
        normalized in {"rank", "ranking", "leaderboard", "leaderboard_position"}
        or normalized.startswith("rank_")
        or normalized.startswith("leaderboard_")
    )


def _reject_ranking_fields(value: Any, *, path: str = "payload") -> None:
    if isinstance(value, Mapping):
        for key, item in value.items():
            if _ranking_field_name(key):
                raise ScoreInputError(f"forbidden ranking field: {path}.{key}")
            _reject_ranking_fields(item, path=f"{path}.{key}")
    elif isinstance(value, list | tuple):
        for index, item in enumerate(value):
            _reject_ranking_fields(item, path=f"{path}[{index}]")


def _exact_fields(
    mapping: Mapping[str, Any], expected: set[str] | frozenset[str], *, context: str
) -> None:
    actual = set(mapping)
    if actual == set(expected):
        return
    details: list[str] = []
    missing = sorted(set(expected) - actual)
    extra = sorted(actual - set(expected))
    if missing:
        details.append("missing " + ", ".join(missing))
    if extra:
        details.append("unexpected " + ", ".join(extra))
    raise ScoreInputError(f"{context} has invalid fields: " + "; ".join(details))


def _validated_tuple(value: Mapping[str, Any], *, context: str) -> dict[str, str]:
    _exact_fields(value, SCORE_TUPLE_FIELDS, context=context)
    normalized: dict[str, str] = {}
    for field in SCORE_TUPLE_FIELDS:
        item = value.get(field)
        if not isinstance(item, str) or not item:
            raise ScoreInputError(f"{context}.{field} must be a non-empty string")
        normalized[field] = item
    return normalized


def _nonnegative_number(value: Any, name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, int | float):
        raise ScoreInputError(f"{name} must be a non-negative finite number")
    numeric = float(value)
    if not math.isfinite(numeric) or numeric < 0:
        raise ScoreInputError(f"{name} must be a non-negative finite number")
    return numeric


def _required_f1_records(
    payload: Mapping[str, Any], profile: Mapping[str, Any]
) -> dict[str, Mapping[str, Any]]:
    raw = _required_mapping(payload, "required_f1_metrics", context="payload")
    configured = profile.get("required_f1_metrics")
    if not isinstance(configured, list) or not all(
        isinstance(metric_id, str) and metric_id for metric_id in configured
    ):
        raise ScoreInputError("score profile required_f1_metrics is invalid")
    _exact_fields(raw, set(configured), context="required_f1_metrics")
    validated: dict[str, Mapping[str, Any]] = {}
    for metric_id in configured:
        record = raw.get(metric_id)
        if not isinstance(record, Mapping):
            raise ScoreInputError(f"required_f1_metrics.{metric_id} must be a mapping")
        if record.get("metric_id") != metric_id:
            raise ScoreInputError(
                f"required_f1_metrics.{metric_id}.metric_id does not match its key"
            )
        if record.get("value_type") != "ratio":
            raise ScoreInputError(
                f"required_f1_metrics.{metric_id}.value_type must be ratio"
            )
        if record.get("status") not in {
            "defined",
            "completed_with_tool_ineligible_zero",
        }:
            raise ScoreInputError(
                f"required_f1_metrics.{metric_id} has a non-scorable status"
            )
        _unit_interval(record.get("value"), f"required_f1_metrics.{metric_id}.value")
        for field in ("numerator", "denominator", "eligible_count"):
            _nonnegative_number(
                record.get(field), f"required_f1_metrics.{metric_id}.{field}"
            )
        manifest_id = record.get("provenance_manifest_id")
        if not isinstance(manifest_id, str) or not _SHA256_RE.fullmatch(manifest_id):
            raise ScoreInputError(
                f"required_f1_metrics.{metric_id}.provenance_manifest_id must be "
                "a lowercase 64-character SHA-256"
            )
        if record.get("strata") != {}:
            raise ScoreInputError(
                f"required_f1_metrics.{metric_id} must be an empty-strata aggregate"
            )
        validated[metric_id] = deepcopy(dict(record))
    return validated


def resource_fraction(value: float, target: float, limit: float) -> float:
    """Return the fixed log-interpolated resource fraction."""

    value = _positive_number(value, "resource value")
    target = _positive_number(target, "resource target")
    limit = _positive_number(limit, "resource limit")
    if target >= limit:
        raise ScoreInputError("resource target must be less than limit")
    if value <= target:
        return 1.0
    if value >= limit:
        return 0.0
    return (math.log(limit) - math.log(value)) / (math.log(limit) - math.log(target))


def _component_points(
    metrics: Mapping[str, Any],
    components: Mapping[str, float],
    *,
    prefix: str,
) -> tuple[float, dict[str, float]]:
    _exact_fields(metrics, set(components), context=prefix)
    points: dict[str, float] = {}
    for metric_name, maximum in components.items():
        value = _unit_interval(metrics.get(metric_name), f"{prefix}.{metric_name}")
        points[f"{prefix}.{metric_name}"] = maximum * value
    return sum(points.values()), points


def _evaluator_points(
    metrics: Mapping[str, Any],
    evaluator_components: Mapping[str, Mapping[str, float]],
) -> tuple[float, dict[str, float], dict[str, float]]:
    _exact_fields(metrics, set(evaluator_components), context="metrics.evaluators")
    total = 0.0
    breakdown: dict[str, float] = {}
    evaluator_scores: dict[str, float] = {}
    for evaluator, components in evaluator_components.items():
        evaluator_metrics = _required_mapping(
            metrics, evaluator, context="metrics.evaluators"
        )
        evaluator_points, component_breakdown = _component_points(
            evaluator_metrics,
            components,
            prefix=f"evaluator.{evaluator}",
        )
        maximum = sum(components.values())
        evaluator_scores[evaluator] = 100.0 * evaluator_points / maximum
        breakdown.update(component_breakdown)
        total += evaluator_points
    return total, breakdown, evaluator_scores


def _resource_points(
    resources: Mapping[str, Any],
    official_mode: str,
    resource_components: Mapping[str, float],
    resource_profiles: Mapping[str, Mapping[str, tuple[float, float]]],
) -> tuple[float, dict[str, float]]:
    if official_mode not in resource_profiles:
        raise ScoreInputError(f"unsupported official mode: {official_mode}")
    _exact_fields(resources, set(resource_components), context="resources")
    profile = resource_profiles[official_mode]
    breakdown: dict[str, float] = {}
    for metric_name, maximum in resource_components.items():
        value = _positive_number(resources.get(metric_name), f"resources.{metric_name}")
        target, limit = profile[metric_name]
        fraction = resource_fraction(value, target, limit)
        breakdown[f"resource.{metric_name}"] = maximum * fraction
    return sum(breakdown.values()), breakdown


def _traceability_points(
    traceability: Mapping[str, Any],
    traceability_components: Mapping[str, float],
) -> tuple[float, dict[str, float], bool]:
    _exact_fields(
        traceability,
        set(traceability_components) | {"core_provenance_valid"},
        context="traceability",
    )
    manifest_completeness = _unit_interval(
        traceability.get("manifest_completeness"),
        "traceability.manifest_completeness",
    )
    gates: dict[str, bool] = {}
    for gate_name in (
        "hash_lineage_complete",
        "environment_complete",
        "run_context_complete",
    ):
        gate_value = traceability.get(gate_name)
        if not isinstance(gate_value, bool):
            raise ScoreInputError(f"traceability.{gate_name} must be boolean")
        gates[gate_name] = gate_value
    core_valid = traceability.get("core_provenance_valid")
    if not isinstance(core_valid, bool):
        raise ScoreInputError("traceability.core_provenance_valid must be boolean")

    breakdown = {
        "traceability.manifest_completeness": traceability_components[
            "manifest_completeness"
        ]
        * manifest_completeness,
        "traceability.hash_lineage_complete": traceability_components[
            "hash_lineage_complete"
        ]
        * float(gates["hash_lineage_complete"]),
        "traceability.environment_complete": traceability_components[
            "environment_complete"
        ]
        * float(gates["environment_complete"]),
        "traceability.run_context_complete": traceability_components[
            "run_context_complete"
        ]
        * float(gates["run_context_complete"]),
    }
    return sum(breakdown.values()), breakdown, core_valid


def _validate_payload_fields(payload: Mapping[str, Any]) -> None:
    actual = set(payload)
    missing = sorted(SCORE_PAYLOAD_FIELDS - actual)
    extra = sorted(actual - (SCORE_PAYLOAD_FIELDS | {"reason"}))
    if "resource_profile" in actual:
        raise ScoreInputError(
            "resource_profile overrides are forbidden; resource budgets come only "
            "from the versioned score profile"
        )
    if not missing and not extra:
        return
    details: list[str] = []
    if missing:
        details.append("missing " + ", ".join(missing))
    if extra:
        details.append("unexpected " + ", ".join(extra))
    raise ScoreInputError("payload has invalid fields: " + "; ".join(details))


def _validate_f1_component_binding(
    required_f1_metrics: Mapping[str, Mapping[str, Any]],
    payload: Mapping[str, Any],
    profile: Mapping[str, Any],
) -> None:
    bound_values: dict[str, float] = {}
    evaluators = _required_mapping(payload, "evaluators", context="payload")
    pangenome = _required_mapping(payload, "pangenome", context="payload")
    layers = profile["layers"]
    for evaluator, evaluator_profile in layers["evaluator_accuracy"][
        "evaluators"
    ].items():
        evaluator_values = _required_mapping(
            evaluators, evaluator, context="payload.evaluators"
        )
        for component_name, component in evaluator_profile["components"].items():
            metric_id = component["metric"]
            if metric_id in required_f1_metrics:
                bound_values[metric_id] = _unit_interval(
                    evaluator_values.get(component_name),
                    f"evaluators.{evaluator}.{component_name}",
                )
    for component_name, component in layers["pangenome_robustness"][
        "components"
    ].items():
        metric_id = component["metric"]
        if metric_id in required_f1_metrics:
            bound_values[metric_id] = _unit_interval(
                pangenome.get(component_name), f"pangenome.{component_name}"
            )
    missing = sorted(set(required_f1_metrics) - set(bound_values))
    if missing:
        raise ScoreInputError(
            "required F1 metrics are not bound to score components: "
            + ", ".join(missing)
        )
    for metric_id, record in required_f1_metrics.items():
        raw_value = _unit_interval(
            record.get("value"), f"required_f1_metrics.{metric_id}.value"
        )
        if not math.isclose(raw_value, bound_values[metric_id], abs_tol=1e-12):
            raise ScoreInputError(
                f"required F1 metric {metric_id} does not match its score component"
            )


def calculate_pgbench_score(
    payload: Mapping[str, Any],
    *,
    expected_tuple: Mapping[str, str],
    evaluation_mode: str,
    score_profile_path: Path | None = None,
) -> ScoreResult:
    """Calculate a score bound to a trusted tuple and evaluation mode."""

    _reject_ranking_fields(payload)
    _validate_payload_fields(payload)
    trusted_tuple = _validated_tuple(expected_tuple, context="expected_tuple")
    if evaluation_mode not in {"formal", "synthetic_smoke"}:
        raise ScoreInputError(f"unsupported evaluation_mode: {evaluation_mode}")
    selected_profile = (
        SCORE_PROFILE
        if score_profile_path is None
        else load_score_profile(score_profile_path)
    )
    (
        evaluator_components,
        evaluator_weights,
        pangenome_components,
        resource_components,
        traceability_components,
        resource_profiles,
    ) = _derive_profile_constants(selected_profile)
    validate_profile_arithmetic(
        evaluator_components=evaluator_components,
        evaluator_weights=evaluator_weights,
        pangenome_components=pangenome_components,
        resource_components=resource_components,
        traceability_components=traceability_components,
    )

    tuple_key = _validated_tuple(
        _required_mapping(payload, "tuple", context="payload"),
        context="tuple",
    )
    if tuple_key != trusted_tuple:
        raise ScoreInputError("payload tuple does not match trusted expected_tuple")

    profile = payload.get("score_profile")
    if profile != selected_profile["profile"]["id"]:
        raise ScoreInputError(f"unsupported score profile: {profile}")
    profile_sha256 = str(selected_profile["_sha256"])
    required_f1_metrics = _required_f1_records(payload, selected_profile)
    _validate_f1_component_binding(required_f1_metrics, payload, selected_profile)

    eligibility_status = payload.get("eligibility_status")
    if eligibility_status in NON_SCORABLE_STATUSES:
        return ScoreResult(
            tuple_key=trusted_tuple,
            score_profile=profile,
            score_profile_sha256=profile_sha256,
            evaluation_mode=str(evaluation_mode),
            score_status=str(eligibility_status),
            pgbench_score=None,
            pgbench_score_raw=None,
            evaluator_points=None,
            pangenome_points=None,
            resource_points=None,
            traceability_points=None,
            point_breakdown={},
            evaluator_scores={},
            required_f1_metrics=required_f1_metrics,
            reason=str(payload.get("reason") or eligibility_status),
        )
    if eligibility_status != "eligible":
        raise ScoreInputError(f"unknown eligibility_status: {eligibility_status}")

    infrastructure_valid = payload.get("infrastructure_valid")
    if not isinstance(infrastructure_valid, bool):
        raise ScoreInputError("infrastructure_valid must be boolean")

    evaluators = _required_mapping(payload, "evaluators", context="payload")
    pangenome = _required_mapping(payload, "pangenome", context="payload")
    resources = _required_mapping(payload, "resources", context="payload")
    traceability = _required_mapping(payload, "traceability", context="payload")

    evaluator_points, evaluator_breakdown, evaluator_scores = _evaluator_points(
        evaluators, evaluator_components
    )
    pangenome_points, pangenome_breakdown = _component_points(
        pangenome,
        pangenome_components,
        prefix="pangenome",
    )
    resource_points, resource_breakdown = _resource_points(
        resources,
        str(tuple_key["official_score_mode"]),
        resource_components,
        resource_profiles,
    )
    traceability_points, traceability_breakdown, core_valid = _traceability_points(
        traceability, traceability_components
    )

    breakdown = {
        **evaluator_breakdown,
        **pangenome_breakdown,
        **resource_breakdown,
        **traceability_breakdown,
    }

    if not infrastructure_valid or not core_valid:
        reason = (
            "benchmark infrastructure failure"
            if not infrastructure_valid
            else "core provenance invalid"
        )
        return ScoreResult(
            tuple_key=trusted_tuple,
            score_profile=profile,
            score_profile_sha256=profile_sha256,
            evaluation_mode=str(evaluation_mode),
            score_status="invalid",
            pgbench_score=None,
            pgbench_score_raw=None,
            evaluator_points=None,
            pangenome_points=None,
            resource_points=None,
            traceability_points=None,
            point_breakdown=breakdown,
            evaluator_scores=evaluator_scores,
            required_f1_metrics=required_f1_metrics,
            reason=reason,
        )

    raw_score = (
        evaluator_points + pangenome_points + resource_points + traceability_points
    )
    if raw_score < -1e-9 or raw_score > 100.0 + 1e-9:
        raise ScoreInputError(f"calculated score outside [0, 100]: {raw_score}")

    score_status = (
        "valid"
        if evaluation_mode == "formal" and math.isclose(traceability_points, 5.0)
        else "provisional"
    )
    return ScoreResult(
        tuple_key=trusted_tuple,
        score_profile=profile,
        score_profile_sha256=profile_sha256,
        evaluation_mode=str(evaluation_mode),
        score_status=score_status,
        pgbench_score=round(raw_score, 2),
        pgbench_score_raw=raw_score,
        evaluator_points=evaluator_points,
        pangenome_points=pangenome_points,
        resource_points=resource_points,
        traceability_points=traceability_points,
        point_breakdown=breakdown,
        evaluator_scores=evaluator_scores,
        required_f1_metrics=required_f1_metrics,
    )


def validate_profile_arithmetic(
    *,
    evaluator_components: Mapping[str, Mapping[str, float]] = EVALUATOR_COMPONENTS,
    evaluator_weights: Mapping[str, float] = EVALUATOR_WEIGHTS,
    pangenome_components: Mapping[str, float] = PANGENOME_COMPONENTS,
    resource_components: Mapping[str, float] = RESOURCE_COMPONENTS,
    traceability_components: Mapping[str, float] = TRACEABILITY_COMPONENTS,
) -> None:
    """Fail fast if a future edit changes the confirmed 100-point profile."""

    evaluator_maximums = {
        evaluator: sum(components.values())
        for evaluator, components in evaluator_components.items()
    }
    expected_maximums = {
        "truvari": 28.0,
        "aardvark": 24.5,
        "vcfdist": 17.5,
    }
    if evaluator_maximums != expected_maximums:
        raise RuntimeError(f"evaluator point profile changed: {evaluator_maximums!r}")
    if not math.isclose(sum(evaluator_weights.values()), 1.0):
        raise RuntimeError("evaluator weights must sum to one")
    for evaluator, maximum in evaluator_maximums.items():
        if not math.isclose(maximum / 70.0, evaluator_weights[evaluator]):
            raise RuntimeError(
                f"{evaluator} point maximum is inconsistent with its weight"
            )
    total = (
        sum(evaluator_maximums.values())
        + sum(pangenome_components.values())
        + sum(resource_components.values())
        + sum(traceability_components.values())
    )
    if not math.isclose(total, 100.0):
        raise RuntimeError(f"score profile must sum to 100, found {total}")


validate_profile_arithmetic()
