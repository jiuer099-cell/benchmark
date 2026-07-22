"""Validate standardized metric documents and derive trusted score inputs."""

from __future__ import annotations

import math
import re
from collections.abc import Mapping, Sequence
from copy import deepcopy
from pathlib import Path
from typing import Any

import jsonschema
import yaml  # type: ignore[import-untyped]


class MetricContractError(ValueError):
    """Raised when standardized metrics violate the scoring contract."""


PROJECT_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_METRICS_SCHEMA_PATH = (
    PROJECT_ROOT / "workflow" / "schemas" / "metrics.schema.yaml"
)
DEFAULT_METRIC_DICTIONARY_PATH = PROJECT_ROOT / "config" / "metric_dictionary.yaml"

SCORE_TUPLE_FIELDS = {
    "run_id",
    "sample",
    "tool",
    "official_score_mode",
    "primary_truth_profile",
}
SCORE_TO_METRICS_TUPLE = {
    "run_id": "run_id",
    "sample": "sample_id",
    "tool": "tool_id",
    "official_score_mode": "official_score_mode",
    "primary_truth_profile": "primary_truth_profile",
}
SCORABLE_STATUSES = {"defined", "completed_with_tool_ineligible_zero"}
_SHA256_RE = re.compile(r"^[0-9a-f]{64}$")


def _load_yaml_mapping(path: Path, label: str) -> dict[str, Any]:
    try:
        with path.open("r", encoding="utf-8") as handle:
            loaded = yaml.safe_load(handle)
    except (OSError, yaml.YAMLError) as exc:
        raise MetricContractError(f"cannot load {label} {path}: {exc}") from exc
    if not isinstance(loaded, dict):
        raise MetricContractError(f"{label} must be a YAML mapping: {path}")
    return loaded


def _json_path(parts: Sequence[Any]) -> str:
    return ".".join(str(part) for part in parts) or "<root>"


def _ensure_finite_numbers(value: Any, *, path: str = "metrics") -> None:
    if isinstance(value, bool) or value is None or isinstance(value, str):
        return
    if isinstance(value, (int, float)):
        if not math.isfinite(float(value)):
            raise MetricContractError(f"{path} must be finite")
        return
    if isinstance(value, Mapping):
        for key, item in value.items():
            _ensure_finite_numbers(item, path=f"{path}.{key}")
        return
    if isinstance(value, Sequence):
        for index, item in enumerate(value):
            _ensure_finite_numbers(item, path=f"{path}[{index}]")


def _reject_ranking_fields(value: Any, *, path: str = "metrics") -> None:
    if isinstance(value, Mapping):
        for key, item in value.items():
            if isinstance(key, str):
                normalized = key.casefold()
                if (
                    normalized
                    in {"rank", "ranking", "leaderboard", "leaderboard_position"}
                    or normalized.startswith("rank_")
                    or normalized.startswith("leaderboard_")
                ):
                    raise MetricContractError(f"forbidden ranking field: {path}.{key}")
            _reject_ranking_fields(item, path=f"{path}.{key}")
    elif isinstance(value, Sequence) and not isinstance(value, (str, bytes)):
        for index, item in enumerate(value):
            _reject_ranking_fields(item, path=f"{path}[{index}]")


def _validate_expected_tuple(expected_tuple: Mapping[str, str]) -> dict[str, str]:
    if set(expected_tuple) != SCORE_TUPLE_FIELDS:
        missing = sorted(SCORE_TUPLE_FIELDS - set(expected_tuple))
        extra = sorted(set(expected_tuple) - SCORE_TUPLE_FIELDS)
        details = []
        if missing:
            details.append("missing " + ", ".join(missing))
        if extra:
            details.append("unexpected " + ", ".join(extra))
        raise MetricContractError(
            "expected score tuple has invalid fields: " + "; ".join(details)
        )
    normalized = dict(expected_tuple)
    for field, value in normalized.items():
        if not isinstance(value, str) or not value:
            raise MetricContractError(
                f"expected score tuple field {field} must be a non-empty string"
            )
    return normalized


def _profile_component_metrics(
    score_profile: Mapping[str, Any],
) -> tuple[
    dict[str, dict[str, str]],
    dict[str, str],
    dict[str, str],
    tuple[str, ...],
]:
    try:
        layers = score_profile["layers"]
        evaluator_profiles = layers["evaluator_accuracy"]["evaluators"]
        pangenome_components = layers["pangenome_robustness"]["components"]
        resource_components = layers["resource_efficiency"]["components"]
        raw_required_f1 = score_profile["required_f1_metrics"]
    except (KeyError, TypeError) as exc:
        raise MetricContractError(
            "score profile is missing metric component mappings"
        ) from exc

    evaluator_metrics: dict[str, dict[str, str]] = {}
    for evaluator, evaluator_profile in evaluator_profiles.items():
        if not isinstance(evaluator, str) or not isinstance(evaluator_profile, Mapping):
            raise MetricContractError("invalid evaluator profile mapping")
        components = evaluator_profile.get("components")
        if not isinstance(components, Mapping):
            raise MetricContractError(f"evaluator {evaluator} has no component mapping")
        evaluator_metrics[evaluator] = _component_metric_map(
            components, f"evaluator_accuracy.{evaluator}"
        )

    if not isinstance(pangenome_components, Mapping) or not isinstance(
        resource_components, Mapping
    ):
        raise MetricContractError("score profile component mappings are invalid")
    if not isinstance(raw_required_f1, list) or not all(
        isinstance(metric_id, str) and metric_id for metric_id in raw_required_f1
    ):
        raise MetricContractError(
            "score profile required_f1_metrics must be a list of metric IDs"
        )
    if len(raw_required_f1) != len(set(raw_required_f1)):
        raise MetricContractError("score profile required_f1_metrics has duplicates")
    return (
        evaluator_metrics,
        _component_metric_map(pangenome_components, "pangenome_robustness"),
        _component_metric_map(resource_components, "resource_efficiency"),
        tuple(raw_required_f1),
    )


def _component_metric_map(
    components: Mapping[str, Any], context: str
) -> dict[str, str]:
    result: dict[str, str] = {}
    for component_name, component in components.items():
        if not isinstance(component_name, str) or not isinstance(component, Mapping):
            raise MetricContractError(f"invalid component in {context}")
        metric_id = component.get("metric")
        if not isinstance(metric_id, str) or not metric_id:
            raise MetricContractError(
                f"{context}.{component_name}.metric must be a metric ID"
            )
        result[component_name] = metric_id
    return result


def _validate_dictionary(
    dictionary: Mapping[str, Any], schema: Mapping[str, Any]
) -> Mapping[str, Mapping[str, Any]]:
    if dictionary.get("schema_version") != 1:
        raise MetricContractError("metric dictionary schema_version must be 1")
    if dictionary.get("dictionary_id") != "pgbench_metrics_v1":
        raise MetricContractError("unsupported metric dictionary ID")
    metrics = dictionary.get("metrics")
    if not isinstance(metrics, Mapping) or not all(
        isinstance(metric_id, str) and isinstance(entry, Mapping)
        for metric_id, entry in metrics.items()
    ):
        raise MetricContractError("metric dictionary metrics must be a mapping")
    try:
        schema_ids = set(schema["$defs"]["metricId"]["enum"])
    except (KeyError, TypeError) as exc:
        raise MetricContractError("metrics schema has no metricId enum") from exc
    if set(metrics) != schema_ids:
        raise MetricContractError(
            "metric dictionary IDs do not exactly match metrics schema IDs"
        )
    return metrics


def _ratio_epsilon(dictionary: Mapping[str, Any]) -> float:
    value_rules = dictionary.get("value_rules")
    if not isinstance(value_rules, Mapping):
        raise MetricContractError("metric dictionary value_rules must be a mapping")
    epsilon = value_rules.get("floating_clamp_epsilon")
    if (
        isinstance(epsilon, bool)
        or not isinstance(epsilon, (int, float))
        or not math.isfinite(float(epsilon))
        or float(epsilon) <= 0
    ):
        raise MetricContractError(
            "metric dictionary floating_clamp_epsilon must be positive and finite"
        )
    return float(epsilon)


def _validate_record_semantics(
    record: Mapping[str, Any],
    dictionary_entry: Mapping[str, Any],
    *,
    index: int,
    ratio_epsilon: float,
) -> None:
    metric_id = str(record["metric_id"])
    expected_type = dictionary_entry.get("value_type")
    if record.get("value_type") != expected_type:
        raise MetricContractError(
            f"records[{index}] {metric_id} value_type must be {expected_type!r}"
        )
    expected_evaluator = dictionary_entry.get("source")
    if record.get("evaluator") != expected_evaluator:
        raise MetricContractError(
            f"records[{index}] {metric_id} evaluator must match dictionary "
            f"source {expected_evaluator!r}"
        )
    manifest_id = record.get("provenance_manifest_id")
    if not isinstance(manifest_id, str) or not _SHA256_RE.fullmatch(manifest_id):
        raise MetricContractError(
            f"records[{index}] {metric_id} provenance_manifest_id must be "
            "a lowercase 64-character SHA-256"
        )
    for field in ("numerator", "denominator", "eligible_count"):
        if record.get(field) is None:
            raise MetricContractError(
                f"records[{index}] {metric_id} requires non-null {field}"
            )
    if expected_type == "ratio" and record.get("status") in SCORABLE_STATUSES:
        value = record.get("value")
        numerator = record.get("numerator")
        denominator = record.get("denominator")
        if (
            isinstance(value, bool)
            or not isinstance(value, (int, float))
            or isinstance(numerator, bool)
            or not isinstance(numerator, (int, float))
            or isinstance(denominator, bool)
            or not isinstance(denominator, (int, float))
        ):
            raise MetricContractError(
                f"records[{index}] {metric_id} ratio fields must be numeric"
            )
        numeric_value = float(value)
        numeric_numerator = float(numerator)
        numeric_denominator = float(denominator)
        if numeric_denominator <= 0:
            raise MetricContractError(
                f"records[{index}] {metric_id} ratio denominator must be positive"
            )
        if not -ratio_epsilon <= numeric_value <= 1.0 + ratio_epsilon:
            raise MetricContractError(
                f"records[{index}] {metric_id} ratio value must be in [0, 1]"
            )
        calculated = numeric_numerator / numeric_denominator
        if not math.isclose(
            numeric_value,
            calculated,
            rel_tol=ratio_epsilon,
            abs_tol=ratio_epsilon,
        ):
            raise MetricContractError(
                f"records[{index}] {metric_id} ratio value does not match "
                "numerator/denominator"
            )


def _aggregate_records(
    records: Sequence[Mapping[str, Any]],
) -> dict[str, Mapping[str, Any]]:
    aggregates: dict[str, Mapping[str, Any]] = {}
    for record in records:
        if record.get("strata") != {}:
            continue
        metric_id = str(record["metric_id"])
        if metric_id in aggregates:
            raise MetricContractError(
                f"metric {metric_id} has more than one empty-strata aggregate record"
            )
        aggregates[metric_id] = record
    return aggregates


def _scoreable_value(record: Mapping[str, Any], metric_id: str) -> Any:
    status = record.get("status")
    if status not in SCORABLE_STATUSES:
        raise MetricContractError(
            f"required scoring metric {metric_id} has non-scorable status {status!r}"
        )
    value = record.get("value")
    if isinstance(value, float) and not math.isfinite(value):
        raise MetricContractError(f"required scoring metric {metric_id} is not finite")
    if value is None:
        raise MetricContractError(f"required scoring metric {metric_id} has no value")
    return value


def build_score_payload_from_metrics(
    document: Mapping[str, Any],
    *,
    expected_tuple: Mapping[str, str],
    score_profile: Mapping[str, Any],
    metrics_schema_path: Path = DEFAULT_METRICS_SCHEMA_PATH,
    metric_dictionary_path: Path = DEFAULT_METRIC_DICTIONARY_PATH,
) -> dict[str, Any]:
    """Validate one metrics document and extract only profile-mapped score values."""

    if not isinstance(document, Mapping):
        raise MetricContractError("metrics document must be a JSON object")
    _reject_ranking_fields(document)
    _ensure_finite_numbers(document)
    schema = _load_yaml_mapping(metrics_schema_path, "metrics schema")
    dictionary = _load_yaml_mapping(metric_dictionary_path, "metric dictionary")
    try:
        jsonschema.Draft202012Validator.check_schema(schema)
    except jsonschema.SchemaError as exc:
        raise MetricContractError(f"invalid metrics schema: {exc}") from exc
    validation_errors = sorted(
        jsonschema.Draft202012Validator(schema).iter_errors(document),
        key=lambda error: [str(part) for part in error.absolute_path],
    )
    if validation_errors:
        details = "; ".join(
            f"{_json_path(error.absolute_path)}: {error.message}"
            for error in validation_errors
        )
        raise MetricContractError(f"invalid metrics document: {details}")

    dictionary_metrics = _validate_dictionary(dictionary, schema)
    ratio_epsilon = _ratio_epsilon(dictionary)
    if document.get("dictionary_id") != dictionary.get("dictionary_id"):
        raise MetricContractError(
            "metrics document dictionary_id does not match metric dictionary"
        )

    normalized_expected = _validate_expected_tuple(expected_tuple)
    metrics_tuple = document["tuple"]
    profile_meta = score_profile.get("profile")
    if not isinstance(profile_meta, Mapping):
        raise MetricContractError("score profile has no profile metadata")
    score_profile_id = profile_meta.get("id")
    if not isinstance(score_profile_id, str) or not score_profile_id:
        raise MetricContractError("score profile ID must be a non-empty string")
    for score_field, metrics_field in SCORE_TO_METRICS_TUPLE.items():
        if metrics_tuple[metrics_field] != normalized_expected[score_field]:
            raise MetricContractError(
                f"metrics tuple {metrics_field} does not match expected score tuple"
            )
    if metrics_tuple["score_profile"] != score_profile_id:
        raise MetricContractError(
            "metrics tuple score_profile does not match selected score profile"
        )

    raw_records = document["records"]
    records = [record for record in raw_records if isinstance(record, Mapping)]
    for index, record in enumerate(records):
        metric_id = str(record["metric_id"])
        _validate_record_semantics(
            record,
            dictionary_metrics[metric_id],
            index=index,
            ratio_epsilon=ratio_epsilon,
        )
        if dictionary_metrics[metric_id].get("source") == "provenance_audit":
            raise MetricContractError(
                f"metrics document must not self-report traceability metric {metric_id}"
            )

    aggregates = _aggregate_records(records)
    if score_profile_id == "pgbench_consensus_v2":
        consensus_metrics = {
            "all_three_correct": "consensus.all_three_correct.count",
            "exactly_two_correct": "consensus.exactly_two_correct.count",
            "exactly_one_correct": "consensus.exactly_one_correct.count",
            "none_correct": "consensus.none_correct.count",
        }
        missing = sorted(set(consensus_metrics.values()) - set(aggregates))
        if missing:
            raise MetricContractError(
                "metrics document is missing consensus counts: " + ", ".join(missing)
            )
        consensus: dict[str, int] = {}
        for category, metric_id in consensus_metrics.items():
            value = _scoreable_value(aggregates[metric_id], metric_id)
            if isinstance(value, bool) or not isinstance(value, int) or value < 0:
                raise MetricContractError(f"{metric_id} must be a non-negative integer")
            consensus[category] = value
        return {
            "tuple": normalized_expected,
            "score_profile": score_profile_id,
            "eligibility_status": "eligible",
            "infrastructure_valid": True,
            "consensus": consensus,
        }
    (
        evaluator_metric_map,
        pangenome_metric_map,
        resource_metric_map,
        required_f1_ids,
    ) = _profile_component_metrics(score_profile)
    point_metric_ids = (
        {
            metric_id
            for components in evaluator_metric_map.values()
            for metric_id in components.values()
        }
        | set(pangenome_metric_map.values())
        | set(resource_metric_map.values())
    )
    required_metric_ids = point_metric_ids | set(required_f1_ids)
    unknown_required = sorted(required_metric_ids - set(dictionary_metrics))
    if unknown_required:
        raise MetricContractError(
            "score profile references unknown metric IDs: "
            + ", ".join(unknown_required)
        )
    missing = sorted(required_metric_ids - set(aggregates))
    if missing:
        raise MetricContractError(
            "metrics document is missing required aggregate metrics: "
            + ", ".join(missing)
        )
    for metric_id in sorted(required_metric_ids):
        _scoreable_value(aggregates[metric_id], metric_id)

    evaluators = {
        evaluator: {
            component: _scoreable_value(aggregates[metric_id], metric_id)
            for component, metric_id in components.items()
        }
        for evaluator, components in evaluator_metric_map.items()
    }
    pangenome = {
        component: _scoreable_value(aggregates[metric_id], metric_id)
        for component, metric_id in pangenome_metric_map.items()
    }
    resources = {
        component: _scoreable_value(aggregates[metric_id], metric_id)
        for component, metric_id in resource_metric_map.items()
    }
    required_f1_metrics = {
        metric_id: deepcopy(dict(aggregates[metric_id]))
        for metric_id in required_f1_ids
    }
    return {
        "tuple": normalized_expected,
        "score_profile": score_profile_id,
        "eligibility_status": "eligible",
        "infrastructure_valid": True,
        "evaluators": evaluators,
        "pangenome": pangenome,
        "resources": resources,
        "required_f1_metrics": required_f1_metrics,
    }
