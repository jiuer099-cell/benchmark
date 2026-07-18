from __future__ import annotations

import hashlib
import math
import sys
from copy import deepcopy
from pathlib import Path

import pytest
import yaml

SCRIPTS = Path(__file__).resolve().parents[2] / "workflow" / "scripts"
sys.path.insert(0, str(SCRIPTS))

from pgbench_scoring import (  # noqa: E402
    EVALUATOR_COMPONENTS,
    EVALUATOR_WEIGHTS,
    PANGENOME_COMPONENTS,
    RESOURCE_COMPONENTS,
    SCORE_PROFILE,
    ScoreInputError,
    calculate_pgbench_score,
    resource_fraction,
)
from score_tools import _apply_provenance_audit  # noqa: E402

METRIC_DICTIONARY = yaml.safe_load(
    (
        Path(__file__).resolve().parents[2] / "config" / "metric_dictionary.yaml"
    ).read_text(encoding="utf-8")
)


def _sync_required_f1(payload: dict) -> None:
    metric_values: dict[str, float] = {}
    required_f1 = set(SCORE_PROFILE["required_f1_metrics"])
    evaluator_profiles = SCORE_PROFILE["layers"]["evaluator_accuracy"]["evaluators"]
    for evaluator, evaluator_profile in evaluator_profiles.items():
        for component_name, component in evaluator_profile["components"].items():
            if component["metric"] in required_f1:
                metric_values[component["metric"]] = payload["evaluators"][evaluator][
                    component_name
                ]
    for component_name, component in SCORE_PROFILE["layers"]["pangenome_robustness"][
        "components"
    ].items():
        if component["metric"] in required_f1:
            metric_values[component["metric"]] = payload["pangenome"][component_name]

    records = {}
    for metric_id in SCORE_PROFILE["required_f1_metrics"]:
        value = metric_values[metric_id]
        records[metric_id] = {
            "metric_id": metric_id,
            "value_type": "ratio",
            "value": value,
            "status": "defined",
            "evaluator": METRIC_DICTIONARY["metrics"][metric_id]["source"],
            "numerator": value,
            "denominator": 1.0,
            "eligible_count": 1,
            "universe_id": "synthetic_test_universe",
            "undefined_reason": None,
            "parser_id": "synthetic_test_parser",
            "parser_source_field": "synthetic.value",
            "provenance_manifest_id": "a" * 64,
            "strata": {},
            "confidence_interval": None,
        }
    payload["required_f1_metrics"] = records


def _perfect_payload() -> dict:
    payload = {
        "tuple": {
            "run_id": "synthetic",
            "sample": "HG002",
            "tool": "mock_genotyper",
            "official_score_mode": "end_to_end_from_reads",
            "primary_truth_profile": "giab_hg002_grch37_v5_0q",
        },
        "score_profile": "pgbench_v1",
        "eligibility_status": "eligible",
        "infrastructure_valid": True,
        "evaluators": {
            evaluator: {metric: 1.0 for metric in components}
            for evaluator, components in EVALUATOR_COMPONENTS.items()
        },
        "pangenome": {metric: 1.0 for metric in PANGENOME_COMPONENTS},
        "resources": {
            "wall_time_hours": 12.0,
            "peak_rss_gib": 64.0,
            "disk_peak_gib": 250.0,
            "cpu_hours": 96.0,
        },
        "traceability": {
            "manifest_completeness": 1.0,
            "hash_lineage_complete": True,
            "environment_complete": True,
            "run_context_complete": True,
            "core_provenance_valid": True,
        },
    }
    _sync_required_f1(payload)
    return payload


def _calculate(
    payload: dict,
    *,
    evaluation_mode: str = "formal",
    expected_tuple: dict[str, str] | None = None,
    **kwargs,
):
    candidate = deepcopy(payload)
    _sync_required_f1(candidate)
    return calculate_pgbench_score(
        candidate,
        expected_tuple=expected_tuple or dict(candidate["tuple"]),
        evaluation_mode=evaluation_mode,
        **kwargs,
    )


def test_confirmed_profile_arithmetic() -> None:
    evaluator_total = sum(
        sum(components.values()) for components in EVALUATOR_COMPONENTS.values()
    )
    assert evaluator_total == 70.0
    assert sum(PANGENOME_COMPONENTS.values()) == 15.0
    assert sum(RESOURCE_COMPONENTS.values()) == 10.0
    assert evaluator_total + 15.0 + 10.0 + 5.0 == 100.0
    assert EVALUATOR_WEIGHTS == {
        "truvari": 0.40,
        "aardvark": 0.35,
        "vcfdist": 0.25,
    }


def test_perfect_result_is_one_hundred_without_rank() -> None:
    result = _calculate(_perfect_payload())
    output = result.to_dict()
    assert output["pgbench_score"] == 100.0
    assert output["score_status"] == "valid"
    assert output["evaluator_points"] == 70.0
    assert output["pangenome_points"] == 15.0
    assert output["resource_points"] == 10.0
    assert output["traceability_points"] == 5.0
    score_profile = (
        Path(__file__).resolve().parents[2] / "config" / "score_weights.yaml"
    )
    expected_hash = hashlib.sha256(score_profile.read_bytes()).hexdigest()
    assert output["score_profile_sha256"] == expected_hash
    assert "rank" not in output
    assert "leaderboard" not in output
    assert set(output["required_f1_metrics"]) == set(
        SCORE_PROFILE["required_f1_metrics"]
    )


def test_evaluator_normalization_matches_fixed_weights() -> None:
    payload = _perfect_payload()
    payload["evaluators"]["truvari"] = {
        metric: 0.5 for metric in EVALUATOR_COMPONENTS["truvari"]
    }
    result = _calculate(payload)
    assert result.evaluator_scores["truvari"] == 50.0
    assert result.evaluator_scores["aardvark"] == 100.0
    assert result.evaluator_scores["vcfdist"] == 100.0
    expected = 0.70 * (0.40 * 50.0 + 0.35 * 100.0 + 0.25 * 100.0)
    assert math.isclose(result.evaluator_points or -1.0, expected)


def test_resource_fraction_is_log_interpolated() -> None:
    assert resource_fraction(12.0, 12.0, 72.0) == 1.0
    assert resource_fraction(72.0, 12.0, 72.0) == 0.0
    midpoint = math.sqrt(12.0 * 72.0)
    assert math.isclose(resource_fraction(midpoint, 12.0, 72.0), 0.5)


def test_incomplete_noncore_traceability_is_provisional() -> None:
    payload = _perfect_payload()
    payload["traceability"]["manifest_completeness"] = 0.75
    result = _calculate(payload)
    assert result.score_status == "provisional"
    assert result.traceability_points == 4.5
    assert result.pgbench_score == 99.5


def test_synthetic_evaluators_can_never_emit_valid_score() -> None:
    payload = _perfect_payload()
    result = _calculate(payload, evaluation_mode="synthetic_smoke")
    assert result.pgbench_score == 100.0
    assert result.score_status == "provisional"
    assert result.evaluation_mode == "synthetic_smoke"


def test_core_provenance_failure_has_no_numeric_score() -> None:
    payload = _perfect_payload()
    payload["traceability"]["core_provenance_valid"] = False
    result = _calculate(payload)
    assert result.score_status == "invalid"
    assert result.pgbench_score is None
    assert result.pgbench_score_raw is None


def test_non_scorable_tool_status_has_no_score() -> None:
    payload = _perfect_payload()
    payload["eligibility_status"] = "not_applicable"
    result = _calculate(payload)
    assert result.score_status == "not_applicable"
    assert result.pgbench_score is None


def test_out_of_range_metric_is_rejected() -> None:
    payload = _perfect_payload()
    payload["evaluators"]["truvari"]["overall_event_f1"] = 1.01
    with pytest.raises(ScoreInputError, match="outside"):
        _calculate(payload)


def test_missing_metric_is_rejected() -> None:
    payload = _perfect_payload()
    del payload["pangenome"]["allele_link_precision"]
    with pytest.raises(ScoreInputError, match="invalid fields"):
        _calculate(payload)


def test_changed_profile_arithmetic_is_rejected(tmp_path: Path) -> None:
    source = Path(__file__).resolve().parents[2] / "config" / "score_weights.yaml"
    profile = yaml.safe_load(source.read_text(encoding="utf-8"))
    profile["layers"]["evaluator_accuracy"]["evaluators"]["truvari"]["components"][
        "overall_event_f1"
    ]["max_points"] = 11.0
    changed = tmp_path / "changed.yaml"
    changed.write_text(yaml.safe_dump(profile), encoding="utf-8")
    with pytest.raises(RuntimeError, match="evaluator point profile changed"):
        _calculate(_perfect_payload(), score_profile_path=changed)


def test_measured_provenance_audit_overrides_fixture_traceability() -> None:
    payload = _perfect_payload()
    measured = {
        "audit_schema_version": "pgbench.provenance_audit.v1",
        "status": "provisional",
        "manifest_completeness": 1.0,
        "hash_lineage_complete": True,
        "environment_complete": False,
        "run_context_complete": True,
        "core_provenance_valid": True,
    }
    merged = _apply_provenance_audit(payload, measured)
    result = _calculate(merged)
    assert result.traceability_points == 4.0
    assert result.pgbench_score == 99.0
    assert payload["traceability"]["environment_complete"] is True


def test_measured_provenance_audit_requires_all_score_gates() -> None:
    with pytest.raises(ScoreInputError, match="missing"):
        _apply_provenance_audit(
            _perfect_payload(),
            {"manifest_completeness": 1.0},
        )


def test_payload_cannot_self_declare_evaluation_mode() -> None:
    payload = _perfect_payload()
    payload["evaluation_mode"] = "formal"
    with pytest.raises(ScoreInputError, match="unexpected evaluation_mode"):
        _calculate(payload, evaluation_mode="synthetic_smoke")


def test_resource_profile_override_is_forbidden() -> None:
    payload = _perfect_payload()
    payload["resource_profile"] = {
        metric: {"target": 1e9, "limit": 1e10} for metric in RESOURCE_COMPONENTS
    }
    with pytest.raises(ScoreInputError, match="overrides are forbidden"):
        _calculate(payload)


def test_nested_ranking_field_is_rejected() -> None:
    payload = _perfect_payload()
    payload["tuple"]["rank"] = 1
    with pytest.raises(ScoreInputError, match="forbidden ranking field"):
        _calculate(payload)


def test_nan_metric_is_rejected() -> None:
    payload = _perfect_payload()
    payload["evaluators"]["truvari"]["overall_event_f1"] = float("nan")
    with pytest.raises(ScoreInputError, match="finite"):
        _calculate(payload)


def test_payload_tuple_must_match_trusted_expected_tuple() -> None:
    payload = _perfect_payload()
    expected = dict(payload["tuple"])
    expected["tool"] = "different_tool"
    with pytest.raises(ScoreInputError, match="trusted expected_tuple"):
        _calculate(payload, expected_tuple=expected)


def test_required_f1_raw_value_must_match_scored_component() -> None:
    payload = _perfect_payload()
    metric_id = "truvari.event.overall.f1"
    payload["required_f1_metrics"][metric_id]["value"] = 0.25
    with pytest.raises(ScoreInputError, match="does not match"):
        calculate_pgbench_score(
            payload,
            expected_tuple=payload["tuple"],
            evaluation_mode="formal",
        )
