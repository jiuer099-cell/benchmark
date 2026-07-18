from __future__ import annotations

import hashlib
import json
import sys
from copy import deepcopy
from pathlib import Path

import pytest
import yaml

ROOT = Path(__file__).resolve().parents[2]
SCRIPTS = ROOT / "workflow" / "scripts"
sys.path.insert(0, str(SCRIPTS))

from pgbench_metrics import (  # noqa: E402
    MetricContractError,
    build_score_payload_from_metrics,
)
from pgbench_scoring import load_score_profile  # noqa: E402
from score_tools import main as score_tools_main  # noqa: E402

SCORE_PROFILE_PATH = ROOT / "config" / "score_weights.yaml"
METRIC_DICTIONARY_PATH = ROOT / "config" / "metric_dictionary.yaml"
METRICS_SCHEMA_PATH = ROOT / "workflow" / "schemas" / "metrics.schema.yaml"
SCORE_PROFILE = load_score_profile(SCORE_PROFILE_PATH)
METRIC_DICTIONARY = yaml.safe_load(METRIC_DICTIONARY_PATH.read_text(encoding="utf-8"))


def _expected_tuple() -> dict[str, str]:
    return {
        "run_id": "formal_run",
        "sample": "HG002",
        "tool": "trusted_tool",
        "official_score_mode": "end_to_end_from_reads",
        "primary_truth_profile": "giab_hg002_grch37_v5_0q",
    }


def _required_metric_ids() -> set[str]:
    layers = SCORE_PROFILE["layers"]
    result = {
        component["metric"]
        for evaluator in layers["evaluator_accuracy"]["evaluators"].values()
        for component in evaluator["components"].values()
    }
    for layer_name in ("pangenome_robustness", "resource_efficiency"):
        result.update(
            component["metric"]
            for component in layers[layer_name]["components"].values()
        )
    result.update(SCORE_PROFILE["required_f1_metrics"])
    return result


def _record(metric_id: str) -> dict:
    definition = METRIC_DICTIONARY["metrics"][metric_id]
    value_type = definition["value_type"]
    value: float | bool = True if value_type == "boolean" else 1.0
    return {
        "metric_id": metric_id,
        "value_type": value_type,
        "value": value,
        "status": "defined",
        "evaluator": definition["source"],
        "numerator": 1.0,
        "denominator": 1.0,
        "eligible_count": 1,
        "universe_id": "formal_test_universe",
        "undefined_reason": None,
        "parser_id": "formal_test_parser_v1",
        "parser_source_field": "summary.value",
        "provenance_manifest_id": "b" * 64,
        "strata": {},
        "confidence_interval": None,
    }


def _metrics_document() -> dict:
    expected = _expected_tuple()
    return {
        "schema_version": 1,
        "dictionary_id": "pgbench_metrics_v1",
        "tuple": {
            "run_id": expected["run_id"],
            "sample_id": expected["sample"],
            "tool_id": expected["tool"],
            "official_score_mode": expected["official_score_mode"],
            "primary_truth_profile": expected["primary_truth_profile"],
            "score_profile": "pgbench_v1",
        },
        "records": [_record(metric_id) for metric_id in sorted(_required_metric_ids())],
    }


def _build(document: dict) -> dict:
    return build_score_payload_from_metrics(
        document,
        expected_tuple=_expected_tuple(),
        score_profile=SCORE_PROFILE,
        metrics_schema_path=METRICS_SCHEMA_PATH,
        metric_dictionary_path=METRIC_DICTIONARY_PATH,
    )


def _write_json(path: Path, payload: dict) -> None:
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


def _cli_args(tmp_path: Path, *, metrics_flag: str = "--metrics") -> list[str]:
    metrics = tmp_path / "metrics.json"
    audit = tmp_path / "audit.json"
    context = tmp_path / "run-context.json"
    output = tmp_path / "score.json"
    _write_json(metrics, _metrics_document())
    _write_json(
        audit,
        {
            "audit_schema_version": "pgbench.provenance_audit.v1",
            "status": "valid",
            "manifest_completeness": 1.0,
            "hash_lineage_complete": True,
            "environment_complete": True,
            "run_context_complete": True,
            "core_provenance_valid": True,
        },
    )
    _write_json(
        context,
        {
            "schema_version": 1,
            "score_profile_sha256": hashlib.sha256(
                SCORE_PROFILE_PATH.read_bytes()
            ).hexdigest(),
        },
    )
    expected = _expected_tuple()
    return [
        metrics_flag,
        str(metrics),
        "--output",
        str(output),
        "--score-profile",
        str(SCORE_PROFILE_PATH),
        "--metric-dictionary",
        str(METRIC_DICTIONARY_PATH),
        "--metrics-schema",
        str(METRICS_SCHEMA_PATH),
        "--provenance-audit",
        str(audit),
        "--run-context",
        str(context),
        "--evaluation-mode",
        "formal",
        "--expected-run-id",
        expected["run_id"],
        "--expected-sample-id",
        expected["sample"],
        "--expected-tool-id",
        expected["tool"],
        "--expected-official-score-mode",
        expected["official_score_mode"],
        "--expected-primary-truth-profile",
        expected["primary_truth_profile"],
    ]


def test_standard_metrics_are_mapped_by_profile_metric_ids() -> None:
    payload = _build(_metrics_document())
    assert payload["evaluators"]["truvari"]["overall_event_f1"] == 1.0
    assert payload["pangenome"]["out_of_panel_truth_f1"] == 1.0
    assert payload["resources"]["wall_time_hours"] == 1.0
    assert set(payload["required_f1_metrics"]) == set(
        SCORE_PROFILE["required_f1_metrics"]
    )
    assert "traceability" not in payload


def test_missing_required_aggregate_metric_is_rejected() -> None:
    document = _metrics_document()
    document["records"] = document["records"][1:]
    with pytest.raises(MetricContractError, match="missing required aggregate"):
        _build(document)


def test_duplicate_empty_strata_aggregate_is_rejected() -> None:
    document = _metrics_document()
    document["records"].append(deepcopy(document["records"][0]))
    with pytest.raises(MetricContractError, match="more than one"):
        _build(document)


def test_metric_evaluator_must_match_dictionary_source() -> None:
    document = _metrics_document()
    document["records"][0]["evaluator"] = "fusion"
    with pytest.raises(MetricContractError, match="dictionary source"):
        _build(document)


def test_metric_value_type_must_match_dictionary() -> None:
    document = _metrics_document()
    document["records"][0]["value_type"] = "count"
    with pytest.raises(MetricContractError, match="value_type"):
        _build(document)


def test_metric_manifest_id_must_be_sha256() -> None:
    document = _metrics_document()
    document["records"][0]["provenance_manifest_id"] = "manifest-1"
    with pytest.raises(MetricContractError, match="64-character SHA-256"):
        _build(document)


def test_metric_counts_must_be_non_null() -> None:
    document = _metrics_document()
    document["records"][0]["eligible_count"] = None
    with pytest.raises(MetricContractError, match="non-null eligible_count"):
        _build(document)


def test_ratio_value_must_match_numerator_and_denominator() -> None:
    document = _metrics_document()
    ratio = next(
        record for record in document["records"] if record["value_type"] == "ratio"
    )
    ratio["value"] = 0.5
    ratio["numerator"] = 3.0
    ratio["denominator"] = 4.0
    with pytest.raises(MetricContractError, match="does not match"):
        _build(document)


def test_ratio_denominator_must_be_positive() -> None:
    document = _metrics_document()
    ratio = next(
        record for record in document["records"] if record["value_type"] == "ratio"
    )
    ratio["value"] = 0.0
    ratio["numerator"] = 0.0
    ratio["denominator"] = 0.0
    with pytest.raises(MetricContractError, match="denominator must be positive"):
        _build(document)


def test_ratio_value_must_be_in_unit_interval() -> None:
    document = _metrics_document()
    ratio = next(
        record for record in document["records"] if record["value_type"] == "ratio"
    )
    ratio["value"] = 1.1
    ratio["numerator"] = 1.1
    ratio["denominator"] = 1.0
    with pytest.raises(MetricContractError, match=r"maximum of 1\.0|must be in"):
        _build(document)


def test_nonfinite_metric_is_rejected_before_json_schema() -> None:
    document = _metrics_document()
    document["records"][0]["value"] = float("nan")
    with pytest.raises(MetricContractError, match="finite"):
        _build(document)


def test_nested_ranking_field_in_metrics_is_rejected() -> None:
    document = _metrics_document()
    document["records"][0]["strata"] = {"rank": "1"}
    with pytest.raises(MetricContractError, match="forbidden ranking field"):
        _build(document)


def test_traceability_metrics_cannot_be_self_reported() -> None:
    document = _metrics_document()
    document["records"].append(_record("traceability.hash_lineage.complete"))
    with pytest.raises(MetricContractError, match="must not self-report"):
        _build(document)


def test_metrics_tuple_must_match_trusted_expected_tuple() -> None:
    document = _metrics_document()
    document["tuple"]["tool_id"] = "forged_tool"
    with pytest.raises(MetricContractError, match="does not match expected"):
        _build(document)


def test_required_metric_status_must_be_scorable() -> None:
    document = _metrics_document()
    record = document["records"][0]
    record["status"] = "undefined"
    record["value"] = None
    record["undefined_reason"] = "fixture deliberately undefined"
    with pytest.raises(MetricContractError, match="non-scorable status"):
        _build(document)


def test_cli_scores_only_bound_standard_metrics_and_preserves_raw_f1(
    tmp_path: Path,
) -> None:
    args = _cli_args(tmp_path)
    assert score_tools_main(args) == 0
    output = json.loads((tmp_path / "score.json").read_text(encoding="utf-8"))
    assert output["score_status"] == "valid"
    assert output["tuple_key"] == _expected_tuple()
    assert set(output["required_f1_metrics"]) == set(
        SCORE_PROFILE["required_f1_metrics"]
    )
    assert all(
        record["status"] == "defined"
        for record in output["required_f1_metrics"].values()
    )


def test_cli_input_alias_accepts_only_standard_metrics(tmp_path: Path) -> None:
    args = _cli_args(tmp_path, metrics_flag="--input")
    assert score_tools_main(args) == 0

    metrics_path = tmp_path / "metrics.json"
    _write_json(metrics_path, {"evaluators": {"truvari": {"f1": 1.0}}})
    assert score_tools_main(args) == 2


def test_cli_rejects_score_profile_not_frozen_in_run_context(
    tmp_path: Path,
) -> None:
    args = _cli_args(tmp_path)
    context = tmp_path / "run-context.json"
    _write_json(context, {"schema_version": 1, "score_profile_sha256": "0" * 64})
    assert score_tools_main(args) == 2
