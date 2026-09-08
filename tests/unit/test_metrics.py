from __future__ import annotations

import json
import sys
from copy import deepcopy
from pathlib import Path

import pytest
import yaml

ROOT = Path(__file__).resolve().parents[2]
SCRIPTS = ROOT / "workflow" / "scripts"
sys.path.insert(0, str(SCRIPTS))

from materialize_synthetic_metrics import materialize  # noqa: E402
from pgbench_metrics import MetricContractError, build_score_payload_from_metrics  # noqa: E402
from pgbench_scoring import load_score_profile  # noqa: E402

PROFILE_PATH = ROOT / "config" / "me_f1_scoring.yaml"
EVALUATOR_PROFILE_PATH = ROOT / "config" / "evaluator_profile.yaml"
SCHEMA_PATH = ROOT / "workflow" / "schemas" / "metrics.schema.yaml"
DICTIONARY_PATH = ROOT / "config" / "metric_dictionary.yaml"
FIXTURE_PATH = ROOT / "tests" / "fixtures" / "synthetic" / "evaluators" / "score_input.json"


def _fixture() -> dict:
    return json.loads(FIXTURE_PATH.read_text(encoding="utf-8"))


def _expected(fixture: dict) -> dict[str, str]:
    return dict(fixture["tuple"])


def _document() -> dict:
    dictionary = yaml.safe_load(DICTIONARY_PATH.read_text(encoding="utf-8"))
    dictionary["_source_path"] = str(DICTIONARY_PATH)
    return materialize(
        fixture=_fixture(),
        metric_dictionary=dictionary,
        score_profile_path=PROFILE_PATH,
        evaluator_profile_path=EVALUATOR_PROFILE_PATH,
        evaluator_manifest_id="a" * 64,
        pangenome_manifest_id="b" * 64,
        resource_manifest_id="c" * 64,
        metrics_schema_path=SCHEMA_PATH,
    )


def _build(document: dict) -> dict:
    return build_score_payload_from_metrics(
        document,
        expected_tuple=_expected(_fixture()),
        score_profile=load_score_profile(PROFILE_PATH),
        metrics_schema_path=SCHEMA_PATH,
        metric_dictionary_path=DICTIONARY_PATH,
    )


def test_me_f1_payload_uses_three_genotype_aware_evaluators() -> None:
    payload = _build(_document())
    assert set(payload["evaluator_f1"]) == {"truvari", "aardvark", "vcfdist"}
    assert payload["truth_eligible_count"] == 100
    assert all(payload["quality_gates"].values())


def test_metrics_tuple_is_bound_to_trusted_run_tuple() -> None:
    document = _document()
    document["tuple"]["tool_id"] = "forged"
    with pytest.raises(MetricContractError, match="does not match expected"):
        _build(document)


def test_me_f1_requires_materialized_evaluator_metrics() -> None:
    document = _document()
    document.pop("analysis")
    with pytest.raises(MetricContractError, match="requires formal analysis"):
        _build(document)


def test_evaluator_denominator_is_not_taken_from_per_tool_subset() -> None:
    document = _document()
    metric = document["analysis"]["evaluator_gt_metrics"]["vcfdist"]
    metric["fn"] = 20
    metric["recall"] = metric["tp"] / (metric["tp"] + metric["fn"])
    metric["f1"] = 2 * metric["precision"] * metric["recall"] / (
        metric["precision"] + metric["recall"]
    )
    with pytest.raises(Exception, match="does not use the frozen truth denominator"):
        from pgbench_scoring import calculate_me_f1

        calculate_me_f1(
            _build(document),
            expected_tuple=_expected(_fixture()),
            evaluation_mode="synthetic_smoke",
            score_profile_path=PROFILE_PATH,
        )


def test_traceability_metric_cannot_be_self_reported() -> None:
    document = _document()
    record = deepcopy(document["records"][0])
    record["metric_id"] = "traceability.hash_lineage.complete"
    record["value_type"] = "boolean"
    record["value"] = True
    record["evaluator"] = "provenance_audit"
    document["records"].append(record)
    with pytest.raises(MetricContractError, match="must not self-report"):
        _build(document)
