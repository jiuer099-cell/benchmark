from __future__ import annotations

import json
import sys
from copy import deepcopy
from pathlib import Path
from typing import Any

import jsonschema
import pytest
import yaml

ROOT = Path(__file__).resolve().parents[2]
SCRIPTS = ROOT / "workflow" / "scripts"
sys.path.insert(0, str(SCRIPTS))

from materialize_synthetic_metrics import (  # noqa: E402
    COMPONENT_PATHS,
    SyntheticMetricError,
    materialize,
)
from pgbench_metrics import build_score_payload_from_metrics  # noqa: E402
from pgbench_scoring import load_score_profile  # noqa: E402

FIXTURE_PATH = (
    ROOT / "tests" / "fixtures" / "synthetic" / "evaluators" / "score_input.json"
)
METRIC_DICTIONARY_PATH = ROOT / "config" / "metric_dictionary.yaml"
METRICS_SCHEMA_PATH = ROOT / "workflow" / "schemas" / "metrics.schema.yaml"
SCORE_PROFILE_PATH = ROOT / "config" / "score_weights.yaml"
EVALUATOR_MANIFEST_ID = "a" * 64
PANGENOME_MANIFEST_ID = "b" * 64
RESOURCE_MANIFEST_ID = "c" * 64


def _fixture() -> dict[str, Any]:
    return json.loads(FIXTURE_PATH.read_text(encoding="utf-8"))


def _metric_dictionary() -> dict[str, Any]:
    dictionary = yaml.safe_load(METRIC_DICTIONARY_PATH.read_text(encoding="utf-8"))
    dictionary["_source_path"] = str(METRIC_DICTIONARY_PATH)
    return dictionary


def _materialize(fixture: dict[str, Any]) -> dict[str, Any]:
    return materialize(
        fixture=fixture,
        metric_dictionary=_metric_dictionary(),
        score_profile_path=SCORE_PROFILE_PATH,
        evaluator_manifest_id=EVALUATOR_MANIFEST_ID,
        pangenome_manifest_id=PANGENOME_MANIFEST_ID,
        resource_manifest_id=RESOURCE_MANIFEST_ID,
        metrics_schema_path=METRICS_SCHEMA_PATH,
    )


def _expected_tuple(fixture: dict[str, Any]) -> dict[str, str]:
    return {
        "run_id": fixture["tuple"]["run_id"],
        "sample": fixture["tuple"]["sample"],
        "tool": fixture["tuple"]["tool"],
        "official_score_mode": fixture["tuple"]["official_score_mode"],
        "primary_truth_profile": fixture["tuple"]["primary_truth_profile"],
    }


def test_complete_fixture_materializes_to_valid_metrics_contract() -> None:
    fixture = _fixture()
    document = _materialize(fixture)

    schema = yaml.safe_load(METRICS_SCHEMA_PATH.read_text(encoding="utf-8"))
    jsonschema.Draft202012Validator(schema).validate(document)
    score_payload = build_score_payload_from_metrics(
        document,
        expected_tuple=_expected_tuple(fixture),
        score_profile=load_score_profile(SCORE_PROFILE_PATH),
        metrics_schema_path=METRICS_SCHEMA_PATH,
        metric_dictionary_path=METRIC_DICTIONARY_PATH,
    )

    records = {record["metric_id"]: record for record in document["records"]}
    dictionary_metrics = _metric_dictionary()["metrics"]
    assert set(records) == set(COMPONENT_PATHS)
    assert all(
        record["value_type"] == dictionary_metrics[metric_id]["value_type"]
        and record["evaluator"] == dictionary_metrics[metric_id]["source"]
        for metric_id, record in records.items()
    )
    assert document["tuple"] == {
        "run_id": "synthetic_smoke",
        "sample_id": "HG002",
        "tool_id": "example_genotyper",
        "official_score_mode": "end_to_end_from_reads",
        "primary_truth_profile": "giab_hg002_grch37_v5_0q",
        "score_profile": "pgbench_v1",
    }
    assert score_payload["evaluators"] == fixture["evaluators"]
    assert score_payload["pangenome"] == fixture["pangenome"]
    assert score_payload["resources"] == fixture["resources"]
    assert records["vcfdist.phase.accuracy"]["status"] == (
        "completed_with_tool_ineligible_zero"
    )


def test_non_synthetic_smoke_fixture_is_rejected() -> None:
    fixture = _fixture()
    fixture["evaluation_mode"] = "formal"

    with pytest.raises(SyntheticMetricError, match="only accepts evaluation_mode"):
        _materialize(fixture)


def test_ratio_outside_unit_interval_is_rejected() -> None:
    fixture = _fixture()
    fixture["evaluators"]["truvari"]["overall_event_f1"] = 1.01

    with pytest.raises(SyntheticMetricError, match=r"ratio outside \[0,1\]"):
        _materialize(fixture)


def test_missing_fixture_metric_is_rejected() -> None:
    fixture = deepcopy(_fixture())
    del fixture["pangenome"]["allele_link_coverage"]

    with pytest.raises(
        SyntheticMetricError,
        match=r"missing pangenome\.allele_link_coverage",
    ):
        _materialize(fixture)
