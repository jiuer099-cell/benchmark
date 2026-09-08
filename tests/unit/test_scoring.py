from __future__ import annotations

from copy import deepcopy
import sys
from pathlib import Path

import pytest

SCRIPTS = Path(__file__).resolve().parents[2] / "workflow" / "scripts"
sys.path.insert(0, str(SCRIPTS))

from pgbench_scoring import ScoreInputError, calculate_me_f1  # noqa: E402
from score_tools import _apply_provenance_audit  # noqa: E402


def _metric(tp: int, fp: int, fn: int) -> dict[str, float | int]:
    precision = tp / (tp + fp) if tp + fp else 0.0
    recall = tp / (tp + fn) if tp + fn else 0.0
    f1 = 2 * precision * recall / (precision + recall) if precision + recall else 0.0
    return {"tp": tp, "fp": fp, "fn": fn, "precision": precision, "recall": recall, "f1": f1}


def _payload() -> dict:
    return {
        "tuple": {
            "run_id": "synthetic",
            "sample": "HG002",
            "tool": "mock_genotyper",
            "official_score_mode": "end_to_end_from_reads",
            "primary_truth_profile": "giab_hg002_grch38_v5_0q",
        },
        "score_profile": "pgbench_me_f1_v1",
        "eligibility_status": "eligible",
        "infrastructure_valid": True,
        "truth_eligible_count": 100,
        "evaluator_f1": {
            "truvari": _metric(80, 10, 20),
            "aardvark": _metric(75, 15, 25),
            "vcfdist": _metric(70, 20, 30),
        },
        "quality_gates": {
            "evaluator_semantic_contract_valid": True,
            "addressability_complete": True,
            "allowed_information_complete": True,
            "tuning_frozen": True,
        },
        "analysis": {"candidate_genotype_summary": {"candidate_output_contract": "all_sites"}},
    }


def _calculate(payload: dict, evaluation_mode: str = "formal"):
    return calculate_me_f1(
        deepcopy(payload), expected_tuple=payload["tuple"], evaluation_mode=evaluation_mode
    )


def test_me_f1_is_unweighted_mean_of_three_recomputed_f1_values() -> None:
    result = _calculate(_payload())
    expected = sum(float(v["f1"]) for v in _payload()["evaluator_f1"].values()) / 3 * 100
    assert result.me_f1 == round(expected, 2)
    assert result.benchmark_score == result.me_f1
    assert result.evaluator_scores == {
        name: round(float(value["f1"]) * 100, 2)
        for name, value in _payload()["evaluator_f1"].items()
    }


def test_f1_must_exactly_match_tp_fp_fn() -> None:
    payload = _payload()
    payload["evaluator_f1"]["aardvark"]["f1"] = 0.99
    with pytest.raises(ScoreInputError, match="does not match TP/FP/FN"):
        _calculate(payload)


def test_all_three_evaluators_use_identical_truth_denominator() -> None:
    payload = _payload()
    payload["evaluator_f1"]["vcfdist"] = _metric(70, 20, 29)
    with pytest.raises(ScoreInputError, match="frozen truth denominator"):
        _calculate(payload)


def test_failed_quality_gate_suppresses_primary_score() -> None:
    payload = _payload()
    payload["quality_gates"]["evaluator_semantic_contract_valid"] = False
    result = _calculate(payload)
    assert result.score_status == "invalid"
    assert result.me_f1 is None
    assert "evaluator_semantic_contract_valid" in result.reason


def test_zero_callset_is_valid_zero_f1_not_missing_data() -> None:
    payload = _payload()
    payload["evaluator_f1"] = {name: _metric(0, 0, 100) for name in payload["evaluator_f1"]}
    result = _calculate(payload)
    assert result.score_status == "valid"
    assert result.me_f1 == 0.0


def test_synthetic_result_is_provisional() -> None:
    assert _calculate(_payload(), "synthetic_smoke").score_status == "provisional"


def test_provenance_is_a_gate_not_points() -> None:
    audit = {
        "audit_schema_version": "pgbench.provenance_audit.v1",
        "status": "provisional",
        "manifest_completeness": 1.0,
        "hash_lineage_complete": True,
        "environment_complete": False,
        "run_context_complete": True,
        "core_provenance_valid": True,
    }
    result = _calculate(_apply_provenance_audit(_payload(), audit))
    assert result.score_status == "provisional"
    assert result.me_f1 is not None


def test_non_all_sites_output_is_rejected() -> None:
    payload = _payload()
    payload["analysis"]["candidate_genotype_summary"]["candidate_output_contract"] = "unsupported_contract"
    with pytest.raises(ScoreInputError, match="all-sites"):
        _calculate(payload)
