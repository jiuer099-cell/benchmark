from __future__ import annotations

import math
import sys
from copy import deepcopy
from pathlib import Path

import pytest

SCRIPTS = Path(__file__).resolve().parents[2] / "workflow" / "scripts"
sys.path.insert(0, str(SCRIPTS))

from pgbench_scoring import ScoreInputError, calculate_pgbench_score  # noqa: E402
from score_tools import _apply_provenance_audit  # noqa: E402


def _payload() -> dict:
    return {
        "tuple": {
            "run_id": "synthetic",
            "sample": "HG002",
            "tool": "mock_genotyper",
            "official_score_mode": "end_to_end_from_reads",
            "primary_truth_profile": "giab_hg002_grch38_v5_0q",
        },
        "score_profile": "pgbench_consensus_v2",
        "eligibility_status": "eligible",
        "infrastructure_valid": True,
        "truth_eligible_count": 100,
        "consensus": {
            "all_three_correct": 60,
            "exactly_two_correct": 20,
            "exactly_one_correct": 10,
            "none_correct": 10,
        },
    }


def _calculate(payload: dict, evaluation_mode: str = "formal"):
    return calculate_pgbench_score(
        deepcopy(payload),
        expected_tuple=payload["tuple"],
        evaluation_mode=evaluation_mode,
    )


def test_equal_vote_consensus_formula_and_counts() -> None:
    result = _calculate(_payload())
    assert result.consensus_counts == {
        "all_three_correct": 60,
        "exactly_two_correct": 20,
        "exactly_one_correct": 10,
        "none_correct": 10,
    }
    assert result.total_evaluated == 100
    assert result.consensus_score == round(230 / 300 * 100, 2)
    assert result.comparable_score == round(2 * (230 / 3) / 200 * 100, 2)
    assert result.pgbench_score == result.consensus_score
    assert result.truth_eligible_count == 100
    assert result.comparable_precision == 230 / 3 / 100
    assert result.comparable_recall == 230 / 3 / 100
    assert result.unanimous_correct_rate == 0.6
    assert result.majority_correct_rate == 0.8
    assert result.evaluator_points is None
    assert result.resource_points is None


def test_panel_macro_f1_and_global_recovery_are_secondary() -> None:
    payload = _payload()
    global_score = 2 * (230 / 3) / 200 * 100
    payload["analysis"] = {
        "comparable_score": global_score,
        "pangenome_genotyping_score": 86.10423,
        "non_reference_f1_score": 84.4781,
        "panel_coverage": 0.0626,
    }

    result = _calculate(payload)

    assert result.pgbench_score == result.consensus_score
    assert result.pangenome_genotyping_score == 86.10
    assert result.non_reference_f1_score == 84.48
    assert result.panel_coverage == 0.0626
    assert result.global_end_to_end_sv_recovery_score == result.comparable_score


def test_variant_sites_contract_uses_universal_consensus_not_no_call_gt_zero() -> None:
    payload = _payload()
    global_score = 2 * (230 / 3) / 200 * 100
    payload["analysis"] = {
        "comparable_score": global_score,
        "pangenome_genotyping_score": 0.0,
        "candidate_genotype_summary": {
            "candidate_output_contract": "variant_sites",
        },
    }

    result = _calculate(payload)

    assert result.pgbench_score == result.consensus_score
    assert result.pangenome_genotyping_score is None


def test_perfect_unanimous_consensus_is_one_hundred() -> None:
    payload = _payload()
    payload["consensus"] = {
        "all_three_correct": 5,
        "exactly_two_correct": 0,
        "exactly_one_correct": 0,
        "none_correct": 0,
    }
    assert _calculate(payload).consensus_score == 100.0


def test_synthetic_result_is_provisional() -> None:
    assert _calculate(_payload(), "synthetic_smoke").score_status == "provisional"


def test_counts_must_be_nonnegative_integers_and_empty_callsets_score_zero() -> None:
    payload = _payload()
    payload["consensus"]["all_three_correct"] = 1.5
    with pytest.raises(ScoreInputError, match="non-negative integer"):
        _calculate(payload)
    payload = _payload()
    payload["consensus"] = {name: 0 for name in payload["consensus"]}
    result = _calculate(payload)
    assert result.comparable_score == 0.0
    assert result.consensus_score == 0.0
    assert result.total_evaluated == 0


def test_soft_true_positive_credit_cannot_exceed_truth_universe() -> None:
    payload = _payload()
    payload["truth_eligible_count"] = 10
    payload["consensus"] = {
        "all_three_correct": 11,
        "exactly_two_correct": 0,
        "exactly_one_correct": 0,
        "none_correct": 0,
    }
    with pytest.raises(ScoreInputError, match="exceeds.*truth universe"):
        _calculate(payload)


def test_no_weights_or_ranking_can_be_injected() -> None:
    payload = _payload()
    payload["weights"] = {"truvari": 0.9}
    with pytest.raises(ScoreInputError, match="unexpected weights"):
        _calculate(payload)
    payload = _payload()
    payload["rank"] = 1
    with pytest.raises(ScoreInputError, match="forbidden ranking"):
        _calculate(payload)


def test_provenance_is_gate_not_points() -> None:
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
    assert result.consensus_score == _calculate(_payload()).consensus_score
    assert result.traceability_points is None


def test_core_provenance_failure_suppresses_numeric_score() -> None:
    payload = _payload()
    payload["traceability"] = {
        "manifest_completeness": 1.0,
        "hash_lineage_complete": True,
        "environment_complete": True,
        "run_context_complete": True,
        "core_provenance_valid": False,
    }
    result = _calculate(payload)
    assert result.score_status == "invalid"
    assert result.consensus_score is None
