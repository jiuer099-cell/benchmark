"""The frozen ME-F1 contract must exist in exactly one place.

A private copy is how the finalizer's expected digest silently drifted from the
digest scoring publishes.  These tests fail closed if a second copy reappears,
or if the diagnostic-vs-validity boundary is re-crossed.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest
import yaml

ROOT = Path(__file__).resolve().parents[2]
SCRIPTS = ROOT / "workflow" / "scripts"
sys.path.insert(0, str(SCRIPTS))

import finalize_score_provenance as finalizer  # noqa: E402
import pgbench_scoring  # noqa: E402
from pgbench_score_contract import (  # noqa: E402
    ME_F1_SCORE_CONTRACT,
    ME_F1_SCORE_CONTRACT_SHA256,
    canonical_contract_sha256,
    contract_differences,
)

PROFILE = yaml.safe_load(
    (ROOT / "config" / "me_f1_scoring.yaml").read_text(encoding="utf-8")
)["score"]


def test_profile_is_the_canonical_contract() -> None:
    assert contract_differences(PROFILE) == []
    assert canonical_contract_sha256(PROFILE) == ME_F1_SCORE_CONTRACT_SHA256


def test_canonical_contract_keeps_the_frozen_v1_semantics() -> None:
    assert ME_F1_SCORE_CONTRACT["id"] == "ME-F1"
    assert ME_F1_SCORE_CONTRACT["version"] == "1.0"
    assert ME_F1_SCORE_CONTRACT["score_semantics_version"] == "1.0"
    assert ME_F1_SCORE_CONTRACT["stratification_affects_primary_score"] is False
    assert ME_F1_SCORE_CONTRACT["renormalize_missing_weights"] is False
    assert ME_F1_SCORE_CONTRACT["require_all_evaluators"] is True
    assert ME_F1_SCORE_CONTRACT["weights"] == {
        "truvari": 0.3333333333333333,
        "aardvark_gt": 0.3333333333333333,
        "vcfdist": 0.3333333333333333,
    }


def test_scoring_publishes_the_canonical_contract_digest() -> None:
    profile = pgbench_scoring.load_score_profile()
    assert profile["_score_contract_sha256"] == ME_F1_SCORE_CONTRACT_SHA256


def test_finalizer_consumes_the_canonical_contract_digest() -> None:
    assert finalizer.ME_F1_SCORE_CONTRACT_SHA256 is ME_F1_SCORE_CONTRACT_SHA256


def test_finalizer_no_longer_keeps_a_private_contract_copy() -> None:
    assert not hasattr(finalizer, "FROZEN_SCORE_CONTRACT")
    assert not hasattr(pgbench_scoring, "_contract_sha256")


def test_a_mutated_profile_is_rejected_by_scoring(tmp_path: Path) -> None:
    profile = yaml.safe_load(
        (ROOT / "config" / "me_f1_scoring.yaml").read_text(encoding="utf-8")
    )
    del profile["score"]["score_semantics_version"]
    mutated = tmp_path / "me_f1_scoring.yaml"
    mutated.write_text(
        yaml.safe_dump(profile, sort_keys=False, allow_unicode=True),
        encoding="utf-8",
    )
    with pytest.raises(
        pgbench_scoring.ScoreInputError, match="missing:score_semantics_version"
    ):
        pgbench_scoring.load_score_profile(mutated)


def test_population_af_is_a_diagnostic_not_a_validity_gate() -> None:
    assert "population_af_stratification_complete" not in (
        pgbench_scoring.QUALITY_GATE_FIELDS
    )
    assert pgbench_scoring.QUALITY_GATE_FIELDS == {
        "evaluator_semantic_contract_valid",
        "addressability_complete",
        "evaluation_query_complete",
        "allowed_information_complete",
        "tuning_frozen",
        "panel_provenance_complete",
        "context_stratification_complete",
        "all_evaluator_evidence_complete",
        "family_aware_loo_complete",
        "statistical_uncertainty_complete",
    }
