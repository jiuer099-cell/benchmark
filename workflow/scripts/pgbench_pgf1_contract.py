"""The single immutable PG-F1 formal scoring contract."""
from __future__ import annotations
import hashlib, json
from collections.abc import Mapping
from typing import Any

_JSON = {"sort_keys": True, "separators": (",", ":"), "ensure_ascii": False}
PG_F1_SCORE_CONTRACT: dict[str, Any] = {
    "id": "PG-F1", "version": "1.0", "score_semantics_version": "1.0",
    "canonical_scoring_unit": "canonical_candidate_x_alt",
    "canonical_unit_id": "candidate_id:A<one_based_alt_index>",
    "evaluators": ["truvari", "aardvark", "vcfdist"],
    "evaluator_evidence": {"relations": ["NONE", "LM", "AM"], "binary_mapping": {"NONE": 0, "LM": 0, "AM": 1}, "require_all_evaluators": True, "no_aggregate_to_unit_inference": True, "mapping_contract": "one_to_one_only"},
    "vote": {"kind": "binary_2_of_3", "threshold": 2, "dynamic_denominator": False},
    "genotype": {"phase_sensitive": False, "allele_order_sensitive": False, "canonicalization": "alt_specific_unphased_dosage", "no_call": "NA"},
    "accounting": {"truth_alt_query_alt_match_exact_gt": "TP", "truth_alt_query_alt_gt_mismatch": "FN+FP", "truth_alt_query_ref_or_no_call": "FN", "truth_ref_query_alt": "FP", "truth_ref_query_ref": "IGNORE"},
    "failure_policy": {"missing_evaluator_evidence": "CONSENSUS_EVIDENCE_INCOMPLETE", "complex_mapping": "AMBIGUOUS_COMPLEX_MAPPING", "normalizer_failure": "NORMALIZER_RECONCILIATION_FAILURE", "failure_is_not_score_zero": True},
    "stratification_affects_primary_score": False,
}
def canonical_contract_sha256(value: Mapping[str, Any] | None = None) -> str:
    return hashlib.sha256(json.dumps(PG_F1_SCORE_CONTRACT if value is None else value, **_JSON).encode("utf-8")).hexdigest()
PG_F1_SCORE_CONTRACT_SHA256 = canonical_contract_sha256()
def contract_differences(candidate: Mapping[str, Any] | None) -> list[str]:
    return [key for key in sorted(set(PG_F1_SCORE_CONTRACT) | set(candidate or {})) if (candidate or {}).get(key) != PG_F1_SCORE_CONTRACT.get(key)]
