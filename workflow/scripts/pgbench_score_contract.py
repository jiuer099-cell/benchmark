"""The single canonical ME-F1 v1.0 score contract.

Scoring and finalization must agree on exactly one definition of the frozen
contract.  Keeping a private copy in each module is how a divergent
``FROZEN_SCORE_CONTRACT`` silently redefined what "the frozen contract" meant:
the finalizer's copy was missing ``score_semantics_version``, so the digest it
expected could never equal the digest scoring published.

This module is therefore the only place where the contract is written down.
``config/me_f1_scoring.yaml`` remains the canonical *data* input and is
validated against this definition; every consumer derives the contract
identity from :data:`ME_F1_SCORE_CONTRACT_SHA256`.

The contract content is ME-F1 v1.0 and must not be edited to make a failing
run pass.  Any change here is a new benchmark contract, not a bug fix.
"""

from __future__ import annotations

import hashlib
import json
from collections.abc import Mapping
from typing import Any

# Canonical contract serialization: sorted keys, no insignificant whitespace.
# This is the exact form whose digest is published as
# ``score_contract_sha256`` inside every sealed score.
_CONTRACT_JSON_OPTIONS = {
    "sort_keys": True,
    "separators": (",", ":"),
    "ensure_ascii": False,
}

ME_F1_SCORE_CONTRACT: dict[str, Any] = {
    "id": "ME-F1",
    "version": "1.0",
    "score_semantics_version": "1.0",
    "evaluators": ["truvari", "aardvark_gt", "vcfdist"],
    "formula": "arithmetic_mean",
    "weights": {
        "truvari": 0.3333333333333333,
        "aardvark_gt": 0.3333333333333333,
        "vcfdist": 0.3333333333333333,
    },
    "require_all_evaluators": True,
    "renormalize_missing_weights": False,
    "primary_score": {
        "expression": "(truvari_f1 + aardvark_gt_f1 + vcfdist_f1) / 3"
    },
    # Frozen fairness rule: no stratification may decide whether the primary
    # score is published.  Stratifications explain score differences only.
    "stratification_affects_primary_score": False,
}


def canonical_contract_sha256(value: Mapping[str, Any] | None = None) -> str:
    """Return the canonical digest of a score contract block."""

    contract = ME_F1_SCORE_CONTRACT if value is None else value
    canonical = json.dumps(contract, **_CONTRACT_JSON_OPTIONS).encode("utf-8")
    return hashlib.sha256(canonical).hexdigest()


ME_F1_SCORE_CONTRACT_SHA256 = canonical_contract_sha256()


def contract_differences(
    candidate: Mapping[str, Any] | None,
    expected: Mapping[str, Any] | None = None,
) -> list[str]:
    """Describe how a candidate contract block differs from the canonical one.

    Returns an empty list when the candidate is identical.  Used to fail
    closed with an actionable message instead of a bare hash mismatch.
    """

    reference = ME_F1_SCORE_CONTRACT if expected is None else expected
    if not isinstance(candidate, Mapping):
        return ["<candidate is not a mapping>"]
    differences: list[str] = []
    for key in sorted(set(reference) | set(candidate)):
        if key not in candidate:
            differences.append(f"missing:{key}")
        elif key not in reference:
            differences.append(f"unexpected:{key}")
        elif candidate[key] != reference[key]:
            differences.append(f"changed:{key}")
    return differences
