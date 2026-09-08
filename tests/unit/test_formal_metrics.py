from __future__ import annotations

import sys
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[2]
SCRIPTS = ROOT / "workflow" / "scripts"
sys.path.insert(0, str(SCRIPTS))

from materialize_formal_metrics import (  # noqa: E402
    EVALUATORS,
    evaluator_genotype_f1,
)


def _row(*, detection: bool, scorable: bool, correct: bool) -> dict[str, object]:
    return {
        "detection_correct": detection,
        "genotype_scorable": scorable,
        "genotype_correct": correct,
    }


def test_three_evaluators_use_identical_gt_aware_tp_fp_fn_definition() -> None:
    event_ids = {"exact", "wrong_gt", "unmatched"}
    ledgers = {
        "truvari": {
            "exact": _row(detection=True, scorable=True, correct=True),
            "wrong_gt": _row(detection=True, scorable=True, correct=False),
            "unmatched": _row(detection=False, scorable=False, correct=False),
        },
        "aardvark": {
            "exact": _row(detection=True, scorable=True, correct=True),
            "wrong_gt": _row(detection=True, scorable=True, correct=True),
            "unmatched": _row(detection=False, scorable=False, correct=False),
        },
        "vcfdist": {
            "exact": _row(detection=True, scorable=True, correct=False),
            "wrong_gt": _row(detection=True, scorable=True, correct=False),
            "unmatched": _row(detection=False, scorable=False, correct=False),
        },
    }

    metrics = evaluator_genotype_f1(ledgers, event_ids, panel_truth_positive=4)

    assert tuple(metrics) == EVALUATORS
    assert metrics["truvari"] == {
        "tp": 1,
        "fp": 2,
        "fn": 3,
        "precision": pytest.approx(1 / 3),
        "recall": pytest.approx(1 / 4),
        "f1": pytest.approx(2 / 7),
    }
    assert metrics["aardvark"]["tp"] == 2
    assert metrics["aardvark"]["fp"] == 1
    assert metrics["aardvark"]["fn"] == 2
    assert metrics["aardvark"]["f1"] == pytest.approx(4 / 7)
    assert metrics["vcfdist"]["tp"] == 0
    assert metrics["vcfdist"]["fp"] == 3
    assert metrics["vcfdist"]["fn"] == 4
    assert metrics["vcfdist"]["f1"] == 0.0


def test_gt_true_positive_cannot_exceed_panel_truth_denominator() -> None:
    ledgers = {
        evaluator: {
            "A": _row(detection=True, scorable=True, correct=True),
            "B": _row(detection=True, scorable=True, correct=True),
        }
        for evaluator in EVALUATORS
    }

    with pytest.raises(ValueError, match="exceed panel-addressable truth"):
        evaluator_genotype_f1(ledgers, {"A", "B"}, panel_truth_positive=1)


def test_empty_panel_truth_denominator_is_rejected() -> None:
    ledgers = {evaluator: {} for evaluator in EVALUATORS}

    with pytest.raises(ValueError, match="no positive candidate"):
        evaluator_genotype_f1(ledgers, set(), panel_truth_positive=0)
