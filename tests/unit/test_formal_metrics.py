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
    genotype_stratified_summary,
)
from sv_matching import SvRecord  # noqa: E402


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


def test_context_and_af_strata_publish_three_identical_gt_contracts(
    tmp_path: Path,
) -> None:
    context_beds: dict[str, Path] = {}
    for name in (
        "non_repeat",
        "tandem_repeat",
        "segmental_duplication",
        "low_complexity",
        "low_mappability",
        "other_difficult",
    ):
        path = tmp_path / f"{name}.bed"
        path.write_text("chr1\t0\t500\n", encoding="utf-8")
        context_beds[name] = path
    records = {
        "PGSV1": SvRecord("PGSV1", "chr1", 50, 109, "DEL", -60, "AA", "A", "0/1"),
        "PGSV2": SvRecord("PGSV2", "chr1", 150, 150, "INS", 80, "A", "AC", "0/1"),
    }
    hidden = {
        "PGSV1": {"truth_label": "positive", "truth_event_id": "T1"},
        "PGSV2": {"truth_label": "negative", "truth_event_id": "T2"},
    }
    addressability = {
        "PGSV1": {
            "tool_native_status": "addressable",
            "output_status": "called",
            "failure_reason": "",
        },
        "PGSV2": {
            "tool_native_status": "addressable",
            "output_status": "called",
            "failure_reason": "",
        },
    }
    ledgers = {
        evaluator: {
            "PGSV1": _row(detection=True, scorable=True, correct=True),
            "PGSV2": _row(detection=False, scorable=False, correct=False),
        }
        for evaluator in EVALUATORS
    }

    summary = genotype_stratified_summary(
        query_records=records,
        hidden_truth=hidden,
        ledgers=ledgers,
        event_ids={"PGSV1", "PGSV2"},
        candidate_af={"PGSV1": 0.005, "PGSV2": 0.2},
        context_beds=context_beds,
        addressability=addressability,
    )

    rare = summary["population_af"]["rare"]
    assert rare["truth_positive_count"] == 1
    assert rare["me_f1"] == pytest.approx(100.0)
    assert set(rare["evaluator_metrics"]) == set(EVALUATORS)
    assert rare["evaluator_metrics"]["truvari"]["f1"] == 1.0
    assert summary["population_af"]["common"]["status"] == "excluded_zero_truth"
    assert summary["genome_context"]["low_mappability"]["addressability_rate"] == 1.0
