from __future__ import annotations

import sys
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[2]
SCRIPTS = ROOT / "workflow" / "scripts"
sys.path.insert(0, str(SCRIPTS))

from materialize_formal_metrics import (  # noqa: E402
    EVALUATORS,
    candidate_genotype_summary,
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

    metrics = evaluator_genotype_f1(
        ledgers, event_ids, canonical_panel_truth_positive=4
    )

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

    with pytest.raises(ValueError, match="exceed canonical-panel truth"):
        evaluator_genotype_f1(
            ledgers, {"A", "B"}, canonical_panel_truth_positive=1
        )


def test_empty_panel_truth_denominator_is_rejected() -> None:
    ledgers = {evaluator: {} for evaluator in EVALUATORS}

    with pytest.raises(ValueError, match="no positive candidate"):
        evaluator_genotype_f1(ledgers, set(), canonical_panel_truth_positive=0)


def test_all_site_gt_diagnostics_use_the_frozen_candidate_denominator(
    tmp_path: Path,
) -> None:
    ledger = tmp_path / "hidden.tsv"
    ledger.write_text(
        "candidate_id\tpangenome_allele_id\ttruth_gt\ttruth_label\ttruth_scorable\ttruth_event_id\n"
        "C1\tA1\t0/1\tpositive\t1\tT1\n"
        "C2\tA2\t1/1\tpositive\t1\tT2\n"
        "C3\tA3\t0/0\tnegative\t1\tT3\n",
        encoding="utf-8",
    )
    query = tmp_path / "query.vcf"
    query.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG002\n"
        "chr1\t100\tC1\tA\tT\t.\tPASS\t.\tGT\t0/1\n"
        "chr1\t200\tC2\tA\tT\t.\tPASS\t.\tGT\t./.\n"
        "chr1\t300\tC3\tA\tT\t.\tPASS\t.\tGT\t0/1\n",
        encoding="utf-8",
    )
    manifest = tmp_path / "tool.yaml"
    manifest.write_text(
        "outputs:\n"
        "  candidate_output_contract: all_sites\n"
        "  absence_semantics: no_call\n",
        encoding="utf-8",
    )

    summary = candidate_genotype_summary(
        query_vcf=query,
        hidden_truth_ledger=ledger,
        tool_manifest=manifest,
        require_phase=False,
    )

    assert summary["all_site_call_rate"] == pytest.approx(2 / 3)
    assert summary["exact_gt_accuracy"] == pytest.approx(1 / 3)
    assert summary["no_call_rate"] == pytest.approx(1 / 3)
    assert summary["genotype_confusion_matrix"]["heterozygous"]["heterozygous"] == 1
    assert summary["genotype_confusion_matrix"]["hom_alt"]["no_call"] == 1
    # These diagnostics are intentionally separate from event-level ME-F1.
    assert "me_f1" not in summary


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
    assert rare["metric_id"] == "ME-F1_AF_RARE"
    assert rare["role"] == "explanatory"
    assert rare["affects_primary_score"] is False
    assert rare["me_f1"] == pytest.approx(100.0)
    assert set(rare["evaluator_metrics"]) == set(EVALUATORS)
    assert rare["evaluator_metrics"]["truvari"]["f1"] == 1.0
    assert summary["population_af"]["common"]["status"] == "excluded_zero_truth"
    assert summary["genome_context"]["low_mappability"]["addressability_rate"] == 1.0
