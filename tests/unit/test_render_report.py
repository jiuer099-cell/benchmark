from __future__ import annotations

import csv
import json
import sys
from pathlib import Path
from typing import Any

import pytest

SCRIPTS = Path(__file__).resolve().parents[2] / "workflow" / "scripts"
sys.path.insert(0, str(SCRIPTS))

import render_report as report  # noqa: E402

TUPLE_KEY = {
    "run_id": "synthetic-run",
    "sample": "HG002",
    "tool": "pangenie",
    "official_score_mode": "end_to_end_from_reads",
    "primary_truth_profile": "giab_hg002_grch38_v5_0q",
}

CONFUSION = {
    "heterozygous": {"heterozygous": 386, "hom_alt": 10, "hom_ref": 277, "no_call": 253},
    "hom_alt": {"heterozygous": 5, "hom_alt": 258, "hom_ref": 76, "no_call": 70},
}


def _score(analysis_overrides: dict[str, Any] | None = None) -> dict[str, Any]:
    analysis: dict[str, Any] = {
        "candidate_genotype_summary": {
            "all_site_call_rate": 0.758,
            "exact_gt_accuracy": 0.482,
            "no_call_rate": 0.242,
            "genotype_confusion_matrix": CONFUSION,
        },
        "addressability_audit": {
            "canonical_candidate_count": 18164,
            "called_count": 12231,
            "explicit_no_call_count": 5933,
            "addressability_rate": 1.0,
        },
    }
    if analysis_overrides:
        analysis.update(analysis_overrides)
    return {
        "tuple_key": dict(TUPLE_KEY),
        "score_profile": "pgbench_me_f1_v1",
        "score_profile_sha256": "a" * 64,
        "score_contract_version": "1.0",
        "score_contract_sha256": "b" * 64,
        "evaluation_mode": "formal",
        "score_status": "valid",
        "benchmark_score": 57.39,
        "required_f1_metrics": {
            "truvari": {"precision": 0.876, "recall": 0.482, "f1": 0.622},
            "aardvark": {"precision": 0.709, "recall": 0.390, "f1": 0.503},
            "vcfdist": {"precision": 0.839, "recall": 0.462, "f1": 0.596},
        },
        "formal_analysis": analysis,
    }


def _single_row(path: Path) -> dict[str, str]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(rows) == 1
    return rows[0]


def test_writes_candidate_diagnostics_and_addressability_no_call_rate(
    tmp_path: Path,
) -> None:
    out = tmp_path / "score.tsv"
    report.write_score_tsv(out, _score())
    row = _single_row(out)
    # The companion metrics requested for the report layer.
    assert float(row["AllSiteCallRate"]) == pytest.approx(0.758)
    assert float(row["ExactGTAccuracy"]) == pytest.approx(0.482)
    # NoCallRate carries the all-site genotype completeness semantics from
    # candidate_genotype_summary; the addressability denominator gets its own
    # differently named column instead of silently overwriting it.
    assert float(row["NoCallRate"]) == pytest.approx(0.242)
    assert float(row["AddressabilityNoCallRate"]) == pytest.approx(5933 / 18164)
    assert json.loads(row["GTConfusionMatrix"]) == CONFUSION


def test_degrades_safely_without_candidate_genotype_summary(tmp_path: Path) -> None:
    out = tmp_path / "score.tsv"
    report.write_score_tsv(
        out, _score(analysis_overrides={"candidate_genotype_summary": None})
    )
    row = _single_row(out)
    assert row["AllSiteCallRate"] == ""
    assert row["ExactGTAccuracy"] == ""
    assert row["NoCallRate"] == ""
    assert row["GTConfusionMatrix"] == ""
    # The addressability diagnostic is independent of the candidate summary.
    assert float(row["AddressabilityNoCallRate"]) == pytest.approx(5933 / 18164)


def test_non_mapping_candidate_summary_degrades_to_empty(tmp_path: Path) -> None:
    out = tmp_path / "score.tsv"
    report.write_score_tsv(
        out, _score(analysis_overrides={"candidate_genotype_summary": [1, 2, 3]})
    )
    row = _single_row(out)
    assert row["AllSiteCallRate"] == ""
    assert row["NoCallRate"] == ""


def test_tsv_header_matches_row_keys_and_has_no_duplicate_semantic_field(
    tmp_path: Path,
) -> None:
    out = tmp_path / "score.tsv"
    report.write_score_tsv(out, _score())
    header = out.read_text(encoding="utf-8").splitlines()[0].split("\t")
    assert len(header) == len(set(header)), "duplicate column in score TSV header"
    assert header.count("NoCallRate") == 1
    assert "AddressabilityNoCallRate" in header
    row = _single_row(out)
    # csv.DictReader stores values that have no header column under the None
    # key, so its absence proves every row key is present in the header.
    assert None not in row


def test_gt_confusion_matrix_serialization_is_stable(tmp_path: Path) -> None:
    first = tmp_path / "first.tsv"
    second = tmp_path / "second.tsv"
    report.write_score_tsv(first, _score())
    report.write_score_tsv(second, _score())
    assert first.read_bytes() == second.read_bytes()
    matrix = _single_row(first)["GTConfusionMatrix"]
    assert json.loads(matrix) == CONFUSION
    assert " " not in matrix, "confusion matrix must be compact JSON"
