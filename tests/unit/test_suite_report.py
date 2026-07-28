from __future__ import annotations

import sys
from pathlib import Path

import pytest

SCRIPTS = Path(__file__).resolve().parents[2] / "workflow" / "scripts"
sys.path.insert(0, str(SCRIPTS))

from render_suite_report import SuiteReportError, render_suite  # noqa: E402


def _score(tool: str, mode: str, value: float) -> dict:
    return {
        "tuple_key": {
            "run_id": f"run-{tool}",
            "sample": "HG002",
            "tool": tool,
            "official_score_mode": mode,
            "primary_truth_profile": "giab_hg002_grch38_v5_0q",
        },
        "score_profile_sha256": "a" * 64,
        "comparable_score": value,
        "consensus_score": value + 1,
        "truth_eligible_count": 100,
        "total_evaluated": 90,
        "score_status": "valid",
    }


def _manifest(tool: str, paradigm: str, technology: str) -> dict:
    return {
        "id": tool,
        "paradigm": paradigm,
        "capabilities": {"technology": [technology]},
    }


def test_cross_track_suite_accepts_same_sample_truth_and_profile() -> None:
    html = render_suite(
        [
            (
                _score("pangenie", "end_to_end_from_reads", 75.0),
                _manifest(
                    "pangenie",
                    "short_read_pangenome_genotyping",
                    "illumina_short_read",
                ),
            ),
            (
                _score("kanpig", "caller_only_shared_alignment", 80.0),
                _manifest("kanpig", "genotyping_only", "pacbio_clr"),
            ),
        ]
    )
    assert "PanGenie".casefold() in html.casefold()
    assert "kanpig" in html
    assert "ComparableScore" in html
    assert "illumina_short_read" in html
    assert "pacbio_clr" in html


def test_cross_track_suite_rejects_different_truth_profile() -> None:
    first = _score("a", "end_to_end_from_reads", 70.0)
    second = _score("b", "caller_only_shared_alignment", 80.0)
    second["tuple_key"]["primary_truth_profile"] = "different"
    with pytest.raises(SuiteReportError, match="share sample, truth profile"):
        render_suite(
            [
                (first, _manifest("a", "mapping_based", "illumina_short_read")),
                (second, _manifest("b", "mapping_based", "pacbio_clr")),
            ]
        )
