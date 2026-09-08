from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

SCRIPTS = Path(__file__).resolve().parents[2] / "workflow" / "scripts"
sys.path.insert(0, str(SCRIPTS))

from summarize_coverage_matrix import CoverageSummaryError, summarize  # noqa: E402


def _score(path: Path, tool: str, run_id: str, value: float) -> Path:
    path.write_text(
        json.dumps(
            {
                "tuple_key": {"tool": tool, "run_id": run_id},
                "benchmark_score": value,
            }
        ),
        encoding="utf-8",
    )
    return path


def test_repeated_coverage_summary_requires_same_three_seeds(tmp_path: Path) -> None:
    paths = []
    for coverage in (10, 20, 30):
        for seed, score in ((1701, 80.0), (1702, 82.0), (1703, 84.0)):
            run_id = f"base_{coverage}x_seed{seed}"
            paths.append(_score(tmp_path / f"{run_id}.json", "tool", run_id, score))
    paths.append(_score(tmp_path / "full.json", "tool", "base_full", 90.0))

    result = summarize(paths)

    assert result["downsampling_seeds"] == [1701, 1702, 1703]
    ten_x = next(row for row in result["rows"] if row["coverage"] == 10)
    assert ten_x["mean_me_f1"] == 82.0
    assert ten_x["metric_id"] == "ME-F1_10X"
    assert ten_x["role"] == "explanatory"
    assert ten_x["affects_primary_score"] is False
    assert ten_x["sd_me_f1"] == pytest.approx(2.0)
    assert ten_x["ci95_lower_me_f1"] == pytest.approx(77.0317245)
    assert ten_x["ci95_upper_me_f1"] == pytest.approx(86.9682755)
    assert ten_x["ci95_method"] == "student_t_across_frozen_seeds"
    full = next(row for row in result["rows"] if row["coverage"] == "full")
    assert full["sd_me_f1"] is None
    assert full["ci95_lower_me_f1"] is None


def test_coverage_summary_rejects_incomplete_matrix(tmp_path: Path) -> None:
    path = _score(tmp_path / "one.json", "tool", "base_10x_seed1701", 80.0)
    with pytest.raises(CoverageSummaryError, match="10x, 20x, and 30x"):
        summarize([path])
