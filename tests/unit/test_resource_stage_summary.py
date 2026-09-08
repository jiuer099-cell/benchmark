from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

SCRIPTS = Path(__file__).resolve().parents[2] / "workflow" / "scripts"
sys.path.insert(0, str(SCRIPTS))

from summarize_resource_stages import ResourceStageError, summarize  # noqa: E402


def _measurement(path: Path, wall: float, cpu: float, ram: float) -> Path:
    path.write_text(
        json.dumps({"s": wall, "cpu_time": cpu, "max_rss_mb": ram}) + "\n",
        encoding="utf-8",
    )
    return path


def test_resource_costs_split_build_and_per_sample(tmp_path: Path) -> None:
    index = tmp_path / "index"
    index.mkdir()
    (index / "asset.bin").write_bytes(b"x" * 10)
    payload = summarize(
        [
            ("index", "one_time_build", _measurement(tmp_path / "build.jsonl", 100, 180, 4000)),
            ("genotype", "per_sample", _measurement(tmp_path / "sample.jsonl", 20, 30, 1000)),
        ],
        index_artifacts=[index],
    )
    assert payload["one_time_build"]["wall_seconds"] == 100
    assert payload["one_time_build"]["index_disk_bytes"] == 10
    assert payload["per_sample"]["wall_seconds"] == 20
    n10 = next(row for row in payload["amortized"] if row["cohort_size"] == 10)
    assert n10["wall_seconds_total"] == 300
    assert n10["wall_seconds_per_sample"] == 30


def test_resource_costs_reject_combined_only_measurement(tmp_path: Path) -> None:
    build = _measurement(tmp_path / "build.jsonl", 1, 1, 1)
    with pytest.raises(ResourceStageError, match="both one_time_build and per_sample"):
        summarize([("index", "one_time_build", build)])
