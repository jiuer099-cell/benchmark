from __future__ import annotations

import json
from pathlib import Path

import pytest

from workflow.scripts import pgbench_exec as MODULE


def _fastq(path: Path, mate: int) -> None:
    path.write_text(
        f"@read-1/{mate}\nACGT\n+\n!!!!\n"
        f"@read-2/{mate}\nGG\n+\n##\n",
        encoding="utf-8",
    )


def _manifest() -> dict[str, object]:
    return {
        "id": "opaque-adapter",
        "supported_modes": {
            "end_to_end_from_reads": {
                "billable_stages": ["calling"],
                "required_inputs": ["short_fastq_r1", "short_fastq_r2"],
                "optional_inputs": [],
            }
        },
    }


def test_core_fastq_validation_cache_is_shared_and_content_bound(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    r1, r2 = tmp_path / "reads.R1.fastq", tmp_path / "reads.R2.fastq"
    _fastq(r1, 1)
    _fastq(r2, 2)
    tool_manifest = tmp_path / "tool.yaml"
    tool_manifest.write_text("id: opaque-adapter\n", encoding="utf-8")
    resolved = tmp_path / "results" / "run-a" / "HG002" / "tool" / "meta" / "inputs.json"
    inputs = {"short_fastq_r1": r1, "short_fastq_r2": r2}

    MODULE.resolve_inputs(
        manifest=_manifest(), tool_manifest_path=tool_manifest,
        mode="end_to_end_from_reads", supplied_inputs=inputs,
        resolved_inputs_path=resolved, run_id="run-a",
    )
    first = json.loads(resolved.read_text(encoding="utf-8"))
    assert first["read_validation"]["cache_status"] == "validated"
    dataset_id = first["read_validation"]["dataset_id"]
    cache_root = tmp_path / "results" / "run-a" / "core" / "input-validation"
    entry = json.loads((cache_root / "entries" / f"{dataset_id}.json").read_text())
    assert entry["validator_version"] == MODULE.FASTQ_VALIDATOR_VERSION
    assert entry["validation"]["paired_fastq"] == {
        "read_pairs": 2, "read_count": 4, "read_bases": 12, "mate_names_match": True,
    }
    assert {"file_sha256", "file_size", "mtime_ns", "gzip_integrity"} <= set(
        entry["files"]["short_fastq_r1"]
    )

    def must_not_validate(_: object) -> dict[str, object]:
        raise AssertionError("cache hit must not decode FASTQ again")

    monkeypatch.setattr(MODULE, "validate_read_evidence", must_not_validate)
    MODULE.resolve_inputs(
        manifest=_manifest(), tool_manifest_path=tool_manifest,
        mode="end_to_end_from_reads", supplied_inputs=inputs,
        resolved_inputs_path=resolved, run_id="run-a",
    )
    second = json.loads(resolved.read_text(encoding="utf-8"))
    assert second["read_validation"]["cache_status"] == "reused"
    assert second["read_validation"]["dataset_id"] == dataset_id


def test_production_runs_layout_uses_the_shared_run_cache(tmp_path: Path) -> None:
    resolved = tmp_path / "runs" / "run-a" / "HG002" / "tool" / "meta" / "inputs.json"
    assert MODULE._input_validation_cache_root(resolved) == (
        tmp_path / "runs" / "run-a" / "core" / "input-validation"
    )
