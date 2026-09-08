from __future__ import annotations

import csv
import hashlib
import json
import sys
from pathlib import Path

import pytest

SCRIPTS = Path(__file__).resolve().parents[2] / "workflow" / "scripts"
sys.path.insert(0, str(SCRIPTS))

from aggregate_run_summary import (  # noqa: E402
    SummaryError,
    aggregate,
    main,
    paired_tool_differences,
)


def _write_json(path: Path, value: object) -> None:
    path.write_text(json.dumps(value), encoding="utf-8")


def _score(tool: str, value: float) -> dict[str, object]:
    return {
        "tuple_key": {
            "run_id": "run",
            "sample": "HG002",
            "tool": tool,
            "official_score_mode": "end_to_end_from_reads",
            "primary_truth_profile": "truth",
        },
        "score_profile": "pgbench_v1",
        "score_profile_sha256": "a" * 64,
        "evaluation_mode": "synthetic_smoke",
        "score_status": "provisional",
        "benchmark_score": value,
        "evaluator_scores": {"truvari": value, "aardvark": value, "vcfdist": value},
        "evaluator_range": 0.0,
        "evaluator_sd": 0.0,
        "point_breakdown": {"evaluator.truvari.overall_event_f1": value / 10},
    }


def _metrics(tool: str) -> dict[str, object]:
    return {
        "schema_version": 1,
        "dictionary_id": "pgbench_metrics_v1",
        "tuple": {
            "run_id": "run",
            "sample_id": "HG002",
            "tool_id": tool,
            "official_score_mode": "end_to_end_from_reads",
            "primary_truth_profile": "truth",
            "score_profile": "pgbench_v1",
        },
        "records": [
            {
                "metric_id": "truvari.event.overall.f1",
                "value_type": "ratio",
                "value": 0.5,
                "status": "defined",
                "evaluator": "truvari",
                "numerator": 1,
                "denominator": 2,
                "eligible_count": 2,
                "universe_id": "u",
                "undefined_reason": None,
                "parser_id": "test",
                "parser_source_field": "f1",
                "provenance_manifest_id": "b" * 64,
                "strata": {},
                "confidence_interval": None,
            }
        ],
    }


def _package(score_path: Path, metrics_path: Path) -> Path:
    score = json.loads(score_path.read_text(encoding="utf-8"))
    payload = {
        "package_schema_version": "pgbench.final_score_package.v1",
        "sealed": True,
        "score": score,
        "score_artifact": {
            "path": str(score_path),
            "sha256": hashlib.sha256(score_path.read_bytes()).hexdigest(),
        },
        "metrics_artifact": {
            "path": str(metrics_path),
            "sha256": hashlib.sha256(metrics_path.read_bytes()).hexdigest(),
        },
        "gates": {
            "core_provenance_valid": True,
            "hash_lineage_complete": True,
            "manifest_completeness": True,
        },
    }
    path = score_path.with_suffix(".package.json")
    _write_json(path, payload)
    return path


def _finalizer_manifest(
    score_path: Path,
    metrics_path: Path,
    package_path: Path,
) -> Path:
    score = json.loads(score_path.read_text(encoding="utf-8"))
    tuple_key = score["tuple_key"]
    payload = {
        "rule_name": "finalize_score_provenance",
        "status": "success",
        "run_id": tuple_key["run_id"],
        "job_key": ".".join(
            (
                tuple_key["sample"],
                tuple_key["tool"],
                tuple_key["official_score_mode"],
            )
        ),
        "wildcards": {"tool": tuple_key["tool"]},
        "input_paths": [str(score_path), str(metrics_path)],
        "input_sha256": {
            str(score_path): hashlib.sha256(score_path.read_bytes()).hexdigest(),
            str(metrics_path): hashlib.sha256(metrics_path.read_bytes()).hexdigest(),
        },
        "output_paths": [str(package_path)],
        "output_sha256": {
            str(package_path): hashlib.sha256(package_path.read_bytes()).hexdigest()
        },
    }
    path = package_path.with_suffix(".manifest.json")
    _write_json(path, payload)
    return path


def test_aggregate_preserves_input_order_and_never_ranks(tmp_path: Path) -> None:
    score_paths: list[Path] = []
    metrics_paths: list[Path] = []
    for tool, value in (("zeta", 10.0), ("alpha", 99.0)):
        score_path = tmp_path / f"{tool}.score.json"
        metrics_path = tmp_path / f"{tool}.metrics.json"
        _write_json(score_path, _score(tool, value))
        _write_json(metrics_path, _metrics(tool))
        score_paths.append(score_path)
        metrics_paths.append(metrics_path)
    score_tsv = tmp_path / "score.tsv"
    point_tsv = tmp_path / "points.tsv"
    metrics_tsv = tmp_path / "metrics.tsv"
    metrics_json = tmp_path / "metrics.json"

    aggregate(
        score_paths=score_paths,
        metrics_paths=metrics_paths,
        output_score_tsv=score_tsv,
        output_point_tsv=point_tsv,
        output_metrics_tsv=metrics_tsv,
        output_metrics_json=metrics_json,
    )

    with score_tsv.open(encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert [row["tool"] for row in rows] == ["zeta", "alpha"]
    assert "rank" not in score_tsv.read_text(encoding="utf-8").lower()
    assert len(point_tsv.read_text(encoding="utf-8").splitlines()) == 3
    assert len(metrics_tsv.read_text(encoding="utf-8").splitlines()) == 3
    assert len(json.loads(metrics_json.read_text())["documents"]) == 2


def test_paired_tool_difference_uses_identical_block_draws() -> None:
    documents = []
    for tool, point, samples in (
        ("tool_a", 80.0, [70.0, 80.0, 90.0]),
        ("tool_b", 70.0, [65.0, 70.0, 75.0]),
    ):
        document = _metrics(tool)
        document["analysis"] = {
            "contract_version": "formal_genotype_me_f1_v1",
            "comparison_track_sha256": "c" * 64,
            "me_f1": point,
            "me_f1_confidence_interval": {
                "replicates": 3,
                "seed_sha256": "d" * 64,
                "replicate_me_f1": samples,
            },
        }
        documents.append(document)

    result = paired_tool_differences(documents)

    assert len(result) == 1
    assert result[0]["contrast"] == "tool_a_minus_tool_b"
    assert result[0]["point_difference"] == 10.0
    assert result[0]["lower"] == pytest.approx(5.25)
    assert result[0]["upper"] == pytest.approx(14.75)


def test_aggregate_rejects_nested_ranking_fields(tmp_path: Path) -> None:
    score_path = tmp_path / "score.json"
    metrics_path = tmp_path / "metrics.json"
    score = _score("tool", 50.0)
    score["tuple_key"]["leaderboard"] = 1  # type: ignore[index]
    _write_json(score_path, score)
    _write_json(metrics_path, _metrics("tool"))

    with pytest.raises(SummaryError, match="forbidden ranking"):
        aggregate(
            score_paths=[score_path],
            metrics_paths=[metrics_path],
            output_score_tsv=tmp_path / "score.tsv",
            output_point_tsv=tmp_path / "points.tsv",
            output_metrics_tsv=tmp_path / "metrics.tsv",
            output_metrics_json=tmp_path / "metrics-out.json",
        )


def test_aggregate_rejects_metrics_changed_after_seal(tmp_path: Path) -> None:
    score_path = tmp_path / "score.json"
    metrics_path = tmp_path / "metrics.json"
    _write_json(score_path, _score("tool", 50.0))
    _write_json(metrics_path, _metrics("tool"))
    package_path = _package(score_path, metrics_path)
    manifest_path = _finalizer_manifest(score_path, metrics_path, package_path)
    metrics = _metrics("tool")
    metrics["records"][0]["value"] = 0.9  # type: ignore[index]
    _write_json(metrics_path, metrics)

    with pytest.raises(SummaryError, match="metrics_artifact hash"):
        aggregate(
            score_paths=[score_path],
            metrics_paths=[metrics_path],
            score_package_paths=[package_path],
            finalizer_manifest_paths=[manifest_path],
            output_score_tsv=tmp_path / "score.tsv",
            output_point_tsv=tmp_path / "points.tsv",
            output_metrics_tsv=tmp_path / "metrics.tsv",
            output_metrics_json=tmp_path / "metrics-out.json",
        )


def test_aggregate_rejects_metrics_not_bound_by_finalizer_manifest(
    tmp_path: Path,
) -> None:
    score_path = tmp_path / "score.json"
    metrics_path = tmp_path / "metrics.json"
    _write_json(score_path, _score("tool", 50.0))
    _write_json(metrics_path, _metrics("tool"))
    package_path = _package(score_path, metrics_path)
    manifest_path = _finalizer_manifest(score_path, metrics_path, package_path)
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["input_sha256"][str(metrics_path)] = "0" * 64
    _write_json(manifest_path, manifest)

    with pytest.raises(SummaryError, match="metrics JSON hash"):
        aggregate(
            score_paths=[score_path],
            metrics_paths=[metrics_path],
            score_package_paths=[package_path],
            finalizer_manifest_paths=[manifest_path],
            output_score_tsv=tmp_path / "score.tsv",
            output_point_tsv=tmp_path / "points.tsv",
            output_metrics_tsv=tmp_path / "metrics.tsv",
            output_metrics_json=tmp_path / "metrics-out.json",
        )


def test_aggregate_rejects_unsuccessful_finalizer_manifest(tmp_path: Path) -> None:
    score_path = tmp_path / "score.json"
    metrics_path = tmp_path / "metrics.json"
    _write_json(score_path, _score("tool", 50.0))
    _write_json(metrics_path, _metrics("tool"))
    package_path = _package(score_path, metrics_path)
    manifest_path = _finalizer_manifest(score_path, metrics_path, package_path)
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["status"] = "failed"
    _write_json(manifest_path, manifest)

    with pytest.raises(SummaryError, match="not successful"):
        aggregate(
            score_paths=[score_path],
            metrics_paths=[metrics_path],
            score_package_paths=[package_path],
            finalizer_manifest_paths=[manifest_path],
            output_score_tsv=tmp_path / "score.tsv",
            output_point_tsv=tmp_path / "points.tsv",
            output_metrics_tsv=tmp_path / "metrics.tsv",
            output_metrics_json=tmp_path / "metrics-out.json",
        )


def test_aggregate_cli_requires_finalizer_manifest(tmp_path: Path) -> None:
    score_path = tmp_path / "score.json"
    metrics_path = tmp_path / "metrics.json"
    _write_json(score_path, _score("tool", 50.0))
    _write_json(metrics_path, _metrics("tool"))
    package_path = _package(score_path, metrics_path)

    with pytest.raises(SystemExit) as exc_info:
        main(
            [
                "--score-json",
                str(score_path),
                "--metrics-json",
                str(metrics_path),
                "--score-package",
                str(package_path),
                "--output-score-tsv",
                str(tmp_path / "score.tsv"),
                "--output-point-tsv",
                str(tmp_path / "points.tsv"),
                "--output-metrics-tsv",
                str(tmp_path / "metrics.tsv"),
                "--output-metrics-json",
                str(tmp_path / "metrics-out.json"),
            ]
        )
    assert exc_info.value.code == 2
