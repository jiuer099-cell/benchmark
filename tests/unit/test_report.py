from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import pytest

SCRIPTS = Path(__file__).resolve().parents[2] / "workflow" / "scripts"
sys.path.insert(0, str(SCRIPTS))

from render_report import ReportError, render_report  # noqa: E402
from render_report_index import render_index  # noqa: E402


def _payload() -> dict:
    score_input = (
        Path(__file__).resolve().parents[1]
        / "fixtures"
        / "synthetic"
        / "evaluators"
        / "score_input.json"
    )
    return json.loads(score_input.read_text(encoding="utf-8"))


def _renderable_score() -> dict:
    return {
        "tuple_key": {
            "run_id": "synthetic",
            "sample": "HG002",
            "tool": "example",
            "official_score_mode": "end_to_end_from_reads",
            "primary_truth_profile": "truth",
        },
        "score_profile": "pgbench_v1",
        "score_profile_sha256": "a" * 64,
        "evaluation_mode": "synthetic_smoke",
        "score_status": "provisional",
        "pgbench_score": 50.0,
        "point_breakdown": {
            "evaluator.truvari.overall_event_f1": 5.0,
            "traceability.manifest_completeness": 2.0,
        },
    }


def _write_package(score_path: Path, score: dict) -> Path:
    package = {
        "package_schema_version": "pgbench.final_score_package.v1",
        "sealed": True,
        "seal_status": score["score_status"],
        "evaluation_mode": score["evaluation_mode"],
        "score": score,
        "score_artifact": {
            "path": str(score_path),
            "sha256": hashlib.sha256(score_path.read_bytes()).hexdigest(),
        },
        "gates": {
            "core_provenance_valid": True,
            "hash_lineage_complete": True,
            "manifest_completeness": True,
            "environment_complete": False,
            "run_context_complete": True,
        },
    }
    path = score_path.with_name("score-package.json")
    path.write_text(json.dumps(package), encoding="utf-8")
    return path


def _write_finalizer_manifest(
    score_path: Path,
    package_path: Path,
    score: dict,
) -> Path:
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
        "input_paths": [str(score_path)],
        "input_sha256": {
            str(score_path): hashlib.sha256(score_path.read_bytes()).hexdigest()
        },
        "output_paths": [str(package_path)],
        "output_sha256": {
            str(package_path): hashlib.sha256(package_path.read_bytes()).hexdigest()
        },
    }
    path = score_path.with_name("finalizer-manifest.json")
    path.write_text(json.dumps(payload), encoding="utf-8")
    return path


def test_report_is_no_ranking_and_preserves_provisional_status(
    tmp_path: Path,
) -> None:
    score = _renderable_score()
    score_path = tmp_path / "score.json"
    score_path.write_text(json.dumps(score), encoding="utf-8")
    package_path = _write_package(score_path, score)
    manifest_path = _write_finalizer_manifest(score_path, package_path, score)
    html = tmp_path / "report" / "index.html"
    score_tsv = tmp_path / "summary" / "score.tsv"
    breakdown = tmp_path / "summary" / "points.tsv"
    render_report(
        score_path,
        score_package=package_path,
        finalizer_manifest=manifest_path,
        html_output=html,
        score_tsv=score_tsv,
        point_breakdown_tsv=breakdown,
    )
    html_text = html.read_text(encoding="utf-8")
    assert "ConsensusScore:" in html_text
    assert "not rank tools" in html_text
    assert "provisional" in html_text
    header = score_tsv.read_text(encoding="utf-8").splitlines()[0]
    assert "rank" not in header.lower()
    assert "ConsensusScore" in header
    assert breakdown.read_text().count("\n") == 3


def test_report_rejects_ranking_fields(tmp_path: Path) -> None:
    score = _renderable_score()
    score["rank"] = 1
    score_path = tmp_path / "score.json"
    score_path.write_text(json.dumps(score), encoding="utf-8")
    package_path = _write_package(score_path, score)
    manifest_path = _write_finalizer_manifest(score_path, package_path, score)
    with pytest.raises(ReportError, match="ranking"):
        render_report(
            score_path,
            score_package=package_path,
            finalizer_manifest=manifest_path,
            html_output=tmp_path / "index.html",
            score_tsv=tmp_path / "score.tsv",
            point_breakdown_tsv=tmp_path / "points.tsv",
        )


def test_report_rejects_nested_ranking_fields(tmp_path: Path) -> None:
    score = _renderable_score()
    score["tuple_key"]["leaderboard"] = 1
    score_path = tmp_path / "score.json"
    score_path.write_text(json.dumps(score), encoding="utf-8")
    package_path = _write_package(score_path, score)
    manifest_path = _write_finalizer_manifest(score_path, package_path, score)

    with pytest.raises(ReportError, match="ranking"):
        render_report(
            score_path,
            score_package=package_path,
            finalizer_manifest=manifest_path,
            html_output=tmp_path / "report.html",
            score_tsv=tmp_path / "score.tsv",
            point_breakdown_tsv=tmp_path / "points.tsv",
        )


def test_report_rejects_score_changed_after_seal(tmp_path: Path) -> None:
    score = _renderable_score()
    score_path = tmp_path / "score.json"
    score_path.write_text(json.dumps(score), encoding="utf-8")
    package_path = _write_package(score_path, score)
    manifest_path = _write_finalizer_manifest(score_path, package_path, score)
    score["pgbench_score"] = 99.0
    score_path.write_text(json.dumps(score), encoding="utf-8")

    with pytest.raises(ReportError, match="does not match"):
        render_report(
            score_path,
            score_package=package_path,
            finalizer_manifest=manifest_path,
            html_output=tmp_path / "report.html",
            score_tsv=tmp_path / "score.tsv",
            point_breakdown_tsv=tmp_path / "points.tsv",
        )


def test_report_rejects_unsuccessful_finalizer_manifest(tmp_path: Path) -> None:
    score = _renderable_score()
    score_path = tmp_path / "score.json"
    score_path.write_text(json.dumps(score), encoding="utf-8")
    package_path = _write_package(score_path, score)
    manifest_path = _write_finalizer_manifest(score_path, package_path, score)
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["status"] = "failed"
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")

    with pytest.raises(ReportError, match="not successful"):
        render_report(
            score_path,
            score_package=package_path,
            finalizer_manifest=manifest_path,
            html_output=tmp_path / "report.html",
            score_tsv=tmp_path / "score.tsv",
            point_breakdown_tsv=tmp_path / "points.tsv",
        )


def test_report_rejects_package_not_bound_by_finalizer_manifest(
    tmp_path: Path,
) -> None:
    score = _renderable_score()
    score_path = tmp_path / "score.json"
    score_path.write_text(json.dumps(score), encoding="utf-8")
    package_path = _write_package(score_path, score)
    manifest_path = _write_finalizer_manifest(score_path, package_path, score)
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["output_sha256"][str(package_path)] = "0" * 64
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")

    with pytest.raises(ReportError, match="score package hash"):
        render_report(
            score_path,
            score_package=package_path,
            finalizer_manifest=manifest_path,
            html_output=tmp_path / "report.html",
            score_tsv=tmp_path / "score.tsv",
            point_breakdown_tsv=tmp_path / "points.tsv",
        )


def test_report_index_lists_tool_cards_without_ranking_language() -> None:
    html = render_index([Path("zeta.html"), Path("alpha.html")])
    assert html.index("alpha") < html.index("zeta")
    assert "不生成排名" in html
    assert "tool_cards/alpha.html" in html
