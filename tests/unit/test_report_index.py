from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SCRIPTS = ROOT / "workflow" / "scripts"
sys.path.insert(0, str(SCRIPTS))

from render_report_index import main, render_index  # noqa: E402


def _score() -> dict:
    return {
        "tuple_key": {
            "run_id": "run1",
            "sample": "HG002",
            "tool": "tool_a",
            "official_score_mode": "end_to_end_from_reads",
            "primary_truth_profile": "truth1",
        },
        "score_profile": "pgbench_v1",
        "score_profile_sha256": "a" * 64,
        "evaluation_mode": "synthetic_smoke",
        "score_status": "provisional",
        "benchmark_score": 75.0,
        "evaluator_scores": {"truvari": 74.0, "aardvark": 75.0, "vcfdist": 76.0},
        "evaluator_range": 2.0,
        "evaluator_sd": 0.82,
        "point_breakdown": {"evaluator_accuracy": 50.0},
    }


def _metrics() -> dict:
    return {
        "schema_version": 1,
        "dictionary_id": "pgbench_metrics_v1",
        "tuple": {
            "run_id": "run1",
            "sample_id": "HG002",
            "tool_id": "tool_a",
            "official_score_mode": "end_to_end_from_reads",
            "primary_truth_profile": "truth1",
            "score_profile": "pgbench_v1",
        },
        "records": [
            {
                "metric_id": "truvari.event.overall.f1",
                "value_type": "ratio",
                "value": 0.8,
                "status": "defined",
                "evaluator": "truvari",
                "numerator": 8,
                "denominator": 10,
                "eligible_count": 10,
                "universe_id": "fixture",
                "undefined_reason": None,
                "parser_id": "fixture",
                "parser_source_field": "f1",
                "provenance_manifest_id": "b" * 64,
                "strata": {},
                "confidence_interval": None,
            }
        ],
    }


def _package(score: Path, metrics: Path) -> Path:
    score_value = json.loads(score.read_text(encoding="utf-8"))
    payload = {
        "package_schema_version": "pgbench.final_score_package.v1",
        "sealed": True,
        "score": score_value,
        "score_artifact": {
            "path": str(score),
            "sha256": hashlib.sha256(score.read_bytes()).hexdigest(),
        },
        "metrics_artifact": {
            "path": str(metrics),
            "sha256": hashlib.sha256(metrics.read_bytes()).hexdigest(),
        },
        "gates": {
            "core_provenance_valid": True,
            "hash_lineage_complete": True,
            "manifest_completeness": True,
        },
    }
    path = score.with_name("score-package.json")
    path.write_text(json.dumps(payload), encoding="utf-8")
    return path


def _finalizer_manifest(score: Path, metrics: Path, package: Path) -> Path:
    score_value = json.loads(score.read_text(encoding="utf-8"))
    tuple_key = score_value["tuple_key"]
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
        "input_paths": [str(score), str(metrics)],
        "input_sha256": {
            str(score): hashlib.sha256(score.read_bytes()).hexdigest(),
            str(metrics): hashlib.sha256(metrics.read_bytes()).hexdigest(),
        },
        "output_paths": [str(package)],
        "output_sha256": {
            str(package): hashlib.sha256(package.read_bytes()).hexdigest()
        },
    }
    path = score.with_name("finalizer-manifest.json")
    path.write_text(json.dumps(payload), encoding="utf-8")
    return path


def test_index_and_run_summaries_are_written_without_ranking(tmp_path: Path) -> None:
    card = tmp_path / "tool_a.html"
    card.write_text("<p>tool</p>\n", encoding="utf-8")
    score = tmp_path / "score.json"
    score.write_text(json.dumps(_score()), encoding="utf-8")
    metrics = tmp_path / "metrics.json"
    metrics.write_text(json.dumps(_metrics()), encoding="utf-8")
    package = _package(score, metrics)
    finalizer_manifest = _finalizer_manifest(score, metrics, package)
    index = tmp_path / "index.html"
    score_tsv = tmp_path / "score.tsv"
    points_tsv = tmp_path / "point_breakdown.tsv"
    metrics_tsv = tmp_path / "metrics.long.tsv"
    metrics_json = tmp_path / "metrics.all.json"

    assert (
        main(
            [
                "--card",
                str(card),
                "--output",
                str(index),
                "--score-json",
                str(score),
                "--metrics-json",
                str(metrics),
                "--score-package",
                str(package),
                "--finalizer-manifest",
                str(finalizer_manifest),
                "--output-score-tsv",
                str(score_tsv),
                "--output-point-tsv",
                str(points_tsv),
                "--output-metrics-tsv",
                str(metrics_tsv),
                "--output-metrics-json",
                str(metrics_json),
            ]
        )
        == 0
    )
    combined = "\n".join(
        path.read_text(encoding="utf-8")
        for path in (index, score_tsv, points_tsv, metrics_tsv, metrics_json)
    ).lower()
    assert "leaderboard" not in combined
    assert "\trank\t" not in combined
    assert "tool_a" in combined


def test_summary_arguments_are_all_or_none(tmp_path: Path) -> None:
    card = tmp_path / "tool.html"
    card.write_text("ok", encoding="utf-8")
    assert (
        main(
            [
                "--card",
                str(card),
                "--output",
                str(tmp_path / "index.html"),
                "--score-json",
                str(tmp_path / "score.json"),
            ]
        )
        == 2
    )


def test_render_index_uses_filename_order_not_score_order() -> None:
    html = render_index([Path("z_tool.html"), Path("a_tool.html")])
    assert html.index("a_tool") < html.index("z_tool")
    assert "no ranking" in html


def test_render_index_displays_comparison_track_and_actual_technology() -> None:
    score = _score()
    score["formal_analysis"] = {
        "evidence_profile": {"actual_technology": "hifi"},
        "comparison_track": {"id": "panel-genotyping-hifi-track"},
        "comparison_track_sha256": "c" * 64,
    }
    html = render_index(
        [Path("tool_a.html")],
        [score],
        [
            {
                "id": "tool_a",
                "paradigm": "graph",
                "capabilities": {"technology": ["short_read"]},
            }
        ],
    )
    assert "panel-genotyping-hifi-track" in html
    assert "c" * 64 in html
    assert ">hifi<" in html
    assert "short_read" not in html
    assert "不生成跨轨道排名" in html
