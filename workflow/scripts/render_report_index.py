#!/usr/bin/env python3
"""Render a cross-tool pangenome-SV comparison from sealed score cards."""

from __future__ import annotations

import argparse
import json
import sys
from collections.abc import Mapping, Sequence
from html import escape
from pathlib import Path

import yaml  # type: ignore[import-untyped]

try:
    from aggregate_run_summary import aggregate
except ModuleNotFoundError:  # pragma: no cover - package-style invocation
    from .aggregate_run_summary import aggregate


def _score_value(score: Mapping[str, object], field: str) -> str:
    value = score.get(field)
    return "N/A" if value is None else f"{float(value):.2f}"


def _track_fields(score: Mapping[str, object]) -> tuple[str, str, str]:
    analysis = score.get("formal_analysis")
    if not isinstance(analysis, Mapping):
        return "untracked", "untracked", "unknown"
    track = analysis.get("comparison_track")
    evidence = analysis.get("evidence_profile")
    technology = (
        str(evidence.get("actual_technology", "unknown"))
        if isinstance(evidence, Mapping)
        else "unknown"
    )
    if not isinstance(track, Mapping):
        return "untracked", "untracked", technology
    track_id = str(track.get("id", "untracked"))
    track_sha256 = str(analysis.get("comparison_track_sha256", "untracked"))
    return track_id, track_sha256, technology


def render_index(
    card_paths: Sequence[Path],
    scores: Sequence[Mapping[str, object]] | None = None,
    manifests: Sequence[Mapping[str, object]] | None = None,
) -> str:
    """Render links only, or a table of sealed ME-F1 results."""

    ordered_cards = sorted(card_paths, key=lambda path: path.name)
    if scores is None:
        links = "".join(
            '<li><a href="tool_cards/'
            + escape(card.name, quote=True)
            + '">'
            + escape(card.stem)
            + "</a></li>"
            for card in ordered_cards
        )
        return (
            "<!doctype html>\n"
            '<html lang="en"><head><meta charset="utf-8">'
            "<title>PGBench report</title></head><body>"
            "<h1>PGBench tool score cards</h1>"
            "<p>Each link is one tool run; no ranking is produced.</p>"
            "<ul>" + links + "</ul></body></html>\n"
        )
    if len(scores) != len(card_paths):
        raise ValueError("score and card inputs must be paired")
    if manifests is not None and len(manifests) != len(scores):
        raise ValueError("tool manifests and scores must be paired")

    cards_by_tool = {card.stem: card for card in card_paths}
    manifests_by_tool = {
        str(manifest.get("id")): manifest for manifest in (manifests or [])
    }
    rows: list[str] = []
    for score in scores:
        tuple_key = score.get("tuple_key")
        if not isinstance(tuple_key, Mapping):
            raise ValueError("score is missing tuple_key")
        tool = str(tuple_key.get("tool"))
        card = cards_by_tool.get(tool)
        if card is None:
            raise ValueError(f"no score card is paired with tool {tool}")
        manifest = manifests_by_tool.get(tool, {})
        capabilities = manifest.get("capabilities", {})
        technologies = (
            capabilities.get("technology", [])
            if isinstance(capabilities, Mapping)
            else []
        )
        track_id, track_sha256, actual_technology = _track_fields(score)
        technology_text = (
            actual_technology
            if actual_technology != "unknown"
            else ", ".join(map(str, technologies))
        )
        rows.append(
            "<tr>"
            f'<td><a href="tool_cards/{escape(card.name, quote=True)}">'
            f"{escape(tool)}</a></td>"
            f"<td>{escape(str(manifest.get('paradigm', 'unknown')))}</td>"
            f"<td>{escape(technology_text)}</td>"
            f"<td>{escape(track_id)}</td>"
            f"<td><code>{escape(track_sha256)}</code></td>"
            f"<td>{_score_value(score, 'benchmark_score')}</td>"
            f"<td>{_score_value(score.get('evaluator_scores', {}), 'truvari')}</td>"
            f"<td>{_score_value(score.get('evaluator_scores', {}), 'aardvark')}</td>"
            f"<td>{_score_value(score.get('evaluator_scores', {}), 'vcfdist')}</td>"
            f"<td>{_score_value(score, 'evaluator_range')}</td>"
            f"<td>{_score_value(score, 'evaluator_sd')}</td>"
            f"<td>{escape(str(score.get('truth_eligible_count', '')))}</td>"
            f"<td>{escape(str(score.get('total_evaluated', '')))}</td>"
            f"<td>{escape(str(score.get('score_status', '')))}</td>"
            "</tr>"
        )
    return """<!doctype html>
<html lang="zh-CN"><head><meta charset="utf-8">
<title>PGBench pangenome SV comparison</title>
<style>
body { font-family: system-ui, sans-serif; margin: 2rem; max-width: 96rem; }
table { border-collapse: collapse; width: 100%; }
th, td { border: 1px solid #ccd; padding: .5rem; text-align: left; }
th { background: #eef2ff; }
.notice { border-left: 4px solid #536dfe; background: #f3f5ff; padding: .8rem; }
</style></head><body>
<h1>PGBench 泛基因组结构变异工具比较</h1>
<p class="notice">ME-F1 是唯一主分：Truvari、Aardvark-GT 与 vcfdist
在同一固定 truth denominator 上重算的 genotype-aware F1 算术平均。评估器差异
仅作诊断，资源消耗不参与得分。
比较只能在完全相同的运行轨道内进行；本表不生成跨轨道排名。</p>
<table><thead><tr>
<th>工具</th><th>范式</th><th>实际测序技术</th><th>运行轨道</th><th>轨道 SHA-256</th>
<th>ME-F1</th><th>Truvari F1</th><th>Aardvark-GT F1</th><th>vcfdist F1</th>
<th>Evaluator range</th><th>Evaluator SD</th>
<th>Truth 数量</th><th>输出数量</th><th>状态</th>
</tr></thead><tbody>""" + "".join(rows) + """</tbody></table>
</body></html>
"""


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--card", action="append", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--score-json", action="append", type=Path)
    parser.add_argument("--metrics-json", action="append", type=Path)
    parser.add_argument("--score-package", action="append", type=Path)
    parser.add_argument("--finalizer-manifest", action="append", type=Path)
    parser.add_argument("--tool-manifest", action="append", type=Path)
    parser.add_argument("--output-score-tsv", type=Path)
    parser.add_argument("--output-point-tsv", type=Path)
    parser.add_argument("--output-metrics-tsv", type=Path)
    parser.add_argument("--output-metrics-json", type=Path)
    return parser.parse_args(argv)


def _load_objects(paths: Sequence[Path], *, yaml_mode: bool) -> list[dict]:
    values = []
    for path in paths:
        text = path.read_text(encoding="utf-8")
        value = yaml.safe_load(text) if yaml_mode else json.loads(text)
        if not isinstance(value, dict):
            raise ValueError(f"expected an object in {path}")
        values.append(value)
    return values


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        for card in args.card:
            if not card.is_file():
                raise OSError(f"report card does not exist: {card}")
        summary_values = (
            args.score_json,
            args.metrics_json,
            args.score_package,
            args.finalizer_manifest,
            args.output_score_tsv,
            args.output_point_tsv,
            args.output_metrics_tsv,
            args.output_metrics_json,
        )
        if any(value is not None for value in summary_values):
            if not all(value is not None for value in summary_values):
                raise OSError(
                    "summary inputs and outputs must be supplied as one complete set"
                )
            aggregate(
                score_paths=args.score_json,
                metrics_paths=args.metrics_json,
                score_package_paths=args.score_package,
                finalizer_manifest_paths=args.finalizer_manifest,
                output_score_tsv=args.output_score_tsv,
                output_point_tsv=args.output_point_tsv,
                output_metrics_tsv=args.output_metrics_tsv,
                output_metrics_json=args.output_metrics_json,
            )
        scores = (
            _load_objects(args.score_json, yaml_mode=False)
            if args.score_json is not None
            else None
        )
        manifests = (
            _load_objects(args.tool_manifest, yaml_mode=True)
            if args.tool_manifest is not None
            else None
        )
        args.output.parent.mkdir(parents=True, exist_ok=True)
        temporary = args.output.with_name(f".{args.output.name}.tmp")
        temporary.write_text(
            render_index(args.card, scores, manifests),
            encoding="utf-8",
        )
        temporary.replace(args.output)
    except (OSError, ValueError, json.JSONDecodeError, yaml.YAMLError) as exc:
        print(f"render_report_index: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
