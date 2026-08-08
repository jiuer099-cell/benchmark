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


def render_index(
    card_paths: Sequence[Path],
    scores: Sequence[Mapping[str, object]] | None = None,
    manifests: Sequence[Mapping[str, object]] | None = None,
) -> str:
    """Render links only for compatibility, or a comparable-score table."""

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
        rows.append(
            "<tr>"
            f'<td><a href="tool_cards/{escape(card.name, quote=True)}">'
            f"{escape(tool)}</a></td>"
            f"<td>{escape(str(manifest.get('paradigm', 'unknown')))}</td>"
            f"<td>{escape(', '.join(map(str, technologies)))}</td>"
            f"<td>{escape(str(tuple_key.get('official_score_mode', '')))}</td>"
            f"<td>{_score_value(score, 'pangenome_genotyping_score')}</td>"
            f"<td>{_score_value(score, 'non_reference_f1_score')}</td>"
            f"<td>{_score_value(score, 'panel_coverage')}</td>"
            f"<td>{_score_value(score, 'global_end_to_end_sv_recovery_score')}</td>"
            f"<td>{_score_value(score, 'consensus_score')}</td>"
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
<p class="notice">ComparableScore 在固定 HG002/GRCh38/GIAB truth universe
上同时惩罚 FP 与 FN。跨长短读长时，该分数比较完整 pipeline 的实际效果；
同一测序轨道内的比较才代表更严格的算法公平性。资源消耗不参与得分。</p>
<table><thead><tr>
<th>工具</th><th>范式</th><th>测序技术</th><th>运行轨道</th>
<th>Pangenome Genotyping Score</th><th>Non-reference F1</th>
<th>Panel coverage (fraction)</th><th>Global End-to-End SV Recovery</th>
<th>ConsensusScore</th>
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
