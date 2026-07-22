#!/usr/bin/env python3
"""Render a no-ranking report index from individual tool score cards."""

from __future__ import annotations

import argparse
import sys
from collections.abc import Sequence
from html import escape
from pathlib import Path

try:
    from aggregate_run_summary import aggregate
except ModuleNotFoundError:  # pragma: no cover - package-style invocation
    from .aggregate_run_summary import aggregate


def render_index(card_paths: Sequence[Path]) -> str:
    cards = []
    for card_path in sorted(card_paths, key=lambda path: path.name):
        cards.append(
            '<li><a href="tool_cards/'
            + escape(card_path.name, quote=True)
            + '">'
            + escape(card_path.stem)
            + "</a></li>"
        )
    return (
        "<!doctype html>\n"
        '<html lang="zh-CN"><head><meta charset="utf-8">'
        "<title>PGBench 报告</title></head><body>"
        "<h1>PGBench 单工具得分报告</h1>"
        "<p>每个链接对应一个工具/运行元组；本页面不生成排名。</p>"
        "<ul>" + "".join(cards) + "</ul></body></html>\n"
    )


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--card", action="append", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--score-json", action="append", type=Path)
    parser.add_argument("--metrics-json", action="append", type=Path)
    parser.add_argument("--score-package", action="append", type=Path)
    parser.add_argument("--finalizer-manifest", action="append", type=Path)
    parser.add_argument("--output-score-tsv", type=Path)
    parser.add_argument("--output-point-tsv", type=Path)
    parser.add_argument("--output-metrics-tsv", type=Path)
    parser.add_argument("--output-metrics-json", type=Path)
    return parser.parse_args(argv)


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
        args.output.parent.mkdir(parents=True, exist_ok=True)
        temporary = args.output.with_name(f".{args.output.name}.tmp")
        temporary.write_text(render_index(args.card), encoding="utf-8")
        temporary.replace(args.output)
    except (OSError, ValueError) as exc:
        print(f"render_report_index: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
