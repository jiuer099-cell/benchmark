#!/usr/bin/env python3
"""Combine finalized scores from multiple input tracks into one HTML comparison."""

from __future__ import annotations

import argparse
import json
import sys
from collections.abc import Mapping, Sequence
from html import escape
from pathlib import Path

import yaml  # type: ignore[import-untyped]

class SuiteReportError(ValueError):
    """Raised when score cards are not cross-track comparable."""


def _manifest(path: Path) -> dict:
    value = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise SuiteReportError(f"tool manifest must be a mapping: {path}")
    return value


def _number(value: object) -> str:
    return "N/A" if value is None else f"{float(value):.2f}"


def render_suite(entries: Sequence[tuple[Mapping, Mapping]]) -> str:
    if not entries:
        raise SuiteReportError("at least one suite entry is required")
    first_score = entries[0][0]
    first_tuple = first_score["tuple_key"]
    comparison_key = (
        first_tuple["sample"],
        first_tuple["primary_truth_profile"],
        first_score["score_profile_sha256"],
    )
    rows = []
    for score, manifest in entries:
        tuple_key = score.get("tuple_key")
        if not isinstance(tuple_key, Mapping):
            raise SuiteReportError("score is missing tuple_key")
        key = (
            tuple_key.get("sample"),
            tuple_key.get("primary_truth_profile"),
            score.get("score_profile_sha256"),
        )
        if key != comparison_key:
            raise SuiteReportError(
                "suite scores must share sample, truth profile, and score profile"
            )
        capabilities = manifest.get("capabilities", {})
        technology = (
            capabilities.get("technology", [])
            if isinstance(capabilities, Mapping)
            else []
        )
        rows.append(
            "<tr>"
            f"<td>{escape(str(tuple_key.get('tool', '')))}</td>"
            f"<td>{escape(str(manifest.get('paradigm', 'unknown')))}</td>"
            f"<td>{escape(', '.join(map(str, technology)))}</td>"
            f"<td>{escape(str(tuple_key.get('official_score_mode', '')))}</td>"
            f"<td>{_number(score.get('comparable_score'))}</td>"
            f"<td>{_number(score.get('consensus_score'))}</td>"
            f"<td>{escape(str(score.get('truth_eligible_count', '')))}</td>"
            f"<td>{escape(str(score.get('total_evaluated', '')))}</td>"
            f"<td>{escape(str(score.get('score_status', '')))}</td>"
            "</tr>"
        )
    return """<!doctype html>
<html lang="zh-CN"><head><meta charset="utf-8">
<title>PGBench cross-track pangenome SV report</title>
<style>
body { font-family: system-ui, sans-serif; margin: 2rem; max-width: 96rem; }
table { border-collapse: collapse; width: 100%; }
th, td { border: 1px solid #ccd; padding: .5rem; text-align: left; }
th { background: #eef2ff; }
.notice { border-left: 4px solid #536dfe; background: #f3f5ff; padding: .8rem; }
</style></head><body>
<h1>HG002/GRCh38 泛基因组结构变异综合比较</h1>
<p class="notice">所有分数使用同一样本、truth profile 和冻结评分合同。
ComparableScore 可比较完整 pipeline 的实际结果；跨测序技术的差异同时包含
测序证据与算法影响，因此必须结合“测序技术”和“工具范式”解释。</p>
<table><thead><tr><th>工具</th><th>范式</th><th>测序技术</th><th>轨道</th>
<th>ComparableScore</th><th>ConsensusScore</th><th>Truth 数</th>
<th>输出数</th><th>状态</th></tr></thead><tbody>""" + "".join(rows) + """
</tbody></table></body></html>
"""


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--entry",
        action="append",
        nargs=2,
        metavar=("SCORE_JSON", "TOOL_YAML"),
        required=True,
    )
    parser.add_argument("--output", required=True, type=Path)
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        entries = []
        for raw_score, raw_manifest in args.entry:
            score = json.loads(Path(raw_score).read_text(encoding="utf-8"))
            if not isinstance(score, dict):
                raise SuiteReportError(f"score must be an object: {raw_score}")
            entries.append((score, _manifest(Path(raw_manifest))))
        args.output.parent.mkdir(parents=True, exist_ok=True)
        temporary = args.output.with_name(f".{args.output.name}.tmp")
        temporary.write_text(render_suite(entries), encoding="utf-8")
        temporary.replace(args.output)
    except (
        OSError,
        json.JSONDecodeError,
        SuiteReportError,
        yaml.YAMLError,
    ) as exc:
        print(f"render_suite_report: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
