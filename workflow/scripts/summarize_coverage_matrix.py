#!/usr/bin/env python3
"""Summarize repeated 10x/20x/30x and full-depth PGBench scores."""

from __future__ import annotations

import argparse
import json
import re
import statistics
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any, Sequence


RUN_PATTERN = re.compile(r"^(?P<base>.+)_(?P<coverage>10|20|30)x_seed(?P<seed>\d+)$")


class CoverageSummaryError(ValueError):
    """Raised when coverage results do not form the frozen repeated matrix."""


def _load(path: Path) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise CoverageSummaryError(f"cannot load score {path}") from exc
    if not isinstance(value, dict):
        raise CoverageSummaryError(f"score must be an object: {path}")
    return value


def summarize(paths: Sequence[Path]) -> dict[str, Any]:
    observations: dict[tuple[str, int | str], list[dict[str, Any]]] = defaultdict(list)
    seeds_by_coverage: dict[int, set[int]] = defaultdict(set)
    base_run: str | None = None
    for path in paths:
        score = _load(path)
        tuple_key = score.get("tuple_key")
        if not isinstance(tuple_key, dict):
            raise CoverageSummaryError(f"score has no tuple_key: {path}")
        tool = tuple_key.get("tool")
        run_id = tuple_key.get("run_id")
        value = score.get("benchmark_score")
        if not isinstance(tool, str) or not isinstance(run_id, str):
            raise CoverageSummaryError(f"score tuple is incomplete: {path}")
        if isinstance(value, bool) or not isinstance(value, int | float):
            raise CoverageSummaryError(f"score is not publishable: {path}")
        match = RUN_PATTERN.fullmatch(run_id)
        if match:
            current_base = match.group("base")
            coverage = int(match.group("coverage"))
            seed = int(match.group("seed"))
            seeds_by_coverage[coverage].add(seed)
            key: int | str = coverage
        elif run_id.endswith("_full"):
            current_base = run_id.removesuffix("_full")
            coverage = "full"
            seed = None
            key = "full"
        else:
            raise CoverageSummaryError(f"run_id does not encode coverage: {run_id}")
        if base_run is None:
            base_run = current_base
        elif base_run != current_base:
            raise CoverageSummaryError("scores come from different base runs")
        observations[(tool, key)].append(
            {"run_id": run_id, "seed": seed, "me_f1": float(value)}
        )

    if set(seeds_by_coverage) != {10, 20, 30}:
        raise CoverageSummaryError("matrix must contain 10x, 20x, and 30x")
    frozen_seeds = next(iter(seeds_by_coverage.values()))
    if len(frozen_seeds) < 3 or any(
        seeds != frozen_seeds for seeds in seeds_by_coverage.values()
    ):
        raise CoverageSummaryError("each downsampled depth must use the same >=3 seeds")
    tools = sorted({tool for tool, _coverage in observations})
    rows = []
    for tool in tools:
        if (tool, "full") not in observations:
            raise CoverageSummaryError(f"{tool} has no full-depth score")
        for coverage in (10, 20, 30, "full"):
            values = observations.get((tool, coverage), [])
            expected = 1 if coverage == "full" else len(frozen_seeds)
            if len(values) != expected:
                raise CoverageSummaryError(
                    f"{tool} {coverage}x has {len(values)} results, expected {expected}"
                )
            scores = [row["me_f1"] for row in values]
            rows.append(
                {
                    "tool": tool,
                    "coverage": coverage,
                    "replicates": len(scores),
                    "mean_me_f1": statistics.fmean(scores),
                    "sd_me_f1": statistics.pstdev(scores) if len(scores) > 1 else 0.0,
                    "min_me_f1": min(scores),
                    "max_me_f1": max(scores),
                    "runs": values,
                }
            )
    return {
        "schema_version": 1,
        "contract": "pgbench_repeated_coverage_summary_v1",
        "base_run_id": base_run,
        "coverages": [10, 20, 30, "full"],
        "downsampling_seeds": sorted(frozen_seeds),
        "tools": tools,
        "rows": rows,
    }


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--score", action="append", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        payload = summarize(args.score)
        args.output.parent.mkdir(parents=True, exist_ok=True)
        temporary = args.output.with_name(f".{args.output.name}.tmp")
        temporary.write_text(
            json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        temporary.replace(args.output)
    except (OSError, CoverageSummaryError, statistics.StatisticsError) as exc:
        print(f"summarize_coverage_matrix: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
