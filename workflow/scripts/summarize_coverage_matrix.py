#!/usr/bin/env python3
"""Summarize repeated 10x/20x/30x and full-depth PGBench scores."""

from __future__ import annotations

import argparse
import json
import math
import re
import statistics
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any, Sequence


RUN_PATTERN = re.compile(r"^(?P<base>.+)_(?P<coverage>10|20|30)x_seed(?P<seed>\d+)$")


class CoverageSummaryError(ValueError):
    """Raised when coverage results do not form the frozen repeated matrix."""


# Two-sided 95% Student-t critical values. Formal downsampling uses three
# frozen seeds (df=2), while the remaining entries keep the helper correct if
# later releases increase rather than decrease replication.
T_CRITICAL_975 = {
    1: 12.7062047364, 2: 4.3026527297, 3: 3.1824463053,
    4: 2.7764451052, 5: 2.5705818356, 6: 2.4469118488,
    7: 2.3646242516, 8: 2.3060041352, 9: 2.2621571629,
    10: 2.2281388520, 11: 2.2009851601, 12: 2.1788128297,
    13: 2.1603686565, 14: 2.1447866879, 15: 2.1314495456,
    16: 2.1199052992, 17: 2.1098155778, 18: 2.1009220402,
    19: 2.0930240544, 20: 2.0859634473, 21: 2.0796138447,
    22: 2.0738730679, 23: 2.0686576104, 24: 2.0638985616,
    25: 2.0595385528, 26: 2.0555294386, 27: 2.0518305165,
    28: 2.0484071418, 29: 2.0452296421, 30: 2.0422724563,
}


def mean_ci95(values: Sequence[float]) -> tuple[float | None, float | None]:
    if len(values) < 2:
        return None, None
    mean = statistics.fmean(values)
    standard_error = statistics.stdev(values) / math.sqrt(len(values))
    critical = T_CRITICAL_975.get(len(values) - 1, 1.9599639845)
    margin = critical * standard_error
    return max(0.0, mean - margin), min(100.0, mean + margin)


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
            ci_lower, ci_upper = mean_ci95(scores)
            rows.append(
                {
                    "tool": tool,
                    "coverage": coverage,
                    "metric_id": f"ME-F1_{str(coverage).upper()}X" if coverage != "full" else "ME-F1_FULL",
                    "role": "explanatory",
                    "affects_primary_score": False,
                    "replicates": len(scores),
                    "mean_me_f1": statistics.fmean(scores),
                    "sd_me_f1": statistics.stdev(scores) if len(scores) > 1 else None,
                    "ci95_lower_me_f1": ci_lower,
                    "ci95_upper_me_f1": ci_upper,
                    "ci95_method": (
                        "student_t_across_frozen_seeds"
                        if len(scores) > 1
                        else "undefined_single_full_depth_run"
                    ),
                    "min_me_f1": min(scores),
                    "max_me_f1": max(scores),
                    "runs": values,
                }
            )
    return {
        "schema_version": 2,
        "contract": "pgbench_repeated_coverage_summary_v2",
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
