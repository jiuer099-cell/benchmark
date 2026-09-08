#!/usr/bin/env python3
"""Summarize separately measured one-time and per-sample resource stages."""

from __future__ import annotations

import argparse
import csv
import json
import statistics
import sys
from collections import defaultdict
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any


SCOPES = {"one_time_build", "per_sample"}
METRIC_ALIASES = {
    "wall_seconds": ("wall_seconds", "s", "seconds", "walltime", "wall_time"),
    "cpu_seconds": ("cpu_seconds", "cpu_time"),
    "peak_rss_mb": ("peak_rss_mb", "max_rss_mb", "max_rss"),
}


class ResourceStageError(ValueError):
    """Raised when runtime measurements cannot support a fair split."""


def _load_rows(path: Path) -> list[Mapping[str, Any]]:
    text = path.read_text(encoding="utf-8").strip()
    if not text:
        raise ResourceStageError(f"resource measurement is empty: {path}")
    if text.lstrip().startswith("{"):
        rows = [json.loads(line) for line in text.splitlines()]
    else:
        dialect = csv.Sniffer().sniff(text[:4096], delimiters="\t,")
        rows = list(csv.DictReader(text.splitlines(), dialect=dialect))
    if not rows or any(not isinstance(row, Mapping) for row in rows):
        raise ResourceStageError(f"resource measurement has invalid rows: {path}")
    return rows


def _metric(row: Mapping[str, Any], metric: str) -> float:
    for name in METRIC_ALIASES[metric]:
        value = row.get(name)
        if value not in {None, ""}:
            numeric = float(value)
            if numeric < 0:
                raise ResourceStageError(f"{metric} must be non-negative")
            return numeric
    raise ResourceStageError(f"measurement row lacks {metric}")


def _disk_bytes(paths: Sequence[Path]) -> int:
    total = 0
    for path in paths:
        if not path.exists():
            raise ResourceStageError(f"index artifact does not exist: {path}")
        if path.is_file():
            total += path.stat().st_size
        else:
            total += sum(item.stat().st_size for item in path.rglob("*") if item.is_file())
    return total


def summarize(
    measurements: Sequence[tuple[str, str, Path]],
    *,
    index_artifacts: Sequence[Path] = (),
    cohort_sizes: Sequence[int] = (1, 10, 100, 1000),
) -> dict[str, Any]:
    if not measurements:
        raise ResourceStageError("at least one stage measurement is required")
    observed: dict[tuple[str, str], list[dict[str, float]]] = defaultdict(list)
    for stage, scope, path in measurements:
        if not stage or scope not in SCOPES:
            raise ResourceStageError(f"invalid stage/scope: {stage!r}/{scope!r}")
        for row in _load_rows(path):
            observed[(stage, scope)].append(
                {metric: _metric(row, metric) for metric in METRIC_ALIASES}
            )
    if {scope for _stage, scope in observed} != SCOPES:
        raise ResourceStageError(
            "formal resource report requires both one_time_build and per_sample measurements"
        )
    stages: list[dict[str, Any]] = []
    scope_medians: dict[str, dict[str, float]] = {
        scope: {metric: 0.0 for metric in METRIC_ALIASES} for scope in SCOPES
    }
    for (stage, scope), rows in sorted(observed.items()):
        median = {
            metric: statistics.median(row[metric] for row in rows)
            for metric in METRIC_ALIASES
        }
        for metric, value in median.items():
            if metric == "peak_rss_mb":
                scope_medians[scope][metric] = max(
                    scope_medians[scope][metric], value
                )
            else:
                scope_medians[scope][metric] += value
        stages.append(
            {
                "stage": stage,
                "scope": scope,
                "repeat_count": len(rows),
                "median": median,
            }
        )
    index_disk_bytes = _disk_bytes(index_artifacts)
    build = scope_medians["one_time_build"]
    sample = scope_medians["per_sample"]
    amortized = [
        {
            "cohort_size": size,
            "wall_seconds_total": build["wall_seconds"] + size * sample["wall_seconds"],
            "cpu_seconds_total": build["cpu_seconds"] + size * sample["cpu_seconds"],
            "wall_seconds_per_sample": build["wall_seconds"] / size + sample["wall_seconds"],
            "cpu_seconds_per_sample": build["cpu_seconds"] / size + sample["cpu_seconds"],
        }
        for size in cohort_sizes
        if size > 0
    ]
    return {
        "schema_version": 1,
        "contract": "pgbench_resource_stage_summary_v1",
        "status": "valid",
        "stage_split_complete": True,
        "aggregation": "sum_stage_medians_peak_ram_max",
        "stages": stages,
        "one_time_build": {
            **build,
            "index_disk_bytes": index_disk_bytes,
        },
        "per_sample": sample,
        "amortized": amortized,
    }


def _measurement(value: str) -> tuple[str, str, Path]:
    parts = value.split("=", 2)
    if len(parts) != 3:
        raise argparse.ArgumentTypeError("measurement must be STAGE=SCOPE=PATH")
    return parts[0], parts[1], Path(parts[2])


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--measurement", action="append", required=True, type=_measurement)
    parser.add_argument("--index-artifact", action="append", default=[], type=Path)
    parser.add_argument("--cohort-size", action="append", type=int)
    parser.add_argument("--output", required=True, type=Path)
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        payload = summarize(
            args.measurement,
            index_artifacts=args.index_artifact,
            cohort_sizes=args.cohort_size or (1, 10, 100, 1000),
        )
        args.output.parent.mkdir(parents=True, exist_ok=True)
        temporary = args.output.with_name(f".{args.output.name}.tmp")
        temporary.write_text(
            json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        temporary.replace(args.output)
    except (OSError, ValueError, json.JSONDecodeError, ResourceStageError) as exc:
        print(f"summarize_resource_stages: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
