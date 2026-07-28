#!/usr/bin/env python3
"""Fuse three standardized formal evaluator ledgers into consensus counts."""

from __future__ import annotations

import argparse
import bisect
import csv
import gzip
import hashlib
import json
import sys
from collections import defaultdict
from pathlib import Path
from typing import TextIO

import yaml  # type: ignore[import-untyped]


class ConsensusMetricError(ValueError):
    """Raised when evaluator ledgers do not describe one identical universe."""


def open_text(path: Path) -> TextIO:
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("r", encoding="utf-8")


def load_benchmark_regions(path: Path) -> dict[str, tuple[list[int], list[int]]]:
    raw: dict[str, list[tuple[int, int]]] = defaultdict(list)
    with path.open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                raise ConsensusMetricError(
                    f"malformed BED record at {path}:{line_number}"
                )
            start, end = int(fields[1]), int(fields[2])
            if start < 0 or end <= start:
                raise ConsensusMetricError(
                    f"invalid BED interval at {path}:{line_number}"
                )
            raw[fields[0]].append((start, end))
    if not raw:
        raise ConsensusMetricError(f"benchmark BED contains no regions: {path}")

    merged: dict[str, tuple[list[int], list[int]]] = {}
    for contig, intervals in raw.items():
        compact: list[list[int]] = []
        for start, end in sorted(intervals):
            if compact and start <= compact[-1][1]:
                compact[-1][1] = max(compact[-1][1], end)
            else:
                compact.append([start, end])
        merged[contig] = (
            [interval[0] for interval in compact],
            [interval[1] for interval in compact],
        )
    return merged


def overlaps_regions(
    regions: dict[str, tuple[list[int], list[int]]],
    contig: str,
    start: int,
    end: int,
) -> bool:
    interval_index = regions.get(contig)
    if interval_index is None:
        return False
    starts, ends = interval_index
    index = bisect.bisect_right(starts, end - 1) - 1
    return index >= 0 and ends[index] > start


def info_end(position_1_based: int, ref: str, info: str) -> int:
    for field in info.split(";"):
        if field.startswith("END="):
            try:
                return max(position_1_based, int(field[4:]))
            except ValueError as exc:
                raise ConsensusMetricError(f"invalid INFO/END value: {field}") from exc
    return position_1_based + max(1, len(ref)) - 1


def eligible_truth_count(truth_vcf: Path, benchmark_bed: Path) -> int:
    regions = load_benchmark_regions(benchmark_bed)
    count = 0
    with open_text(truth_vcf) as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                raise ConsensusMetricError(
                    f"malformed truth VCF record at {truth_vcf}:{line_number}"
                )
            position = int(fields[1])
            end = info_end(position, fields[3], fields[7])
            if overlaps_regions(regions, fields[0], position - 1, end):
                count += 1
    if count == 0:
        raise ConsensusMetricError(
            "primary truth VCF has no records overlapping the benchmark BED"
        )
    return count


def load_ledger(path: Path, expected_evaluator: str) -> dict[str, bool]:
    votes: dict[str, bool] = {}
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"result_id", "evaluator", "correct"}
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise ConsensusMetricError(f"unsupported vote ledger header: {path}")
        for row in reader:
            if row["evaluator"] != expected_evaluator:
                raise ConsensusMetricError(
                    f"{path} contains evaluator {row['evaluator']!r}"
                )
            result_id = row["result_id"]
            if result_id in votes:
                raise ConsensusMetricError(f"duplicate result ID {result_id} in {path}")
            if row["correct"] not in {"0", "1"}:
                raise ConsensusMetricError(f"invalid binary vote in {path}")
            votes[result_id] = row["correct"] == "1"
    if not votes:
        raise ConsensusMetricError(f"empty vote ledger: {path}")
    return votes


def provenance_bundle(paths: list[Path]) -> str:
    digest = hashlib.sha256()
    for path in paths:
        digest.update(path.read_bytes())
    return digest.hexdigest()


def materialize(args: argparse.Namespace) -> dict:
    ledgers = {
        "truvari": load_ledger(args.truvari_ledger, "truvari"),
        "aardvark": load_ledger(args.aardvark_ledger, "aardvark"),
        "vcfdist": load_ledger(args.vcfdist_ledger, "vcfdist"),
    }
    universes = {name: set(votes) for name, votes in ledgers.items()}
    first = universes["truvari"]
    for name, universe in universes.items():
        if universe != first:
            missing = sorted(first - universe)[:5]
            extra = sorted(universe - first)[:5]
            raise ConsensusMetricError(
                f"{name} result universe differs; missing={missing}, extra={extra}"
            )

    counts = [0, 0, 0, 0]
    for result_id in sorted(first):
        vote_count = sum(ledgers[name][result_id] for name in ledgers)
        counts[vote_count] += 1
    total = len(first)
    truth_total = eligible_truth_count(args.truth_vcf, args.benchmark_bed)
    profile = yaml.safe_load(args.score_profile.read_text(encoding="utf-8"))
    profile_id = profile["profile"]["id"]
    provenance_id = provenance_bundle(args.evaluator_manifests)
    categories = (
        ("benchmark.truth.eligible.count", truth_total),
        ("consensus.all_three_correct.count", counts[3]),
        ("consensus.exactly_two_correct.count", counts[2]),
        ("consensus.exactly_one_correct.count", counts[1]),
        ("consensus.none_correct.count", counts[0]),
    )
    records = []
    for metric_id, value in categories:
        records.append(
            {
                "metric_id": metric_id,
                "value_type": "count",
                "value": value,
                "status": "defined",
                "evaluator": "fusion",
                "numerator": value,
                "denominator": (
                    truth_total
                    if metric_id == "benchmark.truth.eligible.count"
                    else total
                ),
                "eligible_count": (
                    truth_total
                    if metric_id == "benchmark.truth.eligible.count"
                    else total
                ),
                "universe_id": (
                    "primary_truth_benchmark_universe_v1"
                    if metric_id == "benchmark.truth.eligible.count"
                    else "formal_query_result_universe_v1"
                ),
                "undefined_reason": None,
                "parser_id": "three_evaluator_binary_vote_fusion_v1",
                "parser_source_field": "correct",
                "provenance_manifest_id": provenance_id,
                "strata": {},
            }
        )
    return {
        "schema_version": 1,
        "dictionary_id": "pgbench_metrics_v1",
        "tuple": {
            "run_id": args.run_id,
            "sample_id": args.sample_id,
            "tool_id": args.tool_id,
            "official_score_mode": args.official_score_mode,
            "primary_truth_profile": args.primary_truth_profile,
            "score_profile": profile_id,
        },
        "records": records,
    }


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--truvari-ledger", required=True, type=Path)
    parser.add_argument("--aardvark-ledger", required=True, type=Path)
    parser.add_argument("--vcfdist-ledger", required=True, type=Path)
    parser.add_argument(
        "--evaluator-manifest",
        action="append",
        dest="evaluator_manifests",
        required=True,
        type=Path,
    )
    parser.add_argument("--score-profile", required=True, type=Path)
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--tool-id", required=True)
    parser.add_argument("--official-score-mode", required=True)
    parser.add_argument("--primary-truth-profile", required=True)
    parser.add_argument("--truth-vcf", required=True, type=Path)
    parser.add_argument("--benchmark-bed", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        payload = materialize(args)
        args.output.parent.mkdir(parents=True, exist_ok=True)
        temporary = args.output.with_name(f".{args.output.name}.tmp")
        temporary.write_text(
            json.dumps(payload, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        temporary.replace(args.output)
    except (OSError, KeyError, TypeError, ConsensusMetricError) as exc:
        print(f"materialize_formal_consensus_metrics: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
