#!/usr/bin/env python3
"""Fuse three standardized formal evaluator ledgers into consensus counts."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
from pathlib import Path

import yaml  # type: ignore[import-untyped]


class ConsensusMetricError(ValueError):
    """Raised when evaluator ledgers do not describe one identical universe."""


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
    profile = yaml.safe_load(args.score_profile.read_text(encoding="utf-8"))
    profile_id = profile["profile"]["id"]
    provenance_id = provenance_bundle(args.evaluator_manifests)
    categories = (
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
                "denominator": total,
                "eligible_count": total,
                "universe_id": "formal_query_result_universe_v1",
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
