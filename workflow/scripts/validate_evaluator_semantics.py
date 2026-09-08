#!/usr/bin/env python3
"""Validate that all evaluator ledgers obey the common GT-F1 contract."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
from pathlib import Path

import yaml  # type: ignore[import-untyped]


EVALUATORS = ("truvari", "aardvark", "vcfdist")


class SemanticValidationError(ValueError):
    pass


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def validate(cases_path: Path, ledgers: dict[str, Path]) -> dict[str, object]:
    suite = yaml.safe_load(cases_path.read_text(encoding="utf-8"))
    if not isinstance(suite, dict) or suite.get("contract") != (
        "pgbench_genotype_evaluator_contract_v1"
    ):
        raise SemanticValidationError("unsupported semantic suite contract")
    cases = suite.get("cases")
    if not isinstance(cases, list) or not cases:
        raise SemanticValidationError("semantic suite has no cases")
    expected = {
        str(case["id"]): {
            key: int(case["expected"][key]) for key in ("tp", "fp", "fn")
        }
        for case in cases
    }
    evaluator_results: dict[str, object] = {}
    for evaluator in EVALUATORS:
        path = ledgers[evaluator]
        with path.open("r", encoding="utf-8", newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        observed = {
            str(row["case_id"]): {
                key: int(row[key]) for key in ("tp", "fp", "fn")
            }
            for row in rows
        }
        if observed != expected:
            missing = sorted(set(expected) - set(observed))
            wrong = sorted(
                case_id
                for case_id in set(expected) & set(observed)
                if expected[case_id] != observed[case_id]
            )
            raise SemanticValidationError(
                f"{evaluator} semantic mismatch; missing={missing}; wrong={wrong}"
            )
        evaluator_results[evaluator] = {
            "status": "valid",
            "ledger_sha256": sha256_file(path),
            "case_count": len(rows),
        }
    return {
        "schema_version": 1,
        "contract": suite["contract"],
        "status": "valid",
        "suite_sha256": sha256_file(cases_path),
        "evaluators": evaluator_results,
    }


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cases", required=True, type=Path)
    for evaluator in EVALUATORS:
        parser.add_argument(f"--{evaluator}-ledger", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        payload = validate(
            args.cases,
            {name: getattr(args, f"{name}_ledger") for name in EVALUATORS},
        )
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(
            json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
    except (OSError, KeyError, TypeError, ValueError, SemanticValidationError) as exc:
        print(f"validate_evaluator_semantics: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
