#!/usr/bin/env python3
"""CLI for calculating one tuple-bound unweighted consensus score."""

from __future__ import annotations

import argparse
import json
import math
import re
import sys
from copy import deepcopy
from pathlib import Path
from typing import Any

from pgbench_metrics import (
    DEFAULT_METRIC_DICTIONARY_PATH,
    DEFAULT_METRICS_SCHEMA_PATH,
    MetricContractError,
    build_score_payload_from_metrics,
)
from pgbench_scoring import (
    ScoreInputError,
    calculate_pgbench_score,
    load_score_profile,
)

_SHA256_RE = re.compile(r"^[0-9a-f]{64}$")
AUDIT_SCHEMA_VERSION = "pgbench.provenance_audit.v1"


def _invalid_json_constant(value: str) -> None:
    raise ScoreInputError(f"non-standard JSON numeric constant is forbidden: {value}")


def _load_json(path: Path, label: str) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        data = json.load(handle, parse_constant=_invalid_json_constant)
    if not isinstance(data, dict):
        raise ScoreInputError(f"{label} must be a JSON object")
    return data


def _atomic_json_dump(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp")
    with temporary.open("w", encoding="utf-8") as handle:
        json.dump(
            payload,
            handle,
            allow_nan=False,
            indent=2,
            sort_keys=True,
            ensure_ascii=False,
        )
        handle.write("\n")
    temporary.replace(path)


def _audit_unit_interval(value: Any, field: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ScoreInputError(f"provenance audit {field} must be numeric")
    numeric = float(value)
    if not math.isfinite(numeric) or numeric < 0 or numeric > 1:
        raise ScoreInputError(f"provenance audit {field} must be finite in [0, 1]")
    return numeric


def _apply_provenance_audit(
    payload: dict[str, Any], audit: dict[str, Any]
) -> dict[str, Any]:
    required = (
        "audit_schema_version",
        "status",
        "manifest_completeness",
        "hash_lineage_complete",
        "environment_complete",
        "run_context_complete",
        "core_provenance_valid",
    )
    missing = [key for key in required if key not in audit]
    if missing:
        raise ScoreInputError("provenance audit is missing: " + ", ".join(missing))
    if audit["audit_schema_version"] != AUDIT_SCHEMA_VERSION:
        raise ScoreInputError("unsupported provenance audit schema version")
    manifest_completeness = _audit_unit_interval(
        audit["manifest_completeness"], "manifest_completeness"
    )
    boolean_fields = (
        "hash_lineage_complete",
        "environment_complete",
        "run_context_complete",
        "core_provenance_valid",
    )
    for field in boolean_fields:
        if not isinstance(audit[field], bool):
            raise ScoreInputError(f"provenance audit {field} must be boolean")
    expected_status = (
        "invalid"
        if not audit["core_provenance_valid"]
        else (
            "valid"
            if manifest_completeness == 1.0
            and audit["hash_lineage_complete"]
            and audit["environment_complete"]
            and audit["run_context_complete"]
            else "provisional"
        )
    )
    if audit["status"] != expected_status:
        raise ScoreInputError(
            "provenance audit status is inconsistent with its scoring gates"
        )
    merged = deepcopy(payload)
    merged["traceability"] = {
        "manifest_completeness": manifest_completeness,
        "hash_lineage_complete": audit["hash_lineage_complete"],
        "environment_complete": audit["environment_complete"],
        "run_context_complete": audit["run_context_complete"],
        "core_provenance_valid": audit["core_provenance_valid"],
    }
    return merged


def _validate_run_context(
    run_context: dict[str, Any], score_profile: dict[str, Any]
) -> None:
    if run_context.get("schema_version") != 1:
        raise ScoreInputError("run context schema_version must be 1")
    actual_hash = run_context.get("score_profile_sha256")
    if not isinstance(actual_hash, str) or not _SHA256_RE.fullmatch(actual_hash):
        raise ScoreInputError(
            "run context score_profile_sha256 must be a lowercase SHA-256"
        )
    if actual_hash != score_profile.get("_sha256"):
        raise ScoreInputError(
            "selected score profile hash does not match the frozen run context"
        )


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Calculate one tuple-bound three-evaluator consensus score."
    )
    parser.add_argument(
        "--metrics",
        "--input",
        dest="metrics",
        required=True,
        type=Path,
        help="Standard metrics document; --input is a compatibility alias",
    )
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--score-profile", required=True, type=Path)
    parser.add_argument(
        "--metric-dictionary",
        type=Path,
        default=DEFAULT_METRIC_DICTIONARY_PATH,
    )
    parser.add_argument(
        "--metrics-schema",
        type=Path,
        default=DEFAULT_METRICS_SCHEMA_PATH,
    )
    parser.add_argument("--provenance-audit", required=True, type=Path)
    parser.add_argument("--run-context", required=True, type=Path)
    parser.add_argument(
        "--evaluation-mode",
        required=True,
        choices=("formal", "synthetic_smoke"),
    )
    parser.add_argument("--expected-run-id", required=True)
    parser.add_argument("--expected-sample-id", required=True)
    parser.add_argument("--expected-tool-id", required=True)
    parser.add_argument(
        "--expected-official-score-mode",
        required=True,
        choices=("caller_only_shared_alignment", "end_to_end_from_reads"),
    )
    parser.add_argument("--expected-primary-truth-profile", required=True)
    return parser.parse_args(argv)


def _expected_tuple(args: argparse.Namespace) -> dict[str, str]:
    return {
        "run_id": args.expected_run_id,
        "sample": args.expected_sample_id,
        "tool": args.expected_tool_id,
        "official_score_mode": args.expected_official_score_mode,
        "primary_truth_profile": args.expected_primary_truth_profile,
    }


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        score_profile = load_score_profile(args.score_profile)
        run_context = _load_json(args.run_context, "run context")
        _validate_run_context(run_context, score_profile)
        expected_tuple = _expected_tuple(args)
        payload = build_score_payload_from_metrics(
            _load_json(args.metrics, "metrics document"),
            expected_tuple=expected_tuple,
            score_profile=score_profile,
            metrics_schema_path=args.metrics_schema,
            metric_dictionary_path=args.metric_dictionary,
        )
        payload = _apply_provenance_audit(
            payload, _load_json(args.provenance_audit, "provenance audit")
        )
        result = calculate_pgbench_score(
            payload,
            expected_tuple=expected_tuple,
            evaluation_mode=args.evaluation_mode,
            score_profile_path=args.score_profile,
        )
        _atomic_json_dump(args.output, result.to_dict())
    except (
        OSError,
        json.JSONDecodeError,
        MetricContractError,
        ScoreInputError,
    ) as exc:
        print(f"score_tools: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
