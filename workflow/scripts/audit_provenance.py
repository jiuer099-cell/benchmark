"""Audit rule manifests and emit the traceability inputs used by scoring."""

from __future__ import annotations

import argparse
from collections.abc import Sequence
from pathlib import Path
from typing import Any

try:
    from pgbench_provenance import (
        ProvenanceError,
        atomic_write_json,
        audit_manifests,
        discover_manifest_paths,
        load_json,
        load_manifests,
    )
except ModuleNotFoundError:  # pragma: no cover - package-style invocation
    from .pgbench_provenance import (
        ProvenanceError,
        atomic_write_json,
        audit_manifests,
        discover_manifest_paths,
        load_json,
        load_manifests,
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", action="append", default=[], type=Path)
    parser.add_argument("--manifest-dir", action="append", default=[], type=Path)
    parser.add_argument(
        "--expected-jobs",
        type=Path,
        help="JSON list, or an object containing a jobs list",
    )
    parser.add_argument(
        "--expected-job",
        action="append",
        default=[],
        metavar="RULE=JOB_KEY",
        help=(
            "Expected executed job; repeat as needed. Jobs supplied this way "
            "are treated as core score dependencies."
        ),
    )
    parser.add_argument("--target-manifest-id", action="append", default=[])
    parser.add_argument("--core-rule", action="append", default=[])
    parser.add_argument("--workspace-root", type=Path, default=Path.cwd())
    parser.add_argument(
        "--skip-companion-check",
        action="store_true",
        help="Do not require conventional rule log and benchmark files",
    )
    parser.add_argument(
        "--verify-paths",
        action="store_true",
        help="Re-hash current declared input and output paths",
    )
    parser.add_argument("--output", required=True, type=Path)
    return parser


def _load_expected_jobs(path: Path | None) -> list[dict[str, Any]] | None:
    if path is None:
        return None
    payload = load_json(path)
    if isinstance(payload, dict):
        payload = payload.get("jobs")
    if not isinstance(payload, list) or not all(
        isinstance(item, dict) for item in payload
    ):
        raise ProvenanceError("expected-jobs JSON must contain a list of objects")
    return payload


def _parse_expected_jobs(
    path: Path | None, inline_values: Sequence[str]
) -> list[dict[str, Any]] | None:
    from_file = _load_expected_jobs(path)
    if from_file is not None and inline_values:
        raise ProvenanceError(
            "--expected-jobs and --expected-job are mutually exclusive"
        )
    if from_file is not None:
        return from_file
    if not inline_values:
        return None
    jobs: list[dict[str, Any]] = []
    for value in inline_values:
        rule_name, separator, job_key = value.partition("=")
        if not separator or not rule_name or not job_key:
            raise ProvenanceError("--expected-job must use RULE=JOB_KEY syntax")
        jobs.append(
            {
                "rule_name": rule_name,
                "job_key": job_key,
                "core": True,
            }
        )
    return jobs


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        paths = discover_manifest_paths(args.manifest, args.manifest_dir)
        if not paths:
            raise ProvenanceError("at least one manifest is required")
        manifests = load_manifests(paths)
        audit = audit_manifests(
            manifests,
            expected_jobs=_parse_expected_jobs(args.expected_jobs, args.expected_job),
            target_manifest_ids=args.target_manifest_id or None,
            core_rule_patterns=args.core_rule or None,
            workspace_root=args.workspace_root,
            require_companions=not args.skip_companion_check,
            verify_paths=args.verify_paths,
        )
        atomic_write_json(args.output, audit)
    except ProvenanceError as exc:
        raise SystemExit(f"audit_provenance: {exc}") from exc
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
