#!/usr/bin/env python3
"""Execute one trusted workflow command and write its formal rule manifest.

This wrapper is for benchmark-owned rules only. External/user-supplied tools
remain behind ``pgbench_exec.py`` and its separate isolation contract.
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import uuid
from collections.abc import Sequence
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

try:
    from pgbench_provenance import (
        MANIFEST_SCHEMA_VERSION,
        ProvenanceError,
        fingerprint_paths,
        load_json,
        sha256_file,
        write_rule_manifest,
    )
    from snapshot_run_context import RunContextError, capture_run_context
except ModuleNotFoundError:  # pragma: no cover - package-style invocation
    from .pgbench_provenance import (
        MANIFEST_SCHEMA_VERSION,
        ProvenanceError,
        fingerprint_paths,
        load_json,
        sha256_file,
        write_rule_manifest,
    )
    from .snapshot_run_context import RunContextError, capture_run_context


class RuleExecutionError(RuntimeError):
    """Raised when the trusted rule execution contract is invalid."""


def _utc_now() -> str:
    return datetime.now(UTC).isoformat().replace("+00:00", "Z")


def _parse_key_value(values: Sequence[str], label: str) -> dict[str, Any]:
    parsed: dict[str, Any] = {}
    for value in values:
        key, separator, raw = value.partition("=")
        if not separator or not key:
            raise RuleExecutionError(
                f"{label} must use KEY=VALUE syntax, found {value!r}"
            )
        if key in parsed:
            raise RuleExecutionError(f"duplicate {label} key: {key}")
        try:
            parsed[key] = json.loads(raw)
        except json.JSONDecodeError:
            parsed[key] = raw
    return parsed


def _unique_paths(values: Sequence[Path]) -> list[str]:
    result: list[str] = []
    seen: set[str] = set()
    for value in values:
        text = str(value)
        if text not in seen:
            result.append(text)
            seen.add(text)
    return result


def _load_mapping(path: Path, label: str) -> dict[str, Any]:
    payload = load_json(path)
    if not isinstance(payload, dict):
        raise RuleExecutionError(f"{label} must be a JSON object: {path}")
    return payload


def _resolve_path(repo_root: Path, path: Path) -> Path:
    return path if path.is_absolute() else repo_root / path


def _load_or_capture_context(
    *,
    repo_root: Path,
    run_context: Path | None,
    score_profile: Path,
    random_seed: int,
    snakemake_version: str,
    execution_profile: str,
) -> dict[str, Any]:
    if run_context is not None:
        resolved = _resolve_path(repo_root, run_context)
        if resolved.is_file():
            return _load_mapping(resolved, "run context")
    return capture_run_context(
        repo_root=repo_root,
        score_profile=_resolve_path(repo_root, score_profile),
        random_seed=random_seed,
        snakemake_version=snakemake_version,
        execution_profile=execution_profile,
    )


def _upstream_manifest_ids(repo_root: Path, paths: Sequence[Path]) -> list[str]:
    manifest_ids: list[str] = []
    seen: set[str] = set()
    for path in paths:
        payload = _load_mapping(
            _resolve_path(repo_root, path), "upstream rule manifest"
        )
        manifest_id = payload.get("manifest_id")
        if not isinstance(manifest_id, str) or not manifest_id:
            raise RuleExecutionError(
                f"upstream rule manifest has no manifest_id: {path}"
            )
        if manifest_id not in seen:
            manifest_ids.append(manifest_id)
            seen.add(manifest_id)
    return manifest_ids


def _optional_hash(repo_root: Path, path: Path | None) -> str | None:
    if path is None:
        return None
    candidate = _resolve_path(repo_root, path)
    return sha256_file(candidate)


def _path_is_present(repo_root: Path, path: Path | None) -> bool:
    if path is None:
        return False
    candidate = _resolve_path(repo_root, path)
    return candidate.exists() or candidate.is_symlink()


def _path_is_declared(
    repo_root: Path,
    path: Path | None,
    declared_paths: Sequence[str],
) -> bool:
    if path is None:
        return False
    candidate = _resolve_path(repo_root, path).absolute()
    return any(
        _resolve_path(repo_root, Path(value)).absolute() == candidate
        for value in declared_paths
    )


def _fingerprint_snapshot(
    paths: Sequence[str],
    *,
    repo_root: Path,
) -> tuple[dict[str, str], dict[str, int], dict[str, int]]:
    return fingerprint_paths(paths, base_dir=repo_root)


def _verify_fingerprint_snapshot(
    paths: Sequence[str],
    *,
    repo_root: Path,
    expected_hashes: dict[str, str],
    expected_sizes: dict[str, int],
    expected_mtimes: dict[str, int],
) -> None:
    try:
        hashes, sizes, mtimes = _fingerprint_snapshot(
            paths,
            repo_root=repo_root,
        )
    except ProvenanceError as exc:
        raise RuleExecutionError(
            f"frozen input/code path became unavailable: {exc}"
        ) from exc
    changed = sorted(
        path
        for path in paths
        if hashes.get(path) != expected_hashes.get(path)
        or sizes.get(path) != expected_sizes.get(path)
        or mtimes.get(path) != expected_mtimes.get(path)
    )
    if changed:
        raise RuleExecutionError(
            "frozen input/code changed during rule execution: " + ", ".join(changed)
        )


def _archive_existing_manifest(
    *,
    repo_root: Path,
    manifest_output: Path,
    attempt_id: str,
) -> Path | None:
    """Move a prior canonical manifest to a recoverable per-attempt archive."""

    manifest = _resolve_path(repo_root, manifest_output)
    if manifest.is_symlink():
        raise RuleExecutionError("existing manifest output must not be a symlink")
    if not manifest.exists():
        return None
    if not manifest.is_file():
        raise RuleExecutionError("existing manifest output is not a regular file")

    history_root = manifest.parent / ".previous"
    if history_root.is_symlink():
        raise RuleExecutionError("manifest archive directory must not be a symlink")
    if history_root.exists() and not history_root.is_dir():
        raise RuleExecutionError("manifest archive path exists but is not a directory")
    history_root.mkdir(mode=0o700, exist_ok=True)
    destination_dir = history_root / attempt_id
    if destination_dir.exists() or destination_dir.is_symlink():
        raise RuleExecutionError(
            f"manifest archive already exists; refusing overwrite: {destination_dir}"
        )
    destination_dir.mkdir(mode=0o700)
    destination = destination_dir / manifest.name
    if destination.exists() or destination.is_symlink():
        raise RuleExecutionError(
            f"manifest archive destination already exists: {destination}"
        )
    try:
        os.replace(manifest, destination)
    except OSError as exc:
        try:
            destination_dir.rmdir()
        except OSError:
            pass
        raise RuleExecutionError(
            f"cannot archive existing manifest safely: {exc}"
        ) from exc
    return destination


def _reject_command_manifest_output(
    *,
    repo_root: Path,
    manifest_output: Path,
) -> None:
    """Ensure only this wrapper can publish the canonical success manifest."""

    manifest = _resolve_path(repo_root, manifest_output)
    if not manifest.exists() and not manifest.is_symlink():
        return
    if manifest.is_symlink() or manifest.is_file():
        manifest.unlink()
        raise RuleExecutionError(
            "rule command created the reserved manifest output; removed it"
        )
    raise RuleExecutionError(
        "rule command created a non-file at the reserved manifest output"
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rule-name", required=True)
    parser.add_argument("--job-key", required=True)
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--module-or-tool-id", required=True)
    parser.add_argument("--snakefile-path", required=True)
    parser.add_argument("--rule-source-path", required=True, type=Path)
    parser.add_argument("--script-or-wrapper-path", required=True, type=Path)
    parser.add_argument("--config-snapshot", required=True, type=Path)
    parser.add_argument("--score-profile", required=True, type=Path)
    parser.add_argument("--run-context", type=Path)
    parser.add_argument("--pangenome-manifest", type=Path)
    parser.add_argument("--reference", type=Path)
    parser.add_argument("--truth-profile", required=True)
    parser.add_argument("--conda-lock", type=Path)
    parser.add_argument("--container-uri")
    parser.add_argument("--container-digest")
    parser.add_argument("--snakemake-version", required=True)
    parser.add_argument("--execution-profile", required=True)
    parser.add_argument("--random-seed", required=True, type=int)
    parser.add_argument("--threads", required=True, type=int)
    parser.add_argument("--resource", action="append", default=[])
    parser.add_argument("--param", action="append", default=[])
    parser.add_argument("--wildcard", action="append", default=[])
    parser.add_argument("--input", action="append", default=[], type=Path)
    parser.add_argument("--output", action="append", default=[], type=Path)
    parser.add_argument("--upstream-manifest", action="append", default=[], type=Path)
    parser.add_argument("--manifest-output", required=True, type=Path)
    parser.add_argument("--repo-root", type=Path, default=Path.cwd())
    parser.add_argument(
        "command",
        nargs=argparse.REMAINDER,
        help="Command argv after a literal -- separator",
    )
    return parser


def _normalize_command(raw: Sequence[str]) -> list[str]:
    command = list(raw)
    if command and command[0] == "--":
        command = command[1:]
    if not command:
        raise RuleExecutionError("a command must follow the -- separator")
    return command


def execute_and_manifest(args: argparse.Namespace) -> int:
    """Execute the parsed command and write a manifest only after success."""

    if args.threads < 1:
        raise RuleExecutionError("--threads must be at least 1")
    repo_root = args.repo_root.resolve(strict=True)
    command = _normalize_command(args.command)
    upstream_paths = list(args.upstream_manifest)
    declared_inputs = _unique_paths([*args.input, *upstream_paths])
    declared_outputs = _unique_paths(args.output)

    score_profile = _resolve_path(repo_root, args.score_profile)
    pangenome_is_output = _path_is_declared(
        repo_root,
        args.pangenome_manifest,
        declared_outputs,
    )
    pangenome_was_frozen = not pangenome_is_output and _path_is_present(
        repo_root, args.pangenome_manifest
    )
    run_context_is_output = _path_is_declared(
        repo_root,
        args.run_context,
        declared_outputs,
    )

    guard_values = [
        *(Path(path) for path in declared_inputs),
        args.rule_source_path,
        args.script_or_wrapper_path,
        Path(__file__).resolve(),
        args.config_snapshot,
        args.score_profile,
    ]
    if args.reference is not None:
        guard_values.append(args.reference)
    if args.conda_lock is not None:
        guard_values.append(args.conda_lock)
    if pangenome_was_frozen and args.pangenome_manifest is not None:
        guard_values.append(args.pangenome_manifest)
    if not run_context_is_output and _path_is_present(repo_root, args.run_context):
        assert args.run_context is not None
        guard_values.append(args.run_context)
    guard_paths = _unique_paths(guard_values)

    guard_hashes, guard_sizes, guard_mtimes = _fingerprint_snapshot(
        guard_paths,
        repo_root=repo_root,
    )
    input_hashes = {path: guard_hashes[path] for path in declared_inputs}
    input_sizes = {path: guard_sizes[path] for path in declared_inputs}
    input_mtimes = {path: guard_mtimes[path] for path in declared_inputs}
    rule_source_sha = guard_hashes[str(args.rule_source_path)]
    script_or_wrapper_sha = guard_hashes[str(args.script_or_wrapper_path)]
    config_snapshot_sha = guard_hashes[str(args.config_snapshot)]
    score_profile_sha = guard_hashes[str(args.score_profile)]
    reference_sha = (
        guard_hashes[str(args.reference)] if args.reference is not None else None
    )
    conda_lock_sha = (
        guard_hashes[str(args.conda_lock)] if args.conda_lock is not None else None
    )
    pangenome_sha = (
        guard_hashes[str(args.pangenome_manifest)]
        if pangenome_was_frozen and args.pangenome_manifest is not None
        else None
    )
    upstream_manifest_ids = _upstream_manifest_ids(repo_root, upstream_paths)
    context = _load_or_capture_context(
        repo_root=repo_root,
        run_context=None if run_context_is_output else args.run_context,
        score_profile=score_profile,
        random_seed=args.random_seed,
        snakemake_version=args.snakemake_version,
        execution_profile=args.execution_profile,
    )
    wildcards = _parse_key_value(args.wildcard, "--wildcard")
    params = _parse_key_value(args.param, "--param")
    requested_resources = _parse_key_value(args.resource, "--resource")
    _verify_fingerprint_snapshot(
        guard_paths,
        repo_root=repo_root,
        expected_hashes=guard_hashes,
        expected_sizes=guard_sizes,
        expected_mtimes=guard_mtimes,
    )

    attempt_id = uuid.uuid4().hex
    previous_manifest_archive = _archive_existing_manifest(
        repo_root=repo_root,
        manifest_output=args.manifest_output,
        attempt_id=attempt_id,
    )
    started_at = _utc_now()
    completed = subprocess.run(
        command,
        cwd=repo_root,
        check=False,
        close_fds=True,
    )
    _reject_command_manifest_output(
        repo_root=repo_root,
        manifest_output=args.manifest_output,
    )
    _verify_fingerprint_snapshot(
        guard_paths,
        repo_root=repo_root,
        expected_hashes=guard_hashes,
        expected_sizes=guard_sizes,
        expected_mtimes=guard_mtimes,
    )
    if completed.returncode != 0:
        print(
            f"pgbench_rule_exec: {args.rule_name}/{args.job_key} "
            f"failed with exit code {completed.returncode}",
            file=sys.stderr,
        )
        return completed.returncode

    if args.pangenome_manifest is not None and not pangenome_was_frozen:
        pangenome_sha = _optional_hash(repo_root, args.pangenome_manifest)
    output_hashes, _, _ = _fingerprint_snapshot(
        declared_outputs,
        repo_root=repo_root,
    )
    not_applicable: dict[str, str] = {}
    if pangenome_sha is None:
        not_applicable["pangenome_manifest_sha256"] = (
            "The pangenome manifest is not available at this bootstrap step"
        )
    if reference_sha is None:
        not_applicable["reference_sha256"] = (
            "The reference is not available at this bootstrap step"
        )
    if conda_lock_sha is None:
        not_applicable["conda_lock_sha256"] = (
            "Phase 1 has an environment specification but no solved lock file"
        )
    if args.container_uri is None:
        not_applicable["container_uri"] = "No container was selected for this rule"
    if args.container_digest is None:
        not_applicable["container_digest"] = "No container was selected for this rule"

    payload: dict[str, Any] = {
        "manifest_schema_version": MANIFEST_SCHEMA_VERSION,
        "attempt_id": attempt_id,
        "run_id": args.run_id,
        "rule_name": args.rule_name,
        "job_key": args.job_key,
        "wildcards": wildcards,
        "module_or_tool_id": args.module_or_tool_id,
        "snakefile_path": args.snakefile_path,
        "rule_source_sha256": rule_source_sha,
        "script_or_wrapper_path": str(args.script_or_wrapper_path),
        "script_or_wrapper_sha256": script_or_wrapper_sha,
        "command": command,
        "params": params,
        "threads": args.threads,
        "requested_resources": requested_resources,
        "input_paths": declared_inputs,
        "input_sha256": input_hashes,
        "input_size": input_sizes,
        "input_mtime": input_mtimes,
        "output_paths": declared_outputs,
        "output_sha256": output_hashes,
        "config_snapshot_sha256": config_snapshot_sha,
        "score_profile_sha256": score_profile_sha,
        "pangenome_manifest_sha256": pangenome_sha,
        "reference_sha256": reference_sha,
        "truth_profile": args.truth_profile,
        "conda_lock_sha256": conda_lock_sha,
        "container_uri": args.container_uri,
        "container_digest": args.container_digest,
        "git_head": context.get("git_head"),
        "git_dirty": context.get("git_dirty"),
        "git_diff_sha256": context.get("git_diff_sha256"),
        "snakemake_version": context.get("snakemake_version", args.snakemake_version),
        "execution_profile": context.get("execution_profile", args.execution_profile),
        "hardware_fingerprint_sha256": context.get("hardware_fingerprint_sha256"),
        "random_seed": context.get("random_seed", args.random_seed),
        "upstream_manifest_ids": upstream_manifest_ids,
        "started_at": started_at,
        "finished_at": _utc_now(),
        "exit_code": 0,
        "status": "success",
        "not_applicable_reason": not_applicable,
    }
    if previous_manifest_archive is not None:
        payload["previous_manifest_archive"] = str(previous_manifest_archive)
    write_rule_manifest(
        payload,
        _resolve_path(repo_root, args.manifest_output),
        base_dir=repo_root,
        calculate_path_metadata=False,
        require_success=True,
    )
    return 0


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        return execute_and_manifest(args)
    except (
        OSError,
        ProvenanceError,
        RuleExecutionError,
        RunContextError,
        subprocess.SubprocessError,
    ) as exc:
        print(f"pgbench_rule_exec: {exc}", file=sys.stderr)
        return 3


if __name__ == "__main__":
    raise SystemExit(main())
