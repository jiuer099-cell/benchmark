#!/usr/bin/env python3
"""Freeze Git, hardware, runtime, profile, and random-seed context."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import subprocess
import sys
from pathlib import Path
from typing import Any


class RunContextError(RuntimeError):
    """Raised when repository state cannot be captured."""


def _git(repo_root: Path, *args: str) -> bytes:
    try:
        completed = subprocess.run(
            ["git", *args],
            cwd=repo_root,
            check=True,
            capture_output=True,
        )
    except (OSError, subprocess.CalledProcessError) as exc:
        raise RunContextError(f"git {' '.join(args)} failed: {exc}") from exc
    return completed.stdout


def _sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def capture_run_context(
    *,
    repo_root: Path,
    score_profile: Path,
    evaluator_profile: Path | None = None,
    random_seed: int,
    snakemake_version: str,
    execution_profile: str,
) -> dict[str, Any]:
    git_head = _git(repo_root, "rev-parse", "HEAD").strip().decode("ascii")
    status_bytes = _git(
        repo_root,
        "status",
        "--porcelain=v1",
        "-z",
        "--untracked-files=all",
    )
    diff_bytes = _git(repo_root, "diff", "--binary", "--no-ext-diff", "HEAD", "--")
    tracked_tree_bytes = _git(repo_root, "ls-files", "-s", "-z")
    untracked_paths = _git(
        repo_root,
        "ls-files",
        "--others",
        "--exclude-standard",
        "-z",
    ).split(b"\0")
    untracked_records: list[bytes] = []
    for raw_path in sorted(item for item in untracked_paths if item):
        path = repo_root / os.fsdecode(raw_path)
        if not path.is_file() or path.is_symlink():
            raise RunContextError(
                f"untracked implementation path must be a regular file: {path}"
            )
        untracked_records.append(
            b"\0" + raw_path + b"\0" + _sha256(path.read_bytes()).encode("ascii")
        )
    untracked_bytes = b"".join(untracked_records)

    status_sha = _sha256(status_bytes)
    diff_sha = _sha256(diff_bytes)
    tracked_sha = _sha256(tracked_tree_bytes)
    untracked_sha = _sha256(untracked_bytes)
    repo_state_sha = _sha256(
        b"\0".join(
            [
                git_head.encode("ascii"),
                status_sha.encode("ascii"),
                diff_sha.encode("ascii"),
                tracked_sha.encode("ascii"),
                untracked_sha.encode("ascii"),
            ]
        )
    )

    hardware = {
        "machine": platform.machine(),
        "platform": platform.platform(),
        "processor": platform.processor(),
        "python": platform.python_version(),
        "cpu_count": os.cpu_count(),
    }
    hardware_bytes = json.dumps(hardware, sort_keys=True, separators=(",", ":")).encode(
        "utf-8"
    )
    payload = {
        "schema_version": 1,
        "git_head": git_head,
        "git_dirty": bool(status_bytes),
        "git_status_sha256": status_sha,
        "git_diff_sha256": diff_sha,
        "tracked_tree_sha256": tracked_sha,
        "untracked_files_sha256": untracked_sha,
        "repo_state_sha256": repo_state_sha,
        "score_profile_path": str(score_profile),
        "score_profile_sha256": _sha256(score_profile.read_bytes()),
        "random_seed": random_seed,
        "snakemake_version": snakemake_version,
        "execution_profile": execution_profile,
        "hardware": hardware,
        "hardware_fingerprint_sha256": _sha256(hardware_bytes),
    }
    if evaluator_profile is not None:
        payload["evaluator_profile_path"] = str(evaluator_profile)
        payload["evaluator_profile_sha256"] = _sha256(
            evaluator_profile.read_bytes()
        )
    return payload


def _atomic_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp")
    temporary.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    temporary.replace(path)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Snapshot PGBench run context.")
    parser.add_argument("--repo-root", type=Path, default=Path.cwd())
    parser.add_argument("--score-profile", required=True, type=Path)
    parser.add_argument("--evaluator-profile", type=Path)
    parser.add_argument("--random-seed", required=True, type=int)
    parser.add_argument("--snakemake-version", required=True)
    parser.add_argument("--execution-profile", required=True)
    parser.add_argument("--output", required=True, type=Path)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        payload = capture_run_context(
            repo_root=args.repo_root.resolve(),
            score_profile=args.score_profile,
            evaluator_profile=args.evaluator_profile,
            random_seed=args.random_seed,
            snakemake_version=args.snakemake_version,
            execution_profile=args.execution_profile,
        )
        _atomic_json(args.output, payload)
    except (OSError, RunContextError) as exc:
        print(f"snapshot_run_context: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
