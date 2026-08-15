#!/usr/bin/env python3
"""Repair a completed GraphTyper2 projection without rerunning genotyping.

This utility is intentionally narrow: it reprojects immutable GraphTyper2
shard VCFs, refreshes the successful attempt's public output, and updates the
rule manifest hashes while preserving backups and a recovery record.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import os
import shutil
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def atomic_json(path: Path, payload: Any) -> None:
    temporary = path.with_name(f".{path.name}.projection-recovery.tmp")
    temporary.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    os.replace(temporary, path)


def backup_once(path: Path) -> Path:
    backup = path.with_name(path.name + ".before-old-variant-id-fix")
    if not backup.exists():
        shutil.copy2(path, backup)
    return backup


def load_runner(path: Path):
    specification = importlib.util.spec_from_file_location(
        "pgbench_graphtyper2_recovery_runner", path
    )
    if specification is None or specification.loader is None:
        raise RuntimeError(f"cannot load GraphTyper2 runner: {path}")
    module = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(module)
    return module


def genotype_counts(path: Path) -> dict[str, int]:
    counts: Counter[str] = Counter()
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            counts[fields[9]] += 1
    return dict(sorted(counts.items()))


def refresh_path_fingerprint(manifest: dict[str, Any], declared: str) -> None:
    path = Path(declared)
    stat = path.stat()
    manifest["input_sha256"][declared] = sha256_file(path)
    manifest["input_size"][declared] = stat.st_size
    manifest["input_mtime"][declared] = stat.st_mtime_ns


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--candidate", required=True, type=Path)
    parser.add_argument("--graphtyper-root", required=True, type=Path)
    parser.add_argument("--runner", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--attempt-record", required=True, type=Path)
    parser.add_argument("--resolved-inputs", required=True, type=Path)
    parser.add_argument("--rule-manifest", required=True, type=Path)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--recovery-record", required=True, type=Path)
    args = parser.parse_args()

    runner = load_runner(args.runner)
    sources = runner.generated_vcfs(args.graphtyper_root)
    temporary_output = args.output.with_name(args.output.name + ".recovery.tmp")
    runner.project_calls(args.candidate, sources, temporary_output, args.sample)

    backup_output = backup_once(args.output)
    backup_attempt = backup_once(args.attempt_record)
    backup_manifest = backup_once(args.rule_manifest)
    os.replace(temporary_output, args.output)

    attempt = json.loads(args.attempt_record.read_text(encoding="utf-8"))
    attempt_output = Path(attempt["attempt_output_vcf"])
    attempt_output.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(args.output, attempt_output)
    output_stat = args.output.stat()
    output_sha = sha256_file(args.output)
    runner_sha = sha256_file(args.runner)
    recovered_at = datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")

    recovery = {
        "schema_version": 1,
        "reason": "GraphTyper2 records identify candidates via INFO/OLD_VARIANT_ID",
        "recovered_at": recovered_at,
        "candidate_sha256": sha256_file(args.candidate),
        "runner_sha256": runner_sha,
        "source_vcf_count": len(sources),
        "previous_output_sha256": sha256_file(backup_output),
        "recovered_output_sha256": output_sha,
        "genotype_counts": genotype_counts(args.output),
        "backups": {
            "output": str(backup_output),
            "attempt_record": str(backup_attempt),
            "rule_manifest": str(backup_manifest),
        },
    }
    args.recovery_record.parent.mkdir(parents=True, exist_ok=True)
    atomic_json(args.recovery_record, recovery)

    attempt["output"].update(
        {
            "mtime_ns": output_stat.st_mtime_ns,
            "sha256": output_sha,
            "size_bytes": output_stat.st_size,
        }
    )
    attempt["runner_sha256"] = runner_sha
    attempt["postprocessing_recovery"] = {
        "record": str(args.recovery_record),
        "record_sha256": sha256_file(args.recovery_record),
        "recovered_at": recovered_at,
    }
    atomic_json(args.attempt_record, attempt)

    # Refresh this declared output's mtime without changing its content/hash so
    # Snakemake does not rerun the expensive tool rule after the runner update.
    resolved = json.loads(args.resolved_inputs.read_text(encoding="utf-8"))
    atomic_json(args.resolved_inputs, resolved)

    manifest = json.loads(args.rule_manifest.read_text(encoding="utf-8"))
    runner_resolved = args.runner.resolve()
    for declared in manifest["input_paths"]:
        path = Path(declared)
        try:
            same_runner = path.resolve() == runner_resolved
        except OSError:
            same_runner = False
        if same_runner:
            refresh_path_fingerprint(manifest, declared)
    manifest["output_sha256"][str(args.output)] = output_sha
    manifest["output_sha256"][str(args.attempt_record)] = sha256_file(
        args.attempt_record
    )
    manifest["output_sha256"][str(args.resolved_inputs)] = sha256_file(
        args.resolved_inputs
    )
    manifest["params"]["projection_recovery"] = {
        "record": str(args.recovery_record),
        "record_sha256": sha256_file(args.recovery_record),
        "recovered_at": recovered_at,
    }
    atomic_json(args.rule_manifest, manifest)

    print(json.dumps(recovery, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
