#!/usr/bin/env python3
"""Import a frozen external caller output for an isolated scoring replay.

This is intentionally a Core-owned transport rule.  It never invokes an
adapter or caller: it verifies a successful source execution manifest, copies
the exact three immutable caller artifacts into the new run, and leaves the
source manifest itself as a hashed input of the new import manifest.
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
import sys
import tempfile
from pathlib import Path
from typing import Any, Sequence

try:
    from pgbench_provenance import atomic_write_json, sha256_file
except ModuleNotFoundError:  # pragma: no cover - package-style invocation
    from .pgbench_provenance import atomic_write_json, sha256_file


class FrozenNativeImportError(ValueError):
    """Raised when a claimed frozen caller output cannot be authenticated."""


def _load_mapping(path: Path, label: str) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise FrozenNativeImportError(f"cannot load {label}: {path}: {error}") from error
    if not isinstance(value, dict):
        raise FrozenNativeImportError(f"{label} must contain a JSON object: {path}")
    return value


def _assert_source_output(
    manifest: dict[str, Any], source: Path, suffix: str, label: str
) -> str:
    entries = manifest.get("output_sha256")
    if not isinstance(entries, dict):
        raise FrozenNativeImportError("source manifest has no output_sha256 mapping")
    expected = [
        digest for path, digest in entries.items()
        if isinstance(path, str) and path.replace("\\", "/").endswith(suffix)
        and isinstance(digest, str) and len(digest) == 64
    ]
    if len(expected) != 1:
        raise FrozenNativeImportError(
            f"source manifest has no unambiguous {label} output: {suffix}"
        )
    if not source.is_file() or source.is_symlink():
        raise FrozenNativeImportError(f"frozen {label} is not a regular file: {source}")
    actual = sha256_file(source)
    if actual != expected[0]:
        raise FrozenNativeImportError(
            f"frozen {label} SHA-256 does not match source manifest"
        )
    return actual


def _copy_immutable(source: Path, destination: Path) -> None:
    if destination.exists() or destination.is_symlink():
        raise FrozenNativeImportError(f"refusing to overwrite replay artifact: {destination}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        dir=destination.parent, prefix=f".{destination.name}.", suffix=".tmp"
    )
    try:
        with os.fdopen(descriptor, "wb") as output, source.open("rb") as input_handle:
            shutil.copyfileobj(input_handle, output, length=1024 * 1024)
            output.flush()
            os.fsync(output.fileno())
        os.replace(temporary_name, destination)
    except BaseException:
        try:
            os.unlink(temporary_name)
        except FileNotFoundError:
            pass
        raise


def import_frozen_native_output(args: argparse.Namespace) -> None:
    source_manifest = Path(args.source_tool_manifest).resolve(strict=True)
    source_vcf = Path(args.source_vcf).resolve(strict=True)
    source_resolved = Path(args.source_resolved_inputs).resolve(strict=True)
    source_attempt = Path(args.source_attempt_record).resolve(strict=True)
    manifest = _load_mapping(source_manifest, "source tool manifest")
    expected_rule = f"tool__{args.tool_id}__execute"
    if manifest.get("status") != "success" or manifest.get("rule_name") != expected_rule:
        raise FrozenNativeImportError(
            "source manifest must be a successful native adapter execution"
        )
    source_run_id = manifest.get("run_id")
    if not isinstance(source_run_id, str) or not source_run_id:
        raise FrozenNativeImportError("source manifest has no run_id")
    vcf_sha = _assert_source_output(manifest, source_vcf, "/raw/calls.vcf", "VCF")
    resolved_sha = _assert_source_output(
        manifest, source_resolved, "/meta/resolved_inputs.json", "resolved-inputs"
    )
    attempt_sha = _assert_source_output(
        manifest, source_attempt, "/meta/attempt.json", "attempt record"
    )
    # Parse the JSON inputs before writing anything, so malformed provenance
    # cannot be carried into an otherwise clean replay root.
    _load_mapping(source_resolved, "source resolved-inputs")
    _load_mapping(source_attempt, "source attempt record")

    outputs = (
        (source_vcf, Path(args.output_vcf), vcf_sha),
        (source_resolved, Path(args.output_resolved_inputs), resolved_sha),
        (source_attempt, Path(args.output_attempt_record), attempt_sha),
    )
    for source, destination, expected_sha in outputs:
        _copy_immutable(source, destination)
        if sha256_file(destination) != expected_sha:
            raise FrozenNativeImportError("copy verification failed for frozen native artifact")
    atomic_write_json(
        Path(args.import_audit),
        {
            "contract": "pgbench_frozen_native_output_import_v1",
            "status": "valid",
            "tool_id": args.tool_id,
            "source_run_id": source_run_id,
            "source_tool_manifest": source_manifest.as_posix(),
            "source_tool_manifest_sha256": sha256_file(source_manifest),
            "source_tool_manifest_id": manifest.get("manifest_id"),
            "source_tool_git_head": manifest.get("git_head"),
            "artifacts": {
                "vcf": {"source": source_vcf.as_posix(), "sha256": vcf_sha},
                "resolved_inputs": {"source": source_resolved.as_posix(), "sha256": resolved_sha},
                "attempt_record": {"source": source_attempt.as_posix(), "sha256": attempt_sha},
            },
        },
    )


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--tool-id", required=True)
    result.add_argument("--source-vcf", required=True)
    result.add_argument("--source-resolved-inputs", required=True)
    result.add_argument("--source-attempt-record", required=True)
    result.add_argument("--source-tool-manifest", required=True)
    result.add_argument("--output-vcf", required=True)
    result.add_argument("--output-resolved-inputs", required=True)
    result.add_argument("--output-attempt-record", required=True)
    result.add_argument("--import-audit", required=True)
    return result


def main(argv: Sequence[str] | None = None) -> int:
    try:
        import_frozen_native_output(parser().parse_args(argv))
    except (FrozenNativeImportError, OSError) as error:
        print(f"import_frozen_native_output: {error}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
