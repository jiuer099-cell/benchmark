#!/usr/bin/env python3
"""Execute reviewed adapter gates and verify content-bound eligibility receipts.

A receipt proves successful recorded commands and unchanged evidence, not that
arbitrary tool flags are biologically fair. trusted_reviewed still requires
human review of the official release, complete mini-panel assertions and argv.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import subprocess
import sys
import time
import tempfile
from pathlib import Path
from typing import Any

import jsonschema
import yaml

GATES = ("schema", "environment", "fixed_mini_panel", "canonical_mapping", "pgf1_evidence", "real_input")
CONTRACT = "pgbench-adapter-conformance-v2"
ADAPTER_CONTRACT = "pgbench-reviewed-adapter-v1"
SCHEMA = Path(__file__).resolve().parents[1] / "schemas/tool.schema.yaml"
VALIDATOR = Path(__file__).resolve()
PRIMARY_POLICY = {
    "configuration_frozen_before_execution": True,
    "allowed_configuration_sources": ["official_default", "official_documentation"],
    "allowed_adapter_exceptions": ["input_compatibility", "resource_limits", "necessary_format_conversion"],
    "extra_genotype_filtering_or_rejection": "forbidden",
    "primary_scores_per_tool_and_comparison_context": 1,
    "select_highest_score_across_configurations": False,
    "superseded_runs": "historical_archive_only_not_official_leaderboard",
}


class AdapterConformanceError(ValueError):
    pass


def sha(path: Path) -> str:
    result = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            result.update(chunk)
    return result.hexdigest()


def artifact(path: Path) -> dict[str, Any]:
    if path.is_symlink() or not path.is_file() or path.stat().st_size == 0:
        raise AdapterConformanceError("missing/empty/non-owned gate evidence: " + str(path))
    return {"path": str(path.resolve()), "sha256": sha(path), "size_bytes": path.stat().st_size}


def source_inventory(root: Path, mutable_artifacts: set[Path] | None = None) -> dict[str, str]:
    excluded = {"__pycache__", ".git", ".pytest_cache"}
    result = {}
    for path in sorted(root.rglob("*")):
        relative = path.relative_to(root)
        if any(part in excluded for part in relative.parts):
            continue
        # Exclude only plan-declared generated artifacts, not whole file types.
        # JSON fixtures and JS/native gate helpers are reviewed source too.
        if path in (mutable_artifacts or set()):
            continue
        if path.is_symlink():
            raise AdapterConformanceError("reviewed plugin source must not contain symlinks")
        if path.is_file():
            result[relative.as_posix()] = sha(path)
    return result


def validate_primary_configuration(manifest: dict[str, Any]) -> None:
    information = manifest["information_contract"]
    mode, source = information["tuning_mode"], information["parameter_source"]
    if (mode, source) not in {("official_default", "official_default"),
                               ("official_recommended_sv_config", "official_documentation")}:
        raise AdapterConformanceError("formal primary configuration must use frozen official defaults/recommendations")
    if manifest.get("formal_primary_policy") != PRIMARY_POLICY:
        raise AdapterConformanceError("formal primary policy is incomplete or permits extra filtering/score selection")
    parameters = manifest.get("parameter_contract")
    if not isinstance(parameters, dict) or not parameters.get("explicit") or not parameters.get("implicit"):
        raise AdapterConformanceError("nonempty explicit/implicit parameter contract is required")
    release = parameters.get("source_release", {})
    if (not release.get("tag") or "REPLACE" in str(release)
            or release.get("commit") == "0" * 40 or release.get("archive_sha256") == "0" * 64):
        raise AdapterConformanceError("replace every native release placeholder before conformance")
    for exception in parameters.get("adapter_exceptions", []):
        if (exception.get("category") not in PRIMARY_POLICY["allowed_adapter_exceptions"]
                or not exception.get("rationale") or not exception.get("parameters")):
            raise AdapterConformanceError("adapter exceptions need a permitted category, rationale and exact parameters")


def _manifest(path: Path) -> dict[str, Any]:
    manifest = yaml.safe_load(path.read_text(encoding="utf-8"))
    jsonschema.Draft202012Validator(yaml.safe_load(SCHEMA.read_text())).validate(manifest)
    validate_primary_configuration(manifest)
    if manifest["execution"].get("adapter_contract") != ADAPTER_CONTRACT:
        raise AdapterConformanceError("reviewed-adapter contract must be declared")
    return manifest


def conformance_path(path: Path, root: Path) -> Path:
    """Reject escapes and symlink parents before resolving an evidence path."""
    path = path.absolute()
    resolved = path.resolve()
    try:
        resolved.relative_to(root / "conformance")
    except ValueError as exc:
        raise AdapterConformanceError("conformance artifact escapes adapter: " + str(path)) from exc
    for component in (path, *path.parents):
        if component.is_symlink():
            raise AdapterConformanceError("conformance artifacts must not use symlinks: " + str(component))
        if component == root:
            break
    return resolved


def load_plan(checks: Path, root: Path, out: Path) -> tuple[dict, dict]:
    checks = conformance_path(checks, root)
    frozen = artifact(checks)
    spec = json.loads(checks.read_text(encoding="utf-8"))
    if artifact(checks) != frozen:
        raise AdapterConformanceError("conformance execution plan changed while being read")
    gates = spec.get("gates", {}) if isinstance(spec, dict) else {}
    if (not isinstance(spec, dict) or spec.get("contract") != "pgbench-adapter-gate-plan-v2"
            or not isinstance(gates, dict) or set(gates) != set(GATES)):
        raise AdapterConformanceError("a reviewed six-gate execution plan is required")
    reserved = {checks, out, *(out.parent / (name + ".log") for name in GATES)}
    proofs = set()
    for name in GATES:
        gate = gates[name]
        argv = gate.get("command") if isinstance(gate, dict) else None
        evidence = gate.get("evidence") if isinstance(gate, dict) else None
        timeout = gate.get("timeout_seconds", 600) if isinstance(gate, dict) else None
        if (not isinstance(argv, list) or not argv or not all(isinstance(arg, str) and arg for arg in argv)
                or not isinstance(evidence, list) or not evidence
                or not all(isinstance(item, str) and item for item in evidence)
                or type(timeout) is not int or not 1 <= timeout <= 3600):
            raise AdapterConformanceError("malformed gate plan: " + name)
        if "REPLACE" in json.dumps(argv):
            raise AdapterConformanceError("replace the gate command placeholder: " + name)
        for item in evidence:
            proof = conformance_path(root / item, root)
            if proof in reserved or proof in proofs:
                raise AdapterConformanceError("gate evidence aliases a plan, receipt, log or another proof: " + name)
            proofs.add(proof)
    return gates, frozen


def mutable_gate_artifacts(gates: dict, root: Path, checks: Path, out: Path) -> set[Path]:
    return {checks, out, *(out.parent / (name + ".log") for name in GATES),
            *(conformance_path(root / item, root) for gate in gates.values() for item in gate["evidence"])}


def verify_receipt(path: Path, tool_manifest: Path) -> dict[str, Any]:
    manifest = _manifest(tool_manifest)
    root = tool_manifest.resolve().parent
    path = conformance_path(path, root)
    value = json.loads(path.read_text(encoding="utf-8"))
    plan = value.get("gate_plan") if isinstance(value, dict) else None
    if not isinstance(plan, dict) or not isinstance(plan.get("path"), str):
        raise AdapterConformanceError("changed/unbound conformance execution plan")
    checks = conformance_path(Path(plan["path"]), root)
    expected_gates, frozen_plan = load_plan(checks, root, path)
    if frozen_plan != plan:
        raise AdapterConformanceError("changed/unbound conformance execution plan")
    if (not isinstance(value, dict) or value.get("contract") != CONTRACT
            or value.get("production_eligibility") != "PASS"
            or value.get("tool_id") != manifest["id"]
            or value.get("tool_manifest_sha256") != sha(tool_manifest)
            or value.get("schema_sha256") != sha(SCHEMA)
            or value.get("validator_sha256") != sha(VALIDATOR)
            or value.get("source_inventory") != source_inventory(root, mutable_gate_artifacts(expected_gates, root, checks, path))
            or value.get("parameter_contract") != manifest["parameter_contract"]
            or value.get("formal_primary_policy") != PRIMARY_POLICY):
        raise AdapterConformanceError("stale/unbound adapter conformance receipt")
    gates = value.get("gates", {})
    if not isinstance(gates, dict) or set(gates) != set(GATES):
        raise AdapterConformanceError("all six conformance gates must be recorded")
    for name in GATES:
        row = gates[name]
        if (not isinstance(row, dict) or row.get("status") != "PASS" or row.get("returncode") != 0
                or not isinstance(row.get("command"), list) or not row["command"]
                or not isinstance(row.get("evidence"), list) or not row["evidence"]
                or not isinstance(row.get("log"), dict)):
            raise AdapterConformanceError("incomplete executed gate: " + name)
        expected = expected_gates[name]
        expected_proofs = [str(conformance_path(root / item, root)) for item in expected["evidence"]]
        if (row["command"] != expected["command"]
                or row.get("timeout_seconds") != expected.get("timeout_seconds", 600)
                or [entry.get("path") if isinstance(entry, dict) else None for entry in row["evidence"]] != expected_proofs
                or row["log"].get("path") != str(path.parent / (name + ".log"))):
            raise AdapterConformanceError("executed gate differs from frozen plan: " + name)
        for entry in [row["log"], *row["evidence"]]:
            if not isinstance(entry, dict) or not isinstance(entry.get("path"), str):
                raise AdapterConformanceError("invalid conformance artifact: " + name)
            evidence_path = conformance_path(Path(entry["path"]), root)
            current = artifact(evidence_path)
            if current != entry:
                raise AdapterConformanceError("changed conformance evidence: " + name)
        with Path(row["log"]["path"]).open("rb") as stream:
            if json.loads(stream.readline()).get("command") != row["command"]:
                raise AdapterConformanceError("gate command differs from execution log: " + name)
    return value


def validate(checks: Path, out: Path, *, tool_manifest: Path, execute: bool = False) -> dict[str, Any]:
    """Execute gate commands, never promote a hand-written list of PASS strings."""
    if not execute:
        raise AdapterConformanceError("PRODUCTION_NOT_ELIGIBLE: generate a receipt with --execute")
    manifest = _manifest(tool_manifest)
    root = tool_manifest.resolve().parent
    checks = conformance_path(checks, root)
    out = conformance_path(out, root)
    if out.parent != root / "conformance":
        raise AdapterConformanceError("receipt must live under the adapter conformance/ directory")
    if out.exists():
        raise AdapterConformanceError("refusing to overwrite an existing frozen conformance receipt")
    gates, frozen_plan = load_plan(checks, root, out)
    frozen_sources = {"tool_manifest_sha256": (tool_manifest, sha(tool_manifest)),
                      "schema_sha256": (SCHEMA, sha(SCHEMA)), "validator_sha256": (VALIDATOR, sha(VALIDATOR))}
    mutable = mutable_gate_artifacts(gates, root, checks, out)
    before = source_inventory(root, mutable)
    out.parent.mkdir(parents=True, exist_ok=True)
    results = {}
    for name in GATES:
        gate = gates[name]
        argv, evidence, timeout = gate["command"], gate["evidence"], gate.get("timeout_seconds", 600)
        if artifact(checks) != frozen_plan:
            raise AdapterConformanceError("conformance execution plan changed while running gates")
        proof_paths = [conformance_path(root / item, root) for item in evidence]
        log = conformance_path(out.parent / (name + ".log"), root)
        started_ns = time.time_ns()
        with log.open("wb") as stream:
            stream.write((json.dumps({"command": argv}) + "\n").encode())
            stream.flush()
            result = subprocess.run(argv, cwd=root, stdout=stream, stderr=subprocess.STDOUT, timeout=timeout, check=False)
        if artifact(conformance_path(checks, root)) != frozen_plan:
            raise AdapterConformanceError("conformance execution plan changed while running gates")
        if result.returncode != 0:
            raise AdapterConformanceError("PRODUCTION_NOT_ELIGIBLE: failed executed gate " + name)
        if any(not path.is_file() or path.stat().st_mtime_ns < started_ns for path in proof_paths):
            raise AdapterConformanceError("gate must produce fresh evidence, not reuse stale artifacts: " + name)
        proof_paths = [conformance_path(root / item, root) for item in evidence]
        log = conformance_path(log, root)
        results[name] = {"status": "PASS", "command": argv, "returncode": result.returncode, "timeout_seconds": timeout,
            "log": artifact(log), "evidence": [artifact(path) for path in proof_paths]}
    if source_inventory(root, mutable) != before:
        raise AdapterConformanceError("adapter source/parameters changed while running gates")
    if any(sha(path) != expected for path, expected in frozen_sources.values()):
        raise AdapterConformanceError("conformance validator/schema/manifest changed while running gates")
    receipt = {"contract": CONTRACT, "tool_id": manifest["id"], "production_eligibility": "PASS",
        **{name: expected for name, (_, expected) in frozen_sources.items()},
        "source_inventory": before, "parameter_contract": manifest["parameter_contract"],
        "formal_primary_policy": PRIMARY_POLICY, "gates": results, "gate_plan": frozen_plan,
        "review_scope": "executed commands and immutable evidence; semantic completeness requires trusted review"}
    # Verify before publishing any PASS receipt. A failed postcondition must not
    # leave a superficially successful report for a downstream consumer.
    with tempfile.NamedTemporaryFile(mode="w", encoding="utf-8", dir=out.parent, prefix=".receipt-", delete=False) as stream:
        temporary = Path(stream.name)
        stream.write(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    try:
        verify_receipt(temporary, tool_manifest)
        # Atomic publication with no overwrite: readers cannot observe a partial
        # PASS report, and a concurrent writer cannot replace an existing one.
        os.link(temporary, out)
    finally:
        temporary.unlink()
    return receipt


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--checks", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--tool-manifest", required=True, type=Path)
    parser.add_argument("--execute", action="store_true")
    parser.add_argument("--verify", type=Path)
    args = parser.parse_args(argv)
    try:
        if args.verify:
            verify_receipt(args.verify, args.tool_manifest)
        elif args.checks and args.output:
            validate(args.checks, args.output, tool_manifest=args.tool_manifest, execute=args.execute)
        else:
            parser.error("use --verify, or --checks/--output/--execute")
    except (OSError, ValueError, jsonschema.ValidationError, subprocess.SubprocessError) as exc:
        print("validate_adapter_conformance: " + str(exc), file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
