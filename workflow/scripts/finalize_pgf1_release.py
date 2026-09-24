#!/usr/bin/env python3
"""Create an immutable, artifact-complete seal for a PG-F1 release."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import tempfile
from pathlib import Path
from typing import Any, Sequence


class PGF1SealError(ValueError):
    """Raised when a formal PG-F1 release cannot be sealed."""


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_json(path: Path, *, label: str) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as error:
        raise PGF1SealError(f"{label} is not valid JSON: {path}") from error
    if not isinstance(value, dict):
        raise PGF1SealError(f"{label} must be a JSON object: {path}")
    return value


def atomic_json(path: Path, value: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, name = tempfile.mkstemp(dir=path.parent, prefix=f".{path.name}.", suffix=".tmp")
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8", newline="\n") as handle:
            json.dump(value, handle, indent=2, sort_keys=True)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(name, path)
    finally:
        Path(name).unlink(missing_ok=True)


def artifact(path: Path) -> dict[str, Any]:
    if not path.is_file():
        raise PGF1SealError(f"required release artifact is missing: {path}")
    return {"path": str(path), "sha256": sha256_file(path), "size_bytes": path.stat().st_size}


def _expected_hash(value: object, *, label: str) -> str:
    if not isinstance(value, str) or len(value) != 64:
        raise PGF1SealError(f"{label} is missing or not a SHA256")
    return value


def seal(
    *,
    score: Path,
    score_audit: Path,
    ledger: Path,
    canonical_units: Path,
    unit_sha_file: Path,
    all_sites: Path,
    query: Path,
    query_audit: Path,
    evidence: dict[str, Path],
    raw_ledgers: dict[str, Path],
    invariants: Path,
    admission: Path,
    run_context: Path,
    config: Path,
    output_package: Path,
    output_audit: Path,
    output_seal: Path,
) -> dict[str, Any]:
    score_data = load_json(score, label="PG-F1 score")
    score_audit_data = load_json(score_audit, label="PG-F1 score audit")
    admission_data = load_json(admission, label="leaderboard admission")
    invariants_data = load_json(invariants, label="leaderboard invariants")

    if score_data.get("score_name") != "PG-F1" or score_data.get("score_status") != "valid":
        raise PGF1SealError("only a valid PG-F1 score can be sealed")
    if admission_data.get("status") != "PASS":
        raise PGF1SealError("leaderboard admission did not pass")
    if admission_data.get("score_sha256") != sha256_file(score):
        raise PGF1SealError("leaderboard admission does not bind this score")

    declared_unit_sha = unit_sha_file.read_text(encoding="utf-8").strip().split()[0]
    units_sha = sha256_file(canonical_units)
    if declared_unit_sha != units_sha:
        raise PGF1SealError("canonical-unit SHA256 file does not match canonical units")
    if score_data.get("canonical_unit_set_sha256") != units_sha:
        raise PGF1SealError("score does not bind the canonical-unit set")
    if invariants_data.get("canonical_unit_set_sha256") != units_sha:
        raise PGF1SealError("leaderboard invariants do not bind the canonical-unit set")
    if score_data.get("pgf1_evidence_sha256") != sha256_file(ledger):
        raise PGF1SealError("score does not bind the PG-F1 evidence ledger")
    if score_audit_data.get("canonical_unit_set_sha256") != units_sha:
        raise PGF1SealError("PG-F1 audit does not bind the canonical-unit set")
    if score_audit_data.get("all_sites_sha256") != sha256_file(all_sites):
        raise PGF1SealError("PG-F1 audit does not bind all-sites input")
    if score_audit_data.get("pgf1_evidence_sha256") != sha256_file(ledger):
        raise PGF1SealError("PG-F1 audit does not bind the evidence ledger")

    for evaluator in ("truvari", "aardvark", "vcfdist"):
        expected = _expected_hash(
            score_audit_data.get("evidence_sha256", {}).get(evaluator)
            if isinstance(score_audit_data.get("evidence_sha256"), dict)
            else None,
            label=f"PG-F1 audit evidence_sha256.{evaluator}",
        )
        actual = sha256_file(evidence[evaluator])
        if expected != actual:
            raise PGF1SealError(f"PG-F1 audit does not bind {evaluator} evidence")

    artifacts = {
        "score": artifact(score),
        "score_audit": artifact(score_audit),
        "pgf1_evidence_ledger": artifact(ledger),
        "canonical_units": artifact(canonical_units),
        "canonical_unit_sha256": artifact(unit_sha_file),
        "all_sites": artifact(all_sites),
        "evaluator_query": artifact(query),
        "evaluator_query_audit": artifact(query_audit),
        "leaderboard_invariants": artifact(invariants),
        "leaderboard_admission": artifact(admission),
        "run_context": artifact(run_context),
        "config": artifact(config),
        "normalized_evidence": {name: artifact(evidence[name]) for name in sorted(evidence)},
        "raw_evaluator_ledgers": {name: artifact(raw_ledgers[name]) for name in sorted(raw_ledgers)},
    }
    tuple_data = score_data.get("tuple")
    if not isinstance(tuple_data, dict):
        raise PGF1SealError("PG-F1 score lacks a run/sample/tool/track tuple")
    audit = {
        "schema_version": "pgbench_pgf1_release_provenance_v1",
        "seal_status": "valid",
        "hash_lineage_complete": True,
        "missing_artifacts": [],
        "canonical_unit_set_sha256": units_sha,
        "score_contract_sha256": score_data.get("score_contract_sha256"),
        "tuple": tuple_data,
        "checks": {
            "score_is_valid_pgf1": True,
            "admission_passed_and_binds_score": True,
            "canonical_unit_identity_consistent": True,
            "score_and_audit_bind_evidence": True,
            "three_raw_evaluator_ledgers_hashed": True,
            "three_normalized_evidence_files_hashed": True,
        },
        "artifacts": artifacts,
    }
    package = {
        "package_schema_version": "pgbench_pgf1_release_package_v1",
        "sealed": True,
        "seal_status": "valid",
        "score": score_data,
        "admission": admission_data,
        "provenance_audit": {"path": str(output_audit), "sha256": None},
        "artifact_hashes": artifacts,
    }
    atomic_json(output_audit, audit)
    package["provenance_audit"]["sha256"] = sha256_file(output_audit)
    atomic_json(output_package, package)
    seal_record = {
        "schema_version": "pgbench_pgf1_release_seal_v1",
        "sealed": True,
        "score_package_sha256": sha256_file(output_package),
        "provenance_audit_sha256": sha256_file(output_audit),
        "score_sha256": sha256_file(score),
        "admission_sha256": sha256_file(admission),
    }
    atomic_json(output_seal, seal_record)
    return package


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("score", "score-audit", "ledger", "canonical-units", "unit-sha-file", "all-sites", "query", "query-audit", "invariants", "admission", "run-context", "config"):
        parser.add_argument(f"--{name}", required=True, type=Path)
    for evaluator in ("truvari", "aardvark", "vcfdist"):
        parser.add_argument(f"--{evaluator}-evidence", required=True, type=Path)
        parser.add_argument(f"--raw-{evaluator}-ledger", required=True, type=Path)
    parser.add_argument("--output-package", required=True, type=Path)
    parser.add_argument("--output-audit", required=True, type=Path)
    parser.add_argument("--output-seal", required=True, type=Path)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        seal(
            score=args.score, score_audit=args.score_audit, ledger=args.ledger,
            canonical_units=args.canonical_units, unit_sha_file=args.unit_sha_file,
            all_sites=args.all_sites, query=args.query, query_audit=args.query_audit,
            evidence={name: getattr(args, f"{name}_evidence") for name in ("truvari", "aardvark", "vcfdist")},
            raw_ledgers={name: getattr(args, f"raw_{name}_ledger") for name in ("truvari", "aardvark", "vcfdist")},
            invariants=args.invariants, admission=args.admission,
            run_context=args.run_context, config=args.config,
            output_package=args.output_package, output_audit=args.output_audit,
            output_seal=args.output_seal,
        )
    except (OSError, PGF1SealError) as error:
        print(f"finalize_pgf1_release: {error}", file=__import__("sys").stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
