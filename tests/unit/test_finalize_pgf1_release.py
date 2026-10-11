from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import pytest

SCRIPTS = Path(__file__).resolve().parents[2] / "workflow" / "scripts"
sys.path.insert(0, str(SCRIPTS))

import finalize_pgf1_release as finalizer  # noqa: E402


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _json(path: Path, value: dict) -> None:
    path.write_text(json.dumps(value), encoding="utf-8")


def _fixture(tmp_path: Path) -> dict[str, Path]:
    paths: dict[str, Path] = {}
    for name in ("units", "all_sites", "query", "query_audit", "ledger", "context", "config"):
        path = tmp_path / name
        path.write_text(f"{name}\n", encoding="utf-8")
        paths[name] = path
    unit_sha = tmp_path / "units.sha256"
    unit_sha.write_text(_sha(paths["units"]) + "  units\n", encoding="utf-8")
    paths["unit_sha"] = unit_sha
    invariants = tmp_path / "invariants"
    _json(invariants, {"canonical_unit_set_sha256": _sha(paths["units"])})
    paths["invariants"] = invariants
    for evaluator in ("truvari", "aardvark", "vcfdist"):
        for suffix in ("evidence", "raw"):
            path = tmp_path / f"{evaluator}.{suffix}"
            path.write_text(f"{evaluator} {suffix}\n", encoding="utf-8")
            paths[f"{evaluator}_{suffix}"] = path
    ledger_sha = _sha(paths["ledger"])
    score_audit = tmp_path / "score-audit.json"
    _json(score_audit, {
        "canonical_unit_set_sha256": _sha(paths["units"]),
        "all_sites_sha256": _sha(paths["all_sites"]),
        "pgf1_evidence_sha256": ledger_sha,
        "evidence_sha256": {e: _sha(paths[f"{e}_evidence"]) for e in ("truvari", "aardvark", "vcfdist")},
    })
    paths["score_audit"] = score_audit
    score = tmp_path / "score.json"
    _json(score, {
        "score_name": "PG-F1", "score_status": "valid",
        "canonical_unit_set_sha256": _sha(paths["units"]),
        "pgf1_evidence_sha256": ledger_sha,
        "score_contract_sha256": "a" * 64,
        "tuple": {"run_id": "run", "sample": "HG002", "tool": "tool", "track": "short_read"},
    })
    paths["score"] = score
    admission = tmp_path / "admission.json"
    _json(admission, {"status": "PASS", "score_sha256": _sha(score)})
    paths["admission"] = admission
    for name in ("package", "audit", "seal"):
        paths[f"output_{name}"] = tmp_path / "release" / f"{name}.json"
    return paths


def _seal(paths: dict[str, Path], *, require_projection=False, **kwargs) -> dict:
    return finalizer.seal(
        score=paths["score"], score_audit=paths["score_audit"], ledger=paths["ledger"],
        canonical_units=paths["units"], unit_sha_file=paths["unit_sha"],
        all_sites=paths["all_sites"], query=paths["query"], query_audit=paths["query_audit"],
        evidence={e: paths[f"{e}_evidence"] for e in ("truvari", "aardvark", "vcfdist")},
        raw_ledgers={e: paths[f"{e}_raw"] for e in ("truvari", "aardvark", "vcfdist")},
        invariants=paths["invariants"], admission=paths["admission"],
        run_context=paths["context"], config=paths["config"],
        output_package=paths["output_package"], output_audit=paths["output_audit"], output_seal=paths["output_seal"],
        require_projection_audit=require_projection,
        **kwargs,
    )


def test_seal_writes_complete_hashed_pgf1_release(tmp_path: Path) -> None:
    paths = _fixture(tmp_path)
    package = _seal(paths)
    assert package["sealed"] is True
    audit = json.loads(paths["output_audit"].read_text(encoding="utf-8"))
    assert audit["hash_lineage_complete"] is True
    assert set(audit["artifacts"]["raw_evaluator_ledgers"]) == {"truvari", "aardvark", "vcfdist"}
    assert paths["output_seal"].is_file()


def test_vg_cannot_seal_without_native_projection_lineage(tmp_path):
    paths = _fixture(tmp_path)
    score = json.loads(paths["score"].read_text())
    score["tuple"]["tool"] = "vg"
    _json(paths["score"], score)
    _json(paths["admission"], {"status": "PASS", "score_sha256": _sha(paths["score"])})
    with pytest.raises(finalizer.PGF1SealError, match="native projection audit"):
        _seal(paths, require_projection=True)


def test_seal_rejects_admission_for_a_different_score(tmp_path: Path) -> None:
    paths = _fixture(tmp_path)
    _json(paths["admission"], {"status": "PASS", "score_sha256": "0" * 64})
    with pytest.raises(finalizer.PGF1SealError, match="does not bind"):
        _seal(paths)


def test_native_audit_binds_policy_and_retained_artifacts(tmp_path):
    paths = _fixture(tmp_path)
    native = tmp_path / "native.vcf"
    native.write_text("native\n")
    audit = tmp_path / "native-audit.json"
    _json(audit, {"contract": "pgbench-native-output-audit-v1", "tool_id": "tool", "truth_used": False,
        "execution_policy": "default", "output_sha256": _sha(paths["all_sites"]),
        "artifacts": {"native": {"path": str(native), "sha256": _sha(native)}}})
    kwargs = dict(native_output_audit=audit, native_output_vcf=paths["all_sites"], require_native_output_audit=True)
    package = _seal(paths, **kwargs)
    assert package["release_qualifiers"]["native_execution_policy"] == "default"
    native.write_text("changed\n")
    with pytest.raises(finalizer.PGF1SealError, match="lineage hash mismatch"):
        _seal(paths, **kwargs)


def test_required_native_audit_cannot_be_omitted(tmp_path):
    with pytest.raises(finalizer.PGF1SealError, match="required"):
        _seal(_fixture(tmp_path), require_native_output_audit=True)
