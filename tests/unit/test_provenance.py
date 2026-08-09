from __future__ import annotations

import json
import os
import sys
from pathlib import Path
from typing import Any

import pytest

SCRIPTS = Path(__file__).resolve().parents[2] / "workflow" / "scripts"
sys.path.insert(0, str(SCRIPTS))

from audit_provenance import _parse_expected_jobs  # noqa: E402
from pgbench_provenance import (  # noqa: E402
    MANIFEST_SCHEMA_VERSION,
    ManifestValidationError,
    PathHashError,
    atomic_write_json,
    audit_manifests,
    build_lineage,
    build_manifest_id,
    fingerprint_path,
    lineage_tsv,
    manifest_validation_errors,
    prepare_rule_manifest,
    sha256_bytes,
    sha256_directory,
    sha256_path,
    validate_manifest,
    write_rule_manifest,
)

SHA_A = sha256_bytes(b"a")
SHA_B = sha256_bytes(b"b")


def _symlink_or_skip(link: Path, target: str) -> None:
    try:
        link.symlink_to(target)
    except OSError as exc:
        if os.name == "nt" and getattr(exc, "winerror", None) == 1314:
            pytest.skip("Windows symlink privilege is not enabled")
        raise


def _manifest_payload(
    *,
    rule_name: str,
    job_key: str,
    attempt_id: str,
    inputs: list[Path],
    outputs: list[Path],
    upstream_manifest_ids: list[str] | None = None,
    run_id: str = "synthetic-run",
    score_profile_sha256: str = SHA_B,
) -> dict[str, Any]:
    return {
        "manifest_schema_version": MANIFEST_SCHEMA_VERSION,
        "attempt_id": attempt_id,
        "run_id": run_id,
        "rule_name": rule_name,
        "job_key": job_key,
        "wildcards": {"sample": "HG002"},
        "module_or_tool_id": "pgbench-core",
        "snakefile_path": "workflow/rules/provenance.smk",
        "rule_source_sha256": SHA_A,
        "script_or_wrapper_path": "workflow/scripts/mock.py",
        "script_or_wrapper_sha256": SHA_B,
        "command": ["python", "mock.py"],
        "params": {"score_profile": "pgbench_v1"},
        "threads": 1,
        "requested_resources": {"mem_mb": 128},
        "input_paths": [str(path) for path in inputs],
        "input_sha256": {},
        "input_size": {},
        "input_mtime": {},
        "output_paths": [str(path) for path in outputs],
        "output_sha256": {},
        "config_snapshot_sha256": SHA_A,
        "score_profile_sha256": score_profile_sha256,
        "pangenome_manifest_sha256": SHA_B,
        "reference_sha256": SHA_A,
        "truth_profile": "giab_hg002_grch38_v5_0q",
        "conda_lock_sha256": SHA_B,
        "container_uri": None,
        "container_digest": None,
        "git_head": "b" * 40,
        "git_dirty": False,
        "git_diff_sha256": SHA_A,
        "snakemake_version": "9.0.0",
        "execution_profile": "pgbench_v1",
        "hardware_fingerprint_sha256": SHA_B,
        "random_seed": 20260717,
        "upstream_manifest_ids": upstream_manifest_ids or [],
        "started_at": "2026-07-17T00:00:00+00:00",
        "finished_at": "2026-07-17T00:00:01+00:00",
        "exit_code": 0,
        "status": "success",
        "not_applicable_reason": {
            "container_uri": "Conda lock is the selected environment identity",
            "container_digest": "Conda lock is the selected environment identity",
        },
    }


def _prepared_manifest(
    *,
    rule_name: str,
    job_key: str,
    attempt_id: str,
    inputs: list[Path],
    outputs: list[Path],
    upstream_manifest_ids: list[str] | None = None,
    run_id: str = "synthetic-run",
    score_profile_sha256: str = SHA_B,
) -> dict[str, Any]:
    return prepare_rule_manifest(
        _manifest_payload(
            rule_name=rule_name,
            job_key=job_key,
            attempt_id=attempt_id,
            inputs=inputs,
            outputs=outputs,
            upstream_manifest_ids=upstream_manifest_ids,
            run_id=run_id,
            score_profile_sha256=score_profile_sha256,
        )
    )


def _write_companions(
    root: Path,
    manifest: dict[str, Any],
    *,
    run_scoped: bool = False,
) -> None:
    scope = [manifest["run_id"]] if run_scoped else []
    log = (
        root
        / "logs"
        / Path(*scope)
        / "rules"
        / manifest["rule_name"]
        / f"{manifest['job_key']}.log"
    )
    benchmark = (
        root
        / "benchmarks"
        / Path(*scope)
        / "rules"
        / manifest["rule_name"]
        / f"{manifest['job_key']}.jsonl"
    )
    log.parent.mkdir(parents=True, exist_ok=True)
    benchmark.parent.mkdir(parents=True, exist_ok=True)
    log.write_text("completed\n", encoding="utf-8")
    benchmark.write_text('{"s": 1.0}\n', encoding="utf-8")


def test_inline_expected_jobs_are_fixed_core_dependencies() -> None:
    jobs = _parse_expected_jobs(
        None,
        [
            "validate_config=config",
            "tool__custom__execute=HG002.end_to_end_from_reads",
        ],
    )
    assert jobs == [
        {
            "rule_name": "validate_config",
            "job_key": "config",
            "core": True,
        },
        {
            "rule_name": "tool__custom__execute",
            "job_key": "HG002.end_to_end_from_reads",
            "core": True,
        },
    ]


def test_file_directory_and_symlink_hashes_are_deterministic(tmp_path: Path) -> None:
    tree = tmp_path / "tree"
    tree.mkdir()
    target = tree / "target.txt"
    target.write_text("alpha\n", encoding="utf-8")
    nested = tree / "nested"
    nested.mkdir()
    (nested / "z.txt").write_text("z\n", encoding="utf-8")
    link = tree / "target-link"
    _symlink_or_skip(link, "target.txt")

    first = sha256_path(tree)
    assert first == sha256_path(tree)
    assert sha256_directory(tree) == first
    link_fingerprint = fingerprint_path(link)
    assert link_fingerprint.path_type == "symlink"
    assert link_fingerprint.link_target == "target.txt"
    assert link_fingerprint.target_sha256 == sha256_path(target)

    target.write_text("beta\n", encoding="utf-8")
    assert sha256_path(tree) != first
    assert fingerprint_path(link).target_sha256 == sha256_path(target)


def test_empty_directory_and_broken_symlink_have_explicit_behavior(
    tmp_path: Path,
) -> None:
    empty_a = tmp_path / "empty-a"
    empty_b = tmp_path / "empty-b"
    empty_a.mkdir()
    empty_b.mkdir()
    assert sha256_path(empty_a) == sha256_path(empty_b)

    broken = tmp_path / "broken"
    _symlink_or_skip(broken, "missing")
    with pytest.raises(PathHashError, match="broken symlink"):
        sha256_path(broken)


def test_manifest_id_uses_only_normative_identity_fields(tmp_path: Path) -> None:
    output = tmp_path / "result.txt"
    output.write_text("result\n", encoding="utf-8")
    manifest = _prepared_manifest(
        rule_name="prepare_reference",
        job_key="hg002",
        attempt_id="attempt-1",
        inputs=[],
        outputs=[output],
    )
    changed_params = dict(manifest)
    changed_params["params"] = {"different": True}
    assert build_manifest_id(changed_params) == manifest["manifest_id"]

    changed_source = dict(manifest)
    changed_source["rule_source_sha256"] = SHA_B
    assert build_manifest_id(changed_source) != manifest["manifest_id"]
    validate_manifest(manifest, verify_paths=True, require_success=True)


def test_manifest_validation_reports_identity_and_null_contract(
    tmp_path: Path,
) -> None:
    output = tmp_path / "result.txt"
    output.write_text("result\n", encoding="utf-8")
    manifest = _prepared_manifest(
        rule_name="prepare_reference",
        job_key="hg002",
        attempt_id="attempt-1",
        inputs=[],
        outputs=[output],
    )
    invalid = dict(manifest)
    invalid["manifest_id"] = "0" * 64
    invalid["container_uri"] = None
    invalid["not_applicable_reason"] = {
        "container_digest": "Conda environment selected"
    }
    errors = manifest_validation_errors(invalid)
    assert any("manifest_id does not match" in error for error in errors)
    assert any("container_uri is null" in error for error in errors)
    with pytest.raises(ManifestValidationError):
        validate_manifest(invalid)


def test_manifest_requires_score_profile_sha256(tmp_path: Path) -> None:
    output = tmp_path / "result.txt"
    output.write_text("result\n", encoding="utf-8")
    manifest = _prepared_manifest(
        rule_name="prepare_reference",
        job_key="hg002",
        attempt_id="attempt-1",
        inputs=[],
        outputs=[output],
    )

    missing = dict(manifest)
    missing.pop("score_profile_sha256")
    assert "missing required field: score_profile_sha256" in manifest_validation_errors(
        missing
    )

    invalid = dict(manifest)
    invalid["score_profile_sha256"] = "not-a-sha256"
    assert any(
        "score_profile_sha256 must be a lowercase SHA-256" in error
        for error in manifest_validation_errors(invalid)
    )


def test_write_rule_manifest_is_atomic_and_excludes_self_hash(
    tmp_path: Path,
) -> None:
    result = tmp_path / "result.txt"
    result.write_text("result\n", encoding="utf-8")
    destination = tmp_path / "manifest.json"
    payload = _manifest_payload(
        rule_name="prepare_reference",
        job_key="hg002",
        attempt_id="attempt-1",
        inputs=[],
        outputs=[result],
    )
    written = write_rule_manifest(payload, destination)
    assert json.loads(destination.read_text(encoding="utf-8")) == written
    assert not list(tmp_path.glob(".manifest.json.*.tmp"))

    recursive_payload = dict(payload)
    recursive_payload["output_paths"] = [str(destination)]
    with pytest.raises(ManifestValidationError, match="must not include itself"):
        write_rule_manifest(recursive_payload, destination)


def test_atomic_json_write_replaces_existing_document(tmp_path: Path) -> None:
    destination = tmp_path / "document.json"
    destination.write_text('{"old": true}\n', encoding="utf-8")
    atomic_write_json(destination, {"new": ["值", 1]})
    assert json.loads(destination.read_text(encoding="utf-8")) == {"new": ["值", 1]}


def test_lineage_links_shared_hashes_and_has_stable_tsv(tmp_path: Path) -> None:
    intermediate = tmp_path / "canonical.vcf"
    intermediate.write_text("canonical\n", encoding="utf-8")
    final = tmp_path / "score.json"
    final.write_text("{}\n", encoding="utf-8")
    upstream = _prepared_manifest(
        rule_name="canonicalize_vcf",
        job_key="hg002",
        attempt_id="attempt-1",
        inputs=[],
        outputs=[intermediate],
    )
    downstream = _prepared_manifest(
        rule_name="compute_pgbench_score",
        job_key="hg002",
        attempt_id="attempt-1",
        inputs=[intermediate],
        outputs=[final],
        upstream_manifest_ids=[upstream["manifest_id"]],
    )

    lineage = build_lineage(
        [downstream, upstream],
        target_manifest_ids=[downstream["manifest_id"]],
    )
    assert lineage["status"] == "complete"
    assert lineage["hash_lineage_complete"] is True
    assert lineage["topological_order"] == [
        upstream["manifest_id"],
        downstream["manifest_id"],
    ]
    assert lineage["edges"][0]["shared_artifacts"][0]["sha256"] == sha256_path(
        intermediate
    )
    tsv = lineage_tsv(lineage)
    assert "canonicalize_vcf" in tsv
    assert "compute_pgbench_score" in tsv
    assert "true" in tsv


def test_lineage_accepts_hashed_upstream_manifest_attestation(tmp_path: Path) -> None:
    upstream_output = tmp_path / "calls.vcf"
    upstream_output.write_text("calls\n", encoding="utf-8")
    upstream = _prepared_manifest(
        rule_name="evaluate_truvari",
        job_key="hg002",
        attempt_id="attempt-1",
        inputs=[],
        outputs=[upstream_output],
    )
    upstream_manifest = tmp_path / "evaluate_truvari.json"
    atomic_write_json(upstream_manifest, upstream)
    downstream_output = tmp_path / "lineage.json"
    downstream_output.write_text("{}\n", encoding="utf-8")
    downstream = _prepared_manifest(
        rule_name="build_rule_lineage",
        job_key="hg002",
        attempt_id="attempt-1",
        inputs=[upstream_manifest],
        outputs=[downstream_output],
        upstream_manifest_ids=[upstream["manifest_id"]],
    )

    lineage = build_lineage(
        [upstream, downstream],
        target_manifest_ids=[downstream["manifest_id"]],
    )

    assert lineage["hash_lineage_complete"] is True
    shared = lineage["edges"][0]["shared_artifacts"]
    assert shared == [
        {
            "sha256": sha256_path(upstream_manifest),
            "upstream_output_path": f"manifest:{upstream['manifest_id']}",
            "downstream_input_path": str(upstream_manifest),
        }
    ]


def test_lineage_reports_missing_upstream_and_cycle(tmp_path: Path) -> None:
    output_a = tmp_path / "a"
    output_b = tmp_path / "b"
    output_a.write_text("a", encoding="utf-8")
    output_b.write_text("b", encoding="utf-8")
    first = _prepared_manifest(
        rule_name="first",
        job_key="job",
        attempt_id="a",
        inputs=[],
        outputs=[output_a],
        upstream_manifest_ids=["f" * 64],
    )
    missing = build_lineage([first], target_manifest_ids=[first["manifest_id"]])
    assert missing["status"] == "incomplete"
    assert missing["missing_manifest_ids"] == ["f" * 64]

    second = _prepared_manifest(
        rule_name="second",
        job_key="job",
        attempt_id="b",
        inputs=[output_a],
        outputs=[output_b],
    )
    first["upstream_manifest_ids"] = [second["manifest_id"]]
    second["upstream_manifest_ids"] = [first["manifest_id"]]
    cyclic = build_lineage([first, second], target_manifest_ids=[first["manifest_id"]])
    assert cyclic["status"] == "incomplete"
    assert cyclic["cycles"]
    assert build_lineage([first, second])["cycles"]


def test_disconnected_multi_node_lineage_is_incomplete(tmp_path: Path) -> None:
    output_a = tmp_path / "a"
    output_b = tmp_path / "b"
    output_a.write_text("a", encoding="utf-8")
    output_b.write_text("b", encoding="utf-8")
    first = _prepared_manifest(
        rule_name="first",
        job_key="job-a",
        attempt_id="a",
        inputs=[],
        outputs=[output_a],
    )
    second = _prepared_manifest(
        rule_name="second",
        job_key="job-b",
        attempt_id="b",
        inputs=[],
        outputs=[output_b],
    )

    lineage = build_lineage([first, second])

    assert lineage["status"] == "incomplete"
    assert lineage["hash_lineage_complete"] is False
    assert lineage["edges"] == []


def test_audit_calculates_exact_traceability_fields(tmp_path: Path) -> None:
    intermediate = tmp_path / "canonical.vcf"
    intermediate.write_text("canonical\n", encoding="utf-8")
    final = tmp_path / "score.json"
    final.write_text("{}\n", encoding="utf-8")
    upstream = _prepared_manifest(
        rule_name="canonicalize_vcf",
        job_key="hg002",
        attempt_id="attempt-1",
        inputs=[],
        outputs=[intermediate],
    )
    downstream = _prepared_manifest(
        rule_name="compute_pgbench_score",
        job_key="hg002",
        attempt_id="attempt-1",
        inputs=[intermediate],
        outputs=[final],
        upstream_manifest_ids=[upstream["manifest_id"]],
    )
    _write_companions(tmp_path, upstream)
    _write_companions(tmp_path, downstream)

    expected = [
        {"rule_name": "canonicalize_vcf", "job_key": "hg002", "core": True},
        {
            "rule_name": "compute_pgbench_score",
            "job_key": "hg002",
            "core": True,
        },
    ]
    audit = audit_manifests(
        [upstream, downstream],
        expected_jobs=expected,
        target_manifest_ids=[downstream["manifest_id"]],
        workspace_root=tmp_path,
        require_companions=True,
        verify_paths=True,
    )
    assert audit["status"] == "valid"
    assert audit["manifest_completeness"] == 1.0
    assert audit["manifest_points"] == 2.0
    assert audit["hash_lineage_complete"] is True
    assert audit["environment_complete"] is True
    assert audit["run_context_complete"] is True
    assert audit["core_provenance_valid"] is True


def test_audit_discovers_run_id_scoped_companions(tmp_path: Path) -> None:
    result = tmp_path / "score.json"
    result.write_text("{}\n", encoding="utf-8")
    manifest = _prepared_manifest(
        rule_name="compute_pgbench_score",
        job_key="hg002",
        attempt_id="attempt-1",
        inputs=[],
        outputs=[result],
    )
    manifest["params"]["score_profile"] = "pgbench_v1"
    _write_companions(tmp_path, manifest, run_scoped=True)

    audit = audit_manifests(
        [manifest],
        expected_jobs=[
            {
                "rule_name": manifest["rule_name"],
                "job_key": manifest["job_key"],
                "manifest_id": manifest["manifest_id"],
                "core": True,
            }
        ],
        target_manifest_ids=[manifest["manifest_id"]],
        workspace_root=tmp_path,
        require_companions=True,
        verify_paths=True,
    )

    assert audit["status"] == "valid"
    assert audit["valid_manifest_log_benchmark_count"] == 1
    assert not any(
        issue["code"] in {"missing_log", "missing_benchmark"}
        for issue in audit["issues"]
    )


def test_missing_noncore_package_is_provisional_but_missing_core_is_invalid(
    tmp_path: Path,
) -> None:
    result = tmp_path / "result"
    result.write_text("ok", encoding="utf-8")
    manifest = _prepared_manifest(
        rule_name="prepare_reference",
        job_key="hg002",
        attempt_id="attempt-1",
        inputs=[],
        outputs=[result],
    )
    _write_companions(tmp_path, manifest)

    noncore = audit_manifests(
        [manifest],
        expected_jobs=[
            {"rule_name": "prepare_reference", "job_key": "hg002"},
            {"rule_name": "optional_qc", "job_key": "hg002"},
        ],
        target_manifest_ids=[manifest["manifest_id"]],
        workspace_root=tmp_path,
        require_companions=True,
    )
    assert noncore["status"] == "provisional"
    assert noncore["manifest_completeness"] == 0.5
    assert noncore["core_provenance_valid"] is True

    core = audit_manifests(
        [manifest],
        expected_jobs=[
            {"rule_name": "prepare_reference", "job_key": "hg002"},
            {"rule_name": "compute_pgbench_score", "job_key": "hg002", "core": True},
        ],
        target_manifest_ids=[manifest["manifest_id"]],
        workspace_root=tmp_path,
        require_companions=True,
    )
    assert core["status"] == "invalid"
    assert core["core_provenance_valid"] is False


def test_legacy_validate_profile_drift_is_superseded_by_frozen_context(
    tmp_path: Path,
) -> None:
    profile = tmp_path / "config" / "consensus_scoring.yaml"
    profile.parent.mkdir(parents=True)
    profile.write_text("version: old\n", encoding="utf-8")
    validated = tmp_path / "validated.json"
    validated.write_text("{}\n", encoding="utf-8")
    validate = _prepared_manifest(
        rule_name="validate_config",
        job_key="config",
        attempt_id="validate",
        inputs=[profile],
        outputs=[validated],
    )

    profile.write_text("version: frozen\n", encoding="utf-8")
    context = tmp_path / "run-context.json"
    context.write_text("{}\n", encoding="utf-8")
    snapshot = _prepared_manifest(
        rule_name="snapshot_run_context",
        job_key="context",
        attempt_id="snapshot",
        inputs=[validated, profile],
        outputs=[context],
        upstream_manifest_ids=[validate["manifest_id"]],
    )

    audit = audit_manifests(
        [validate, snapshot],
        expected_jobs=[
            {"rule_name": "validate_config", "job_key": "config", "core": True},
            {
                "rule_name": "snapshot_run_context",
                "job_key": "context",
                "core": True,
            },
        ],
        target_manifest_ids=[snapshot["manifest_id"]],
        workspace_root=tmp_path,
        verify_paths=True,
    )

    assert audit["core_provenance_valid"] is True
    assert any(
        issue["code"] == "superseded_validation_profile"
        for issue in audit["issues"]
    )


def test_duplicate_manifest_identity_invalidates_audit(tmp_path: Path) -> None:
    result = tmp_path / "result"
    result.write_text("ok", encoding="utf-8")
    manifest = _prepared_manifest(
        rule_name="prepare_reference",
        job_key="hg002",
        attempt_id="attempt-1",
        inputs=[],
        outputs=[result],
    )
    _write_companions(tmp_path, manifest)
    audit = audit_manifests(
        [manifest, dict(manifest)],
        expected_jobs=[
            {
                "rule_name": "prepare_reference",
                "job_key": "hg002",
                "manifest_id": manifest["manifest_id"],
            }
        ],
        target_manifest_ids=[manifest["manifest_id"]],
        workspace_root=tmp_path,
        require_companions=True,
    )
    assert audit["status"] == "invalid"
    assert audit["core_provenance_valid"] is False
    assert any(issue["code"] == "duplicate_manifest_id" for issue in audit["issues"])


def test_audit_requires_one_consistent_frozen_run_context(tmp_path: Path) -> None:
    intermediate = tmp_path / "intermediate"
    final = tmp_path / "final"
    intermediate.write_text("one", encoding="utf-8")
    final.write_text("two", encoding="utf-8")
    upstream = _prepared_manifest(
        rule_name="canonicalize_vcf",
        job_key="hg002",
        attempt_id="one",
        inputs=[],
        outputs=[intermediate],
    )
    downstream = _prepared_manifest(
        rule_name="compute_pgbench_score",
        job_key="hg002",
        attempt_id="two",
        inputs=[intermediate],
        outputs=[final],
        upstream_manifest_ids=[upstream["manifest_id"]],
    )
    downstream["hardware_fingerprint_sha256"] = SHA_A
    _write_companions(tmp_path, upstream)
    _write_companions(tmp_path, downstream)

    audit = audit_manifests(
        [upstream, downstream],
        expected_jobs=[
            {"rule_name": "canonicalize_vcf", "job_key": "hg002"},
            {"rule_name": "compute_pgbench_score", "job_key": "hg002"},
        ],
        target_manifest_ids=[downstream["manifest_id"]],
        workspace_root=tmp_path,
        require_companions=True,
    )

    assert audit["run_context_complete"] is False
    assert audit["status"] == "invalid"
    assert audit["core_provenance_valid"] is False
    assert any(issue["code"] == "run_context_mismatch" for issue in audit["issues"])


def test_audit_requires_nonempty_snakemake_version(tmp_path: Path) -> None:
    result = tmp_path / "result"
    result.write_text("ok", encoding="utf-8")
    manifest = _prepared_manifest(
        rule_name="canonicalize_vcf",
        job_key="hg002",
        attempt_id="attempt-1",
        inputs=[],
        outputs=[result],
    )
    manifest["snakemake_version"] = None
    manifest["not_applicable_reason"]["snakemake_version"] = "test fixture"
    _write_companions(tmp_path, manifest)

    audit = audit_manifests(
        [manifest],
        expected_jobs=[
            {"rule_name": "canonicalize_vcf", "job_key": "hg002", "core": True}
        ],
        target_manifest_ids=[manifest["manifest_id"]],
        workspace_root=tmp_path,
        require_companions=True,
    )

    assert audit["run_context_complete"] is False
    assert audit["status"] == "provisional"


@pytest.mark.parametrize("mismatch", ["run_id", "score_profile_sha256"])
def test_audit_rejects_mixed_run_or_score_profile(
    tmp_path: Path, mismatch: str
) -> None:
    intermediate = tmp_path / "intermediate"
    final = tmp_path / "final"
    intermediate.write_text("one", encoding="utf-8")
    final.write_text("two", encoding="utf-8")
    upstream = _prepared_manifest(
        rule_name="canonicalize_vcf",
        job_key="hg002",
        attempt_id="one",
        inputs=[],
        outputs=[intermediate],
    )
    downstream_kwargs: dict[str, str] = {}
    if mismatch == "run_id":
        downstream_kwargs["run_id"] = "another-run"
    else:
        downstream_kwargs["score_profile_sha256"] = SHA_A
    downstream = _prepared_manifest(
        rule_name="compute_pgbench_score",
        job_key="hg002",
        attempt_id="two",
        inputs=[intermediate],
        outputs=[final],
        upstream_manifest_ids=[upstream["manifest_id"]],
        **downstream_kwargs,
    )
    _write_companions(tmp_path, upstream)
    _write_companions(tmp_path, downstream)

    audit = audit_manifests(
        [upstream, downstream],
        expected_jobs=[
            {"rule_name": "canonicalize_vcf", "job_key": "hg002"},
            {"rule_name": "compute_pgbench_score", "job_key": "hg002"},
        ],
        target_manifest_ids=[downstream["manifest_id"]],
        workspace_root=tmp_path,
        require_companions=True,
    )

    assert audit["run_context_complete"] is False
    assert audit["core_provenance_valid"] is False
    assert audit["status"] == "invalid"
    assert any(
        issue.get("code") == "run_context_mismatch" and issue.get("field") == mismatch
        for issue in audit["issues"]
    )
