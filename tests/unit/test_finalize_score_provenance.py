from __future__ import annotations

import json
import os
import sys
from pathlib import Path
from typing import Any

import pytest

SCRIPTS = Path(__file__).resolve().parents[2] / "workflow" / "scripts"
sys.path.insert(0, str(SCRIPTS))

import finalize_score_provenance as finalizer  # noqa: E402
from pgbench_provenance import (  # noqa: E402
    MANIFEST_SCHEMA_VERSION,
    atomic_write_json,
    audit_manifests,
    prepare_rule_manifest,
    sha256_bytes,
)

SHA_A = sha256_bytes(b"a")
SHA_B = sha256_bytes(b"b")
F1_METRIC_ID = "truvari.event.overall.f1"


def _manifest(
    *,
    rule_name: str,
    job_key: str,
    attempt_id: str,
    inputs: list[Path],
    outputs: list[Path],
    upstream_manifest_ids: list[str] | None = None,
    environment_complete: bool = True,
    run_id: str = "synthetic-run",
    truth_profile: str = "giab_hg002_grch38_v5_0q",
    score_profile_sha256: str = SHA_B,
    params: dict[str, Any] | None = None,
) -> dict[str, Any]:
    reasons = {
        "container_uri": "Conda environment selected",
        "container_digest": "Conda environment selected",
    }
    if not environment_complete:
        reasons["conda_lock_sha256"] = "Synthetic environment is not locked"
    payload = {
        "manifest_schema_version": MANIFEST_SCHEMA_VERSION,
        "attempt_id": attempt_id,
        "run_id": run_id,
        "rule_name": rule_name,
        "job_key": job_key,
        "wildcards": {"sample": "HG002", "tool": "example_genotyper"},
        "module_or_tool_id": "pgbench-core",
        "snakefile_path": "workflow/rules/provenance.smk",
        "rule_source_sha256": SHA_A,
        "script_or_wrapper_path": "workflow/scripts/mock.py",
        "script_or_wrapper_sha256": SHA_B,
        "command": ["python", "mock.py"],
        "params": params or {"score_profile": "pgbench_v1"},
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
        "truth_profile": truth_profile,
        "conda_lock_sha256": SHA_B if environment_complete else None,
        "container_uri": None,
        "container_digest": None,
        "git_head": "b" * 40,
        "git_dirty": False,
        "git_diff_sha256": SHA_A,
        "snakemake_version": "9.0.0",
        "execution_profile": "local",
        "hardware_fingerprint_sha256": SHA_B,
        "random_seed": 20260717,
        "upstream_manifest_ids": upstream_manifest_ids or [],
        "started_at": "2026-07-17T00:00:00+00:00",
        "finished_at": "2026-07-17T00:00:01+00:00",
        "exit_code": 0,
        "status": "success",
        "not_applicable_reason": reasons,
    }
    return prepare_rule_manifest(payload)


def _write_manifest(path: Path, manifest: dict[str, Any]) -> None:
    atomic_write_json(path, manifest)


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


def _metric_record(provenance_manifest_id: str) -> dict[str, Any]:
    return {
        "metric_id": F1_METRIC_ID,
        "value_type": "ratio",
        "value": 0.82,
        "status": "defined",
        "evaluator": "truvari",
        "numerator": 82,
        "denominator": 100,
        "eligible_count": 100,
        "universe_id": "unit-test",
        "undefined_reason": None,
        "parser_id": "unit-test-parser",
        "parser_source_field": "f1",
        "provenance_manifest_id": provenance_manifest_id,
        "strata": {},
        "confidence_interval": None,
    }


def _metrics_payload(record: dict[str, Any]) -> dict[str, Any]:
    return {
        "schema_version": 1,
        "dictionary_id": "pgbench_metrics_v1",
        "tuple": {
            "run_id": "synthetic-run",
            "sample_id": "HG002",
            "tool_id": "example_genotyper",
            "official_score_mode": "end_to_end_from_reads",
            "primary_truth_profile": "giab_hg002_grch38_v5_0q",
            "score_profile": "pgbench_v1",
        },
        "records": [record],
    }


def _score_payload(
    *,
    status: str,
    evaluation_mode: str,
    required_f1_record: dict[str, Any],
) -> dict[str, Any]:
    return {
        "tuple_key": {
            "run_id": "synthetic-run",
            "sample": "HG002",
            "tool": "example_genotyper",
            "official_score_mode": "end_to_end_from_reads",
            "primary_truth_profile": "giab_hg002_grch38_v5_0q",
        },
        "score_profile": "pgbench_v1",
        "score_profile_sha256": SHA_B,
        "evaluation_mode": evaluation_mode,
        "score_status": status,
        "pgbench_score": 78.06,
        "point_breakdown": {"traceability.environment_complete": 1.0},
        "required_f1_metrics": {
            F1_METRIC_ID: dict(required_f1_record),
        },
    }


def _seal_fixture(
    root: Path,
    *,
    score_status: str = "valid",
    evaluation_mode: str = "formal",
    environment_complete: bool = True,
    score_mutator: Any | None = None,
    metrics_mutator: Any | None = None,
    audit_mutator: Any | None = None,
    link_score_to_audit: bool = True,
    score_manifest_profile_sha256: str = SHA_B,
    audit_manifest_profile_sha256: str = SHA_B,
    run_scoped_companions: bool = False,
) -> dict[str, Any]:
    root.mkdir(parents=True, exist_ok=True)
    artifacts = root / "artifacts"
    manifests = root / "manifests"
    artifacts.mkdir()
    manifests.mkdir()

    linked = artifacts / "linked.vcf"
    linked.write_text("linked\n", encoding="utf-8")
    tuple_job_key = "HG002.example_genotyper.end_to_end_from_reads"
    ancestor = _manifest(
        rule_name="link_pangenome_alleles",
        job_key=tuple_job_key,
        attempt_id="ancestor-1",
        inputs=[],
        outputs=[linked],
        environment_complete=environment_complete,
    )
    ancestor_path = manifests / "ancestor.json"
    _write_manifest(ancestor_path, ancestor)
    _write_companions(root, ancestor, run_scoped=run_scoped_companions)

    metric_record = _metric_record(ancestor["manifest_id"])
    metrics_payload = _metrics_payload(metric_record)
    if metrics_mutator is not None:
        metrics_mutator(metrics_payload)
    metrics = artifacts / "metrics.json"
    atomic_write_json(metrics, metrics_payload)
    metrics_manifest = _manifest(
        rule_name="fuse_evaluator_metrics",
        job_key=tuple_job_key,
        attempt_id="metrics-1",
        inputs=[linked],
        outputs=[metrics],
        upstream_manifest_ids=[ancestor["manifest_id"]],
        environment_complete=environment_complete,
    )
    metrics_manifest_path = manifests / "metrics.json"
    _write_manifest(metrics_manifest_path, metrics_manifest)
    _write_companions(root, metrics_manifest, run_scoped=run_scoped_companions)

    pre_score_audit_payload = audit_manifests(
        [ancestor, metrics_manifest],
        expected_jobs=[
            {
                "rule_name": ancestor["rule_name"],
                "job_key": ancestor["job_key"],
                "manifest_id": ancestor["manifest_id"],
                "core": True,
            },
            {
                "rule_name": metrics_manifest["rule_name"],
                "job_key": metrics_manifest["job_key"],
                "manifest_id": metrics_manifest["manifest_id"],
                "core": True,
            },
        ],
        target_manifest_ids=[metrics_manifest["manifest_id"]],
        workspace_root=root,
        require_companions=True,
        verify_paths=True,
    )
    if audit_mutator is not None:
        audit_mutator(pre_score_audit_payload)
    pre_score_audit = artifacts / "pre-score-audit.json"
    atomic_write_json(pre_score_audit, pre_score_audit_payload)

    audit_manifest = _manifest(
        rule_name="audit_score_inputs",
        job_key=tuple_job_key,
        attempt_id="audit-1",
        inputs=[linked, metrics],
        outputs=[pre_score_audit],
        upstream_manifest_ids=[
            ancestor["manifest_id"],
            metrics_manifest["manifest_id"],
        ],
        environment_complete=environment_complete,
        score_profile_sha256=audit_manifest_profile_sha256,
    )
    audit_manifest_path = manifests / "audit.json"
    _write_manifest(audit_manifest_path, audit_manifest)
    _write_companions(root, audit_manifest, run_scoped=run_scoped_companions)

    score_payload = _score_payload(
        status=score_status,
        evaluation_mode=evaluation_mode,
        required_f1_record=dict(metrics_payload["records"][0]),
    )
    if score_mutator is not None:
        score_mutator(score_payload)
    score = artifacts / "score.json"
    atomic_write_json(score, score_payload)
    score_manifest = _manifest(
        rule_name="compute_pgbench_score",
        job_key=tuple_job_key,
        attempt_id="score-1",
        inputs=[pre_score_audit, metrics],
        outputs=[score],
        upstream_manifest_ids=(
            [audit_manifest["manifest_id"], metrics_manifest["manifest_id"]]
            if link_score_to_audit
            else [metrics_manifest["manifest_id"]]
        ),
        environment_complete=environment_complete,
        score_profile_sha256=score_manifest_profile_sha256,
        params={
            "score_profile": "config/score_weights.yaml",
            "evaluation_mode": evaluation_mode,
        },
    )
    score_manifest_path = manifests / "score.json"
    _write_manifest(score_manifest_path, score_manifest)
    _write_companions(root, score_manifest, run_scoped=run_scoped_companions)

    outputs = {
        "package": root / "published" / "package" / "score-package.json",
        "lineage_json": root / "published" / "lineage" / "lineage.json",
        "lineage_tsv": root / "published" / "tables" / "lineage.tsv",
        "audit": root / "published" / "audit" / "audit.json",
    }
    return {
        "score": score,
        "score_manifest": score_manifest_path,
        "pre_score_audit": pre_score_audit,
        "audit_manifest": audit_manifest_path,
        "metrics": metrics,
        "metrics_manifest": metrics_manifest_path,
        "manifests": [ancestor_path, metrics_manifest_path],
        "outputs": outputs,
    }


def _run_seal(root: Path, fixture: dict[str, Any]) -> int:
    arguments = [
        "--score-json",
        str(fixture["score"]),
        "--score-rule-manifest",
        str(fixture["score_manifest"]),
        "--pre-score-audit",
        str(fixture["pre_score_audit"]),
        "--audit-rule-manifest",
        str(fixture["audit_manifest"]),
        "--workspace-root",
        str(root),
        "--output-package",
        str(fixture["outputs"]["package"]),
        "--output-lineage-json",
        str(fixture["outputs"]["lineage_json"]),
        "--output-lineage-tsv",
        str(fixture["outputs"]["lineage_tsv"]),
        "--output-audit",
        str(fixture["outputs"]["audit"]),
    ]
    for manifest in fixture["manifests"]:
        arguments.extend(["--manifest", str(manifest)])
    return finalizer.main(arguments)


def _assert_no_outputs(fixture: dict[str, Any]) -> None:
    assert all(not path.exists() for path in fixture["outputs"].values())


def test_seals_valid_score_with_complete_final_lineage(tmp_path: Path) -> None:
    root = tmp_path / "workspace"
    fixture = _seal_fixture(root)

    assert _run_seal(root, fixture) == 0

    package = json.loads(fixture["outputs"]["package"].read_text(encoding="utf-8"))
    lineage = json.loads(fixture["outputs"]["lineage_json"].read_text(encoding="utf-8"))
    audit = json.loads(fixture["outputs"]["audit"].read_text(encoding="utf-8"))
    assert package["sealed"] is True
    assert package["seal_status"] == "valid"
    assert package["score"]["pgbench_score"] == 78.06
    assert package["metrics_artifact"]["path"] == str(fixture["metrics"])
    assert len(package["metrics_artifact"]["sha256"]) == 64
    assert all(package["gates"].values())
    assert lineage["status"] == "complete"
    assert lineage["hash_lineage_complete"] is True
    assert len(lineage["nodes"]) == 4
    assert audit["seal_status"] == "valid"
    assert audit["score_artifact_hash_verified"] is True
    assert audit["metrics_artifact_hash_verified"] is True
    assert audit["metrics_provenance_ids_verified"] is True
    assert audit["required_f1_metrics_verified"] is True
    assert "compute_pgbench_score" in fixture["outputs"]["lineage_tsv"].read_text(
        encoding="utf-8"
    )


def test_seals_score_with_run_id_scoped_companions(tmp_path: Path) -> None:
    root = tmp_path / "workspace"
    fixture = _seal_fixture(root, run_scoped_companions=True)

    assert _run_seal(root, fixture) == 0

    audit = json.loads(fixture["outputs"]["audit"].read_text(encoding="utf-8"))
    assert audit["core_provenance_valid"] is True
    assert audit["valid_manifest_log_benchmark_count"] == 4


def test_seals_synthetic_score_as_provisional_without_environment_lock(
    tmp_path: Path,
) -> None:
    root = tmp_path / "workspace"
    fixture = _seal_fixture(
        root,
        score_status="provisional",
        evaluation_mode="synthetic_smoke",
        environment_complete=False,
    )

    assert _run_seal(root, fixture) == 0

    package = json.loads(fixture["outputs"]["package"].read_text(encoding="utf-8"))
    assert package["seal_status"] == "provisional"
    assert package["gates"]["core_provenance_valid"] is True
    assert package["gates"]["hash_lineage_complete"] is True
    assert package["gates"]["environment_complete"] is False


def test_rejects_valid_score_when_environment_gate_is_incomplete(
    tmp_path: Path,
) -> None:
    root = tmp_path / "workspace"
    fixture = _seal_fixture(root, environment_complete=False)

    assert _run_seal(root, fixture) == 2
    _assert_no_outputs(fixture)


def test_rejects_nested_ranking_field_without_partial_outputs(
    tmp_path: Path,
) -> None:
    root = tmp_path / "workspace"
    fixture = _seal_fixture(
        root,
        score_mutator=lambda score: score.update(
            {"metadata": {"leaderboard_position": 1}}
        ),
    )

    assert _run_seal(root, fixture) == 2
    _assert_no_outputs(fixture)


@pytest.mark.parametrize(
    "field_name",
    ["ranking", "rank_percentile", "leaderboard_url"],
)
def test_rejects_all_scorer_ranking_field_variants(
    tmp_path: Path,
    field_name: str,
) -> None:
    root = tmp_path / field_name
    fixture = _seal_fixture(
        root,
        score_mutator=lambda score: score.update(
            {"metadata": {field_name: "forbidden"}}
        ),
    )

    assert _run_seal(root, fixture) == 2
    _assert_no_outputs(fixture)


@pytest.mark.parametrize(
    "tuple_mutator",
    [
        lambda value: value.pop("tool"),
        lambda value: value.update({"extra": "not-normative"}),
        lambda value: value.update({"sample": 2}),
        lambda value: value.update({"official_score_mode": "unsupported"}),
    ],
    ids=("missing-field", "extra-field", "non-string", "unsupported-mode"),
)
def test_rejects_non_normative_score_tuple(
    tmp_path: Path,
    tuple_mutator: Any,
) -> None:
    root = tmp_path / "workspace"

    def mutate_score(score: dict[str, Any]) -> None:
        tuple_mutator(score["tuple_key"])

    fixture = _seal_fixture(root, score_mutator=mutate_score)

    assert _run_seal(root, fixture) == 2
    _assert_no_outputs(fixture)


@pytest.mark.parametrize(
    ("field", "value"),
    [
        ("run_id", "another-run"),
        ("tool", "another-tool"),
        ("primary_truth_profile", "another-truth"),
    ],
)
def test_rejects_tuple_not_bound_to_score_and_audit_manifests(
    tmp_path: Path,
    field: str,
    value: str,
) -> None:
    root = tmp_path / "workspace"
    fixture = _seal_fixture(
        root,
        score_mutator=lambda score: score["tuple_key"].update({field: value}),
    )

    assert _run_seal(root, fixture) == 2
    _assert_no_outputs(fixture)


def test_rejects_score_profile_hash_not_bound_to_both_manifests(
    tmp_path: Path,
) -> None:
    score_root = tmp_path / "score-profile"
    score_fixture = _seal_fixture(
        score_root,
        score_mutator=lambda score: score.update({"score_profile_sha256": SHA_A}),
    )
    assert _run_seal(score_root, score_fixture) == 2
    _assert_no_outputs(score_fixture)

    audit_root = tmp_path / "audit-profile"
    audit_fixture = _seal_fixture(
        audit_root,
        audit_manifest_profile_sha256=SHA_A,
    )
    assert _run_seal(audit_root, audit_fixture) == 2
    _assert_no_outputs(audit_fixture)


def test_rejects_metric_provenance_outside_pre_score_lineage(
    tmp_path: Path,
) -> None:
    root = tmp_path / "workspace"
    fixture = _seal_fixture(
        root,
        metrics_mutator=lambda metrics: metrics["records"][0].update(
            {"provenance_manifest_id": "f" * 64}
        ),
    )

    assert _run_seal(root, fixture) == 2
    _assert_no_outputs(fixture)


def test_accepts_attested_evaluator_bundle_provenance(tmp_path: Path) -> None:
    root = tmp_path / "workspace"
    bundle_id = "f" * 64

    def attach_bundle(metrics: dict[str, Any]) -> None:
        metrics["analysis"] = {"evaluator_bundle_sha256": bundle_id}
        metrics["records"][0]["provenance_manifest_id"] = bundle_id

    fixture = _seal_fixture(root, metrics_mutator=attach_bundle)

    assert _run_seal(root, fixture) == 0


def test_rejects_required_f1_record_not_equal_to_metrics_aggregate(
    tmp_path: Path,
) -> None:
    root = tmp_path / "workspace"

    def mutate_score(score: dict[str, Any]) -> None:
        score["required_f1_metrics"][F1_METRIC_ID]["value"] = 0.81

    fixture = _seal_fixture(root, score_mutator=mutate_score)

    assert _run_seal(root, fixture) == 2
    _assert_no_outputs(fixture)


def test_rejects_metrics_tuple_mismatch_and_post_manifest_tampering(
    tmp_path: Path,
) -> None:
    tuple_root = tmp_path / "tuple-workspace"
    tuple_fixture = _seal_fixture(
        tuple_root,
        metrics_mutator=lambda metrics: metrics["tuple"].update(
            {"tool_id": "another-tool"}
        ),
    )
    assert _run_seal(tuple_root, tuple_fixture) == 2
    _assert_no_outputs(tuple_fixture)

    hash_root = tmp_path / "hash-workspace"
    hash_fixture = _seal_fixture(hash_root)
    metrics = json.loads(hash_fixture["metrics"].read_text(encoding="utf-8"))
    metrics["records"][0]["value"] = 0.81
    atomic_write_json(hash_fixture["metrics"], metrics)
    assert _run_seal(hash_root, hash_fixture) == 2
    _assert_no_outputs(hash_fixture)


def test_rejects_unreproduced_audit_and_missing_upstream_relationship(
    tmp_path: Path,
) -> None:
    audit_root = tmp_path / "audit-workspace"
    audit_fixture = _seal_fixture(
        audit_root,
        audit_mutator=lambda audit: audit.update(environment_complete=False),
    )
    assert _run_seal(audit_root, audit_fixture) == 2
    _assert_no_outputs(audit_fixture)

    upstream_root = tmp_path / "upstream-workspace"
    upstream_fixture = _seal_fixture(
        upstream_root,
        link_score_to_audit=False,
    )
    assert _run_seal(upstream_root, upstream_fixture) == 2
    _assert_no_outputs(upstream_fixture)


def test_rejects_score_hash_mismatch_and_missing_companion(tmp_path: Path) -> None:
    hash_root = tmp_path / "hash-workspace"
    hash_fixture = _seal_fixture(hash_root)
    score = json.loads(hash_fixture["score"].read_text(encoding="utf-8"))
    score["tampered_after_manifest"] = True
    atomic_write_json(hash_fixture["score"], score)
    assert _run_seal(hash_root, hash_fixture) == 2
    _assert_no_outputs(hash_fixture)

    companion_root = tmp_path / "companion-workspace"
    companion_fixture = _seal_fixture(companion_root)
    score_manifest = json.loads(
        companion_fixture["score_manifest"].read_text(encoding="utf-8")
    )
    score_log = (
        companion_root
        / "logs"
        / "rules"
        / score_manifest["rule_name"]
        / f"{score_manifest['job_key']}.log"
    )
    score_log.unlink()
    assert _run_seal(companion_root, companion_fixture) == 2
    _assert_no_outputs(companion_fixture)


def test_transaction_rolls_back_outputs_in_different_directories(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    outputs = {
        tmp_path / "a" / "one.json": "new-one\n",
        tmp_path / "b" / "two.json": "new-two\n",
        tmp_path / "c" / "three.tsv": "new-three\n",
        tmp_path / "d" / "four.json": "new-four\n",
    }
    for index, destination in enumerate(outputs, start=1):
        destination.parent.mkdir(parents=True)
        destination.write_text(f"old-{index}\n", encoding="utf-8")
    original = os.replace
    failed_destination = list(outputs)[1].resolve()

    def fail_second_install(source: Any, destination: Any) -> None:
        source_path = Path(source)
        destination_path = Path(destination).resolve(strict=False)
        if (
            source_path.name.endswith(".stage")
            and destination_path == failed_destination
        ):
            raise OSError("simulated multi-output commit failure")
        original(source, destination)

    monkeypatch.setattr(finalizer.os, "replace", fail_second_install)
    with pytest.raises(OSError, match="simulated"):
        finalizer._transactional_write(outputs)

    assert [path.read_text(encoding="utf-8") for path in outputs] == [
        "old-1\n",
        "old-2\n",
        "old-3\n",
        "old-4\n",
    ]
    assert not list(tmp_path.rglob("*.stage"))
    assert not list(tmp_path.rglob("*.backup"))
