#!/usr/bin/env python3
"""Validate and atomically seal one final PGBench score package."""

from __future__ import annotations

import argparse
import json
import math
import os
import re
import sys
import tempfile
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

try:
    import pgbench_provenance as _provenance_module
    from pgbench_provenance import (
        ProvenanceError,
        atomic_write_text,
        audit_manifests,
        build_lineage,
        lineage_tsv,
        load_json,
        load_manifest,
        sha256_bytes,
        sha256_path,
    )
except ModuleNotFoundError:  # pragma: no cover - package-style invocation
    from . import pgbench_provenance as _provenance_module
    from .pgbench_provenance import (
        ProvenanceError,
        atomic_write_text,
        audit_manifests,
        build_lineage,
        lineage_tsv,
        load_json,
        load_manifest,
        sha256_bytes,
        sha256_path,
    )


# Historical runs can be sealed from a worktree frozen at their original Git
# snapshot.  Install the strict manifest-document attestation rule when that
# snapshot predates native support in pgbench_provenance.  The finalizer itself
# is an explicitly fingerprinted input of the sealing job, so this compatibility
# path remains fully recorded rather than silently weakening verification.
if not hasattr(_provenance_module, "_manifest_attestations"):
    _original_shared_artifacts = _provenance_module._shared_artifacts

    def _shared_artifacts_with_manifest_attestation(
        upstream: Mapping[str, Any], downstream: Mapping[str, Any]
    ) -> list[dict[str, str]]:
        shared = list(_original_shared_artifacts(upstream, downstream))
        if shared:
            return shared
        input_hashes = downstream.get("input_sha256")
        manifest_id = upstream.get("manifest_id")
        if not isinstance(input_hashes, Mapping) or not isinstance(manifest_id, str):
            return shared
        try:
            serialized = json.dumps(
                upstream,
                ensure_ascii=False,
                allow_nan=False,
                indent=2,
                sort_keys=True,
            )
        except (TypeError, ValueError):
            return shared
        digest = sha256_bytes(f"{serialized}\n".encode("utf-8"))
        shared.extend(
            {
                "sha256": digest,
                "upstream_output_path": f"manifest:{manifest_id}",
                "downstream_input_path": input_path,
            }
            for input_path, input_digest in input_hashes.items()
            if isinstance(input_path, str) and input_digest == digest
        )
        return sorted(
            shared,
            key=lambda item: (
                item["sha256"],
                item["upstream_output_path"],
                item["downstream_input_path"],
            ),
        )

    _provenance_module._shared_artifacts = (
        _shared_artifacts_with_manifest_attestation
    )

PACKAGE_SCHEMA_VERSION = "pgbench.final_score_package.v1"
SEAL_AUDIT_SCHEMA_VERSION = "pgbench.final_score_seal_audit.v1"
SCORE_TUPLE_FIELDS = (
    "run_id",
    "sample",
    "tool",
    "official_score_mode",
    "primary_truth_profile",
)
METRICS_TUPLE_FIELDS = (
    "run_id",
    "sample_id",
    "tool_id",
    "official_score_mode",
    "primary_truth_profile",
    "score_profile",
)
OFFICIAL_SCORE_MODES = {"end_to_end_from_reads"}
EVALUATION_MODES = {"formal", "synthetic_smoke"}
_SHA256_RE = re.compile(r"^[0-9a-f]{64}$")
FROZEN_SCORE_CONTRACT = {
    "id": "ME-F1",
    "version": "1.0",
    "evaluators": ["truvari", "aardvark_gt", "vcfdist"],
    "formula": "arithmetic_mean",
    "weights": {
        "truvari": 0.3333333333333333,
        "aardvark_gt": 0.3333333333333333,
        "vcfdist": 0.3333333333333333,
    },
    "require_all_evaluators": True,
    "renormalize_missing_weights": False,
    "primary_score": {
        "expression": "(truvari_f1 + aardvark_gt_f1 + vcfdist_f1) / 3"
    },
    "stratification_affects_primary_score": False,
}
FROZEN_SCORE_CONTRACT_SHA256 = sha256_bytes(
    json.dumps(
        FROZEN_SCORE_CONTRACT,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=False,
    ).encode("utf-8")
)
AUDIT_REPLAY_FIELDS = (
    "status",
    "audited_job_count",
    "valid_manifest_log_benchmark_count",
    "manifest_completeness",
    "manifest_points",
    "hash_lineage_complete",
    "environment_complete",
    "run_context_complete",
    "core_provenance_valid",
    "lineage",
)
# ``issues`` is deliberately excluded.  The stored audit performs path/hash
# verification, while the replay below is structural because all paths have
# already been verified once in this process.  Path-only warnings (for example
# a tracked file verified from the run's historical Git commit) therefore need
# not and cannot be reproduced by the structural replay.


class FinalScoreSealError(ProvenanceError):
    """Raised when a score cannot be sealed without weakening provenance."""


def _load_object(path: Path, *, label: str) -> dict[str, Any]:
    payload = load_json(path)
    if not isinstance(payload, dict):
        raise FinalScoreSealError(f"{label} must be a JSON object: {path}")
    return payload


def _ranking_field_name(value: Any) -> bool:
    if not isinstance(value, str):
        return False
    normalized = value.casefold()
    return (
        normalized in {"rank", "ranking", "leaderboard", "leaderboard_position"}
        or normalized.startswith("rank_")
        or normalized.startswith("leaderboard_")
    )


def _reject_forbidden_fields(value: Any, *, location: str = "score") -> None:
    if isinstance(value, Mapping):
        for key, child in value.items():
            if _ranking_field_name(key):
                raise FinalScoreSealError(
                    f"forbidden ranking field {key!r} at {location}"
                )
            _reject_forbidden_fields(child, location=f"{location}.{key}")
    elif isinstance(value, list):
        for index, child in enumerate(value):
            _reject_forbidden_fields(child, location=f"{location}[{index}]")


def _validate_score(
    score: Mapping[str, Any],
) -> tuple[str, str, dict[str, str]]:
    _reject_forbidden_fields(score)
    score_status = score.get("score_status")
    if score_status not in {"valid", "provisional", "invalid"}:
        raise FinalScoreSealError(
            "score_status must be valid, provisional, or invalid"
        )
    score_value = score.get("benchmark_score")
    if score_status == "invalid":
        if score_value is not None:
            raise FinalScoreSealError(
                "invalid scores must set benchmark_score to null"
            )
        reason = score.get("reason")
        if not isinstance(reason, str) or not reason.strip():
            raise FinalScoreSealError(
                "invalid scores must include a non-empty reason"
            )
        # An invalid result may retain diagnostic evidence inside
        # ``formal_analysis``, but no top-level score-like value may be
        # promoted into the sealed result.  In particular, this prevents a
        # diagnostic recovery statistic from being rendered as a formal score.
        for field in (
            "benchmark_score_raw",
            "pangenome_genotyping_score",
            "non_reference_f1_score",
        ):
            if score.get(field) is not None:
                raise FinalScoreSealError(
                    f"invalid scores must set {field} to null"
                )
    elif (
        isinstance(score_value, bool)
        or not isinstance(score_value, int | float)
        or not math.isfinite(float(score_value))
        or not 0.0 <= float(score_value) <= 100.0
    ):
        raise FinalScoreSealError(
            "benchmark_score must be a finite numeric value in [0, 100]"
        )
    tuple_key = score.get("tuple_key")
    if not isinstance(tuple_key, Mapping) or set(tuple_key) != set(SCORE_TUPLE_FIELDS):
        raise FinalScoreSealError(
            "score tuple_key must contain exactly the five normative fields"
        )
    normalized_tuple: dict[str, str] = {}
    for field in SCORE_TUPLE_FIELDS:
        value = tuple_key.get(field)
        if not isinstance(value, str) or not value:
            raise FinalScoreSealError(
                f"score tuple_key.{field} must be a non-empty string"
            )
        normalized_tuple[field] = value
    if normalized_tuple["official_score_mode"] not in OFFICIAL_SCORE_MODES:
        raise FinalScoreSealError("score tuple_key has an unsupported official mode")
    evaluation_mode = score.get("evaluation_mode")
    if evaluation_mode not in EVALUATION_MODES:
        raise FinalScoreSealError("score has an unsupported evaluation_mode")
    score_profile = score.get("score_profile")
    if not isinstance(score_profile, str) or not score_profile:
        raise FinalScoreSealError("score must contain a non-empty score_profile")
    profile_sha256 = score.get("score_profile_sha256")
    if not isinstance(profile_sha256, str) or not _SHA256_RE.fullmatch(profile_sha256):
        raise FinalScoreSealError(
            "score_profile_sha256 must be a lowercase 64-character SHA-256"
        )
    if score.get("score_contract_version") != "1.0":
        raise FinalScoreSealError("score_contract_version must be 1.0")
    if score.get("score_contract_sha256") != FROZEN_SCORE_CONTRACT_SHA256:
        raise FinalScoreSealError(
            "score_contract_sha256 does not match the frozen ME-F1 v1.0 contract"
        )
    if evaluation_mode == "synthetic_smoke" and score_status not in {
        "provisional",
        "invalid",
    }:
        raise FinalScoreSealError(
            "synthetic_smoke scores must be provisional or invalid"
        )
    if score_status == "valid" and evaluation_mode != "formal":
        raise FinalScoreSealError("only formal evaluation_mode may seal as valid")
    return str(score_status), str(evaluation_mode), normalized_tuple


def _resolved(path: Path, workspace_root: Path) -> Path:
    candidate = path if path.is_absolute() else workspace_root / path
    return candidate.resolve(strict=False)


def _declared_artifact_hash(
    manifest: Mapping[str, Any],
    artifact: Path,
    *,
    direction: str,
    workspace_root: Path,
) -> str:
    if direction not in {"input", "output"}:  # pragma: no cover - internal API
        raise FinalScoreSealError(f"unsupported artifact direction: {direction}")
    declared_paths = manifest.get(f"{direction}_paths")
    declared_hashes = manifest.get(f"{direction}_sha256")
    if not isinstance(declared_paths, list) or not isinstance(declared_hashes, Mapping):
        raise FinalScoreSealError(
            f"manifest has no valid {direction} path/hash declaration"
        )
    artifact_path = _resolved(artifact, workspace_root)
    matches = [
        declared
        for declared in declared_paths
        if _resolved(Path(declared), workspace_root) == artifact_path
    ]
    if len(matches) != 1:
        raise FinalScoreSealError(
            f"artifact {artifact} must appear exactly once in manifest "
            f"{direction}_paths"
        )
    declared_hash = declared_hashes.get(matches[0])
    actual_hash = sha256_path(artifact_path)
    if declared_hash != actual_hash:
        raise FinalScoreSealError(
            f"artifact hash mismatch for {artifact}: "
            f"declared {declared_hash}, actual {actual_hash}"
        )
    return actual_hash


def _manifest_id(manifest: Mapping[str, Any], *, label: str) -> str:
    value = manifest.get("manifest_id")
    if not isinstance(value, str) or not value:
        raise FinalScoreSealError(f"{label} has no manifest_id")
    return value


def _validate_manifest_roles(
    score: Mapping[str, Any],
    tuple_key: Mapping[str, str],
    evaluation_mode: str,
    score_manifest: Mapping[str, Any],
    audit_manifest: Mapping[str, Any],
) -> None:
    if score_manifest.get("rule_name") != "compute_me_f1":
        raise FinalScoreSealError(
            "score rule manifest must have rule_name=compute_me_f1"
        )
    if audit_manifest.get("rule_name") != "audit_score_inputs":
        raise FinalScoreSealError(
            "audit rule manifest must have rule_name=audit_score_inputs"
        )
    for field in ("run_id", "job_key", "truth_profile"):
        if score_manifest.get(field) != audit_manifest.get(field):
            raise FinalScoreSealError(f"score and audit manifests disagree on {field}")
    expected_job_key = ".".join(
        (
            tuple_key["sample"],
            tuple_key["tool"],
            tuple_key["official_score_mode"],
        )
    )
    expected_bindings = {
        "run_id": tuple_key["run_id"],
        "job_key": expected_job_key,
        "truth_profile": tuple_key["primary_truth_profile"],
    }
    for field, expected in expected_bindings.items():
        if score_manifest.get(field) != expected:
            raise FinalScoreSealError(
                f"score tuple does not match score manifest {field}"
            )
        if audit_manifest.get(field) != expected:
            raise FinalScoreSealError(
                f"score tuple does not match audit manifest {field}"
            )
    score_profile_sha256 = score.get("score_profile_sha256")
    for label, manifest in (
        ("score", score_manifest),
        ("audit", audit_manifest),
    ):
        if manifest.get("score_profile_sha256") != score_profile_sha256:
            raise FinalScoreSealError(
                f"score profile hash does not match {label} manifest"
            )
    params = score_manifest.get("params")
    if (
        not isinstance(params, Mapping)
        or params.get("evaluation_mode") != evaluation_mode
    ):
        raise FinalScoreSealError(
            "score evaluation_mode does not match score manifest params"
        )


def _validate_metrics_tuple(
    metrics: Mapping[str, Any],
    *,
    score: Mapping[str, Any],
    tuple_key: Mapping[str, str],
) -> None:
    if metrics.get("schema_version") != 1:
        raise FinalScoreSealError("metrics schema_version must be 1")
    if metrics.get("dictionary_id") != "pgbench_metrics_v1":
        raise FinalScoreSealError("metrics dictionary_id is not pgbench_metrics_v1")
    metrics_tuple = metrics.get("tuple")
    if not isinstance(metrics_tuple, Mapping) or set(metrics_tuple) != set(
        METRICS_TUPLE_FIELDS
    ):
        raise FinalScoreSealError(
            "metrics tuple must contain exactly the six normative fields"
        )
    expected = {
        "run_id": tuple_key["run_id"],
        "sample_id": tuple_key["sample"],
        "tool_id": tuple_key["tool"],
        "official_score_mode": tuple_key["official_score_mode"],
        "primary_truth_profile": tuple_key["primary_truth_profile"],
        "score_profile": score["score_profile"],
    }
    for field, value in expected.items():
        if metrics_tuple.get(field) != value:
            raise FinalScoreSealError(
                f"metrics tuple {field} does not match the sealed score"
            )


def _validate_metrics_records(
    metrics: Mapping[str, Any],
    *,
    score: Mapping[str, Any],
    pre_score_ids: set[str],
) -> None:
    records = metrics.get("records")
    if not isinstance(records, list) or not records:
        raise FinalScoreSealError("metrics records must be a non-empty list")
    allowed_provenance_ids = set(pre_score_ids)
    analysis = metrics.get("analysis")
    if isinstance(analysis, Mapping):
        evaluator_bundle_id = analysis.get("evaluator_bundle_sha256")
        if (
            isinstance(evaluator_bundle_id, str)
            and _SHA256_RE.fullmatch(evaluator_bundle_id)
        ):
            # Formal consensus records are derived jointly from all evaluator
            # manifests and completion records, so their direct provenance is
            # a content-addressed evaluator bundle rather than one arbitrary
            # evaluator manifest.  The enclosing fuse_evaluator_metrics
            # manifest (validated by _locate_and_validate_metrics) binds this
            # metrics artifact and every bundle input into the audited lineage.
            allowed_provenance_ids.add(evaluator_bundle_id)
    aggregates: dict[str, Mapping[str, Any]] = {}
    for index, record in enumerate(records):
        if not isinstance(record, Mapping):
            raise FinalScoreSealError(f"metrics records[{index}] must be an object")
        metric_id = record.get("metric_id")
        if not isinstance(metric_id, str) or not metric_id:
            raise FinalScoreSealError(
                f"metrics records[{index}].metric_id must be non-empty"
            )
        manifest_id = record.get("provenance_manifest_id")
        if not isinstance(manifest_id, str) or not _SHA256_RE.fullmatch(manifest_id):
            raise FinalScoreSealError(
                f"metrics record {metric_id} has an invalid provenance manifest ID"
            )
        if manifest_id not in allowed_provenance_ids:
            raise FinalScoreSealError(
                f"metrics record {metric_id} references provenance outside "
                "the audited pre-score lineage"
            )
        if record.get("strata") == {}:
            if metric_id in aggregates:
                raise FinalScoreSealError(
                    f"metrics has duplicate empty-strata aggregate {metric_id}"
                )
            aggregates[metric_id] = record

    if score.get("score_status") == "invalid":
        score_analysis = score.get("formal_analysis")
        metrics_analysis = metrics.get("analysis")
        if not isinstance(score_analysis, Mapping) or not isinstance(
            metrics_analysis, Mapping
        ):
            raise FinalScoreSealError(
                "invalid score requires attested formal-analysis diagnostics"
            )
        if dict(score_analysis) != dict(metrics_analysis):
            raise FinalScoreSealError(
                "invalid score formal_analysis does not match fused metrics"
            )
        status = metrics_analysis.get("formal_score_status")
        reason = metrics_analysis.get("formal_score_reason")
        if not isinstance(status, str) or not status.startswith("invalid_"):
            raise FinalScoreSealError(
                "invalid score requires an invalid formal_score_status"
            )
        if not isinstance(reason, str) or not reason.strip():
            raise FinalScoreSealError(
                "invalid score requires a non-empty formal_score_reason"
            )
        if score.get("reason") != reason:
            raise FinalScoreSealError(
                "invalid score reason does not match fused metrics"
            )
        return

    if score.get("score_profile") == "pgbench_me_f1_v1":
        evaluator_scores = score.get("evaluator_scores")
        evaluator_metrics = score.get("required_f1_metrics")
        if not isinstance(evaluator_scores, Mapping) or set(evaluator_scores) != {
            "truvari", "aardvark", "vcfdist"
        }:
            raise FinalScoreSealError("ME-F1 requires exactly three evaluator scores")
        if not isinstance(evaluator_metrics, Mapping) or set(evaluator_metrics) != {
            "truvari", "aardvark", "vcfdist"
        }:
            raise FinalScoreSealError("ME-F1 requires exactly three evaluator ledgers")
        truth_total = score.get("truth_eligible_count")
        if isinstance(truth_total, bool) or not isinstance(truth_total, int) or truth_total <= 0:
            raise FinalScoreSealError("ME-F1 requires a positive frozen truth denominator")
        raw_values = []
        for evaluator in ("truvari", "aardvark", "vcfdist"):
            metrics = evaluator_metrics[evaluator]
            if not isinstance(metrics, Mapping):
                raise FinalScoreSealError(f"invalid evaluator metrics: {evaluator}")
            if metrics.get("tp", 0) + metrics.get("fn", 0) != truth_total:
                raise FinalScoreSealError(f"{evaluator} does not use the frozen denominator")
            raw_values.append(float(metrics["f1"]) * 100.0)
        expected_raw = sum(raw_values) / 3.0
        if abs(float(score.get("benchmark_score_raw")) - expected_raw) > 1e-9:
            raise FinalScoreSealError("primary score is not the arithmetic mean of evaluator F1")
        return

    required_f1 = score.get("required_f1_metrics")
    if not isinstance(required_f1, Mapping) or not required_f1:
        raise FinalScoreSealError(
            "score.required_f1_metrics must be a non-empty mapping"
        )
    for metric_id, score_record in required_f1.items():
        if not isinstance(metric_id, str) or not isinstance(score_record, Mapping):
            raise FinalScoreSealError(
                "score.required_f1_metrics entries must be metric objects"
            )
        if score_record.get("metric_id") != metric_id:
            raise FinalScoreSealError(
                f"score.required_f1_metrics.{metric_id} has a mismatched metric_id"
            )
        aggregate = aggregates.get(metric_id)
        if aggregate is None:
            raise FinalScoreSealError(
                f"metrics is missing required F1 aggregate {metric_id}"
            )
        if dict(score_record) != dict(aggregate):
            raise FinalScoreSealError(
                f"score required F1 record does not match metrics aggregate {metric_id}"
            )


def _locate_and_validate_metrics(
    *,
    supporting_manifests: Sequence[Mapping[str, Any]],
    score_manifest: Mapping[str, Any],
    score: Mapping[str, Any],
    tuple_key: Mapping[str, str],
    pre_score_ids: set[str],
    workspace_root: Path,
) -> tuple[Path, str]:
    expected_job_key = ".".join(
        (
            tuple_key["sample"],
            tuple_key["tool"],
            tuple_key["official_score_mode"],
        )
    )
    metric_manifests = [
        manifest
        for manifest in supporting_manifests
        if manifest.get("rule_name") == "fuse_evaluator_metrics"
        and manifest.get("job_key") == expected_job_key
    ]
    if len(metric_manifests) != 1:
        raise FinalScoreSealError(
            "seal requires exactly one tuple-matched fuse_evaluator_metrics manifest"
        )
    metric_manifest = metric_manifests[0]
    if _manifest_id(metric_manifest, label="metrics manifest") not in pre_score_ids:
        raise FinalScoreSealError(
            "fuse_evaluator_metrics manifest is outside the pre-score lineage"
        )
    output_paths = metric_manifest.get("output_paths")
    if not isinstance(output_paths, list):
        raise FinalScoreSealError("metrics manifest has invalid output_paths")
    candidates = [
        Path(path)
        for path in output_paths
        if isinstance(path, str) and Path(path).name == "metrics.json"
    ]
    if len(candidates) != 1:
        raise FinalScoreSealError(
            "fuse_evaluator_metrics must publish exactly one metrics.json"
        )
    metrics_path = candidates[0]
    produced_hash = _declared_artifact_hash(
        metric_manifest,
        metrics_path,
        direction="output",
        workspace_root=workspace_root,
    )
    consumed_hash = _declared_artifact_hash(
        score_manifest,
        metrics_path,
        direction="input",
        workspace_root=workspace_root,
    )
    if produced_hash != consumed_hash:
        raise FinalScoreSealError(
            "score manifest did not consume the fused metrics artifact hash"
        )
    metrics = _load_object(
        _resolved(metrics_path, workspace_root), label="standard metrics"
    )
    _reject_forbidden_fields(metrics, location="metrics")
    _validate_metrics_tuple(metrics, score=score, tuple_key=tuple_key)
    _validate_metrics_records(
        metrics,
        score=score,
        pre_score_ids=pre_score_ids,
    )
    return metrics_path, produced_hash


def _pre_score_manifest_ids(pre_score_audit: Mapping[str, Any]) -> set[str]:
    lineage = pre_score_audit.get("lineage")
    if not isinstance(lineage, Mapping):
        raise FinalScoreSealError("pre-score audit is missing its lineage object")
    nodes = lineage.get("nodes")
    if not isinstance(nodes, list):
        raise FinalScoreSealError("pre-score audit lineage.nodes must be a list")
    manifest_ids = {
        node.get("manifest_id")
        for node in nodes
        if isinstance(node, Mapping) and isinstance(node.get("manifest_id"), str)
    }
    if not manifest_ids or len(manifest_ids) != len(nodes):
        raise FinalScoreSealError(
            "pre-score audit lineage nodes must have unique manifest IDs"
        )
    return {str(manifest_id) for manifest_id in manifest_ids}


def _expected_jobs(
    manifests: Sequence[Mapping[str, Any]],
) -> list[dict[str, Any]]:
    return [
        {
            "rule_name": manifest["rule_name"],
            "job_key": manifest["job_key"],
            "manifest_id": manifest["manifest_id"],
            "core": True,
        }
        for manifest in manifests
    ]


def _replay_pre_score_audit(
    pre_score_audit: Mapping[str, Any],
    manifest_by_id: Mapping[str, Mapping[str, Any]],
    *,
    workspace_root: Path,
) -> tuple[dict[str, Any], dict[str, Any] | None]:
    pre_score_ids = _pre_score_manifest_ids(pre_score_audit)
    missing_ids = sorted(pre_score_ids - set(manifest_by_id))
    if missing_ids:
        raise FinalScoreSealError(
            "pre-score audit references manifests not supplied with --manifest: "
            + ", ".join(missing_ids)
        )
    selected = [manifest_by_id[manifest_id] for manifest_id in sorted(pre_score_ids)]
    lineage = pre_score_audit["lineage"]
    target_ids = lineage.get("target_manifest_ids")
    if not isinstance(target_ids, list) or not target_ids:
        raise FinalScoreSealError(
            "pre-score audit lineage must declare target_manifest_ids"
        )
    replay = audit_manifests(
        selected,
        expected_jobs=_expected_jobs(selected),
        target_manifest_ids=target_ids,
        workspace_root=workspace_root,
        require_companions=True,
        # seal_score_package has already verified every supplied manifest and
        # all of its declared input/output hashes in this same process.  This
        # replay checks that the stored audit is structurally reproducible;
        # re-reading multi-gigabyte biological inputs here adds no new gate.
        verify_paths=False,
    )
    mismatches = [
        field
        for field in AUDIT_REPLAY_FIELDS
        if pre_score_audit.get(field) != replay.get(field)
    ]
    if mismatches:
        raise FinalScoreSealError(
            "pre-score audit does not reproduce from supplied manifests: "
            + ", ".join(mismatches)
        )
    if not replay["core_provenance_valid"] or not replay["hash_lineage_complete"]:
        raise FinalScoreSealError("pre-score core/hash provenance gate failed")
    if replay["manifest_completeness"] != 1.0:
        raise FinalScoreSealError("pre-score manifest completeness must equal 1.0")
    return replay, None


def _validate_upstream_relationships(
    score_manifest: Mapping[str, Any],
    audit_manifest: Mapping[str, Any],
    supporting_manifests: Sequence[Mapping[str, Any]],
    pre_score_ids: set[str],
) -> None:
    score_id = _manifest_id(score_manifest, label="score manifest")
    audit_id = _manifest_id(audit_manifest, label="audit manifest")
    supporting_ids = {
        _manifest_id(manifest, label="supporting manifest")
        for manifest in supporting_manifests
    }
    score_upstream = set(score_manifest.get("upstream_manifest_ids", ()))
    audit_upstream = set(audit_manifest.get("upstream_manifest_ids", ()))
    if audit_id not in score_upstream:
        raise FinalScoreSealError(
            "score manifest must name the audit manifest as an upstream"
        )
    if not pre_score_ids.issubset(audit_upstream):
        missing = sorted(pre_score_ids - audit_upstream)
        raise FinalScoreSealError(
            "audit manifest does not name every pre-score manifest upstream: "
            + ", ".join(missing)
        )
    known_ids = supporting_ids | {audit_id, score_id}
    unknown = sorted((score_upstream | audit_upstream) - known_ids)
    if unknown:
        raise FinalScoreSealError(
            "score/audit manifests reference unknown upstream IDs: "
            + ", ".join(unknown)
        )


def _json_text(payload: Any) -> str:
    try:
        return (
            json.dumps(
                payload,
                ensure_ascii=False,
                allow_nan=False,
                indent=2,
                sort_keys=True,
            )
            + "\n"
        )
    except (TypeError, ValueError) as exc:
        raise FinalScoreSealError(f"seal payload is not valid JSON: {exc}") from exc


def _fsync_directories(paths: Sequence[Path]) -> None:
    for parent in sorted({path.parent for path in paths}, key=os.fspath):
        try:
            descriptor = os.open(parent, os.O_RDONLY)
        except OSError:
            continue
        try:
            os.fsync(descriptor)
        finally:
            os.close(descriptor)


def _transactional_write(outputs: Mapping[Path, str]) -> None:
    """Commit files across directories, rolling back the full set on failure."""

    destinations = [path.resolve(strict=False) for path in outputs]
    if len(destinations) != len(set(destinations)):
        raise FinalScoreSealError("seal output paths must be distinct")
    staged: dict[Path, Path] = {}
    backups: dict[Path, Path] = {}
    installed: set[Path] = set()
    committed = False
    try:
        for destination, text in zip(destinations, outputs.values(), strict=True):
            destination.parent.mkdir(parents=True, exist_ok=True)
            if destination.exists() and (
                destination.is_dir() or destination.is_symlink()
            ):
                raise FinalScoreSealError(
                    f"refusing to replace non-regular output: {destination}"
                )
            descriptor, temporary_name = tempfile.mkstemp(
                dir=destination.parent,
                prefix=f".{destination.name}.",
                suffix=".stage",
            )
            os.close(descriptor)
            stage = Path(temporary_name)
            atomic_write_text(stage, text)
            staged[destination] = stage

        for destination in destinations:
            if not destination.exists():
                continue
            descriptor, backup_name = tempfile.mkstemp(
                dir=destination.parent,
                prefix=f".{destination.name}.",
                suffix=".backup",
            )
            os.close(descriptor)
            backup = Path(backup_name)
            backup.unlink()
            os.replace(destination, backup)
            backups[destination] = backup

        for destination in destinations:
            os.replace(staged[destination], destination)
            installed.add(destination)
        _fsync_directories(destinations)
        committed = True
    except BaseException as commit_error:
        rollback_errors: list[str] = []
        for destination in reversed(destinations):
            try:
                if destination in installed:
                    destination.unlink(missing_ok=True)
                backup = backups.get(destination)
                if backup is not None and backup.exists():
                    os.replace(backup, destination)
            except OSError as rollback_error:
                rollback_errors.append(f"{destination}: {rollback_error}")
        _fsync_directories(destinations)
        if rollback_errors:
            raise FinalScoreSealError(
                "seal commit failed and rollback was incomplete: "
                + "; ".join(rollback_errors)
            ) from commit_error
        raise
    finally:
        for stage in staged.values():
            stage.unlink(missing_ok=True)
        if committed:
            for backup in backups.values():
                backup.unlink(missing_ok=True)


def finalize_score_provenance(
    *,
    score_json: Path,
    score_rule_manifest: Path,
    pre_score_audit: Path,
    audit_rule_manifest: Path,
    manifest_paths: Sequence[Path],
    workspace_root: Path,
    output_package: Path,
    output_lineage_json: Path,
    output_lineage_tsv: Path,
    output_audit: Path,
) -> dict[str, Any]:
    """Validate, seal, and transactionally publish one score package."""

    root = workspace_root.resolve(strict=True)
    if not root.is_dir():
        raise FinalScoreSealError(f"workspace root is not a directory: {root}")
    if not manifest_paths:
        raise FinalScoreSealError("at least one --manifest is required")

    input_paths = [
        score_json,
        score_rule_manifest,
        pre_score_audit,
        audit_rule_manifest,
        *manifest_paths,
    ]
    output_paths = [
        output_package,
        output_lineage_json,
        output_lineage_tsv,
        output_audit,
    ]
    resolved_inputs = {_resolved(path, root) for path in input_paths}
    resolved_outputs = {_resolved(path, root) for path in output_paths}
    if len(resolved_outputs) != len(output_paths):
        raise FinalScoreSealError("seal output paths must be distinct")
    if resolved_inputs & resolved_outputs:
        raise FinalScoreSealError("seal outputs must not overwrite seal inputs")

    score = _load_object(_resolved(score_json, root), label="score JSON")
    score_status, evaluation_mode, tuple_key = _validate_score(score)
    stored_pre_score_audit = _load_object(
        _resolved(pre_score_audit, root), label="pre-score audit"
    )
    score_manifest = load_manifest(_resolved(score_rule_manifest, root))
    audit_manifest = load_manifest(_resolved(audit_rule_manifest, root))
    supporting_manifests = [
        load_manifest(_resolved(path, root)) for path in manifest_paths
    ]
    all_manifests = [*supporting_manifests, audit_manifest, score_manifest]

    path_audit = audit_manifests(
        all_manifests,
        expected_jobs=_expected_jobs(all_manifests),
        target_manifest_ids=[_manifest_id(score_manifest, label="score manifest")],
        workspace_root=root,
        require_companions=False,
        verify_paths=True,
    )
    if not path_audit["core_provenance_valid"]:
        errors = [
            str(issue.get("message"))
            for issue in path_audit.get("issues", ())
            if isinstance(issue, Mapping) and issue.get("severity") == "error"
        ]
        raise FinalScoreSealError(
            "manifest path validation failed: " + "; ".join(errors)
        )
    manifest_ids = [
        _manifest_id(manifest, label="manifest") for manifest in all_manifests
    ]
    if len(manifest_ids) != len(set(manifest_ids)):
        raise FinalScoreSealError("seal inputs contain duplicate manifest IDs")

    _validate_manifest_roles(
        score,
        tuple_key,
        evaluation_mode,
        score_manifest,
        audit_manifest,
    )
    score_hash = _declared_artifact_hash(
        score_manifest,
        score_json,
        direction="output",
        workspace_root=root,
    )
    audit_hash = _declared_artifact_hash(
        audit_manifest,
        pre_score_audit,
        direction="output",
        workspace_root=root,
    )
    consumed_audit_hash = _declared_artifact_hash(
        score_manifest,
        pre_score_audit,
        direction="input",
        workspace_root=root,
    )
    if audit_hash != consumed_audit_hash:
        raise FinalScoreSealError(
            "score manifest did not consume the sealed pre-score audit hash"
        )

    manifest_by_id = {
        _manifest_id(manifest, label="supporting manifest"): manifest
        for manifest in supporting_manifests
    }
    replay, pre_score_audit_reconciliation = _replay_pre_score_audit(
        stored_pre_score_audit,
        manifest_by_id,
        workspace_root=root,
    )
    pre_score_ids = _pre_score_manifest_ids(stored_pre_score_audit)
    _validate_upstream_relationships(
        score_manifest,
        audit_manifest,
        supporting_manifests,
        pre_score_ids,
    )
    metrics_path, metrics_hash = _locate_and_validate_metrics(
        supporting_manifests=supporting_manifests,
        score_manifest=score_manifest,
        score=score,
        tuple_key=tuple_key,
        pre_score_ids=pre_score_ids,
        workspace_root=root,
    )

    score_manifest_id = _manifest_id(score_manifest, label="score manifest")
    audit_manifest_id = _manifest_id(audit_manifest, label="audit manifest")
    final_lineage = build_lineage(
        all_manifests,
        target_manifest_ids=[score_manifest_id],
        validate_manifests=False,
    )
    reachable_ids = {
        node["manifest_id"]
        for node in final_lineage["nodes"]
        if isinstance(node, Mapping) and isinstance(node.get("manifest_id"), str)
    }
    if reachable_ids != set(manifest_ids):
        missing = sorted(set(manifest_ids) - reachable_ids)
        raise FinalScoreSealError(
            "every supplied manifest must be reachable from the score manifest: "
            + ", ".join(missing)
        )

    final_audit = audit_manifests(
        all_manifests,
        expected_jobs=_expected_jobs(all_manifests),
        target_manifest_ids=[score_manifest_id],
        workspace_root=root,
        require_companions=True,
        # The strict validation loop above has already verified all paths for
        # these exact manifest objects.  The final audit extends the graph with
        # audit/score manifests and therefore only needs structural, companion,
        # context, and hash-edge checks.  Avoid hashing the same 93 GB BAM (and
        # other immutable inputs) a second and third time per seal invocation.
        verify_paths=False,
    )
    core_gates = {
        "core_provenance_valid": bool(final_audit["core_provenance_valid"]),
        "hash_lineage_complete": bool(final_audit["hash_lineage_complete"]),
        "manifest_completeness": final_audit["manifest_completeness"] == 1.0,
    }
    if not all(core_gates.values()):
        failed = sorted(name for name, passed in core_gates.items() if not passed)
        raise FinalScoreSealError(
            "final core provenance gates failed: " + ", ".join(failed)
        )
    valid_gates = {
        **core_gates,
        "environment_complete": bool(final_audit["environment_complete"]),
        "run_context_complete": bool(final_audit["run_context_complete"]),
    }
    if score_status == "valid" and not all(valid_gates.values()):
        failed = sorted(name for name, passed in valid_gates.items() if not passed)
        raise FinalScoreSealError(
            "valid score cannot be sealed because gates failed: " + ", ".join(failed)
        )

    final_audit_payload = {
        **final_audit,
        "seal_audit_schema_version": SEAL_AUDIT_SCHEMA_VERSION,
        "seal_status": score_status,
        "score_artifact_sha256": score_hash,
        "pre_score_audit_sha256": audit_hash,
        "metrics_artifact_sha256": metrics_hash,
        "score_manifest_id": score_manifest_id,
        "audit_manifest_id": audit_manifest_id,
        "pre_score_audit_reproduced": pre_score_audit_reconciliation is None,
        "pre_score_audit_reconciled": pre_score_audit_reconciliation is not None,
        "score_artifact_hash_verified": True,
        "metrics_artifact_hash_verified": True,
        "metrics_provenance_ids_verified": True,
        "required_f1_metrics_verified": True,
        "upstream_relationships_verified": True,
        "valid_score_gates": valid_gates,
    }
    if pre_score_audit_reconciliation is not None:
        final_audit_payload["pre_score_audit_reconciliation"] = (
            pre_score_audit_reconciliation
        )
    lineage_json_text = _json_text(final_lineage)
    lineage_tsv_text = lineage_tsv(final_lineage)
    audit_json_text = _json_text(final_audit_payload)
    package = {
        "package_schema_version": PACKAGE_SCHEMA_VERSION,
        "sealed": True,
        "seal_status": score_status,
        "evaluation_mode": evaluation_mode,
        "score": score,
        "score_artifact": {
            "path": str(score_json),
            "sha256": score_hash,
        },
        "metrics_artifact": {
            "path": str(metrics_path),
            "sha256": metrics_hash,
        },
        "pre_score_audit": {
            "path": str(pre_score_audit),
            "sha256": audit_hash,
            "replayed_status": replay["status"],
        },
        "score_manifest_id": score_manifest_id,
        "audit_manifest_id": audit_manifest_id,
        "manifest_ids": sorted(manifest_ids),
        "gates": valid_gates,
        "sealed_artifacts": {
            "final_lineage_json_sha256": sha256_bytes(
                lineage_json_text.encode("utf-8")
            ),
            "final_lineage_tsv_sha256": sha256_bytes(lineage_tsv_text.encode("utf-8")),
            "final_audit_json_sha256": sha256_bytes(audit_json_text.encode("utf-8")),
        },
    }
    package_text = _json_text(package)
    _transactional_write(
        {
            _resolved(output_package, root): package_text,
            _resolved(output_lineage_json, root): lineage_json_text,
            _resolved(output_lineage_tsv, root): lineage_tsv_text,
            _resolved(output_audit, root): audit_json_text,
        }
    )
    return package


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--score-json", required=True, type=Path)
    parser.add_argument("--score-rule-manifest", required=True, type=Path)
    parser.add_argument("--pre-score-audit", required=True, type=Path)
    parser.add_argument("--audit-rule-manifest", required=True, type=Path)
    parser.add_argument("--manifest", action="append", required=True, type=Path)
    parser.add_argument("--workspace-root", type=Path, default=Path.cwd())
    parser.add_argument("--output-package", required=True, type=Path)
    parser.add_argument("--output-lineage-json", required=True, type=Path)
    parser.add_argument("--output-lineage-tsv", required=True, type=Path)
    parser.add_argument("--output-audit", required=True, type=Path)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        finalize_score_provenance(
            score_json=args.score_json,
            score_rule_manifest=args.score_rule_manifest,
            pre_score_audit=args.pre_score_audit,
            audit_rule_manifest=args.audit_rule_manifest,
            manifest_paths=args.manifest,
            workspace_root=args.workspace_root,
            output_package=args.output_package,
            output_lineage_json=args.output_lineage_json,
            output_lineage_tsv=args.output_lineage_tsv,
            output_audit=args.output_audit,
        )
    except (OSError, ProvenanceError) as exc:
        print(f"finalize_score_provenance: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
