#!/usr/bin/env python3
"""Render a minimal no-ranking PGBench score report."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
from collections.abc import Mapping
from html import escape
from pathlib import Path
from typing import Any


class ReportError(ValueError):
    """Raised when a score artifact cannot be rendered safely."""


PACKAGE_SCHEMA_VERSION = "pgbench.final_score_package.v1"
FINALIZER_RULE_NAME = "finalize_score_provenance"


def _is_ranking_field(value: Any) -> bool:
    if not isinstance(value, str):
        return False
    normalized = value.strip().lower()
    return (
        normalized in {"rank", "ranking", "leaderboard", "leaderboard_position"}
        or normalized.startswith("rank_")
        or normalized.startswith("leaderboard_")
    )


def _reject_forbidden_fields(value: Any, path: str = "score") -> None:
    if isinstance(value, Mapping):
        for key, nested in value.items():
            if _is_ranking_field(key):
                raise ReportError(f"forbidden ranking field at {path}.{key}")
            _reject_forbidden_fields(nested, f"{path}.{key}")
    elif isinstance(value, list):
        for index, nested in enumerate(value):
            _reject_forbidden_fields(nested, f"{path}[{index}]")


def _load_score(path: Path) -> dict[str, Any]:
    try:
        with path.open("r", encoding="utf-8") as handle:
            score = json.load(handle)
    except (OSError, json.JSONDecodeError) as exc:
        raise ReportError(f"cannot load score JSON {path}: {exc}") from exc
    if not isinstance(score, dict):
        raise ReportError("score JSON must be an object")
    _reject_forbidden_fields(score)
    tuple_key = score.get("tuple_key")
    expected_tuple_fields = {
        "run_id",
        "sample",
        "tool",
        "official_score_mode",
        "primary_truth_profile",
    }
    if not isinstance(tuple_key, dict) or set(tuple_key) != expected_tuple_fields:
        raise ReportError("score JSON is missing tuple_key")
    return score


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _load_finalizer_manifest(path: Path) -> dict[str, Any]:
    try:
        with path.open("r", encoding="utf-8") as handle:
            manifest = json.load(handle)
    except (OSError, json.JSONDecodeError) as exc:
        raise ReportError(f"cannot load finalizer manifest {path}: {exc}") from exc
    if not isinstance(manifest, dict):
        raise ReportError("finalizer manifest must be an object")
    if manifest.get("rule_name") != FINALIZER_RULE_NAME:
        raise ReportError("manifest is not a finalize_score_provenance manifest")
    if manifest.get("status") != "success":
        raise ReportError("finalizer manifest is not successful")
    return manifest


def _manifest_artifact_hash(
    manifest: Mapping[str, Any],
    artifact_path: Path,
    *,
    direction: str,
    label: str,
) -> str:
    paths_field = f"{direction}_paths"
    hashes_field = f"{direction}_sha256"
    declared_paths = manifest.get(paths_field)
    declared_hashes = manifest.get(hashes_field)
    if not isinstance(declared_paths, list) or not all(
        isinstance(path, str) and path for path in declared_paths
    ):
        raise ReportError(f"finalizer manifest has invalid {paths_field}")
    if (
        len(declared_paths) != len(set(declared_paths))
        or not isinstance(declared_hashes, Mapping)
        or set(declared_hashes) != set(declared_paths)
    ):
        raise ReportError(
            f"finalizer manifest {hashes_field} keys must exactly match {paths_field}"
        )

    expected = artifact_path.resolve(strict=False)
    matching_paths = [
        path for path in declared_paths if Path(path).resolve(strict=False) == expected
    ]
    if len(matching_paths) != 1:
        raise ReportError(
            f"finalizer manifest must bind exactly one {label} {direction} path"
        )
    digest = declared_hashes.get(matching_paths[0])
    if (
        not isinstance(digest, str)
        or len(digest) != 64
        or any(character not in "0123456789abcdef" for character in digest)
    ):
        raise ReportError(f"finalizer manifest has invalid {label} {direction} hash")
    if _sha256_file(artifact_path) != digest:
        raise ReportError(
            f"{label} hash does not match the successful finalizer manifest"
        )
    return digest


def _validate_finalizer_tuple(
    manifest: Mapping[str, Any], score: Mapping[str, Any]
) -> None:
    tuple_key = score["tuple_key"]
    if manifest.get("run_id") != tuple_key["run_id"]:
        raise ReportError("finalizer manifest run_id does not match score tuple")
    wildcards = manifest.get("wildcards")
    if not isinstance(wildcards, Mapping) or wildcards.get("tool") != tuple_key["tool"]:
        raise ReportError("finalizer manifest tool does not match score tuple")
    expected_job_key = ".".join(
        (
            tuple_key["sample"],
            tuple_key["tool"],
            tuple_key["official_score_mode"],
        )
    )
    if manifest.get("job_key") != expected_job_key:
        raise ReportError("finalizer manifest job_key does not match score tuple")


def _load_sealed_score(
    score_path: Path,
    score_package_path: Path,
    finalizer_manifest_path: Path,
) -> dict[str, Any]:
    score = _load_score(score_path)
    try:
        with score_package_path.open("r", encoding="utf-8") as handle:
            package = json.load(handle)
    except (OSError, json.JSONDecodeError) as exc:
        raise ReportError(
            f"cannot load score package {score_package_path}: {exc}"
        ) from exc
    if not isinstance(package, dict):
        raise ReportError("score package must be an object")
    _reject_forbidden_fields(package, "score_package")
    if package.get("package_schema_version") != PACKAGE_SCHEMA_VERSION:
        raise ReportError("unsupported score package schema")
    if package.get("sealed") is not True:
        raise ReportError("score package is not sealed")
    if package.get("score") != score:
        raise ReportError("score package payload does not match score JSON")
    if package.get("seal_status") != score.get("score_status"):
        raise ReportError("score package status does not match score JSON")
    if package.get("evaluation_mode") != score.get("evaluation_mode"):
        raise ReportError("score package evaluation mode does not match score JSON")
    gates = package.get("gates")
    if not isinstance(gates, Mapping) or not all(
        gates.get(name) is True
        for name in (
            "core_provenance_valid",
            "hash_lineage_complete",
            "manifest_completeness",
        )
    ):
        raise ReportError("score package core provenance gates are not sealed")
    artifact = package.get("score_artifact")
    if not isinstance(artifact, Mapping):
        raise ReportError("score package is missing score_artifact")
    artifact_path = artifact.get("path")
    artifact_hash = artifact.get("sha256")
    if not isinstance(artifact_path, str) or not isinstance(artifact_hash, str):
        raise ReportError("score package has an invalid score_artifact")
    if Path(artifact_path).resolve(strict=False) != score_path.resolve(strict=False):
        raise ReportError("score package points to a different score artifact")
    if _sha256_file(score_path) != artifact_hash:
        raise ReportError("score JSON hash does not match the sealed package")

    finalizer_manifest = _load_finalizer_manifest(finalizer_manifest_path)
    _validate_finalizer_tuple(finalizer_manifest, score)
    manifest_package_hash = _manifest_artifact_hash(
        finalizer_manifest,
        score_package_path,
        direction="output",
        label="score package",
    )
    manifest_score_hash = _manifest_artifact_hash(
        finalizer_manifest,
        score_path,
        direction="input",
        label="score JSON",
    )
    if manifest_score_hash != artifact_hash:
        raise ReportError(
            "score package and finalizer manifest bind different score hashes"
        )
    if manifest_package_hash != _sha256_file(score_package_path):
        raise ReportError("score package hash is not bound by finalizer manifest")
    return score


def _atomic_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp")
    temporary.write_text(text, encoding="utf-8")
    temporary.replace(path)


def _write_score_tsv(path: Path, score: Mapping[str, Any]) -> None:
    tuple_key = score["tuple_key"]
    fields = [
        "run_id",
        "sample",
        "tool",
        "official_score_mode",
        "primary_truth_profile",
        "score_profile",
        "score_profile_sha256",
        "evaluation_mode",
        "score_status",
        "PGBenchConsensusScore",
        "PangenomeGenotypingScore",
        "NonReferenceF1",
        "PanelCoverage",
        "GlobalEndToEndSVRecoveryScore",
        "ComparableScore",
        "ComparableScore_CI95_lower",
        "ComparableScore_CI95_upper",
        "ConsensusScore",
        "truth_eligible_count",
        "query_result_count",
        "comparable_precision",
        "comparable_recall",
        "DiagnosticGlobalEndToEndSVRecoveryScore",
        "DiagnosticEligibleTruthCount",
        "DiagnosticQueryEventCount",
        "DiagnosticCreditedTruthEvents",
        "DiagnosticGlobalRecoveryStatus",
    ]
    analysis = score.get("formal_analysis")
    analysis = analysis if isinstance(analysis, Mapping) else {}
    interval = analysis.get("comparable_score_confidence_interval")
    interval = interval if isinstance(interval, Mapping) else {}
    matching = analysis.get("matching")
    matching = matching if isinstance(matching, Mapping) else {}
    invalid_formal_score = score.get("score_status") == "invalid"
    row = {
        **tuple_key,
        "score_profile": score["score_profile"],
        "score_profile_sha256": score["score_profile_sha256"],
        "evaluation_mode": score["evaluation_mode"],
        "score_status": score["score_status"],
        "PGBenchConsensusScore": score.get("pgbench_score"),
        "PangenomeGenotypingScore": score.get("pangenome_genotyping_score"),
        "NonReferenceF1": score.get("non_reference_f1_score"),
        "PanelCoverage": score.get("panel_coverage"),
        "GlobalEndToEndSVRecoveryScore": score.get(
            "global_end_to_end_sv_recovery_score"
        ),
        "ComparableScore": score.get("comparable_score", score.get("pgbench_score")),
        "ComparableScore_CI95_lower": interval.get("lower"),
        "ComparableScore_CI95_upper": interval.get("upper"),
        "ConsensusScore": score.get("consensus_score"),
        # A fail-closed score intentionally has no official denominator or
        # precision/recall.  Leave those TSV cells empty instead of showing
        # structural zero placeholders from the score object; diagnostics are
        # supplied in separately named, non-ranking fields below.
        "truth_eligible_count": (
            None if invalid_formal_score else score.get("truth_eligible_count")
        ),
        "query_result_count": (
            None if invalid_formal_score else score.get("total_evaluated")
        ),
        "comparable_precision": (
            None if invalid_formal_score else score.get("comparable_precision")
        ),
        "comparable_recall": (
            None if invalid_formal_score else score.get("comparable_recall")
        ),
        "DiagnosticGlobalEndToEndSVRecoveryScore": analysis.get(
            "global_end_to_end_sv_recovery_score"
        ),
        "DiagnosticEligibleTruthCount": matching.get("eligible_truth_event_count"),
        "DiagnosticQueryEventCount": matching.get("query_event_count"),
        "DiagnosticCreditedTruthEvents": matching.get("credited_truth_events"),
        "DiagnosticGlobalRecoveryStatus": analysis.get("global_recovery_status"),
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp")
    with temporary.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerow(row)
    temporary.replace(path)


def _write_breakdown_tsv(path: Path, score: Mapping[str, Any]) -> None:
    breakdown = score.get("point_breakdown")
    if not isinstance(breakdown, dict):
        raise ReportError("score JSON is missing point_breakdown")
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp")
    with temporary.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["category", "count"])
        for component, points in sorted(breakdown.items()):
            writer.writerow([component, int(points)])
    temporary.replace(path)


def _html(score: Mapping[str, Any]) -> str:
    tuple_key = score["tuple_key"]
    analysis = score.get("formal_analysis")
    analysis = analysis if isinstance(analysis, Mapping) else {}
    candidate = analysis.get("candidate_genotype_summary")
    candidate = candidate if isinstance(candidate, Mapping) else {}
    detection_only_contract = (
        candidate.get("candidate_output_contract") == "variant_sites"
    )
    not_applicable = "N/A — detection-only output contract"
    pgbench_value = (
        "not available"
        if score.get("pgbench_score") is None
        else f"{float(score['pgbench_score']):.2f}"
    )
    panel_value = (
        not_applicable
        if detection_only_contract
        else "not available"
        if score.get("pangenome_genotyping_score") is None
        else f"{float(score['pangenome_genotyping_score']):.2f}"
    )
    comparable_value = (
        "not available"
        if score.get("comparable_score", score.get("pgbench_score")) is None
        else f"{float(score.get('comparable_score', score.get('pgbench_score'))):.2f}"
    )
    matching = analysis.get("matching")
    matching = matching if isinstance(matching, Mapping) else {}
    diagnostic_global_score = analysis.get("global_end_to_end_sv_recovery_score")
    diagnostic_global_text = "not available"
    if isinstance(diagnostic_global_score, int | float) and not isinstance(
        diagnostic_global_score, bool
    ):
        diagnostic_global_text = f"{float(diagnostic_global_score):.2f}"
    diagnostic_truth_count = matching.get("eligible_truth_event_count")
    diagnostic_query_count = matching.get("query_event_count")
    diagnostic_credited_count = matching.get("credited_truth_events")
    formal_score_status = analysis.get("formal_score_status", "not available")
    invalid_formal_score = score.get("score_status") == "invalid"
    if invalid_formal_score:
        formal_count_text = (
            "正式计分计数不可用，因为结果已被 fail-closed 门控。"
            f"诊断性全局回收分母：eligible truth={diagnostic_truth_count if diagnostic_truth_count is not None else 'not available'}；"
            f"submitted events={diagnostic_query_count if diagnostic_query_count is not None else 'not available'}；"
            f"credited truth={diagnostic_credited_count if diagnostic_credited_count is not None else 'not available'}。"
        )
    else:
        formal_count_text = (
            f"固定 truth 数量：{score.get('truth_eligible_count', '')}；"
            f"检测事件数量：{score.get('total_evaluated', '')}；"
            f"comparable precision：{score.get('comparable_precision', '')}；"
            f"comparable recall：{score.get('comparable_recall', '')}。"
        )
    consensus_value = (
        "not available"
        if score.get("consensus_score") is None
        else f"{float(score['consensus_score']):.2f}"
    )
    point_rows = "\n".join(
        f"<tr><td>{escape(component)}</td><td>{float(points):.4f}</td></tr>"
        for component, points in sorted(score["point_breakdown"].items())
    )
    interval = analysis.get("comparable_score_confidence_interval")
    interval_text = "not available"
    if isinstance(interval, Mapping):
        interval_text = (
            f"95% CI {float(interval['lower']):.2f}–"
            f"{float(interval['upper']):.2f} "
            f"({int(interval['replicates'])} bootstrap replicates)"
        )
    semantics = analysis.get("semantic_summary")
    semantic_rows = ""
    if isinstance(semantics, Mapping):
        semantic_rows = "\n".join(
            "<tr>"
            f"<td>{escape(str(evaluator))}</td>"
            f"<td>{int(values.get('detection_events', 0))}</td>"
            f"<td>{int(values.get('detection_correct', 0))}</td>"
            f"<td>{int(values.get('genotype_scorable', 0))}</td>"
            f"<td>{int(values.get('genotype_correct', 0))}</td>"
            f"<td>{int(values.get('no_call', 0))}</td>"
            "</tr>"
            for evaluator, values in sorted(semantics.items())
            if isinstance(values, Mapping)
        )
    stratified = analysis.get("stratified_summary")
    strata_rows = ""
    if isinstance(stratified, Mapping):
        strata_rows = "\n".join(
            "<tr>"
            f"<td>{escape(str(dimension))}</td>"
            f"<td>{escape(str(entry.get('stratum', '')))}</td>"
            f"<td>{int(entry.get('truth_events', 0))}</td>"
            f"<td>{int(entry.get('query_events', 0))}</td>"
            f"<td>{float(entry.get('comparable_score', 0.0)):.2f}</td>"
            "</tr>"
            for dimension, entries in sorted(stratified.items())
            if isinstance(entries, list)
            for entry in entries
            if isinstance(entry, Mapping)
        )
    resources = analysis.get("resource_summary")
    resource_text = "not available"
    if isinstance(resources, Mapping):
        median = resources.get("median")
        median = median if isinstance(median, Mapping) else {}
        resource_text = (
            f"policy={escape(str(resources.get('cache_policy', '')))}, "
            f"repeats={int(resources.get('repeat_count', 0))}, "
            f"median wall={median.get('wall_seconds', 'n/a')} s, "
            f"median RSS={median.get('max_rss_mb', 'n/a')} MB"
        )
    candidate_text = "not available"
    if isinstance(candidate, Mapping):
        accuracy = candidate.get("genotype_accuracy")
        accuracy_text = (
            "undefined" if accuracy is None else f"{float(accuracy):.4f}"
        )
        called_accuracy = candidate.get("called_only_genotype_accuracy")
        called_accuracy_text = (
            "undefined"
            if called_accuracy is None
            else f"{float(called_accuracy):.4f}"
        )
        balanced_accuracy = candidate.get("balanced_accuracy")
        balanced_accuracy_text = (
            "undefined"
            if balanced_accuracy is None
            else f"{float(balanced_accuracy):.4f}"
        )
        macro_f1 = candidate.get("genotype_macro_f1")
        macro_f1_text = (
            "undefined" if macro_f1 is None else f"{float(macro_f1):.4f}"
        )
        candidate_text = (
            f"contract={escape(str(candidate.get('candidate_output_contract', '')))}; "
            f"absence={escape(str(candidate.get('absence_semantics', '')))}; "
            f"sites={int(candidate.get('candidate_count', 0))}; "
            f"truth-scorable={int(candidate.get('truth_scorable', 0))}; "
            f"observed={int(candidate.get('observed_candidates', 0))}; "
            f"GT scorable={int(candidate.get('genotype_scorable', 0))}; "
            f"called={int(candidate.get('called_candidates', 0))}; "
            f"GT correct={int(candidate.get('genotype_correct', 0))}; "
            f"GT accuracy={accuracy_text}; "
            f"called-only accuracy={called_accuracy_text}; "
            f"balanced accuracy={balanced_accuracy_text}; "
            f"GT macro-F1={macro_f1_text}; "
            f"no-call={int(candidate.get('no_call', 0))}; "
            f"direct/linked={int(candidate.get('direct_candidate_records', 0))}/"
            f"{int(candidate.get('linked_candidate_records', 0))}; "
            f"link conflicts={int(candidate.get('candidate_link_conflicts', 0))}"
        )
    nonref_value = (
        not_applicable
        if detection_only_contract
        else "not available"
        if score.get("non_reference_f1_score") is None
        else f"{float(score['non_reference_f1_score']):.2f}"
    )
    panel_coverage = score.get("panel_coverage")
    panel_coverage_text = (
        "not available"
        if panel_coverage is None
        else f"{float(panel_coverage) * 100.0:.2f}%"
    )
    no_call_count = (
        int(candidate.get("no_call", 0)) if isinstance(candidate, Mapping) else 0
    )
    mapping = analysis.get("evaluator_event_mapping")
    mapping_rows = ""
    if isinstance(mapping, Mapping):
        mapping_rows = "\n".join(
            "<tr>"
            f"<td>{escape(str(name))}</td>"
            f"<td>{int(values.get('resolved_events', 0))}</td>"
            f"<td>{int(values.get('unresolved_events', 0))}</td>"
            f"<td>{float(values.get('mapping_coverage', 0.0)) * 100.0:.2f}%</td>"
            "</tr>"
            for name, values in sorted(mapping.items())
            if isinstance(values, Mapping)
        )
    native = analysis.get("evaluator_native_metrics")
    native_rows = ""
    if isinstance(native, Mapping):
        native_rows = "\n".join(
            "<tr>"
            f"<td>{escape(str(name))}</td>"
            f"<td>{escape(str(values.get('unit', '')))}</td>"
            f"<td>{float(values.get('precision', 0.0)):.4f}</td>"
            f"<td>{float(values.get('recall', 0.0)):.4f}</td>"
            f"<td>{float(values.get('f1', 0.0)):.4f}</td>"
            "</tr>"
            for name, values in sorted(native.items())
            if isinstance(values, Mapping)
        )
    comparison_track = analysis.get("comparison_track")
    comparison_track_id = (
        comparison_track.get("id", "not available")
        if isinstance(comparison_track, Mapping)
        else "not available"
    )
    asset_hashes = analysis.get("asset_hashes")
    asset_rows = ""
    if isinstance(asset_hashes, Mapping):
        asset_rows = "\n".join(
            "<tr>"
            f"<td>{escape(str(name))}</td>"
            f"<td><code>{escape(str(value))}</code></td>"
            "</tr>"
            for name, value in sorted(asset_hashes.items())
        )
    return f"""<!doctype html>
<html lang="zh-CN">
<head>
  <meta charset="utf-8">
  <title>PGBench pangenome SV score card: {escape(str(tuple_key["tool"]))}</title>
  <style>
    body {{ font-family: system-ui, sans-serif; margin: 2rem; max-width: 72rem; }}
    table {{ border-collapse: collapse; width: 100%; }}
    th, td {{ border: 1px solid #ccd; padding: 0.45rem; text-align: left; }}
    .score {{ font-size: 2rem; font-weight: 700; }}
    .notice {{ background: #f3f5ff; padding: 0.8rem; border-left: 4px solid #536dfe; }}
  </style>
</head>
<body>
  <h1>PGBench 泛基因组结构变异工具得分卡</h1>
  <p class="notice">PGBench Consensus Score 对所有工具采用同一公式：三套冻结
  评测器各投一张等权的检测票。每个 truth 事件最多只记一次；全基因组召回、
  基因型与 no-call 分开报告，不使用资源权重，也不生成排名。</p>
  <dl>
    <dt>Run</dt><dd>{escape(str(tuple_key["run_id"]))}</dd>
    <dt>Sample</dt><dd>{escape(str(tuple_key["sample"]))}</dd>
    <dt>Tool</dt><dd>{escape(str(tuple_key["tool"]))}</dd>
    <dt>Official mode</dt><dd>{escape(str(tuple_key["official_score_mode"]))}</dd>
    <dt>Truth profile</dt><dd>{escape(str(tuple_key["primary_truth_profile"]))}</dd>
    <dt>Score profile</dt><dd>{escape(str(score.get("score_profile", "not available")))}</dd>
    <dt>Score profile hash</dt><dd><code>{escape(str(score.get("score_profile_sha256", "not available")))}</code></dd>
    <dt>Comparison track</dt><dd>{escape(str(comparison_track_id))}</dd>
    <dt>Comparison track hash</dt><dd><code>{escape(str(analysis.get("comparison_track_sha256", "not available")))}</code></dd>
    <dt>Evaluation mode</dt><dd>{escape(str(score["evaluation_mode"]))}</dd>
    <dt>Status</dt><dd>{escape(str(score["score_status"]))}</dd>
  </dl>
  <p><strong>Formal score gate:</strong> {escape(str(formal_score_status))};
  {escape(str(analysis.get('formal_score_reason') or 'no blocking reason recorded'))}.</p>
  <p class="score">PGBench Consensus Score: {pgbench_value} / 100</p>
  <p>This universal primary score is the mean of three frozen evaluators'
  binary detection votes on submitted in-scope events.</p>
  <p><strong>Pangenome Genotyping Score:</strong> {panel_value} / 100.
  This secondary score is genotype macro-F1 on the frozen blinded panel;
  it is not applicable to discovery-only output contracts.</p>
  <p><strong>Non-reference F1:</strong> {nonref_value} / 100;
  <strong>Panel coverage:</strong> {panel_coverage_text};
  <strong>No-call:</strong> {no_call_count}.</p>
  <p><strong>Global End-to-End SV Recovery Score (ComparableScore):</strong>
  {comparable_value} / 100. This secondary score measures recovery against the
  full eligible HG002 truth set and is not the primary panel-genotyping score.</p>
  <p><strong>Diagnostic Global Recovery:</strong> {diagnostic_global_text} / 100;
  status={escape(str(analysis.get('global_recovery_status', 'not available')))}.
  This diagnostic is never promoted to a formal score when the formal gate is invalid.</p>
  <p>ComparableScore: {comparable_value} / 100 (legacy field name).</p>
  <p>{interval_text}</p>
  <p>ConsensusScore: {consensus_value} / 100</p>
  <p>{escape(formal_count_text)}</p>
  <p>All three correct: {score.get("consensus_counts", {}).get("all_three_correct", "")}; exactly two: {score.get("consensus_counts", {}).get("exactly_two_correct", "")}; exactly one: {score.get("consensus_counts", {}).get("exactly_one_correct", "")}; none: {score.get("consensus_counts", {}).get("none_correct", "")}.</p>
  <h2>Consensus counts</h2>
  <table>
    <thead><tr><th>Category</th><th>Count</th></tr></thead>
    <tbody>{point_rows}</tbody>
  </table>
  <h2>Evaluator-native metrics</h2>
  <p>These are the evaluator's own summary units and are not overwritten by
  event-ledger conversion.</p>
  <table>
    <thead><tr><th>Evaluator</th><th>Native unit</th><th>Precision</th><th>Recall</th><th>F1</th></tr></thead>
    <tbody>{native_rows}</tbody>
  </table>
  <h2>Event mapping audit</h2>
  <p>Global recovery status: {escape(str(analysis.get('global_recovery_status', 'not available')))}.</p>
  <table>
    <thead><tr><th>Evaluator</th><th>Resolved</th><th>Unresolved</th><th>Coverage</th></tr></thead>
    <tbody>{mapping_rows}</tbody>
  </table>
  <h2>Hidden candidate-site genotype diagnostics</h2>
  <p>{candidate_text}</p>
  <h2>检测、GT 与 no-call 分项</h2>
  <table>
    <thead><tr><th>Evaluator</th><th>Detection events</th><th>Detection correct</th><th>GT scorable</th><th>GT correct</th><th>No-call</th></tr></thead>
    <tbody>{semantic_rows}</tbody>
  </table>
  <h2>分层结果</h2>
  <table>
    <thead><tr><th>Dimension</th><th>Stratum</th><th>Truth</th><th>Query</th><th>ComparableScore</th></tr></thead>
    <tbody>{strata_rows}</tbody>
  </table>
  <h2>资源测量</h2>
  <p>{resource_text}</p>
  <h2>冻结 benchmark 资源身份</h2>
  <table>
    <thead><tr><th>Asset</th><th>SHA-256 / null</th></tr></thead>
    <tbody>{asset_rows}</tbody>
  </table>
  <p>Evaluator profile: {escape(str(analysis.get("evaluator_profile_id", "not available")))};
  hash: <code>{escape(str(analysis.get("evaluator_profile_sha256", "not available")))}</code>.</p>
</body>
</html>
"""


def render_report(
    score_json: Path,
    *,
    score_package: Path,
    finalizer_manifest: Path,
    html_output: Path,
    score_tsv: Path,
    point_breakdown_tsv: Path,
) -> None:
    score = _load_sealed_score(score_json, score_package, finalizer_manifest)
    _write_score_tsv(score_tsv, score)
    _write_breakdown_tsv(point_breakdown_tsv, score)
    _atomic_text(html_output, _html(score))


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Render a PGBench score card.")
    parser.add_argument("--score-json", required=True, type=Path)
    parser.add_argument("--score-package", required=True, type=Path)
    parser.add_argument("--finalizer-manifest", required=True, type=Path)
    parser.add_argument("--html", required=True, type=Path)
    parser.add_argument("--score-tsv", required=True, type=Path)
    parser.add_argument("--point-breakdown-tsv", required=True, type=Path)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        render_report(
            args.score_json,
            score_package=args.score_package,
            finalizer_manifest=args.finalizer_manifest,
            html_output=args.html,
            score_tsv=args.score_tsv,
            point_breakdown_tsv=args.point_breakdown_tsv,
        )
    except ReportError as exc:
        print(f"render_report: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
