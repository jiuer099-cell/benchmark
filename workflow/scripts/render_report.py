#!/usr/bin/env python3
"""Render one sealed ME-F1 score card without compatibility score fields."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
from collections.abc import Mapping, Sequence
from html import escape
from pathlib import Path
from typing import Any


class ReportError(ValueError):
    """Raised when a sealed score cannot be rendered safely."""


TUPLE_FIELDS = ("run_id", "sample", "tool", "official_score_mode", "primary_truth_profile")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_object(path: Path, label: str) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ReportError(f"{label} must contain an object")
    return value


def reject_ranking_fields(value: Any, path: str = "score") -> None:
    if isinstance(value, Mapping):
        for key, child in value.items():
            name = str(key).casefold()
            if name in {"rank", "ranking", "leaderboard"} or name.startswith(("rank_", "leaderboard_")):
                raise ReportError(f"forbidden ranking field at {path}.{key}")
            reject_ranking_fields(child, f"{path}.{key}")
    elif isinstance(value, list):
        for index, child in enumerate(value):
            reject_ranking_fields(child, f"{path}[{index}]")


def load_sealed_score(score_path: Path, package_path: Path, finalizer_path: Path) -> dict[str, Any]:
    score = load_object(score_path, "score")
    package = load_object(package_path, "score package")
    finalizer = load_object(finalizer_path, "finalizer manifest")
    reject_ranking_fields(score)
    if package.get("package_schema_version") != "pgbench.final_score_package.v1" or package.get("sealed") is not True:
        raise ReportError("score package is not sealed")
    if package.get("score") != score or package.get("seal_status") != score.get("score_status"):
        raise ReportError("score package does not match score JSON")
    artifact = package.get("score_artifact")
    if not isinstance(artifact, Mapping) or artifact.get("sha256") != sha256(score_path):
        raise ReportError("sealed score hash does not match score JSON")
    if finalizer.get("rule_name") != "finalize_score_provenance" or finalizer.get("status") != "success":
        raise ReportError("finalizer manifest is not successful")
    tuple_key = score.get("tuple_key")
    if not isinstance(tuple_key, Mapping) or set(tuple_key) != set(TUPLE_FIELDS):
        raise ReportError("score tuple is incomplete")
    expected_job = ".".join(str(tuple_key[field]) for field in ("sample", "tool", "official_score_mode"))
    if finalizer.get("run_id") != tuple_key["run_id"] or finalizer.get("job_key") != expected_job:
        raise ReportError("finalizer manifest does not match score tuple")
    output_hashes = finalizer.get("output_sha256")
    if not isinstance(output_hashes, Mapping) or sha256(package_path) not in output_hashes.values():
        raise ReportError("finalizer manifest does not bind the score package")
    return score


def atomic_text(path: Path, value: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp")
    temporary.write_text(value, encoding="utf-8")
    temporary.replace(path)


def evaluator_metrics(score: Mapping[str, Any]) -> Mapping[str, Any]:
    value = score.get("required_f1_metrics")
    return value if isinstance(value, Mapping) else {}


def write_score_tsv(path: Path, score: Mapping[str, Any]) -> None:
    metrics = evaluator_metrics(score)
    analysis = score.get("formal_analysis")
    analysis = analysis if isinstance(analysis, Mapping) else {}
    addressability = analysis.get("addressability_audit")
    addressability = addressability if isinstance(addressability, Mapping) else {}
    me_f1_ci = analysis.get("me_f1_confidence_interval")
    me_f1_ci = me_f1_ci if isinstance(me_f1_ci, Mapping) else {}
    fields = [
        *TUPLE_FIELDS, "ME-F1", "TruvariPrecision", "TruvariRecall", "TruvariF1",
        "AardvarkGTPrecision", "AardvarkGTRecall", "AardvarkGTF1",
        "vcfdistPrecision", "vcfdistRecall", "vcfdistF1", "EvaluatorRange",
        "score_profile", "score_profile_sha256", "score_contract_version",
        "score_contract_sha256", "evaluation_mode", "score_status",
        "EvaluatorSD", "truth_eligible_count", "query_result_count",
        "DiagnosticGTMacroF1", "DiagnosticNonReferenceF1", "DiagnosticPanelCoverage",
        "ME-F1_AddressableDiagnostic", "CallRate", "NoCallRate", "AddressabilityRate",
        "PanelTotal", "AddressableCount", "CalledCount", "ExplicitNoCallCount",
        "MissingOutputCount", "LinkingFailureCount",
        "UnsupportedRepresentationCount", "AdapterConversionFailureCount",
        "IndexBuildFailureCount", "AmbiguousMappingCount",
        "ME-F1_CI95_Lower", "ME-F1_CI95_Upper",
    ]

    def value(evaluator: str, field: str) -> Any:
        row = metrics.get(evaluator)
        return row.get(field) if isinstance(row, Mapping) else None

    row = {
        **score["tuple_key"], "score_profile": score.get("score_profile"),
        "score_profile_sha256": score.get("score_profile_sha256"),
        "score_contract_version": score.get("score_contract_version"),
        "score_contract_sha256": score.get("score_contract_sha256"),
        "evaluation_mode": score.get("evaluation_mode"), "score_status": score.get("score_status"),
        "ME-F1": score.get("benchmark_score"),
        "TruvariPrecision": value("truvari", "precision"), "TruvariRecall": value("truvari", "recall"), "TruvariF1": value("truvari", "f1"),
        "AardvarkGTPrecision": value("aardvark", "precision"), "AardvarkGTRecall": value("aardvark", "recall"), "AardvarkGTF1": value("aardvark", "f1"),
        "vcfdistPrecision": value("vcfdist", "precision"), "vcfdistRecall": value("vcfdist", "recall"), "vcfdistF1": value("vcfdist", "f1"),
        "EvaluatorRange": score.get("evaluator_range"), "EvaluatorSD": score.get("evaluator_sd"),
        "truth_eligible_count": score.get("truth_eligible_count"), "query_result_count": score.get("total_evaluated"),
        "DiagnosticGTMacroF1": score.get("pangenome_genotyping_score"),
        "DiagnosticNonReferenceF1": score.get("non_reference_f1_score"),
        "DiagnosticPanelCoverage": score.get("panel_coverage"),
        "AllSiteCallRate": candidate.get("all_site_call_rate"),
        "ExactGTAccuracy": candidate.get("exact_gt_accuracy"),
        "NoCallRate": candidate.get("no_call_rate"),
        "GTConfusionMatrix": json.dumps(
            candidate.get("genotype_confusion_matrix"), sort_keys=True
        ) if candidate.get("genotype_confusion_matrix") is not None else None,
        "ME-F1_AddressableDiagnostic": (
            analysis.get("addressable_subset_diagnostic", {}).get("me_f1")
            if isinstance(analysis.get("addressable_subset_diagnostic"), Mapping)
            else None
        ),
        "CallRate": (
            float(addressability.get("called_count", 0))
            / float(addressability["canonical_candidate_count"])
            if addressability.get("canonical_candidate_count")
            else None
        ),
        "NoCallRate": (
            float(addressability.get("explicit_no_call_count", 0))
            / float(addressability["canonical_candidate_count"])
            if addressability.get("canonical_candidate_count")
            else None
        ),
        "AddressabilityRate": addressability.get("addressability_rate"),
        "PanelTotal": addressability.get("canonical_candidate_count"),
        "AddressableCount": addressability.get("addressable_count"),
        "CalledCount": addressability.get("called_count"),
        "ExplicitNoCallCount": addressability.get("explicit_no_call_count"),
        "MissingOutputCount": addressability.get("missing_output_count"),
        "LinkingFailureCount": addressability.get("linking_failure_count"),
        "UnsupportedRepresentationCount": addressability.get("unsupported_representation_count"),
        "AdapterConversionFailureCount": addressability.get("adapter_conversion_failure_count"),
        "IndexBuildFailureCount": addressability.get("index_build_failure_count"),
        "AmbiguousMappingCount": addressability.get("ambiguous_mapping_count"),
        "ME-F1_CI95_Lower": me_f1_ci.get("lower"),
        "ME-F1_CI95_Upper": me_f1_ci.get("upper"),
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp")
    with temporary.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerow(row)
    temporary.replace(path)


def write_breakdown(path: Path, score: Mapping[str, Any]) -> None:
    values = score.get("point_breakdown")
    if not isinstance(values, Mapping):
        raise ReportError("score has no evaluator breakdown")
    lines = ["metric\tvalue"] + [f"{key}\t{float(value):.10g}" for key, value in sorted(values.items())]
    atomic_text(path, "\n".join(lines) + "\n")


def number(value: Any, digits: int = 2) -> str:
    return "N/A" if value is None else f"{float(value):.{digits}f}"


def genotype_strata_html(analysis: Mapping[str, Any]) -> str:
    summary = analysis.get("genotype_stratified_summary")
    if not isinstance(summary, Mapping):
        return ""
    rows = []
    for dimension, raw_strata in summary.items():
        if not isinstance(raw_strata, Mapping):
            continue
        for stratum, raw in raw_strata.items():
            if not isinstance(raw, Mapping):
                continue
            evaluator_metrics = raw.get("evaluator_metrics")
            evaluator_metrics = (
                evaluator_metrics if isinstance(evaluator_metrics, Mapping) else {}
            )
            evaluator_f1 = []
            for name in ("truvari", "aardvark", "vcfdist"):
                values = evaluator_metrics.get(name)
                evaluator_f1.append(
                    number(values.get("f1"), 4)
                    if isinstance(values, Mapping)
                    else "N/A"
                )
            rows.append(
                "<tr>"
                f"<td>{escape(str(raw.get('metric_id', '')))}</td>"
                f"<td>{escape(str(dimension))}</td>"
                f"<td>{escape(str(stratum))}</td>"
                f"<td>{escape(str(raw.get('status', '')))}</td>"
                f"<td>{raw.get('canonical_candidate_count', '')}</td>"
                f"<td>{raw.get('truth_positive_count', '')}</td>"
                f"<td>{number(raw.get('me_f1'))}</td><td>{evaluator_f1[0]}</td>"
                f"<td>{evaluator_f1[1]}</td><td>{evaluator_f1[2]}</td>"
                f"<td>{number(raw.get('addressability_rate'), 4)}</td>"
                "</tr>"
            )
    if not rows:
        return ""
    return (
        "<h2>Genotype-aware stratified results</h2>"
        "<table><thead><tr><th>Metric</th><th>Dimension</th><th>Stratum</th><th>Status</th>"
        "<th>Candidates</th><th>Truth+</th><th>ME-F1</th><th>Truvari F1</th>"
        "<th>Aardvark-GT F1</th><th>vcfdist F1</th>"
        "<th>Addressability</th></tr></thead><tbody>"
        + "".join(rows)
        + "</tbody></table>"
    )


def html(score: Mapping[str, Any]) -> str:
    tuple_key = score["tuple_key"]
    metrics = evaluator_metrics(score)
    analysis = score.get("formal_analysis")
    analysis = analysis if isinstance(analysis, Mapping) else {}
    candidate = analysis.get("candidate_genotype_summary")
    candidate = candidate if isinstance(candidate, Mapping) else {}
    addressability = analysis.get("addressability_audit")
    addressability = addressability if isinstance(addressability, Mapping) else {}
    me_f1_ci = analysis.get("me_f1_confidence_interval")
    me_f1_ci = me_f1_ci if isinstance(me_f1_ci, Mapping) else {}
    strata_html = genotype_strata_html(analysis)
    rows = []
    for evaluator, label in (("truvari", "Truvari"), ("aardvark", "Aardvark-GT"), ("vcfdist", "vcfdist")):
        values = metrics.get(evaluator)
        values = values if isinstance(values, Mapping) else {}
        rows.append(
            f"<tr><td>{label}</td><td>{values.get('tp', '')}</td><td>{values.get('fp', '')}</td><td>{values.get('fn', '')}</td>"
            f"<td>{number(values.get('precision'), 4)}</td><td>{number(values.get('recall'), 4)}</td><td>{number(values.get('f1'), 4)}</td></tr>"
        )
    reason = escape(str(score.get("reason") or ""))
    return f"""<!doctype html>
<html lang="zh-CN"><head><meta charset="utf-8"><title>PGBench ME-F1 — {escape(str(tuple_key['tool']))}</title>
<style>body{{font-family:system-ui,sans-serif;margin:2rem;max-width:75rem}}table{{border-collapse:collapse;width:100%}}th,td{{border:1px solid #ccd;padding:.45rem;text-align:left}}th{{background:#eef2ff}}.score{{font-size:2rem}}</style></head><body>
<h1>{escape(str(tuple_key['tool']))}</h1>
<p>同一短读长、同一冻结 canonical panel、同一完整 panel truth denominator 的 genotype-aware non-reference SV 评估；unsupported、missing、no-call 与 linking failure 均不会缩小主分分母。</p>
<p class="score"><strong>ME-F1: {number(score.get('benchmark_score'))}</strong></p>
<p>Paired genomic-block bootstrap 95% CI：{number(me_f1_ci.get('lower'))}–{number(me_f1_ci.get('upper'))}。</p>
<p>状态：{escape(str(score.get('score_status')))}；evaluator range：{number(score.get('evaluator_range'))}；evaluator SD：{number(score.get('evaluator_sd'))}。</p>
<p>{reason}</p>
<table><thead><tr><th>Evaluator</th><th>TP</th><th>FP</th><th>FN</th><th>Precision</th><th>Recall</th><th>F1</th></tr></thead><tbody>{''.join(rows)}</tbody></table>
<h2>All-sites / addressability diagnostics</h2>
<p>Truth denominator：{score.get('truth_eligible_count', '')}；candidate count：{candidate.get('candidate_count', '')}；called：{addressability.get('called_count', '')}；explicit no-call：{addressability.get('explicit_no_call_count', '')}；missing output：{addressability.get('missing_output_count', '')}；linking failure：{addressability.get('linking_failure_count', '')}。</p>
<p>Unsupported representation：{addressability.get('unsupported_representation_count', '')}；linking failure：{addressability.get('linking_failure_count', '')}；ambiguous mapping：{addressability.get('ambiguous_mapping_count', '')}；adapter conversion failure：{addressability.get('adapter_conversion_failure_count', '')}；index build failure：{addressability.get('index_build_failure_count', '')}；addressability rate：{number(addressability.get('addressability_rate'), 4)}；panel coverage：{number(score.get('panel_coverage'), 4)}。</p>
<p>All-site Call Rate（诊断）：{number(candidate.get('all_site_call_rate'), 4)}；Exact GT Accuracy（诊断）：{number(candidate.get('exact_gt_accuracy'), 4)}；No-call Rate（诊断）：{number(candidate.get('no_call_rate'), 4)}。这些字段不参与 ME-F1。</p>
<p>Candidate genotype macro-F1（诊断）：{number(score.get('pangenome_genotyping_score'))}；non-reference F1（诊断）：{number(score.get('non_reference_f1_score'))}。</p>
<details><summary>GT confusion matrix（诊断）</summary><pre>{escape(json.dumps(candidate.get('genotype_confusion_matrix', {}), indent=2, sort_keys=True))}</pre></details>
{strata_html}
<p>Score contract：<code>{escape(str(score.get('score_contract_version')))}</code>；SHA-256：<code>{escape(str(score.get('score_contract_sha256')))}</code></p>
</body></html>\n"""


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--score-json", required=True, type=Path)
    parser.add_argument("--score-package", required=True, type=Path)
    parser.add_argument("--finalizer-manifest", required=True, type=Path)
    parser.add_argument("--html", required=True, type=Path)
    parser.add_argument("--score-tsv", required=True, type=Path)
    parser.add_argument("--point-breakdown-tsv", required=True, type=Path)
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        score = load_sealed_score(args.score_json, args.score_package, args.finalizer_manifest)
        write_score_tsv(args.score_tsv, score)
        write_breakdown(args.point_breakdown_tsv, score)
        atomic_text(args.html, html(score))
    except (OSError, json.JSONDecodeError, ReportError, TypeError, ValueError) as exc:
        print(f"render_report: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
