#!/usr/bin/env python3
"""Combine compatible finalized runs into one unified HTML comparison."""

from __future__ import annotations

import argparse
import json
import re
import sys
from collections.abc import Mapping, Sequence
from html import escape
from pathlib import Path

import yaml  # type: ignore[import-untyped]


class SuiteReportError(ValueError):
    """Raised when score cards cannot be compared under one frozen contract."""


SHA256_RE = re.compile(r"^[0-9a-f]{64}$")
SHARED_ASSET_KEYS = (
    "reference",
    "truth_vcf",
    "benchmark_bed",
    "pangenome_manifest",
    "challenge_hidden_ledger",
)
ALL_ASSET_KEYS = (*SHARED_ASSET_KEYS, "graph_asset_lock")


def _manifest(path: Path) -> dict:
    value = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise SuiteReportError(f"tool manifest must be a mapping: {path}")
    return value


def _number(value: object) -> str:
    return "N/A" if value is None else f"{float(value):.2f}"


def _formal_analysis(score: Mapping) -> list[tuple[str, Mapping]]:
    analyses: list[tuple[str, Mapping]] = []
    for field in ("formal_analysis", "analysis"):
        value = score.get(field)
        if value is None:
            continue
        if not isinstance(value, Mapping):
            raise SuiteReportError(f"score field {field} must be a mapping")
        analyses.append((field, value))
    return analyses


def _evaluator_profile_hash(score: Mapping) -> str | None:
    observed: list[str] = []
    for field, analysis in _formal_analysis(score):
        for key in ("evaluator_profile_sha256", "evaluator_profile_hash"):
            value = analysis.get(key)
            if value is None:
                continue
            if not isinstance(value, str) or not SHA256_RE.fullmatch(value):
                raise SuiteReportError(f"{field}.{key} must be a SHA-256")
            observed.append(value)
    if not observed:
        return None
    if len(set(observed)) != 1:
        raise SuiteReportError("score contains conflicting evaluator profile hashes")
    return observed[0]


def _evaluator_version_hashes(
    score: Mapping,
) -> tuple[tuple[str, str], ...] | None:
    observed: list[tuple[tuple[str, str], ...]] = []
    for field, analysis in _formal_analysis(score):
        versions = analysis.get("evaluator_versions")
        if versions is None:
            continue
        if not isinstance(versions, Mapping) or not versions:
            raise SuiteReportError(
                f"{field}.evaluator_versions must be a non-empty mapping"
            )
        normalized: list[tuple[str, str]] = []
        for name, metadata in versions.items():
            if not isinstance(name, str) or not isinstance(metadata, Mapping):
                raise SuiteReportError(
                    f"{field}.evaluator_versions has an invalid entry"
                )
            digest = metadata.get("sha256")
            if not isinstance(digest, str) or not SHA256_RE.fullmatch(digest):
                raise SuiteReportError(
                    f"{field}.evaluator_versions.{name}.sha256 is invalid"
                )
            normalized.append((name, digest))
        observed.append(tuple(sorted(normalized)))
    if not observed:
        return None
    if len(set(observed)) != 1:
        raise SuiteReportError(
            "score contains conflicting evaluator version fingerprints"
        )
    return observed[0]


def _asset_hashes(score: Mapping) -> tuple[tuple[str, str | None], ...] | None:
    observed: list[tuple[tuple[str, str | None], ...]] = []
    for field, analysis in _formal_analysis(score):
        assets = analysis.get("asset_hashes")
        if assets is None:
            if field == "formal_analysis":
                raise SuiteReportError("formal_analysis.asset_hashes is required")
            continue
        if not isinstance(assets, Mapping) or set(assets) != set(ALL_ASSET_KEYS):
            raise SuiteReportError(
                f"{field}.asset_hashes must contain exactly reference, truth_vcf, "
                "benchmark_bed, pangenome_manifest, challenge_hidden_ledger, "
                "and graph_asset_lock"
            )
        normalized: list[tuple[str, str | None]] = []
        for name in ALL_ASSET_KEYS:
            value = assets[name]
            if value is None and name != "graph_asset_lock":
                raise SuiteReportError(
                    f"{field}.asset_hashes.{name} must be a non-empty SHA-256"
                )
            if value is not None and not (
                isinstance(value, str) and SHA256_RE.fullmatch(value)
            ):
                raise SuiteReportError(
                    f"{field}.asset_hashes.{name} is not a SHA-256"
                )
            normalized.append((name, value))
        observed.append(tuple(normalized))
    if not observed:
        return None
    if len(set(observed)) != 1:
        raise SuiteReportError("score contains conflicting benchmark asset hashes")
    return observed[0]


def _evidence_contract(score: Mapping) -> tuple | None:
    observed: list[tuple] = []
    for field, analysis in _formal_analysis(score):
        evidence = analysis.get("evidence_profile")
        if evidence is None:
            continue
        if not isinstance(evidence, Mapping):
            raise SuiteReportError(f"{field}.evidence_profile must be a mapping")
        required = (
            "contract",
            "actual_technology",
            "library_id",
            "source_evidence_id",
            "evidence_kind",
            "official_score_mode",
            "input_assets",
            "resolved_inputs_sha256",
        )
        if any(key not in evidence for key in required):
            raise SuiteReportError(
                f"{field}.evidence_profile is missing required identity fields"
            )
        assets = evidence["input_assets"]
        if not isinstance(assets, Mapping) or not assets:
            raise SuiteReportError(
                f"{field}.evidence_profile.input_assets must be non-empty"
            )
        normalized_assets: list[tuple[str, str, int]] = []
        for name, metadata in assets.items():
            if not isinstance(name, str) or not isinstance(metadata, Mapping):
                raise SuiteReportError("evidence input asset entry is invalid")
            digest = metadata.get("sha256")
            size = metadata.get("size_bytes")
            if not isinstance(digest, str) or not SHA256_RE.fullmatch(digest):
                raise SuiteReportError(f"evidence asset {name}.sha256 is invalid")
            if not isinstance(size, int) or size <= 0:
                raise SuiteReportError(f"evidence asset {name}.size_bytes is invalid")
            normalized_assets.append((name, digest, size))
        resolved_hash = evidence["resolved_inputs_sha256"]
        if not isinstance(resolved_hash, str) or not SHA256_RE.fullmatch(resolved_hash):
            raise SuiteReportError("evidence resolved_inputs_sha256 is invalid")
        observed.append(
            (
                evidence["contract"],
                evidence["actual_technology"],
                evidence["library_id"],
                evidence["source_evidence_id"],
                evidence["evidence_kind"],
                evidence["official_score_mode"],
                tuple(sorted(normalized_assets)),
                resolved_hash,
                evidence.get("coverage_x"),
                evidence.get("read_count"),
                evidence.get("read_bases"),
                evidence.get("downsampling_seed"),
            )
        )
    if not observed:
        return None
    if len(set(observed)) != 1:
        raise SuiteReportError("score contains conflicting evidence profiles")
    return observed[0]


def _require_uniform_optional_contract(
    values: Sequence[object | None],
    *,
    missing_message: str,
    mismatch_message: str,
) -> None:
    if not any(value is not None for value in values):
        return
    if any(value is None for value in values):
        raise SuiteReportError(missing_message)
    if len(set(values)) != 1:
        raise SuiteReportError(mismatch_message)


def render_suite(entries: Sequence[tuple[Mapping, Mapping]]) -> str:
    if not entries:
        raise SuiteReportError("at least one suite entry is required")
    first_score = entries[0][0]
    first_tuple = first_score["tuple_key"]
    comparison_key = (
        first_tuple["sample"],
        first_tuple["primary_truth_profile"],
        first_score["score_profile_sha256"],
    )

    evaluator_hashes = [_evaluator_profile_hash(score) for score, _ in entries]
    _require_uniform_optional_contract(
        evaluator_hashes,
        missing_message=(
            "suite cannot mix formal scores with and without an evaluator profile hash"
        ),
        mismatch_message="suite formal scores must share the same evaluator profile hash",
    )
    evaluator_versions = [_evaluator_version_hashes(score) for score, _ in entries]
    _require_uniform_optional_contract(
        evaluator_versions,
        missing_message=(
            "suite cannot mix formal scores with and without evaluator version fingerprints"
        ),
        mismatch_message="suite formal scores must share evaluator version fingerprints",
    )

    benchmark_assets = [_asset_hashes(score) for score, _ in entries]
    shared_assets = [
        (
            tuple((name, dict(contract)[name]) for name in SHARED_ASSET_KEYS)
            if contract is not None
            else None
        )
        for contract in benchmark_assets
    ]
    _require_uniform_optional_contract(
        shared_assets,
        missing_message=(
            "suite cannot mix formal scores with and without benchmark asset hashes"
        ),
        mismatch_message=(
            "suite formal scores must share the same reference, truth, BED, "
            "pangenome manifest, and hidden challenge hashes"
        ),
    )

    evidence_contracts = [_evidence_contract(score) for score, _ in entries]
    if any(value is not None for value in evidence_contracts) and any(
        value is None for value in evidence_contracts
    ):
        raise SuiteReportError(
            "suite cannot mix scores with and without evidence profiles"
        )

    graph_contracts: dict[tuple[object, object], str | None] = {}
    evidence_by_technology_and_mode: dict[tuple[object, object], tuple] = {}
    rows: list[str] = []
    for index, ((score, manifest), asset_contract) in enumerate(
        zip(entries, benchmark_assets, strict=True)
    ):
        tuple_key = score.get("tuple_key")
        if not isinstance(tuple_key, Mapping):
            raise SuiteReportError("score is missing tuple_key")
        key = (
            tuple_key.get("sample"),
            tuple_key.get("primary_truth_profile"),
            score.get("score_profile_sha256"),
        )
        if key != comparison_key:
            raise SuiteReportError(
                "suite scores must share sample, truth profile, and score profile"
            )
        tool_mode = (tuple_key.get("tool"), tuple_key.get("official_score_mode"))
        if asset_contract is not None:
            graph_hash = dict(asset_contract)["graph_asset_lock"]
            previous = graph_contracts.setdefault(tool_mode, graph_hash)
            if previous != graph_hash:
                raise SuiteReportError(
                    "suite results for the same tool and mode must share the same "
                    "graph asset lock hash"
                )

        evidence = evidence_contracts[index]
        if evidence is not None:
            evidence_key = (evidence[1], evidence[5])
            previous_evidence = evidence_by_technology_and_mode.setdefault(
                evidence_key, evidence
            )
            if previous_evidence != evidence:
                raise SuiteReportError(
                    "suite results for the same actual technology and score mode must "
                    "use identical frozen read/alignment evidence"
                )
            technology = str(evidence[1])
            evidence_label = (
                f"{evidence[4]}; library={evidence[2]}; source={evidence[3]}"
            )
        else:
            capabilities = manifest.get("capabilities", {})
            supported = (
                capabilities.get("technology", [])
                if isinstance(capabilities, Mapping)
                else []
            )
            technology = ", ".join(map(str, supported))
            evidence_label = "旧结果：未记录输入证据"

        rows.append(
            "<tr>"
            f"<td>{escape(str(tuple_key.get('tool', '')))}</td>"
            f"<td>{escape(str(manifest.get('paradigm', 'unknown')))}</td>"
            f"<td>{escape(technology)}</td>"
            f"<td>{escape(evidence_label)}</td>"
            f"<td>{escape(str(tuple_key.get('official_score_mode', '')))}</td>"
            f"<td>{_number(score.get('pangenome_genotyping_score'))}</td>"
            f"<td>{_number(score.get('non_reference_f1_score'))}</td>"
            f"<td>{_number(score.get('panel_coverage'))}</td>"
            f"<td>{_number(score.get('global_end_to_end_sv_recovery_score'))}</td>"
            f"<td>{_number(score.get('consensus_score'))}</td>"
            f"<td>{escape(str(score.get('truth_eligible_count', '')))}</td>"
            f"<td>{escape(str(score.get('total_evaluated', '')))}</td>"
            f"<td>{escape(str(score.get('score_status', '')))}</td>"
            "</tr>"
        )

    evaluator_notice = (
        "正式结果使用同一个冻结评测器 profile（hash："
        f"<code>{escape(str(evaluator_hashes[0]))}</code>），"
        "并已校验评测器版本指纹和 benchmark 资源哈希一致。"
        if evaluator_hashes[0] is not None
        else "这些旧式或合成结果未声明正式评测器 profile。"
    )
    return """<!doctype html>
<html lang="zh-CN"><head><meta charset="utf-8">
<title>PGBench unified pangenome SV report</title>
<style>
body { font-family: system-ui, sans-serif; margin: 2rem; max-width: 96rem; }
table { border-collapse: collapse; width: 100%; }
th, td { border: 1px solid #ccd; padding: .5rem; text-align: left; }
th { background: #eef2ff; }
.notice { border-left: 4px solid #536dfe; background: #f3f5ff; padding: .8rem; }
</style></head><body>
<h1>HG002/GRCh38 泛基因组结构变异综合比较</h1>
<p class="notice">所有分数使用同一样本、truth profile 和冻结评分合同。""" + evaluator_notice + """
ComparableScore 比较完整 pipeline 的实际检测结果；跨测序技术的差异同时包含
测序证据与算法影响，因此必须结合“实际测序证据”和“工具范式”解释。</p>
<table><thead><tr><th>工具</th><th>范式</th><th>实际测序技术</th><th>冻结输入证据</th><th>评分模式</th>
<th>Pangenome Genotyping Score</th><th>Non-reference F1</th>
<th>Panel coverage</th><th>Global End-to-End SV Recovery</th>
<th>ConsensusScore</th><th>Truth 数</th>
<th>输出数</th><th>状态</th></tr></thead><tbody>""" + "".join(rows) + """
</tbody></table></body></html>
"""


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--entry",
        action="append",
        nargs=2,
        metavar=("SCORE_JSON", "TOOL_YAML"),
        required=True,
    )
    parser.add_argument("--output", required=True, type=Path)
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        entries = []
        for raw_score, raw_manifest in args.entry:
            score = json.loads(Path(raw_score).read_text(encoding="utf-8"))
            if not isinstance(score, dict):
                raise SuiteReportError(f"score must be an object: {raw_score}")
            entries.append((score, _manifest(Path(raw_manifest))))
        args.output.parent.mkdir(parents=True, exist_ok=True)
        temporary = args.output.with_name(f".{args.output.name}.tmp")
        temporary.write_text(render_suite(entries), encoding="utf-8")
        temporary.replace(args.output)
    except (OSError, json.JSONDecodeError, SuiteReportError, yaml.YAMLError) as exc:
        print(f"render_suite_report: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
