#!/usr/bin/env python3
"""Convert the fixed smoke fixture into the normative metrics contract.

This adapter is development-only.  It cannot emit a formal evaluation mode;
the score command independently binds ``synthetic_smoke`` from the workflow
configuration.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

import yaml  # type: ignore[import-untyped]
from pgbench_metrics import (
    DEFAULT_METRICS_SCHEMA_PATH,
    MetricContractError,
    build_score_payload_from_metrics,
)
from pgbench_provenance import load_manifest
from pgbench_scoring import load_score_profile


class SyntheticMetricError(ValueError):
    """Raised when the smoke fixture cannot be represented safely."""


COMPONENT_PATHS = {
    "consensus.all_three_correct.count": ("consensus", "all_three_correct"),
    "consensus.exactly_two_correct.count": ("consensus", "exactly_two_correct"),
    "consensus.exactly_one_correct.count": ("consensus", "exactly_one_correct"),
    "consensus.none_correct.count": ("consensus", "none_correct"),
    "truvari.event.overall.f1": ("evaluators", "truvari", "overall_event_f1"),
    "truvari.event.svtype.macro_f1": (
        "evaluators",
        "truvari",
        "svtype_macro_f1",
    ),
    "truvari.event.length.macro_f1": (
        "evaluators",
        "truvari",
        "length_macro_f1",
    ),
    "truvari.gt.exact_concordance": (
        "evaluators",
        "truvari",
        "exact_gt_concordance",
    ),
    "truvari.match.mean_similarity": ("evaluators", "truvari", "fidelity"),
    "truvari.context.difficult.macro_f1": (
        "evaluators",
        "truvari",
        "difficult_context_f1",
    ),
    "aardvark.haplotype.exact_f1": (
        "evaluators",
        "aardvark",
        "exact_haplotype_f1",
    ),
    "aardvark.haplotype.partial_f1": (
        "evaluators",
        "aardvark",
        "partial_credit_f1",
    ),
    "aardvark.gt.exact_concordance": (
        "evaluators",
        "aardvark",
        "exact_gt_concordance",
    ),
    "aardvark.complex.exact_f1": (
        "evaluators",
        "aardvark",
        "complex_cluster_f1",
    ),
    "aardvark.haplotype.sequence_similarity": (
        "evaluators",
        "aardvark",
        "allele_sequence_fidelity",
    ),
    "vcfdist.representation.f1": (
        "evaluators",
        "vcfdist",
        "representation_f1",
    ),
    "vcfdist.gt.exact_concordance": (
        "evaluators",
        "vcfdist",
        "gt_concordance",
    ),
    "vcfdist.phase.accuracy": ("evaluators", "vcfdist", "phase_accuracy"),
    "vcfdist.complex.consistency": (
        "evaluators",
        "vcfdist",
        "complex_consistency",
    ),
    "vcfdist.eligible.coverage": (
        "evaluators",
        "vcfdist",
        "eligible_coverage",
    ),
    "pangenome.in_panel.gt.fused_macro_f1": (
        "pangenome",
        "in_panel_gt_macro_f1",
    ),
    "pangenome.link.precision": ("pangenome", "allele_link_precision"),
    "pangenome.link.coverage": ("pangenome", "allele_link_coverage"),
    "pangenome.af.fused_macro_f1": (
        "pangenome",
        "population_af_macro_f1",
    ),
    "pangenome.graph.fused_macro_f1": (
        "pangenome",
        "graph_complexity_macro_f1",
    ),
    "pangenome.context.fused_macro_f1": (
        "pangenome",
        "difficult_context_macro_f1",
    ),
    "pangenome.novel.fused_f1": (
        "pangenome",
        "out_of_panel_truth_f1",
    ),
    "resource.wall_time.hours": ("resources", "wall_time_hours"),
    "resource.peak_rss.gib": ("resources", "peak_rss_gib"),
    "resource.disk_peak.gib": ("resources", "disk_peak_gib"),
    "resource.cpu.hours": ("resources", "cpu_hours"),
}


def _load_json_object(path: Path, label: str) -> dict[str, Any]:
    try:
        with path.open("r", encoding="utf-8") as handle:
            value = json.load(handle)
    except (OSError, json.JSONDecodeError) as exc:
        raise SyntheticMetricError(f"cannot load {label} {path}: {exc}") from exc
    if not isinstance(value, dict):
        raise SyntheticMetricError(f"{label} must be a JSON object")
    return value


def _load_yaml_object(path: Path, label: str) -> dict[str, Any]:
    try:
        with path.open("r", encoding="utf-8") as handle:
            value = yaml.safe_load(handle)
    except (OSError, yaml.YAMLError) as exc:
        raise SyntheticMetricError(f"cannot load {label} {path}: {exc}") from exc
    if not isinstance(value, dict):
        raise SyntheticMetricError(f"{label} must be a YAML mapping")
    return value


def _lookup(payload: Mapping[str, Any], path: Sequence[str]) -> float | int:
    value: Any = payload
    for part in path:
        if not isinstance(value, Mapping) or part not in value:
            raise SyntheticMetricError("smoke fixture is missing " + ".".join(path))
        value = value[part]
    if (
        isinstance(value, bool)
        or not isinstance(value, (int, float))
        or not math.isfinite(float(value))
    ):
        raise SyntheticMetricError(
            "smoke fixture metric must be finite: " + ".".join(path)
        )
    return value


def _manifest_id(path: Path) -> str:
    manifest = load_manifest(path)
    manifest_id = manifest.get("manifest_id")
    if not isinstance(manifest_id, str) or len(manifest_id) != 64:
        raise SyntheticMetricError(f"invalid provenance manifest ID: {path}")
    return manifest_id


def materialize(
    *,
    fixture: Mapping[str, Any],
    metric_dictionary: Mapping[str, Any],
    score_profile_path: Path,
    evaluator_manifest_id: str,
    pangenome_manifest_id: str,
    resource_manifest_id: str,
    metrics_schema_path: Path = DEFAULT_METRICS_SCHEMA_PATH,
) -> dict[str, Any]:
    if fixture.get("evaluation_mode") != "synthetic_smoke":
        raise SyntheticMetricError(
            "the fixture adapter only accepts evaluation_mode=synthetic_smoke"
        )
    tuple_value = fixture.get("tuple")
    if not isinstance(tuple_value, Mapping):
        raise SyntheticMetricError("smoke fixture is missing tuple")
    expected_tuple = {
        field: tuple_value.get(field)
        for field in (
            "run_id",
            "sample",
            "tool",
            "official_score_mode",
            "primary_truth_profile",
        )
    }
    if not all(isinstance(value, str) and value for value in expected_tuple.values()):
        raise SyntheticMetricError("smoke fixture tuple fields must be non-empty")

    dictionary_id = metric_dictionary.get("dictionary_id")
    dictionary_metrics = metric_dictionary.get("metrics")
    if dictionary_id != "pgbench_metrics_v1" or not isinstance(
        dictionary_metrics, Mapping
    ):
        raise SyntheticMetricError("unsupported metric dictionary")

    records: list[dict[str, Any]] = []
    for metric_id, source_path in COMPONENT_PATHS.items():
        entry = dictionary_metrics.get(metric_id)
        if not isinstance(entry, Mapping):
            raise SyntheticMetricError(f"metric dictionary is missing {metric_id}")
        value = _lookup(fixture, source_path)
        if str(entry.get("value_type")) == "ratio" and not 0.0 <= value <= 1.0:
            raise SyntheticMetricError(f"ratio outside [0,1]: {metric_id}")
        source = str(entry.get("source"))
        if source == "resource_collector":
            provenance_id = resource_manifest_id
        elif source == "pangenome_linker":
            provenance_id = pangenome_manifest_id
        else:
            provenance_id = evaluator_manifest_id
        ineligible_zero = metric_id == "vcfdist.phase.accuracy" and value == 0.0
        denominator = 1_000_000 if entry.get("value_type") == "ratio" else 1
        numerator = value * denominator
        records.append(
            {
                "metric_id": metric_id,
                "value_type": entry.get("value_type"),
                "value": value,
                "status": (
                    "completed_with_tool_ineligible_zero"
                    if ineligible_zero
                    else "defined"
                ),
                "evaluator": source,
                "numerator": numerator,
                "denominator": denominator,
                "eligible_count": denominator,
                "universe_id": "synthetic_smoke_contract_v1",
                "undefined_reason": (
                    "example genotyper declares phasing=false"
                    if ineligible_zero
                    else None
                ),
                "parser_id": "synthetic_fixture_adapter_v1",
                "parser_source_field": ".".join(source_path),
                "provenance_manifest_id": provenance_id,
                "strata": {},
                "confidence_interval": None,
            }
        )

    score_profile = load_score_profile(score_profile_path)
    document = {
        "schema_version": 1,
        "dictionary_id": dictionary_id,
        "tuple": {
            "run_id": expected_tuple["run_id"],
            "sample_id": expected_tuple["sample"],
            "tool_id": expected_tuple["tool"],
            "official_score_mode": expected_tuple["official_score_mode"],
            "primary_truth_profile": expected_tuple["primary_truth_profile"],
            "score_profile": score_profile["profile"]["id"],
        },
        "records": records,
    }
    build_score_payload_from_metrics(
        document,
        expected_tuple=expected_tuple,  # type: ignore[arg-type]
        score_profile=score_profile,
        metrics_schema_path=metrics_schema_path,
        metric_dictionary_path=Path(
            str(metric_dictionary.get("_source_path", "config/metric_dictionary.yaml"))
        ),
    )
    return document


def _atomic_json(path: Path, payload: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp")
    temporary.write_text(
        json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    temporary.replace(path)


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fixture", required=True, type=Path)
    parser.add_argument("--metric-dictionary", required=True, type=Path)
    parser.add_argument("--metrics-schema", required=True, type=Path)
    parser.add_argument("--score-profile", required=True, type=Path)
    parser.add_argument("--evaluator-manifest", required=True, type=Path)
    parser.add_argument("--pangenome-manifest", required=True, type=Path)
    parser.add_argument("--resource-manifest", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        dictionary = _load_yaml_object(args.metric_dictionary, "metric dictionary")
        dictionary["_source_path"] = str(args.metric_dictionary)
        document = materialize(
            fixture=_load_json_object(args.fixture, "smoke fixture"),
            metric_dictionary=dictionary,
            score_profile_path=args.score_profile,
            evaluator_manifest_id=_manifest_id(args.evaluator_manifest),
            pangenome_manifest_id=_manifest_id(args.pangenome_manifest),
            resource_manifest_id=_manifest_id(args.resource_manifest),
            metrics_schema_path=args.metrics_schema,
        )
        _atomic_json(args.output, document)
    except (
        MetricContractError,
        OSError,
        SyntheticMetricError,
        yaml.YAMLError,
    ) as exc:
        print(f"materialize_synthetic_metrics: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
