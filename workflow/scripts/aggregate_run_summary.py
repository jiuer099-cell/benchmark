#!/usr/bin/env python3
"""Aggregate per-tool artifacts without calculating a run-level score.

The input order is the display order.  Scores are never used for sorting and
no ranking field is accepted anywhere in the input documents.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import sys
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

TUPLE_FIELDS = (
    "run_id",
    "sample",
    "tool",
    "official_score_mode",
    "primary_truth_profile",
)
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


class SummaryError(ValueError):
    """Raised when per-tool artifacts cannot be safely aggregated."""


def _load_object(path: Path, label: str) -> dict[str, Any]:
    try:
        with path.open("r", encoding="utf-8") as handle:
            value = json.load(handle)
    except (OSError, json.JSONDecodeError) as exc:
        raise SummaryError(f"cannot load {label} {path}: {exc}") from exc
    if not isinstance(value, dict):
        raise SummaryError(f"{label} must be a JSON object: {path}")
    _reject_forbidden_fields(value, label)
    return value


def _reject_forbidden_fields(value: Any, path: str = "document") -> None:
    if isinstance(value, Mapping):
        for key, nested in value.items():
            if _is_ranking_field(key):
                raise SummaryError(f"forbidden ranking field at {path}.{key}")
            _reject_forbidden_fields(nested, f"{path}.{key}")
    elif isinstance(value, list):
        for index, nested in enumerate(value):
            _reject_forbidden_fields(nested, f"{path}[{index}]")


def _score_tuple(score: Mapping[str, Any]) -> dict[str, str]:
    raw = score.get("tuple_key")
    if not isinstance(raw, Mapping) or set(raw) != set(TUPLE_FIELDS):
        raise SummaryError("score tuple_key must contain exactly the five tuple fields")
    result: dict[str, str] = {}
    for field in TUPLE_FIELDS:
        value = raw.get(field)
        if not isinstance(value, str) or not value:
            raise SummaryError(f"score tuple_key.{field} must be non-empty")
        result[field] = value
    return result


def _metric_tuple(document: Mapping[str, Any]) -> dict[str, str]:
    raw = document.get("tuple")
    if not isinstance(raw, Mapping):
        raise SummaryError("metrics document is missing tuple")
    mapping = {
        "run_id": "run_id",
        "sample": "sample_id",
        "tool": "tool_id",
        "official_score_mode": "official_score_mode",
        "primary_truth_profile": "primary_truth_profile",
    }
    result: dict[str, str] = {}
    for target, source in mapping.items():
        value = raw.get(source)
        if not isinstance(value, str) or not value:
            raise SummaryError(f"metrics tuple.{source} must be non-empty")
        result[target] = value
    return result


def _tsv_value(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, float):
        if not math.isfinite(value):
            raise SummaryError("summary values must be finite")
        return f"{value:.12g}"
    if isinstance(value, dict | list):
        return json.dumps(value, sort_keys=True, separators=(",", ":"))
    return str(value)


def _write_tsv(
    path: Path,
    fields: Sequence[str],
    rows: Sequence[Mapping[str, Any]],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp")
    with temporary.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(fields),
            delimiter="\t",
            extrasaction="raise",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow({field: _tsv_value(row.get(field)) for field in fields})
    temporary.replace(path)


def _atomic_json(path: Path, payload: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp")
    temporary.write_text(
        json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    temporary.replace(path)


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _validate_sealed_package(
    *,
    package_path: Path,
    score_path: Path,
    metrics_path: Path,
    score: Mapping[str, Any],
) -> dict[str, str]:
    package = _load_object(package_path, "score package")
    if package.get("package_schema_version") != "pgbench.final_score_package.v1":
        raise SummaryError("unsupported score package schema")
    if package.get("sealed") is not True or package.get("score") != score:
        raise SummaryError("score package is not a sealed copy of score JSON")
    gates = package.get("gates")
    if not isinstance(gates, Mapping) or not all(
        gates.get(name) is True
        for name in (
            "core_provenance_valid",
            "hash_lineage_complete",
            "manifest_completeness",
        )
    ):
        raise SummaryError("score package core provenance gates are not sealed")
    artifacts = (
        ("score_artifact", score_path),
        ("metrics_artifact", metrics_path),
    )
    artifact_hashes: dict[str, str] = {}
    for label, path in artifacts:
        artifact = package.get(label)
        if not isinstance(artifact, Mapping):
            raise SummaryError(f"score package is missing {label}")
        artifact_path = artifact.get("path")
        artifact_hash = artifact.get("sha256")
        if not isinstance(artifact_path, str) or not isinstance(artifact_hash, str):
            raise SummaryError(f"score package has an invalid {label}")
        if Path(artifact_path).resolve(strict=False) != path.resolve(strict=False):
            raise SummaryError(f"score package {label} points to a different file")
        if _sha256_file(path) != artifact_hash:
            raise SummaryError(f"{label} hash does not match the sealed package")
        artifact_hashes[label] = artifact_hash
    return artifact_hashes


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
        raise SummaryError(f"finalizer manifest has invalid {paths_field}")
    if (
        len(declared_paths) != len(set(declared_paths))
        or not isinstance(declared_hashes, Mapping)
        or set(declared_hashes) != set(declared_paths)
    ):
        raise SummaryError(
            f"finalizer manifest {hashes_field} keys must exactly match {paths_field}"
        )

    expected = artifact_path.resolve(strict=False)
    matching_paths = [
        path for path in declared_paths if Path(path).resolve(strict=False) == expected
    ]
    if len(matching_paths) != 1:
        raise SummaryError(
            f"finalizer manifest must bind exactly one {label} {direction} path"
        )
    digest = declared_hashes.get(matching_paths[0])
    if (
        not isinstance(digest, str)
        or len(digest) != 64
        or any(character not in "0123456789abcdef" for character in digest)
    ):
        raise SummaryError(f"finalizer manifest has invalid {label} {direction} hash")
    if _sha256_file(artifact_path) != digest:
        raise SummaryError(
            f"{label} hash does not match the successful finalizer manifest"
        )
    return digest


def _validate_finalizer_manifest(
    *,
    manifest_path: Path,
    package_path: Path,
    score_path: Path,
    metrics_path: Path,
    score: Mapping[str, Any],
    package_artifact_hashes: Mapping[str, str],
) -> None:
    manifest = _load_object(manifest_path, "finalizer manifest")
    if manifest.get("rule_name") != FINALIZER_RULE_NAME:
        raise SummaryError("manifest is not a finalize_score_provenance manifest")
    if manifest.get("status") != "success":
        raise SummaryError("finalizer manifest is not successful")

    tuple_key = _score_tuple(score)
    if manifest.get("run_id") != tuple_key["run_id"]:
        raise SummaryError("finalizer manifest run_id does not match score tuple")
    wildcards = manifest.get("wildcards")
    if not isinstance(wildcards, Mapping) or wildcards.get("tool") != tuple_key["tool"]:
        raise SummaryError("finalizer manifest tool does not match score tuple")
    expected_job_key = ".".join(
        (
            tuple_key["sample"],
            tuple_key["tool"],
            tuple_key["official_score_mode"],
        )
    )
    if manifest.get("job_key") != expected_job_key:
        raise SummaryError("finalizer manifest job_key does not match score tuple")

    _manifest_artifact_hash(
        manifest,
        package_path,
        direction="output",
        label="score package",
    )
    manifest_score_hash = _manifest_artifact_hash(
        manifest,
        score_path,
        direction="input",
        label="score JSON",
    )
    manifest_metrics_hash = _manifest_artifact_hash(
        manifest,
        metrics_path,
        direction="input",
        label="metrics JSON",
    )
    if manifest_score_hash != package_artifact_hashes["score_artifact"]:
        raise SummaryError(
            "score package and finalizer manifest bind different score hashes"
        )
    if manifest_metrics_hash != package_artifact_hashes["metrics_artifact"]:
        raise SummaryError(
            "score package and finalizer manifest bind different metrics hashes"
        )


def aggregate(
    *,
    score_paths: Sequence[Path],
    metrics_paths: Sequence[Path],
    score_package_paths: Sequence[Path] | None = None,
    finalizer_manifest_paths: Sequence[Path] | None = None,
    output_score_tsv: Path,
    output_point_tsv: Path,
    output_metrics_tsv: Path,
    output_metrics_json: Path,
) -> None:
    if not score_paths or len(score_paths) != len(metrics_paths):
        raise SummaryError("score and metrics inputs must be non-empty paired lists")
    if score_package_paths is not None and len(score_package_paths) != len(score_paths):
        raise SummaryError("score package inputs must pair with score inputs")
    if finalizer_manifest_paths is not None and len(finalizer_manifest_paths) != len(
        score_paths
    ):
        raise SummaryError("finalizer manifest inputs must pair with score inputs")
    if (score_package_paths is None) != (finalizer_manifest_paths is None):
        raise SummaryError(
            "score package and finalizer manifest inputs must be supplied together"
        )

    score_rows: list[dict[str, Any]] = []
    point_rows: list[dict[str, Any]] = []
    metric_rows: list[dict[str, Any]] = []
    metric_documents: list[dict[str, Any]] = []
    seen_tuples: set[tuple[str, ...]] = set()

    package_paths: Sequence[Path | None] = (
        score_package_paths
        if score_package_paths is not None
        else [None] * len(score_paths)
    )
    finalizer_paths: Sequence[Path | None] = (
        finalizer_manifest_paths
        if finalizer_manifest_paths is not None
        else [None] * len(score_paths)
    )
    for score_path, metrics_path, package_path, finalizer_path in zip(
        score_paths,
        metrics_paths,
        package_paths,
        finalizer_paths,
        strict=True,
    ):
        score = _load_object(score_path, "score")
        metrics = _load_object(metrics_path, "metrics")
        if package_path is not None:
            if finalizer_path is None:  # guarded above; narrows the type for mypy
                raise SummaryError("finalizer manifest input is missing")
            artifact_hashes = _validate_sealed_package(
                package_path=package_path,
                score_path=score_path,
                metrics_path=metrics_path,
                score=score,
            )
            _validate_finalizer_manifest(
                manifest_path=finalizer_path,
                package_path=package_path,
                score_path=score_path,
                metrics_path=metrics_path,
                score=score,
                package_artifact_hashes=artifact_hashes,
            )
        tuple_key = _score_tuple(score)
        if _metric_tuple(metrics) != tuple_key:
            raise SummaryError(
                f"score/metrics tuple mismatch: {score_path} vs {metrics_path}"
            )
        tuple_identity = tuple(tuple_key[field] for field in TUPLE_FIELDS)
        if tuple_identity in seen_tuples:
            raise SummaryError(f"duplicate tool tuple: {tuple_identity!r}")
        seen_tuples.add(tuple_identity)

        score_rows.append(
            {
                **tuple_key,
                "score_profile": score.get("score_profile"),
                "score_profile_sha256": score.get("score_profile_sha256"),
                "evaluation_mode": score.get("evaluation_mode"),
                "score_status": score.get("score_status"),
                "ConsensusScore": score.get("consensus_score", score.get("pgbench_score")),
                "total_evaluated": score.get("total_evaluated"),
                "unanimous_correct_rate": score.get("unanimous_correct_rate"),
                "majority_correct_rate": score.get("majority_correct_rate"),
            }
        )
        breakdown = score.get("point_breakdown")
        if not isinstance(breakdown, Mapping):
            raise SummaryError(f"score is missing point_breakdown: {score_path}")
        for component, points in sorted(breakdown.items()):
            if (
                not isinstance(component, str)
                or isinstance(points, bool)
                or not isinstance(points, int | float)
                or not math.isfinite(float(points))
            ):
                raise SummaryError(f"invalid point breakdown in {score_path}")
            point_rows.append(
                {**tuple_key, "category": component, "count": int(points)}
            )

        records = metrics.get("records")
        if not isinstance(records, list) or not records:
            raise SummaryError(f"metrics records must be non-empty: {metrics_path}")
        for record in records:
            if not isinstance(record, Mapping):
                raise SummaryError(f"metric record must be an object: {metrics_path}")
            metric_rows.append({**tuple_key, **dict(record)})
        metric_documents.append(metrics)

    _write_tsv(
        output_score_tsv,
        (
            *TUPLE_FIELDS,
            "score_profile",
            "score_profile_sha256",
            "evaluation_mode",
            "score_status",
            "ConsensusScore",
            "total_evaluated",
            "unanimous_correct_rate",
            "majority_correct_rate",
        ),
        score_rows,
    )
    _write_tsv(
        output_point_tsv,
        (*TUPLE_FIELDS, "category", "count"),
        point_rows,
    )
    metric_fields = (
        *TUPLE_FIELDS,
        "metric_id",
        "value_type",
        "value",
        "status",
        "evaluator",
        "numerator",
        "denominator",
        "eligible_count",
        "universe_id",
        "undefined_reason",
        "parser_id",
        "parser_source_field",
        "provenance_manifest_id",
        "strata",
        "confidence_interval",
    )
    _write_tsv(output_metrics_tsv, metric_fields, metric_rows)
    _atomic_json(
        output_metrics_json,
        {"schema_version": 1, "documents": metric_documents},
    )


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--score-json", action="append", required=True, type=Path)
    parser.add_argument("--metrics-json", action="append", required=True, type=Path)
    parser.add_argument("--score-package", action="append", required=True, type=Path)
    parser.add_argument(
        "--finalizer-manifest", action="append", required=True, type=Path
    )
    parser.add_argument("--output-score-tsv", required=True, type=Path)
    parser.add_argument("--output-point-tsv", required=True, type=Path)
    parser.add_argument("--output-metrics-tsv", required=True, type=Path)
    parser.add_argument("--output-metrics-json", required=True, type=Path)
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        aggregate(
            score_paths=args.score_json,
            metrics_paths=args.metrics_json,
            score_package_paths=args.score_package,
            finalizer_manifest_paths=args.finalizer_manifest,
            output_score_tsv=args.output_score_tsv,
            output_point_tsv=args.output_point_tsv,
            output_metrics_tsv=args.output_metrics_tsv,
            output_metrics_json=args.output_metrics_json,
        )
    except (OSError, SummaryError) as exc:
        print(f"aggregate_run_summary: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
