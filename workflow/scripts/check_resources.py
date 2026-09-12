#!/usr/bin/env python3
"""Read-only preflight for a configured PGBench fixed-panel track."""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

import yaml  # type: ignore[import-untyped]


TRACK = "end_to_end_from_reads"
class ResourceCheckError(ValueError):
    """Raised when the preflight configuration cannot be interpreted."""


def load_config(path: Path) -> dict[str, Any]:
    try:
        value = yaml.safe_load(path.read_text(encoding="utf-8"))
    except (OSError, yaml.YAMLError) as exc:
        raise ResourceCheckError(f"cannot load {path}: {exc}") from exc
    if not isinstance(value, dict):
        raise ResourceCheckError(f"{path} must contain a YAML mapping")
    return value


def selected_tracks(config: Mapping[str, Any], track_selection: str = "configured") -> list[str]:
    if track_selection not in {"configured", "short-read"}:
        raise ResourceCheckError(f"unsupported track selection: {track_selection}")
    configured = config.get("execution", {}).get("tracks", [])
    if configured != [TRACK]:
        raise ResourceCheckError(f"execution.tracks must be exactly [{TRACK}]")
    return [TRACK]


def _resolve(repo_root: Path, value: str) -> Path:
    path = Path(value)
    return (path if path.is_absolute() else repo_root / path).resolve(strict=False)


def _inspect(
    resource_id: str,
    category: str,
    value: Any,
    repo_root: Path,
    *,
    required: bool = True,
) -> dict[str, Any]:
    if not isinstance(value, str) or not value:
        return {
            "id": resource_id, "category": category, "required": required,
            "status": "unconfigured", "path": None, "size_bytes": None,
            "candidate_paths": [], "detail": None,
        }
    path = _resolve(repo_root, value)
    status = "ok" if path.is_file() and path.stat().st_size > 0 else "missing"
    return {
        "id": resource_id, "category": category, "required": required,
        "status": status, "path": str(path),
        "size_bytes": path.stat().st_size if status == "ok" else None,
        "candidate_paths": [str(path)], "detail": None,
    }


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _inspect_frozen_bundle_lock(value: Any, repo_root: Path, *, required: bool) -> dict[str, Any]:
    """Check that the declared lock actually seals every source input by hash."""

    item = _inspect(
        "pangenome.frozen_haplotype_source_bundle.content_lock",
        "pangenome_bundle_lock",
        value,
        repo_root,
        required=required,
    )
    if item["status"] != "ok":
        return item
    try:
        lock_path = Path(str(item["path"]))
        lock = json.loads(lock_path.read_text(encoding="utf-8"))
        required_assets = {
            "gfa_or_gbz", "population_vcf", "sample_roster", "haplotype_roster",
            "family_exclusion_manifest", "reference",
        }
        assets = lock.get("assets") if isinstance(lock, Mapping) else None
        software = lock.get("software") if isinstance(lock, Mapping) else None
        if (
            lock.get("schema_version") != 1
            or lock.get("contract") != "pgbench_frozen_haplotype_source_bundle_v1"
            or not isinstance(assets, Mapping)
            or not required_assets.issubset(assets)
            or not isinstance(software, list)
            or not software
        ):
            raise ResourceCheckError("invalid bundle lock contract")
        for name in required_assets:
            record = assets[name]
            if not isinstance(record, Mapping):
                raise ResourceCheckError(f"bundle lock asset {name} is invalid")
            raw_path, expected = record.get("path"), record.get("sha256")
            if not isinstance(raw_path, str) or not isinstance(expected, str):
                raise ResourceCheckError(f"bundle lock asset {name} lacks path or SHA256")
            asset_path = Path(raw_path)
            asset_path = asset_path if asset_path.is_absolute() else lock_path.parent / asset_path
            if not asset_path.is_file() or _sha256(asset_path) != expected:
                raise ResourceCheckError(f"bundle lock asset {name} is missing or hash-mismatched")
        for tool in software:
            if not isinstance(tool, Mapping) or not all(
                isinstance(tool.get(field), str) and tool[field]
                for field in ("name", "version", "sha256")
            ):
                raise ResourceCheckError("bundle lock software record is invalid")
    except (OSError, ValueError, TypeError, ResourceCheckError) as exc:
        item["status"] = "invalid"
        item["detail"] = str(exc)
    return item


def _development_value(config: Mapping[str, Any], key: str, production: Any) -> Any:
    development = config.get("development", {})
    if development.get("synthetic_mode") and development.get(key) is not None:
        return development[key]
    return production


def build_inventory(
    config: Mapping[str, Any], *, repo_root: Path, track_selection: str = "configured"
) -> dict[str, Any]:
    tracks = selected_tracks(config, track_selection)
    reference = config.get("reference", {})
    sample = config.get("sample", {})
    pangenome = config.get("pangenome", {})
    evaluation = config.get("evaluation", {})
    synthetic = bool(config.get("development", {}).get("synthetic_mode"))
    truth = _development_value(config, "truth_vcf", evaluation.get("truth_vcf"))
    bed = _development_value(config, "benchmark_bed", evaluation.get("benchmark_bed"))
    population = _development_value(config, "population_vcf", pangenome.get("population_vcf"))
    frozen_bundle = pangenome.get("frozen_haplotype_source_bundle", {})
    assets = [
        _inspect("reference.fasta", "reference", reference.get("fasta"), repo_root),
        _inspect("reference.fai", "reference", reference.get("fai"), repo_root, required=not synthetic),
        _inspect("reference.dict", "reference", reference.get("dict"), repo_root, required=not synthetic),
        _inspect("sample.fastq_r1", "short_reads", sample.get("fastq_r1"), repo_root),
        _inspect("sample.fastq_r2", "short_reads", sample.get("fastq_r2"), repo_root),
        _inspect("evaluation.truth_vcf", "truth", truth, repo_root),
        _inspect("evaluation.benchmark_bed", "truth", bed, repo_root),
        _inspect("pangenome.population_vcf", "pangenome", population, repo_root),
        # The lock is mandatory for a formal run.  It is a content-addressed
        # record of the graph/VCF/rosters/exclusion manifest and tool versions;
        # synthetic fixtures intentionally do not need a production bundle.
        _inspect_frozen_bundle_lock(
            frozen_bundle.get("content_lock") if isinstance(frozen_bundle, Mapping) else None,
            repo_root,
            required=not synthetic,
        ),
    ]
    if not synthetic:
        context_beds = evaluation.get("context_beds", {})
        if not isinstance(context_beds, Mapping):
            raise ResourceCheckError("evaluation.context_beds must be a mapping")
        assets.extend(
            _inspect(
                f"evaluation.context_beds.{name}",
                "genome_context",
                value,
                repo_root,
            )
            for name, value in sorted(context_beds.items())
        )
    if isinstance(truth, str) and truth.endswith(".gz"):
        assets.append(_inspect("evaluation.truth_vcf_tbi", "truth", truth + ".tbi", repo_root))
    if isinstance(population, str) and population.endswith(".gz"):
        assets.append(_inspect("pangenome.population_vcf_tbi", "pangenome", population + ".tbi", repo_root))
    graph = pangenome.get("graph_assets", {})
    profile = graph.get("profile", "none")
    if profile != "none":
        for name, value in sorted(graph.items()):
            if name not in {"profile", "reference_path"} and value:
                assets.append(
                    _inspect(
                        f"pangenome.graph_assets.{name}",
                        "graph",
                        value,
                        repo_root,
                    )
                )
    missing = [item["id"] for item in assets if item["required"] and item["status"] != "ok"]
    return {
        "schema_version": 2,
        "status": "ready" if not missing else "not_ready",
        "tracks": tracks,
        "summary": {
            "required": sum(bool(item["required"]) for item in assets),
            "ready_required": sum(item["required"] and item["status"] == "ok" for item in assets),
            "missing_required": len(missing),
        },
        "missing_required": missing,
        "assets": assets,
    }


def render_text(report: Mapping[str, Any]) -> str:
    lines = [f"resource preflight: {report['status']}", f"track: {TRACK}"]
    lines.extend(
        f"{item['status']:12} {item['id']}: {item['path'] or '-'}"
        for item in report["assets"]
    )
    return "\n".join(lines)


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=Path("config/config.example.yaml"))
    parser.add_argument("--repo-root", type=Path, default=Path.cwd())
    parser.add_argument("--track", choices=("configured", "short-read"), default="configured")
    parser.add_argument("--json", action="store_true")
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        report = build_inventory(load_config(args.config), repo_root=args.repo_root, track_selection=args.track)
    except (OSError, ResourceCheckError, KeyError, TypeError) as exc:
        print(f"check_resources: {exc}", file=sys.stderr)
        return 2
    if args.json:
        json.dump(report, sys.stdout, indent=2, sort_keys=True, ensure_ascii=False)
        sys.stdout.write("\n")
    else:
        print(render_text(report))
    return 0 if report["status"] == "ready" else 1


if __name__ == "__main__":
    raise SystemExit(main())
