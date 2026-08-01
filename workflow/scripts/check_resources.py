#!/usr/bin/env python3
"""Read-only preflight for files required by configured benchmark tracks."""

from __future__ import annotations

import argparse
import json
import sys
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

import yaml  # type: ignore[import-untyped]


CALLER_ONLY = "caller_only_shared_alignment"
END_TO_END = "end_to_end_from_reads"
TRACK_ALIASES = {
    "caller-only": [CALLER_ONLY],
    "end-to-end": [END_TO_END],
    "all": [CALLER_ONLY, END_TO_END],
}
GRAPH_PROFILE_ASSETS = {
    "none": (),
    "vg_gbz_min_dist": ("manifest", "gbz", "min", "dist", "sample_list"),
    "vg_legacy_xg": ("manifest", "gbz", "xg", "min", "dist", "sample_list"),
}


class ResourceCheckError(ValueError):
    """Raised when the preflight configuration cannot be interpreted."""


def load_config(path: Path) -> dict[str, Any]:
    """Load a YAML configuration without changing it or any resource."""

    try:
        with path.open("r", encoding="utf-8") as handle:
            loaded = yaml.safe_load(handle)
    except (OSError, yaml.YAMLError) as exc:
        raise ResourceCheckError(f"cannot load {path}: {exc}") from exc
    if not isinstance(loaded, dict):
        raise ResourceCheckError(f"{path} must contain a YAML mapping")
    return loaded


def selected_tracks(
    config: Mapping[str, Any],
    track_selection: str = "configured",
) -> list[str]:
    """Return the canonical tracks covered by this preflight."""

    if track_selection != "configured":
        return TRACK_ALIASES[track_selection]
    configured = config.get("execution", {}).get("tracks", [])
    if not isinstance(configured, list):
        raise ResourceCheckError("execution.tracks must be a list")
    tracks = [track for track in configured if track in (CALLER_ONLY, END_TO_END)]
    if not tracks:
        raise ResourceCheckError("execution.tracks selects no supported track")
    return tracks


def _configured_path(value: Any) -> str | None:
    if value is None:
        return None
    if not isinstance(value, str) or not value.strip():
        return None
    return value


def _resolve(repo_root: Path, configured_path: str) -> Path:
    path = Path(configured_path)
    if not path.is_absolute():
        path = repo_root / path
    return path.resolve(strict=False)


def _inspect_file(
    *,
    resource_id: str,
    category: str,
    configured_path: str | None,
    repo_root: Path,
    required: bool,
    candidates: Sequence[str] = (),
) -> dict[str, Any]:
    candidate_paths = [
        _resolve(repo_root, item)
        for item in ([configured_path] if configured_path else [])
    ]
    candidate_paths.extend(_resolve(repo_root, item) for item in candidates)
    # Preserve order while removing the duplicate primary BAI candidate.
    candidate_paths = list(dict.fromkeys(candidate_paths))
    existing = next((path for path in candidate_paths if path.is_file()), None)
    inspected = existing or (candidate_paths[0] if candidate_paths else None)

    status = "unconfigured"
    size_bytes: int | None = None
    detail: str | None = None
    if inspected is not None:
        try:
            if not inspected.exists():
                status = "missing"
            elif not inspected.is_file():
                status = "not_file"
            else:
                size_bytes = inspected.stat().st_size
                with inspected.open("rb") as handle:
                    handle.read(1)
                status = "ok" if size_bytes > 0 else "empty"
        except OSError as exc:
            status = "unreadable"
            detail = str(exc)

    return {
        "id": resource_id,
        "category": category,
        "required": required,
        "status": status,
        "path": str(inspected) if inspected is not None else None,
        "size_bytes": size_bytes,
        "candidate_paths": [str(path) for path in candidate_paths],
        "detail": detail,
    }


def _development_value(
    config: Mapping[str, Any],
    development_key: str,
    production_value: Any,
) -> Any:
    development = config.get("development", {})
    if development.get("synthetic_mode", False):
        override = development.get(development_key)
        if override is not None:
            return override
    return production_value


def _bam_index_candidates(bam: str | None) -> tuple[str, ...]:
    if bam is None:
        return ()
    path = Path(bam)
    candidates = [f"{bam}.bai"]
    if path.suffix.lower() == ".bam":
        candidates.append(str(path.with_suffix(".bai")))
    return tuple(dict.fromkeys(candidates))


def build_inventory(
    config: Mapping[str, Any],
    *,
    repo_root: Path,
    track_selection: str = "configured",
) -> dict[str, Any]:
    """Inspect configured paths and return a JSON-serializable inventory."""

    tracks = selected_tracks(config, track_selection)
    needs_bam = CALLER_ONLY in tracks
    needs_fastq = END_TO_END in tracks
    reference = config.get("reference", {})
    pangenome = config.get("pangenome", {})
    evaluation = config.get("evaluation", {})
    sample = config.get("sample", {})
    caller_only = config.get("caller_only", {})
    synthetic_mode = bool(
        config.get("development", {}).get("synthetic_mode", False)
    )

    fasta = _configured_path(reference.get("fasta"))
    fai = _configured_path(reference.get("fai"))
    if fai is None and fasta is not None:
        fai = f"{fasta}.fai"
    sequence_dict = _configured_path(reference.get("dict"))
    if sequence_dict is None and fasta is not None:
        sequence_dict = str(Path(fasta).with_suffix(".dict"))

    bam = _configured_path(caller_only.get("shared_bam") or sample.get("bam"))
    fastq = _configured_path(
        _development_value(config, "canonical_fastq", sample.get("fastq"))
    )
    fastq_r1 = _configured_path(sample.get("fastq_r1"))
    fastq_r2 = _configured_path(sample.get("fastq_r2"))
    truth_vcf = _configured_path(
        _development_value(config, "truth_vcf", evaluation.get("truth_vcf"))
    )
    benchmark_bed = _configured_path(
        _development_value(
            config,
            "benchmark_bed",
            evaluation.get("benchmark_bed"),
        )
    )
    population_vcf = _configured_path(
        _development_value(
            config,
            "population_vcf",
            pangenome.get("population_vcf"),
        )
    )
    graph_required = bool(pangenome.get("build_graph_assets", False))

    assets = [
        _inspect_file(
            resource_id="reference.fasta",
            category="reference",
            configured_path=fasta,
            repo_root=repo_root,
            required=True,
        ),
        _inspect_file(
            resource_id="reference.fai",
            category="reference",
            configured_path=fai,
            repo_root=repo_root,
            required=not synthetic_mode,
        ),
        _inspect_file(
            resource_id="reference.dict",
            category="reference",
            configured_path=sequence_dict,
            repo_root=repo_root,
            required=not synthetic_mode,
        ),
        _inspect_file(
            resource_id="sample.bam",
            category="caller_only",
            configured_path=bam,
            repo_root=repo_root,
            required=needs_bam,
        ),
        _inspect_file(
            resource_id="sample.bai",
            category="caller_only",
            configured_path=None,
            candidates=_bam_index_candidates(bam),
            repo_root=repo_root,
            required=needs_bam,
        ),
        _inspect_file(
            resource_id="evaluation.truth_vcf",
            category="truth",
            configured_path=truth_vcf,
            repo_root=repo_root,
            required=True,
        ),
        _inspect_file(
            resource_id="evaluation.truth_vcf_tbi",
            category="truth",
            configured_path=f"{truth_vcf}.tbi" if truth_vcf else None,
            repo_root=repo_root,
            required=bool(truth_vcf and truth_vcf.endswith(".gz")),
        ),
        _inspect_file(
            resource_id="evaluation.benchmark_bed",
            category="truth",
            configured_path=benchmark_bed,
            repo_root=repo_root,
            required=True,
        ),
        _inspect_file(
            resource_id="pangenome.population_vcf",
            category="pangenome",
            configured_path=population_vcf,
            repo_root=repo_root,
            required=True,
        ),
        _inspect_file(
            resource_id="pangenome.population_vcf_tbi",
            category="pangenome",
            configured_path=(
                f"{population_vcf}.tbi"
                if population_vcf and population_vcf.endswith(".gz")
                else None
            ),
            repo_root=repo_root,
            required=bool(population_vcf and population_vcf.endswith(".gz")),
        ),
    ]

    graph_assets = pangenome.get("graph_assets", {})
    graph_profile = graph_assets.get("profile", "none")
    required_graph_assets = set(GRAPH_PROFILE_ASSETS.get(graph_profile, ()))
    for asset_name in ("manifest", "gbz", "xg", "min", "dist", "sample_list"):
        assets.append(
            _inspect_file(
                resource_id=f"pangenome.graph_assets.{asset_name}",
                category="graph",
                configured_path=_configured_path(graph_assets.get(asset_name)),
                repo_root=repo_root,
                required=graph_required and asset_name in required_graph_assets,
            )
        )

    paired_configured = bool(fastq_r1 or fastq_r2)
    assets.append(
        _inspect_file(
            resource_id="sample.fastq",
            category="end_to_end",
            configured_path=fastq,
            repo_root=repo_root,
            required=needs_fastq and not paired_configured,
        )
    )
    for resource_id, configured_path in (
        ("sample.fastq_r1", fastq_r1),
        ("sample.fastq_r2", fastq_r2),
    ):
        assets.append(
            _inspect_file(
                resource_id=resource_id,
                category="end_to_end",
                configured_path=configured_path,
                repo_root=repo_root,
                required=needs_fastq and paired_configured,
            )
        )

    incomplete_statuses = {"unconfigured", "missing", "empty", "not_file", "unreadable"}
    missing_required = [
        asset["id"]
        for asset in assets
        if asset["required"] and asset["status"] in incomplete_statuses
    ]
    return {
        "status": "ready" if not missing_required else "incomplete",
        "tracks": tracks,
        "repo_root": str(repo_root.resolve()),
        "summary": {
            "assets": len(assets),
            "required": sum(asset["required"] for asset in assets),
            "ready_required": sum(
                asset["required"] and asset["status"] == "ok" for asset in assets
            ),
            "missing_required": len(missing_required),
        },
        "missing_required": missing_required,
        "assets": assets,
    }


def _human_size(size: int | None) -> str:
    if size is None:
        return "-"
    value = float(size)
    units = ("B", "KiB", "MiB", "GiB", "TiB")
    for unit in units:
        if value < 1024 or unit == units[-1]:
            return f"{int(value)} B" if unit == "B" else f"{value:.1f} {unit}"
        value /= 1024
    raise AssertionError("unreachable")


def render_text(report: Mapping[str, Any]) -> str:
    """Render a compact, terminal-friendly report."""

    lines = [
        f"resource preflight: {report['status']}",
        f"tracks: {', '.join(report['tracks'])}",
        "",
        f"{'required':8} {'status':12} {'size':10} {'resource':38} path",
    ]
    for asset in report["assets"]:
        lines.append(
            f"{('yes' if asset['required'] else 'no'):8} "
            f"{asset['status']:12} "
            f"{_human_size(asset['size_bytes']):10} "
            f"{asset['id']:38} "
            f"{asset['path'] or '-'}"
        )
    summary = report["summary"]
    lines.extend(
        [
            "",
            (
                f"required ready: {summary['ready_required']}/"
                f"{summary['required']}; missing/invalid: "
                f"{summary['missing_required']}"
            ),
        ]
    )
    return "\n".join(lines)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Read-only preflight for configured PGBench data resources."
    )
    parser.add_argument(
        "--config",
        type=Path,
        default=Path("config/config.example.yaml"),
    )
    parser.add_argument("--repo-root", type=Path, default=Path.cwd())
    parser.add_argument(
        "--track",
        choices=("configured", "caller-only", "end-to-end", "all"),
        default="configured",
        help="check configured tracks (default) or an explicit track",
    )
    parser.add_argument(
        "--json",
        action="store_true",
        help="emit machine-readable JSON to stdout",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        config = load_config(args.config)
        report = build_inventory(
            config,
            repo_root=args.repo_root,
            track_selection=args.track,
        )
    except (ResourceCheckError, KeyError, TypeError) as exc:
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
