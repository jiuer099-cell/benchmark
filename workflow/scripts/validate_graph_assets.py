#!/usr/bin/env python3
"""Validate a vg graph asset bundle and write a content-addressed lock file."""

from __future__ import annotations

import argparse
import hashlib
import sys
from collections.abc import Mapping
from pathlib import Path
from typing import Any

import yaml  # type: ignore[import-untyped]


ASSET_FILENAMES = {
    "gbz": "graph.gbz",
    "xg": "graph.xg",
    "min": "graph.min",
    "zipcodes": "graph.shortread.zipcodes",
    "dist": "graph.dist",
    "sample_list": "samples.txt",
    "gfa": "graph.gfa.gz",
    "variation_calls": "graph.variation.calls.bed.gz",
}
PROFILE_ASSETS = {
    "vg_gbz_min_dist": {"gbz", "min", "dist", "sample_list"},
    "vg_giraffe_shortread": {
        "gbz", "min", "zipcodes", "dist", "sample_list"
    },
}


class GraphAssetError(ValueError):
    """Raised when a graph bundle cannot be frozen safely."""


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _load_manifest(path: Path) -> dict[str, Any]:
    try:
        with path.open("r", encoding="utf-8") as handle:
            loaded = yaml.safe_load(handle)
    except (OSError, yaml.YAMLError) as exc:
        raise GraphAssetError(f"cannot load source manifest {path}: {exc}") from exc
    if not isinstance(loaded, dict) or not loaded:
        raise GraphAssetError("source graph manifest must be a non-empty YAML mapping")
    return loaded


def _validate_regular_file(path: Path, label: str) -> None:
    if path.is_symlink():
        raise GraphAssetError(f"{label} must not be a symlink: {path}")
    if not path.is_file():
        raise GraphAssetError(f"{label} is not a regular file: {path}")
    if path.stat().st_size == 0:
        raise GraphAssetError(f"{label} is empty: {path}")


def _asset_record(path: Path) -> dict[str, str | int]:
    return {
        "path": str(path),
        "sha256": sha256_file(path),
        "size_bytes": path.stat().st_size,
    }


def _validate_sample_list(path: Path, excluded_samples: tuple[str, ...]) -> int:
    try:
        samples = [
            line.strip()
            for line in path.read_text(encoding="utf-8").splitlines()
            if line.strip() and not line.lstrip().startswith("#")
        ]
    except (OSError, UnicodeError) as exc:
        raise GraphAssetError(f"cannot read sample list {path}: {exc}") from exc
    if not samples:
        raise GraphAssetError("sample_list contains no sample names")
    if len(samples) != len(set(samples)):
        raise GraphAssetError("sample_list contains duplicate sample names")
    if any(any(character.isspace() for character in sample) for sample in samples):
        raise GraphAssetError("sample_list entries must not contain whitespace")
    normalized_samples = {sample.split("#", 1)[0] for sample in samples}
    leaked = sorted(normalized_samples.intersection(excluded_samples))
    if leaked:
        raise GraphAssetError(
            "sample_list contains excluded HG002 truth sample aliases: "
            + ", ".join(leaked)
        )
    return len(samples)


def _validate_declared_lock(
    source: Mapping[str, Any],
    *,
    reference_path: str,
    records: Mapping[str, Mapping[str, Any]],
) -> None:
    """Check identities when the source manifest already declares lock fields."""

    declared_reference = source.get("reference_path")
    if declared_reference is not None and declared_reference != reference_path:
        raise GraphAssetError(
            "source manifest reference_path does not match configured reference_path"
        )
    declared_assets = source.get("assets")
    if declared_assets is None:
        return
    if not isinstance(declared_assets, Mapping):
        raise GraphAssetError("source manifest assets must be a mapping")
    for name, actual in records.items():
        declared = declared_assets.get(name)
        if declared is None:
            continue
        if not isinstance(declared, Mapping):
            raise GraphAssetError(f"source manifest assets.{name} must be a mapping")
        expected_sha = declared.get("sha256")
        if expected_sha is not None and expected_sha != actual["sha256"]:
            raise GraphAssetError(
                f"source manifest checksum mismatch for {name}: "
                f"expected {expected_sha}, observed {actual['sha256']}"
            )
        expected_size = declared.get("size_bytes")
        if expected_size is not None and expected_size != actual["size_bytes"]:
            raise GraphAssetError(
                f"source manifest size mismatch for {name}: "
                f"expected {expected_size}, observed {actual['size_bytes']}"
            )
        expected_filename = declared.get("filename")
        if expected_filename is not None and expected_filename != Path(
            str(actual["path"])
        ).name:
            raise GraphAssetError(
                f"source manifest filename mismatch for {name}: "
                f"expected {expected_filename!r}"
            )


def _validate_declared_exclusions(
    source: Mapping[str, Any], excluded_samples: tuple[str, ...]
) -> None:
    """Require a graph without the benchmark sample when no sample list exists.

    Minigraph rGFA releases do not necessarily ship a separate sample list.
    Their immutable source manifest therefore has to make the leave-out
    declaration explicit before it can be used by a formal benchmark.
    """

    declared = source.get("excluded_samples")
    if not isinstance(declared, list) or not all(
        isinstance(item, str) and item for item in declared
    ):
        raise GraphAssetError(
            "source graph manifest must declare excluded_samples for this profile"
        )
    missing = sorted(set(excluded_samples) - set(declared))
    if missing:
        raise GraphAssetError(
            "source graph manifest does not declare excluded benchmark aliases: "
            + ", ".join(missing)
        )


def lock_graph_assets(
    *,
    source_manifest: Path,
    gbz: Path | None = None,
    xg: Path | None = None,
    min_index: Path | None = None,
    dist: Path | None = None,
    sample_list: Path | None = None,
    reference_path: str,
    zipcodes: Path | None = None,
    gfa: Path | None = None,
    variation_calls: Path | None = None,
    profile: str = "vg_giraffe_shortread",
    excluded_samples: tuple[str, ...] = ("HG002", "NA24385"),
) -> dict[str, Any]:
    """Validate one conventional graph directory and return its lock payload."""

    if not reference_path.strip() or "\n" in reference_path or "\r" in reference_path:
        raise GraphAssetError("reference_path must be a non-empty single-line value")

    if profile not in PROFILE_ASSETS:
        raise GraphAssetError(f"unsupported graph asset profile: {profile}")
    supplied_paths = {
        "gbz": gbz,
        "xg": xg,
        "min": min_index,
        "zipcodes": zipcodes,
        "dist": dist,
        "sample_list": sample_list,
        "gfa": gfa,
        "variation_calls": variation_calls,
    }
    missing = sorted(
        name for name in PROFILE_ASSETS[profile] if supplied_paths.get(name) is None
    )
    if missing:
        raise GraphAssetError(f"graph profile {profile} is missing assets: {missing}")
    paths = {
        name: path
        for name, path in supplied_paths.items()
        if name in PROFILE_ASSETS[profile] and path is not None
    }
    _validate_regular_file(source_manifest, "source manifest")
    source = _load_manifest(source_manifest)

    roots = {path.parent.resolve(strict=True) for path in paths.values()}
    if len(roots) != 1:
        raise GraphAssetError("all graph assets must reside in one directory")
    asset_root = next(iter(roots))
    recorded_asset_root = asset_root
    if source_manifest.parent.resolve(strict=True) != asset_root:
        raise GraphAssetError(
            "source manifest must reside in the graph asset directory"
        )
    if source_manifest.name != "graph-assets.lock.yaml":
        raise GraphAssetError(
            "source manifest must use conventional filename "
            "'graph-assets.lock.yaml'"
        )

    for name, path in paths.items():
        expected_name = (
            "graph.shortread.withzip.min"
            if profile == "vg_giraffe_shortread" and name == "min"
            else ASSET_FILENAMES[name]
        )
        if path.name != expected_name:
            raise GraphAssetError(
                f"{name} must use conventional filename {expected_name!r}, "
                f"found {path.name!r}"
            )
        _validate_regular_file(path, name)

    if not excluded_samples or len(set(excluded_samples)) != len(excluded_samples):
        raise GraphAssetError("excluded_samples must be non-empty and unique")
    sample_count = (
        _validate_sample_list(sample_list, excluded_samples)
        if "sample_list" in paths
        else None
    )
    if "sample_list" not in paths:
        _validate_declared_exclusions(source, excluded_samples)
    records = {name: _asset_record(path) for name, path in paths.items()}
    _validate_declared_lock(
        source,
        reference_path=reference_path,
        records=records,
    )
    return {
        "schema_version": 1,
        "profile": profile,
        "asset_root": str(recorded_asset_root),
        "reference_path": reference_path,
        "excluded_samples": list(excluded_samples),
        "source_manifest": _asset_record(source_manifest),
        "assets": records,
        "sample_count": sample_count,
    }


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Validate a conventional vg graph bundle and write its lock."
    )
    parser.add_argument("--source-manifest", required=True, type=Path)
    parser.add_argument("--gbz", type=Path)
    parser.add_argument("--xg", type=Path)
    parser.add_argument("--min", dest="min_index", type=Path)
    parser.add_argument("--zipcodes", type=Path)
    parser.add_argument("--dist", type=Path)
    parser.add_argument("--sample-list", type=Path)
    parser.add_argument("--gfa", type=Path)
    parser.add_argument("--variation-calls", type=Path)
    parser.add_argument("--reference-path", required=True)
    parser.add_argument(
        "--profile",
        required=True,
        choices=sorted(PROFILE_ASSETS),
    )
    parser.add_argument(
        "--exclude-sample",
        action="append",
        dest="excluded_samples",
        default=[],
    )
    parser.add_argument("--output", required=True, type=Path)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        payload = lock_graph_assets(
            source_manifest=args.source_manifest,
            gbz=args.gbz,
            xg=args.xg,
            min_index=args.min_index,
            zipcodes=args.zipcodes,
            gfa=args.gfa,
            variation_calls=args.variation_calls,
            dist=args.dist,
            sample_list=args.sample_list,
            reference_path=args.reference_path,
            profile=args.profile,
            excluded_samples=tuple(
                args.excluded_samples or ("HG002", "NA24385")
            ),
        )
        args.output.parent.mkdir(parents=True, exist_ok=True)
        temporary = args.output.with_name(f".{args.output.name}.tmp")
        with temporary.open("w", encoding="utf-8") as handle:
            yaml.safe_dump(payload, handle, sort_keys=False, allow_unicode=True)
        temporary.replace(args.output)
    except (OSError, GraphAssetError) as exc:
        print(f"validate_graph_assets: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
