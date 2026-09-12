#!/usr/bin/env python3
"""Create a byte-addressed lock for a Frozen Haplotype Source Bundle."""

from __future__ import annotations

import argparse
import hashlib
import json
import re
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_asset(value: str) -> tuple[str, Path]:
    name, separator, raw_path = value.partition("=")
    if not separator or not name or not raw_path:
        raise ValueError("--asset must be NAME=PATH")
    path = Path(raw_path).resolve(strict=True)
    if not path.is_file():
        raise ValueError(f"bundle asset is not a regular file: {path}")
    return name, path


def parse_software(value: str) -> dict[str, str]:
    """Parse a versioned, byte-addressed build-tool declaration.

    A bare version string cannot establish reproducibility: the executable or
    container digest is part of the Frozen Haplotype Source Bundle as well.
    """

    name, separator, remainder = value.partition("=")
    version, separator2, digest = remainder.partition("=")
    if not separator or not separator2 or not name or not version:
        raise ValueError("--software must be NAME=VERSION=SHA256")
    if not re.fullmatch(r"[0-9a-f]{64}", digest):
        raise ValueError("--software SHA256 must be 64 lowercase hexadecimal characters")
    return {"name": name, "version": version, "sha256": digest}


def freeze(*, bundle_id: str, release: str, reference: str, excluded_samples: list[str],
           assets: list[tuple[str, Path]], software: list[dict[str, str]]) -> dict[str, object]:
    names = [name for name, _ in assets]
    if len(names) != len(set(names)):
        raise ValueError("bundle assets must have unique names")
    required = {
        "gfa_or_gbz",
        "population_vcf",
        "sample_roster",
        "haplotype_roster",
        "family_exclusion_manifest",
        "reference",
    }
    missing = sorted(required - set(names))
    if missing:
        raise ValueError("bundle lock is missing required assets: " + ", ".join(missing))
    return {
        "schema_version": 1,
        "contract": "pgbench_frozen_haplotype_source_bundle_v1",
        "bundle_id": bundle_id,
        "release": release,
        "reference": reference,
        "family_exclusion": sorted(excluded_samples),
        "assets": {
            name: {"path": str(path), "sha256": sha256(path), "size_bytes": path.stat().st_size}
            for name, path in sorted(assets)
        },
        "software": sorted(software, key=lambda item: (item["name"], item["version"])),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--bundle-id", required=True)
    parser.add_argument("--release", required=True)
    parser.add_argument("--reference", required=True)
    parser.add_argument("--exclude-sample", action="append", default=[])
    parser.add_argument("--asset", action="append", default=[])
    parser.add_argument("--software", action="append", default=[])
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    try:
        payload = freeze(
            bundle_id=args.bundle_id,
            release=args.release,
            reference=args.reference,
            excluded_samples=args.exclude_sample,
            assets=[parse_asset(value) for value in args.asset],
            software=[parse_software(value) for value in args.software],
        )
        args.output.parent.mkdir(parents=True, exist_ok=True)
        temporary = args.output.with_name(f".{args.output.name}.tmp")
        temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        temporary.replace(args.output)
    except (OSError, ValueError) as exc:
        parser.error(str(exc))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
