#!/usr/bin/env python3
"""Fetch and freeze the exact GIAB context BED assets used by PGBench."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import os
import sys
import tempfile
import urllib.request
from pathlib import Path, PurePosixPath
from typing import Any

import yaml  # type: ignore[import-untyped]


class StratificationError(ValueError):
    """Raised when the frozen stratification catalogue cannot be verified."""


def _digest(path: Path, algorithm: str) -> str:
    digest = hashlib.new(algorithm)
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _safe_relative_path(value: Any) -> Path:
    if not isinstance(value, str) or not value:
        raise StratificationError("stratum relative_path must be non-empty")
    candidate = PurePosixPath(value)
    if candidate.is_absolute() or ".." in candidate.parts:
        raise StratificationError(f"unsafe stratum relative_path: {value!r}")
    return Path(*candidate.parts)


def _validate_bed(path: Path) -> int:
    opener = gzip.open if path.suffix == ".gz" else open
    count = 0
    with opener(path, "rt", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                raise StratificationError(f"malformed BED row {path}:{line_number}")
            try:
                start, end = int(fields[1]), int(fields[2])
            except ValueError as exc:
                raise StratificationError(
                    f"non-integer BED interval {path}:{line_number}"
                ) from exc
            if start < 0 or end <= start:
                raise StratificationError(f"invalid BED interval {path}:{line_number}")
            count += 1
    if count == 0:
        raise StratificationError(f"stratification BED is empty: {path}")
    return count


def _download(url: str, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{destination.name}.", suffix=".download", dir=destination.parent
    )
    os.close(descriptor)
    temporary = Path(temporary_name)
    try:
        request = urllib.request.Request(url, headers={"User-Agent": "PGBench/1"})
        with urllib.request.urlopen(request, timeout=120) as response:  # noqa: S310
            with temporary.open("wb") as output:
                while chunk := response.read(1024 * 1024):
                    output.write(chunk)
        os.replace(temporary, destination)
    finally:
        temporary.unlink(missing_ok=True)


def prepare(catalogue_path: Path, output: Path, *, download: bool) -> dict[str, Any]:
    try:
        catalogue = yaml.safe_load(catalogue_path.read_text(encoding="utf-8"))
    except (OSError, yaml.YAMLError) as exc:
        raise StratificationError(f"cannot load catalogue: {catalogue_path}") from exc
    if not isinstance(catalogue, dict) or catalogue.get("schema_version") != 2:
        raise StratificationError("stratification catalogue schema_version must be 2")
    source_root = catalogue.get("source_root")
    local_root = catalogue.get("local_root")
    strata = catalogue.get("strata")
    if not isinstance(source_root, str) or not source_root.startswith("https://"):
        raise StratificationError("source_root must be an HTTPS URL")
    if not isinstance(local_root, str) or not local_root:
        raise StratificationError("local_root must be a path")
    if not isinstance(strata, dict) or not strata:
        raise StratificationError("catalogue must define strata")

    assets: dict[str, dict[str, Any]] = {}
    for stratum, raw in sorted(strata.items()):
        if not isinstance(stratum, str) or not isinstance(raw, dict):
            raise StratificationError("strata must map names to objects")
        relative = _safe_relative_path(raw.get("relative_path"))
        expected_md5 = raw.get("md5")
        if not isinstance(expected_md5, str) or len(expected_md5) != 32:
            raise StratificationError(f"{stratum} has no frozen MD5")
        path = Path(local_root) / relative
        url = f"{source_root.rstrip('/')}/{relative.as_posix()}"
        if not path.is_file():
            if not download:
                raise StratificationError(f"missing stratification BED: {path}")
            _download(url, path)
        actual_md5 = _digest(path, "md5")
        if actual_md5 != expected_md5:
            raise StratificationError(
                f"{stratum} MD5 mismatch: expected {expected_md5}, got {actual_md5}"
            )
        assets[stratum] = {
            "category": raw.get("category"),
            "path": str(path),
            "source_url": url,
            "size_bytes": path.stat().st_size,
            "interval_count": _validate_bed(path),
            "md5": actual_md5,
            "sha256": _digest(path, "sha256"),
        }

    payload = {
        "schema_version": 1,
        "contract": "pgbench_context_stratifications_v1",
        "profile_id": catalogue.get("profile_id"),
        "reference_id": catalogue.get("reference_id"),
        "release": catalogue.get("release"),
        "catalogue_sha256": _digest(catalogue_path, "sha256"),
        "assets": assets,
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_name(f".{output.name}.tmp")
    temporary.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    temporary.replace(output)
    return payload


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--catalogue", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--download", action="store_true")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        prepare(args.catalogue, args.output, download=args.download)
    except (OSError, StratificationError) as exc:
        print(f"prepare_stratifications: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
