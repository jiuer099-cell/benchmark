#!/usr/bin/env python3
"""Create and optionally execute all frozen coverage/seed benchmark runs."""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
from pathlib import Path
from typing import Any, Sequence

import yaml  # type: ignore[import-untyped]


class MatrixError(ValueError):
    """Raised when the coverage matrix violates the benchmark contract."""


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_mapping(path: Path, label: str) -> dict[str, Any]:
    value = yaml.safe_load(path.read_text(encoding="utf-8")) if path.suffix in {".yaml", ".yml"} else json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise MatrixError(f"{label} must contain a mapping")
    return value


def atomic_yaml(path: Path, value: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp")
    temporary.write_text(yaml.safe_dump(value, sort_keys=False), encoding="utf-8")
    temporary.replace(path)


def build_configs(base: dict[str, Any], manifest: dict[str, Any], output_dir: Path) -> list[Path]:
    contract = base.get("benchmark_contract", {})
    expected_coverages = [value for value in contract.get("coverage_levels", []) if value != "full"]
    expected_seeds = contract.get("downsampling_seeds", [])
    if expected_coverages != [10, 20, 30] or len(expected_seeds) < 3:
        raise MatrixError("base config must freeze 10x/20x/30x and at least three seeds")
    sample = base.get("sample")
    if not isinstance(sample, dict):
        raise MatrixError("base config has no sample mapping")
    source_r1, source_r2 = Path(sample["fastq_r1"]), Path(sample["fastq_r2"])
    if manifest.get("source_fastq_r1_sha256") != sha256(source_r1) or manifest.get("source_fastq_r2_sha256") != sha256(source_r2):
        raise MatrixError("coverage manifest is not bound to the base FASTQs")
    rows = manifest.get("replicates")
    if not isinstance(rows, list):
        raise MatrixError("coverage manifest has no replicate list")
    indexed = {(row.get("coverage_x"), row.get("seed")): row for row in rows if isinstance(row, dict)}
    expected = {(coverage, seed) for coverage in expected_coverages for seed in expected_seeds}
    if set(indexed) != expected:
        raise MatrixError("coverage manifest does not contain the exact frozen matrix")
    run_root = str(base["project"]["run_id"])
    generated: list[Path] = []
    for coverage, seed in sorted(expected):
        row = indexed[(coverage, seed)]
        config = json.loads(json.dumps(base))
        config["project"]["run_id"] = f"{run_root}_{coverage}x_seed{seed}"
        config["sample"].update({
            "coverage_x": coverage,
            "downsampling_seed": seed,
            "read_count": row["read_count"],
            "read_bases": row["read_bases"],
            "fastq_r1": row["fastq_r1"],
            "fastq_r2": row["fastq_r2"],
        })
        destination = output_dir / f"{coverage}x.seed-{seed}.yaml"
        atomic_yaml(destination, config)
        generated.append(destination)
    full = manifest.get("full_depth")
    if not isinstance(full, dict):
        raise MatrixError("coverage manifest has no full-depth entry")
    config = json.loads(json.dumps(base))
    config["project"]["run_id"] = f"{run_root}_full"
    config["sample"].update({
        "coverage_x": full["coverage_x"],
        "downsampling_seed": None,
        "read_count": full["read_count"],
        "read_bases": full["read_bases"],
        "fastq_r1": full["fastq_r1"],
        "fastq_r2": full["fastq_r2"],
    })
    destination = output_dir / "full.yaml"
    atomic_yaml(destination, config)
    generated.append(destination)
    return generated


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base-config", required=True, type=Path)
    parser.add_argument("--coverage-manifest", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--execute", action="store_true")
    parser.add_argument("--cores", type=int, default=1)
    parser.add_argument("--snakemake-arg", action="append", default=[])
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        configs = build_configs(
            load_mapping(args.base_config, "base config"),
            load_mapping(args.coverage_manifest, "coverage manifest"),
            args.output_dir,
        )
        for config in configs:
            print(config)
            if args.execute:
                subprocess.run(
                    [sys.executable, "-m", "snakemake", "--cores", str(args.cores), "--configfile", str(config), *args.snakemake_arg],
                    check=True,
                )
    except (OSError, KeyError, MatrixError, subprocess.CalledProcessError, yaml.YAMLError, json.JSONDecodeError) as exc:
        print(f"run_coverage_matrix: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
