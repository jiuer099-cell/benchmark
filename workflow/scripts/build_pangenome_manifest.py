#!/usr/bin/env python3
"""Assign stable pangenome allele IDs and write a manifest/ledger."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import re
import sys
from collections.abc import Iterable
from datetime import UTC, datetime
from pathlib import Path
from typing import TextIO, cast

import yaml  # type: ignore[import-untyped]

GRAPH_CLASSES = {
    "simple_biallelic",
    "multiallelic",
    "nested_or_overlapping",
    "repeat_ambiguous",
}


class PangenomeManifestError(ValueError):
    """Raised when a population panel cannot be converted safely."""


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_graph_assets_lock(path: Path) -> dict:
    """Load a graph lock and re-verify every frozen asset before embedding it."""

    try:
        with path.open("r", encoding="utf-8") as handle:
            lock = yaml.safe_load(handle)
    except (OSError, yaml.YAMLError) as exc:
        raise PangenomeManifestError(
            f"cannot load graph assets lock {path}: {exc}"
        ) from exc
    if not isinstance(lock, dict):
        raise PangenomeManifestError("graph assets lock must contain a YAML mapping")
    if lock.get("schema_version") != 1:
        raise PangenomeManifestError("unsupported graph assets lock schema_version")
    assets = lock.get("assets")
    required_assets = {"gbz", "xg", "min", "dist", "sample_list"}
    if not isinstance(assets, dict) or set(assets) != required_assets:
        raise PangenomeManifestError(
            "graph assets lock must contain exactly gbz, xg, min, dist, sample_list"
        )
    reference_path = lock.get("reference_path")
    asset_root = lock.get("asset_root")
    sample_count = lock.get("sample_count")
    excluded_samples = lock.get("excluded_samples")
    if not isinstance(reference_path, str) or not reference_path:
        raise PangenomeManifestError("graph assets lock has no reference_path")
    if not isinstance(asset_root, str) or not asset_root:
        raise PangenomeManifestError("graph assets lock has no asset_root")
    if not isinstance(sample_count, int) or sample_count < 1:
        raise PangenomeManifestError("graph assets lock has invalid sample_count")
    if (
        not isinstance(excluded_samples, list)
        or "HG002" not in excluded_samples
        or len(excluded_samples) != len(set(excluded_samples))
    ):
        raise PangenomeManifestError(
            "graph assets lock has invalid excluded_samples"
        )
    source_record = lock.get("source_manifest")
    if not isinstance(source_record, dict):
        raise PangenomeManifestError("graph assets lock has no source_manifest")
    source_path_raw = source_record.get("path")
    if not isinstance(source_path_raw, str) or not source_path_raw:
        raise PangenomeManifestError("graph source manifest has no path")
    source_path = Path(source_path_raw)
    if not source_path.is_file():
        raise PangenomeManifestError(
            f"graph source manifest does not exist: {source_path}"
        )
    source_sha = sha256_file(source_path)
    source_size = source_path.stat().st_size
    if (
        source_record.get("sha256") != source_sha
        or source_record.get("size_bytes") != source_size
    ):
        raise PangenomeManifestError(
            "graph source manifest no longer matches its content identity"
        )
    try:
        source_path.resolve(strict=True).relative_to(
            Path(asset_root).resolve(strict=True)
        )
    except (OSError, ValueError) as exc:
        raise PangenomeManifestError(
            "graph source manifest is outside asset_root"
        ) from exc

    verified: dict[str, dict[str, str | int]] = {}
    for name in sorted(required_assets):
        record = assets[name]
        if not isinstance(record, dict):
            raise PangenomeManifestError(f"graph asset {name} record must be a mapping")
        raw_path = record.get("path")
        expected_sha = record.get("sha256")
        expected_size = record.get("size_bytes")
        if not isinstance(raw_path, str) or not raw_path:
            raise PangenomeManifestError(f"graph asset {name} has no path")
        asset_path = Path(raw_path)
        if not asset_path.is_file():
            raise PangenomeManifestError(
                f"locked graph asset {name} does not exist: {asset_path}"
            )
        observed_sha = sha256_file(asset_path)
        observed_size = asset_path.stat().st_size
        if observed_sha != expected_sha or observed_size != expected_size:
            raise PangenomeManifestError(
                f"locked graph asset {name} no longer matches its content identity"
            )
        try:
            asset_path.resolve(strict=True).relative_to(
                Path(asset_root).resolve(strict=True)
            )
        except (OSError, ValueError) as exc:
            raise PangenomeManifestError(
                f"locked graph asset {name} is outside asset_root"
            ) from exc
        verified[name] = {
            "path": raw_path,
            "sha256": observed_sha,
            "size_bytes": observed_size,
        }

    return {
        "manifest": {
            "path": source_path_raw,
            "sha256": source_sha,
            "size_bytes": source_size,
        },
        "lock_manifest": {
            "path": str(path),
            "sha256": sha256_file(path),
            "size_bytes": path.stat().st_size,
        },
        "asset_root": asset_root,
        "reference_path": reference_path,
        "excluded_samples": excluded_samples,
        "gbz": verified["gbz"],
        "xg": verified["xg"],
        "min": verified["min"],
        "dist": verified["dist"],
        "sample_list": verified["sample_list"],
        "sample_count": sample_count,
    }


def _open_text(path: Path, mode: str) -> TextIO:
    if path.suffix == ".gz":
        return cast(TextIO, gzip.open(path, mode, encoding="utf-8"))
    return cast(TextIO, path.open(mode, encoding="utf-8"))


def parse_info(raw_info: str) -> dict[str, str | bool]:
    info: dict[str, str | bool] = {}
    if raw_info in {"", "."}:
        return info
    for item in raw_info.split(";"):
        if not item:
            continue
        if "=" in item:
            key, value = item.split("=", 1)
            info[key] = value
        else:
            info[item] = True
    return info


def format_info(info: dict[str, str | bool]) -> str:
    fields: list[str] = []
    for key, value in info.items():
        if value is True:
            fields.append(key)
        elif value is not False:
            fields.append(f"{key}={value}")
    return ";".join(fields) if fields else "."


def infer_svtype(alt: str, info: dict[str, str | bool]) -> str:
    declared = info.get("SVTYPE")
    if isinstance(declared, str) and declared:
        return declared.upper()
    symbolic = re.fullmatch(r"<([^>]+)>", alt)
    if symbolic:
        return symbolic.group(1).upper()
    if "[" in alt or "]" in alt:
        return "BND"
    return "INS" if len(alt) > 1 else "OTHER"


def _integer_info(info: dict[str, str | bool], key: str, *, default: int) -> int:
    value = info.get(key)
    if value is None or value is True:
        return default
    first = str(value).split(",", 1)[0]
    try:
        return int(first)
    except ValueError as exc:
        raise PangenomeManifestError(
            f"{key} must be an integer, found {value!r}"
        ) from exc


def stable_allele_id(
    *,
    namespace: str,
    chrom: str,
    pos: int,
    ref: str,
    alt: str,
    svtype: str,
    end: int,
    svlen: int,
) -> str:
    canonical = json.dumps(
        {
            "alt": alt,
            "chrom": chrom,
            "end": end,
            "pos": pos,
            "ref": ref,
            "svlen": svlen,
            "svtype": svtype,
        },
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8")
    return f"{namespace}_{hashlib.sha256(canonical).hexdigest()[:20]}"


def _format_af(info: dict[str, str | bool]) -> str:
    value = info.get("AF")
    if value is None or value is True:
        return "."
    first = str(value).split(",", 1)[0]
    try:
        numeric = float(first)
    except ValueError as exc:
        raise PangenomeManifestError(f"AF must be numeric, found {value!r}") from exc
    if numeric < 0 or numeric > 1:
        raise PangenomeManifestError(f"AF must be in [0, 1], found {numeric}")
    return f"{numeric:.12g}"


def assign_stable_alleles(
    input_vcf: Path,
    output_vcf: Path,
    ledger_path: Path,
    *,
    namespace: str,
    excluded_samples: Iterable[str] = (),
) -> int:
    """Write an ID-annotated panel and allele ledger; return record count."""

    output_vcf.parent.mkdir(parents=True, exist_ok=True)
    ledger_path.parent.mkdir(parents=True, exist_ok=True)
    seen_ids: set[str] = set()
    record_count = 0
    info_header_written = False
    excluded = set(excluded_samples)

    with (
        _open_text(input_vcf, "rt") as source,
        _open_text(output_vcf, "wt") as panel,
        ledger_path.open("w", encoding="utf-8") as ledger,
    ):
        ledger.write(
            "\t".join(
                [
                    "PANGENOME_ALLELE_ID",
                    "CHROM",
                    "POS",
                    "END",
                    "SVTYPE",
                    "SVLEN",
                    "REF",
                    "ALT",
                    "AF",
                    "GRAPH_COMPLEXITY_CLASS",
                    "SOURCE_RECORD_ID",
                ]
            )
            + "\n"
        )
        for raw_line in source:
            if raw_line.startswith("##"):
                if raw_line.startswith("##INFO=<ID=PANGENOME_ALLELE_ID,"):
                    info_header_written = True
                panel.write(raw_line)
                continue
            if raw_line.startswith("#CHROM"):
                header_fields = raw_line.rstrip("\n").split("\t")
                leaked = sorted(excluded.intersection(header_fields[9:]))
                if leaked:
                    raise PangenomeManifestError(
                        "population VCF contains excluded HG002 truth sample aliases: "
                        + ", ".join(leaked)
                    )
                if not info_header_written:
                    panel.write(
                        "##INFO=<ID=PANGENOME_ALLELE_ID,Number=1,Type=String,"
                        'Description="Stable population pangenome allele ID">\n'
                    )
                panel.write(raw_line)
                continue
            if not raw_line.strip():
                continue

            fields = raw_line.rstrip("\n").split("\t")
            if len(fields) < 8:
                raise PangenomeManifestError(
                    f"VCF record must have at least 8 fields: {raw_line!r}"
                )
            chrom, pos_raw, source_id, ref, alt = fields[:5]
            try:
                pos = int(pos_raw)
            except ValueError as exc:
                raise PangenomeManifestError(
                    f"VCF POS must be integer, found {pos_raw!r}"
                ) from exc
            info = parse_info(fields[7])
            svtype = infer_svtype(alt, info)
            end = _integer_info(info, "END", default=pos)
            default_svlen = len(alt) - len(ref)
            svlen = _integer_info(info, "SVLEN", default=default_svlen)
            graph_class_raw = info.get("GRAPH_COMPLEXITY", "simple_biallelic")
            graph_class = str(graph_class_raw)
            if graph_class not in GRAPH_CLASSES:
                raise PangenomeManifestError(
                    f"unsupported GRAPH_COMPLEXITY={graph_class!r}"
                )
            allele_id = stable_allele_id(
                namespace=namespace,
                chrom=chrom,
                pos=pos,
                ref=ref,
                alt=alt,
                svtype=svtype,
                end=end,
                svlen=svlen,
            )
            if allele_id in seen_ids:
                raise PangenomeManifestError(
                    f"duplicate canonical allele generated ID {allele_id}"
                )
            seen_ids.add(allele_id)
            info["PANGENOME_ALLELE_ID"] = allele_id
            fields[7] = format_info(info)
            panel.write("\t".join(fields) + "\n")
            ledger.write(
                "\t".join(
                    [
                        allele_id,
                        chrom,
                        str(pos),
                        str(end),
                        svtype,
                        str(svlen),
                        ref,
                        alt,
                        _format_af(info),
                        graph_class,
                        source_id,
                    ]
                )
                + "\n"
            )
            record_count += 1

    if record_count == 0:
        raise PangenomeManifestError("population panel contains no variant records")
    return record_count


def build_manifest(
    *,
    pangenome_id: str,
    backbone_id: str,
    reference: Path,
    population_source_id: str,
    population_source: Path,
    panel_vcf: Path,
    allele_ledger: Path,
    namespace: str,
    excluded_truth_samples: Iterable[str],
    graph_build_recipe_sha256: str | None,
    graph_assets_lock: Path | None,
    generated_at: str,
) -> dict:
    excluded = list(excluded_truth_samples)
    if "HG002" not in excluded:
        raise PangenomeManifestError("truth_samples_excluded must contain HG002")
    with _open_text(panel_vcf, "rt") as handle:
        record_count = sum(
            1 for line in handle if line.strip() and not line.startswith("#")
        )
    return {
        "schema_version": 1,
        "pangenome_id": pangenome_id,
        "backbone_reference": {
            "id": backbone_id,
            "path": str(reference),
            "sha256": sha256_file(reference),
        },
        "population_sources": [
            {
                "id": population_source_id,
                "path": str(population_source),
                "sha256": sha256_file(population_source),
            }
        ],
        "allele_id_namespace": namespace,
        "truth_samples_excluded": excluded,
        "panel_vcf": {
            "path": str(panel_vcf),
            "sha256": sha256_file(panel_vcf),
            "record_count": record_count,
        },
        "allele_ledger": {
            "path": str(allele_ledger),
            "sha256": sha256_file(allele_ledger),
        },
        "graph_assets": (
            load_graph_assets_lock(graph_assets_lock)
            if graph_assets_lock is not None
            else None
        ),
        "graph_build_recipe_sha256": graph_build_recipe_sha256,
        "generated_at": generated_at,
    }


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Assign stable pangenome allele IDs and write a manifest."
    )
    parser.add_argument("--pangenome-id", required=True)
    parser.add_argument("--backbone-id", required=True)
    parser.add_argument("--reference", required=True, type=Path)
    parser.add_argument("--population-source-id", required=True)
    parser.add_argument("--population-vcf", required=True, type=Path)
    parser.add_argument("--output-panel", required=True, type=Path)
    parser.add_argument("--output-ledger", required=True, type=Path)
    parser.add_argument("--output-manifest", required=True, type=Path)
    parser.add_argument("--allele-namespace", default="PGSV")
    parser.add_argument(
        "--exclude-truth-sample",
        action="append",
        default=[],
    )
    parser.add_argument("--graph-build-recipe-sha256")
    parser.add_argument("--graph-assets-lock", type=Path)
    parser.add_argument("--generated-at")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    generated_at = args.generated_at or datetime.now(UTC).isoformat()
    try:
        excluded_samples = args.exclude_truth_sample or ["HG002", "NA24385"]
        assign_stable_alleles(
            args.population_vcf,
            args.output_panel,
            args.output_ledger,
            namespace=args.allele_namespace,
            excluded_samples=excluded_samples,
        )
        manifest = build_manifest(
            pangenome_id=args.pangenome_id,
            backbone_id=args.backbone_id,
            reference=args.reference,
            population_source_id=args.population_source_id,
            population_source=args.population_vcf,
            panel_vcf=args.output_panel,
            allele_ledger=args.output_ledger,
            namespace=args.allele_namespace,
            excluded_truth_samples=excluded_samples,
            graph_build_recipe_sha256=args.graph_build_recipe_sha256,
            graph_assets_lock=args.graph_assets_lock,
            generated_at=generated_at,
        )
        args.output_manifest.parent.mkdir(parents=True, exist_ok=True)
        temporary = args.output_manifest.with_name(f".{args.output_manifest.name}.tmp")
        with temporary.open("w", encoding="utf-8") as handle:
            yaml.safe_dump(manifest, handle, sort_keys=False, allow_unicode=True)
        temporary.replace(args.output_manifest)
    except (OSError, PangenomeManifestError) as exc:
        print(f"build_pangenome_manifest: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
