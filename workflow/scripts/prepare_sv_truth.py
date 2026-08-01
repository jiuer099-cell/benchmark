#!/usr/bin/env python3
"""Materialize one frozen, indexed SV-only truth universe."""

from __future__ import annotations

import argparse
import bisect
import gzip
import hashlib
import json
import re
import shutil
import subprocess
import sys
import tempfile
from collections import defaultdict
from pathlib import Path
from typing import Any, TextIO

import yaml  # type: ignore[import-untyped]

from build_pangenome_manifest import default_variant_end
from sv_matching import (
    SvMatchError,
    infer_svtype,
    load_evaluator_profile,
    parse_info,
)


class TruthPreparationError(ValueError):
    """Raised when the primary truth cannot be materialized unambiguously."""


SHA256_RE = re.compile(r"^[0-9a-f]{64}$")


def open_text(path: Path, mode: str) -> TextIO:
    if path.suffix.casefold() == ".gz":
        return gzip.open(path, mode, encoding="utf-8")
    return path.open(mode, encoding="utf-8")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_regions(path: Path) -> dict[str, tuple[list[int], list[int]]]:
    raw: dict[str, list[tuple[int, int]]] = defaultdict(list)
    with path.open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                raise TruthPreparationError(
                    f"malformed BED record at {path}:{line_number}"
                )
            start, end = int(fields[1]), int(fields[2])
            if start < 0 or end <= start:
                raise TruthPreparationError(
                    f"invalid BED interval at {path}:{line_number}"
                )
            raw[fields[0]].append((start, end))
    if not raw:
        raise TruthPreparationError(f"benchmark BED is empty: {path}")
    merged: dict[str, tuple[list[int], list[int]]] = {}
    for contig, intervals in raw.items():
        compact: list[list[int]] = []
        for start, end in sorted(intervals):
            if compact and start <= compact[-1][1]:
                compact[-1][1] = max(compact[-1][1], end)
            else:
                compact.append([start, end])
        merged[contig] = (
            [interval[0] for interval in compact],
            [interval[1] for interval in compact],
        )
    return merged


def fully_contained(
    regions: dict[str, tuple[list[int], list[int]]],
    contig: str,
    start: int,
    end: int,
) -> bool:
    index = regions.get(contig)
    if index is None:
        return False
    starts, ends = index
    candidate = bisect.bisect_right(starts, start) - 1
    return candidate >= 0 and ends[candidate] >= end


def _integer_values(value: str | bool | None) -> list[int]:
    if value in {None, True, "", "."}:
        return []
    try:
        return [int(float(token)) for token in str(value).split(",")]
    except ValueError as exc:
        raise TruthPreparationError(f"invalid integer INFO value: {value!r}") from exc


def variant_size(fields: list[str], info: dict[str, str | bool]) -> int:
    values = _integer_values(info.get("SVLEN"))
    if values:
        return max(abs(value) for value in values)
    ref = fields[3]
    alts = fields[4].split(",")
    resolved = [
        abs(len(alt) - len(ref))
        for alt in alts
        if alt not in {"", "."}
        and not alt.startswith("<")
        and "[" not in alt
        and "]" not in alt
    ]
    if resolved:
        return max(resolved)
    end_values = _integer_values(info.get("END"))
    if end_values:
        return max(abs(value - int(fields[1])) for value in end_values)
    return 0


def _stable_truth_id(fields: list[str]) -> str:
    payload = "\0".join((fields[0], fields[1], fields[3], fields[4], fields[7]))
    return f"TRUTH_{hashlib.sha256(payload.encode('utf-8')).hexdigest()[:20]}"


def load_filter_contract(
    catalog_path: Path,
    profile_id: str,
) -> dict[str, Any]:
    try:
        catalog = yaml.safe_load(catalog_path.read_text(encoding="utf-8"))
        profile = catalog["truthsets"][profile_id]
        filters = profile["filters"]
        files = profile["files"]
    except (OSError, yaml.YAMLError, KeyError, TypeError) as exc:
        raise TruthPreparationError(
            f"cannot load truth filter contract {profile_id!r}"
        ) from exc
    allowed = filters.get("allowed_svtypes")
    if not isinstance(allowed, list) or not allowed:
        raise TruthPreparationError("truth filters.allowed_svtypes must be non-empty")
    contract = {
        "pass_only": filters.get("pass_only") is True,
        "inside_benchmark_bed": filters.get("inside_benchmark_bed") is True,
        "region_policy": filters.get("region_policy"),
        "minimum_sv_size": filters.get("minimum_sv_size"),
        "maximum_sv_size": filters.get("maximum_sv_size"),
        "allowed_svtypes": [str(value) for value in allowed],
        "reject_multiallelic_sv": (
            filters.get("reject_multiallelic_sv") is True
        ),
        "multiallelic_policy": filters.get("multiallelic_policy"),
        "required_sha256_before_execution": (
            profile.get("required_sha256_before_execution") is True
        ),
        "expected_sha256": {
            "source_vcf": files["vcf"].get("sha256"),
            "source_index": files["vcf_index"].get("sha256"),
            "benchmark_bed": files["benchmark_bed"].get("sha256"),
        },
    }
    if not isinstance(contract["minimum_sv_size"], int):
        raise TruthPreparationError("truth minimum_sv_size must be an integer")
    maximum_sv_size = contract["maximum_sv_size"]
    if maximum_sv_size is not None and (
        not isinstance(maximum_sv_size, int)
        or maximum_sv_size < contract["minimum_sv_size"]
    ):
        raise TruthPreparationError(
            "truth maximum_sv_size must be an integer greater than or equal "
            "to minimum_sv_size"
        )
    if not contract["inside_benchmark_bed"]:
        raise TruthPreparationError(
            "formal truth must be restricted to the benchmark BED"
        )
    if contract["region_policy"] != "fully_contained":
        raise TruthPreparationError(
            "truth filters.region_policy must be fully_contained"
        )
    if contract["multiallelic_policy"] != "exclude":
        raise TruthPreparationError(
            "truth filters.multiallelic_policy must be exclude"
        )
    if contract["required_sha256_before_execution"]:
        missing = sorted(
            name
            for name, digest in contract["expected_sha256"].items()
            if not isinstance(digest, str) or SHA256_RE.fullmatch(digest) is None
        )
        if missing:
            raise TruthPreparationError(
                "formal truth catalog requires lowercase SHA-256 values for: "
                + ", ".join(missing)
            )
    return contract


def verify_frozen_inputs(
    *,
    source_vcf: Path,
    source_index: Path,
    benchmark_bed: Path,
    contract: dict[str, Any],
) -> dict[str, str]:
    observed = {
        "source_vcf": sha256_file(source_vcf),
        "source_index": sha256_file(source_index),
        "benchmark_bed": sha256_file(benchmark_bed),
    }
    if contract["required_sha256_before_execution"]:
        expected = contract["expected_sha256"]
        mismatches = [
            f"{name}: expected {expected[name]}, observed {observed[name]}"
            for name in sorted(observed)
            if observed[name] != expected[name]
        ]
        if mismatches:
            raise TruthPreparationError(
                "formal truth asset checksum mismatch: " + "; ".join(mismatches)
            )
    return observed


def verify_universe_contract(
    contract: dict[str, Any],
    evaluator_profile: dict[str, Any],
) -> None:
    """Fail closed if truth preparation and evaluation define different events."""

    universe = evaluator_profile["universe"]
    expected = {
        "pass_only": universe["truth_filter_policy"] == "pass_only",
        "inside_benchmark_bed": universe["region_policy"] == "fully_contained",
        "region_policy": universe["region_policy"],
        "minimum_sv_size": universe["minimum_sv_size"],
        "maximum_sv_size": universe["maximum_sv_size"],
        "allowed_svtypes": sorted(universe["allowed_svtypes"]),
        "reject_multiallelic_sv": bool(universe["biallelic_only"]),
        "multiallelic_policy": universe["multiallelic_policy"],
    }
    observed = {
        "pass_only": contract["pass_only"],
        "inside_benchmark_bed": contract["inside_benchmark_bed"],
        "region_policy": contract["region_policy"],
        "minimum_sv_size": contract["minimum_sv_size"],
        "maximum_sv_size": contract["maximum_sv_size"],
        "allowed_svtypes": sorted(contract["allowed_svtypes"]),
        "reject_multiallelic_sv": contract["reject_multiallelic_sv"],
        "multiallelic_policy": contract["multiallelic_policy"],
    }
    if observed != expected:
        raise TruthPreparationError(
            "truth catalog filters differ from the frozen evaluator universe: "
            f"catalog={observed!r}, evaluator={expected!r}"
        )


def materialize_plain_truth(
    *,
    source_vcf: Path,
    benchmark_bed: Path,
    output_vcf: Path,
    contract: dict[str, Any],
) -> dict[str, int]:
    """Filter records and write a plain VCF for later BGZF materialization."""

    regions = load_regions(benchmark_bed)
    allowed = set(contract["allowed_svtypes"])
    counts = {
        "source_records": 0,
        "excluded_non_pass": 0,
        "excluded_outside_bed": 0,
        "excluded_below_minimum_size": 0,
        "excluded_above_maximum_size": 0,
        "excluded_svtype": 0,
        "excluded_multiallelic": 0,
        "excluded_duplicate": 0,
        "eligible_records": 0,
    }
    saw_header = False
    inserted_contract = False
    seen_ids: set[str] = set()
    with open_text(source_vcf, "rt") as source, output_vcf.open(
        "w", encoding="utf-8"
    ) as output:
        for line_number, line in enumerate(source, start=1):
            if line.startswith("##"):
                output.write(line)
                continue
            if line.startswith("#CHROM"):
                output.write(
                    "##pgbench_truth_universe=PASS,"
                    "FULLY_CONTAINED_BED,SVTYPE,"
                    f"{contract['minimum_sv_size']}<=SVLEN<="
                    f"{contract['maximum_sv_size']}\n"
                )
                output.write(line)
                saw_header = True
                inserted_contract = True
                continue
            if not line.strip():
                continue
            if not saw_header:
                raise TruthPreparationError("truth VCF has no #CHROM header")
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                raise TruthPreparationError(
                    f"malformed VCF record at {source_vcf}:{line_number}"
                )
            counts["source_records"] += 1
            # VCF permits both PASS and "." for records without a failing
            # filter. GIAB assembly-based draft truth uses "." for some
            # unflagged records. Accept those source records, but normalize
            # them to PASS in the materialized truth so every downstream
            # evaluator observes the frozen truth_filter_policy: pass_only
            # contract literally and reproducibly.
            if contract["pass_only"] and fields[6] not in {"PASS", "."}:
                counts["excluded_non_pass"] += 1
                continue
            if contract["pass_only"] and fields[6] == ".":
                fields[6] = "PASS"
            info = parse_info(fields[7])
            svtype = infer_svtype(fields[4], info, ref=fields[3])
            if svtype not in allowed:
                counts["excluded_svtype"] += 1
                continue
            size = variant_size(fields, info)
            if size < int(contract["minimum_sv_size"]):
                counts["excluded_below_minimum_size"] += 1
                continue
            if size > int(contract["maximum_sv_size"]):
                counts["excluded_above_maximum_size"] += 1
                continue
            if "," in fields[4] and contract["reject_multiallelic_sv"]:
                counts["excluded_multiallelic"] += 1
                continue
            pos = int(fields[1])
            end_values = _integer_values(info.get("END"))
            end = (
                max(end_values)
                if end_values
                else default_variant_end(
                    pos=pos,
                    ref=fields[3],
                    alt=fields[4],
                    svtype=svtype,
                )
            )
            if not fully_contained(
                regions,
                fields[0],
                pos - 1,
                max(pos, end),
            ):
                counts["excluded_outside_bed"] += 1
                continue
            if fields[2] in {"", "."}:
                fields[2] = _stable_truth_id(fields)
            if fields[2] in seen_ids:
                counts["excluded_duplicate"] += 1
                continue
            seen_ids.add(fields[2])
            output.write("\t".join(fields) + "\n")
            counts["eligible_records"] += 1
    if not inserted_contract:
        raise TruthPreparationError("truth VCF has no column header")
    if counts["eligible_records"] == 0:
        raise TruthPreparationError("SV-only truth universe is empty")
    return counts


def prepare(args: argparse.Namespace) -> dict[str, Any]:
    contract = load_filter_contract(args.truth_catalog, args.profile_id)
    evaluator_profile = load_evaluator_profile(args.evaluator_profile)
    # The evaluator profile is the normative common-support envelope. Older
    # truth catalogs may not yet carry an upper bound, so materialize and audit
    # the effective maximum from the frozen evaluator contract. If a catalog
    # does declare one, verify_universe_contract still fails on disagreement.
    if contract["maximum_sv_size"] is None:
        contract["maximum_sv_size"] = evaluator_profile["universe"][
            "maximum_sv_size"
        ]
    verify_universe_contract(contract, evaluator_profile)
    observed_hashes = verify_frozen_inputs(
        source_vcf=args.source_vcf,
        source_index=args.source_index,
        benchmark_bed=args.benchmark_bed,
        contract=contract,
    )
    args.output_vcf.parent.mkdir(parents=True, exist_ok=True)
    args.audit_json.parent.mkdir(parents=True, exist_ok=True)
    work = Path(tempfile.mkdtemp(prefix=".truth-", dir=args.output_vcf.parent))
    try:
        plain = work / "truth.vcf"
        counts = materialize_plain_truth(
            source_vcf=args.source_vcf,
            benchmark_bed=args.benchmark_bed,
            output_vcf=plain,
            contract=contract,
        )
        staged = work / "truth.vcf.gz"
        subprocess.run(
            ["bcftools", "sort", "-Oz", "-o", str(staged), str(plain)],
            check=True,
        )
        subprocess.run(
            ["bcftools", "index", "--tbi", "--force", str(staged)],
            check=True,
        )
        staged_index = Path(f"{staged}.tbi")
        output_index = Path(f"{args.output_vcf}.tbi")
        shutil.copyfile(staged, args.output_vcf)
        shutil.copyfile(staged_index, output_index)
        audit = {
            "schema_version": 1,
            "profile_id": args.profile_id,
            "contract": contract,
            "counts": counts,
            "source_vcf": str(args.source_vcf),
            "source_vcf_sha256": observed_hashes["source_vcf"],
            "source_index": str(args.source_index),
            "source_index_sha256": observed_hashes["source_index"],
            "benchmark_bed": str(args.benchmark_bed),
            "benchmark_bed_sha256": observed_hashes["benchmark_bed"],
            "truth_catalog": str(args.truth_catalog),
            "truth_catalog_sha256": sha256_file(args.truth_catalog),
            "evaluator_profile": str(args.evaluator_profile),
            "evaluator_profile_sha256": evaluator_profile["_sha256"],
            "output_vcf": str(args.output_vcf),
            "output_vcf_sha256": sha256_file(args.output_vcf),
            "output_index_sha256": sha256_file(output_index),
        }
        temporary = args.audit_json.with_name(f".{args.audit_json.name}.tmp")
        temporary.write_text(
            json.dumps(audit, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        temporary.replace(args.audit_json)
        return audit
    finally:
        shutil.rmtree(work, ignore_errors=True)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-vcf", required=True, type=Path)
    parser.add_argument("--source-index", required=True, type=Path)
    parser.add_argument("--benchmark-bed", required=True, type=Path)
    parser.add_argument("--truth-catalog", required=True, type=Path)
    parser.add_argument("--evaluator-profile", required=True, type=Path)
    parser.add_argument("--profile-id", required=True)
    parser.add_argument("--output-vcf", required=True, type=Path)
    parser.add_argument("--audit-json", required=True, type=Path)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        prepare(args)
    except (
        OSError,
        subprocess.CalledProcessError,
        TruthPreparationError,
        SvMatchError,
    ) as exc:
        print(f"prepare_sv_truth: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
