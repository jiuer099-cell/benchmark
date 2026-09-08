#!/usr/bin/env python3
"""Freeze and gate PanGenie's native, family-excluded haplotype context.

The canonical scoring panel is intentionally *not* an index input for
PanGenie.  This rule accepts only a panel previously rebuilt from the frozen
HPRC graph/haplotypes after family exclusion, validates its provenance and
phased genotypes, and records the projection back to the fixed scoring
universe.  It never imputes, phases, or removes a sample column itself.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import sys
from collections import Counter
from pathlib import Path
from typing import TextIO

import yaml  # type: ignore[import-untyped]


FAMILY = {"HG002", "NA24385", "HG003", "NA24149", "HG004", "NA24143"}
CONTRACT = "pgbench_pangenie_native_context_v1"
METHOD = "family_excluded_hprc_graph_haplotypes"


class PrivatePanelError(ValueError):
    pass


def open_text(path: Path, mode: str) -> TextIO:
    if path.suffix.casefold() == ".gz":
        return gzip.open(path, mode, encoding="utf-8")  # type: ignore[return-value]
    return path.open(mode, encoding="utf-8")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def parse_info(raw: str) -> dict[str, str | bool]:
    result: dict[str, str | bool] = {}
    if raw in {"", "."}:
        return result
    for item in raw.split(";"):
        if not item:
            continue
        key, separator, value = item.partition("=")
        result[key] = value if separator else True
    return result


def format_info(values: dict[str, str | bool]) -> str:
    return ";".join(
        key if value is True else f"{key}={value}"
        for key, value in values.items()
    ) or "."


def load_provenance(
    path: Path,
    *,
    canonical_source: Path,
    source_graph: Path,
    source_haplotype_manifest: Path,
    excluded_samples: set[str],
    cohort_id: str,
) -> dict[str, object]:
    try:
        value = yaml.safe_load(path.read_text(encoding="utf-8"))
    except (OSError, yaml.YAMLError) as exc:
        raise PrivatePanelError(f"cannot read PanGenie context provenance: {exc}") from exc
    if not isinstance(value, dict) or value.get("contract") != CONTRACT:
        raise PrivatePanelError("unsupported PanGenie native-context provenance contract")
    required = {
        "source_cohort_id",
        "canonical_population_source_sha256",
        "source_graph_sha256",
        "source_haplotype_manifest_sha256",
        "reference_build",
        "generation_method",
        "family_exclusion_applied_before_native_panel_generation",
        "excluded_samples",
    }
    missing = sorted(required - set(value))
    if missing:
        raise PrivatePanelError("native-context provenance missing: " + ", ".join(missing))
    checks = {
        "source_cohort_id": cohort_id,
        "canonical_population_source_sha256": sha256(canonical_source),
        "source_graph_sha256": sha256(source_graph),
        "source_haplotype_manifest_sha256": sha256(source_haplotype_manifest),
        "reference_build": "grch38",
        "generation_method": METHOD,
        "family_exclusion_applied_before_native_panel_generation": True,
    }
    for key, expected in checks.items():
        if value.get(key) != expected:
            raise PrivatePanelError(f"native-context provenance {key} does not match frozen input")
    declared = value.get("excluded_samples")
    if not isinstance(declared, list) or set(declared) != excluded_samples:
        raise PrivatePanelError("native-context provenance exclusion set differs from run")
    return value


def candidate_keys(path: Path) -> dict[tuple[str, str, str, str], str]:
    result: dict[tuple[str, str, str, str], str] = {}
    with open_text(path, "rt") as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8 or not fields[2].startswith("CAND_"):
                raise PrivatePanelError(f"canonical scoring panel line {line_number} is invalid")
            key = tuple(fields[index] for index in (0, 1, 3, 4))
            if key in result:
                raise PrivatePanelError("canonical scoring panel has duplicate allele key")
            result[key] = fields[2]
    if not result:
        raise PrivatePanelError("canonical scoring panel has no candidates")
    return result


def prepare(
    *,
    source: Path,
    scoring_panel: Path,
    output: Path,
    projection: Path,
    gate: Path,
    excluded_samples: set[str],
) -> dict[str, int]:
    candidates = candidate_keys(scoring_panel)
    native_matches: dict[str, list[tuple[str, int]]] = {key: [] for key in candidates.values()}
    statistics: Counter[str] = Counter(
        {
            "N_GT_total": 0,
            "N_GT_missing": 0,
            "N_GT_unphased": 0,
            "N_GT_phased": 0,
            "native_record_count": 0,
            "index_addressable": 0,
            "index_unaddressable": 0,
            "ambiguous_projection": 0,
        }
    )
    samples: list[str] | None = None
    output.parent.mkdir(parents=True, exist_ok=True)
    with open_text(source, "rt") as reader, output.open("w", encoding="utf-8") as writer:
        for line_number, line in enumerate(reader, 1):
            if line.startswith("##"):
                writer.write(line)
                continue
            fields = line.rstrip("\n").split("\t")
            if line.startswith("#CHROM"):
                samples = fields[9:]
                if not samples:
                    raise PrivatePanelError("private phased panel has no sample columns")
                leaked = sorted(excluded_samples.intersection(samples))
                if leaked:
                    raise PrivatePanelError("private phased panel contains family samples: " + ", ".join(leaked))
                writer.write(line)
                continue
            if not line.strip():
                continue
            if samples is None or len(fields) != 9 + len(samples):
                raise PrivatePanelError(f"private panel line {line_number} has invalid sample columns")
            format_fields = fields[8].split(":")
            if "GT" not in format_fields:
                raise PrivatePanelError(f"private panel line {line_number} has no GT")
            gt_index = format_fields.index("GT")
            alt_count = len(fields[4].split(","))
            ac = [0] * alt_count
            an = 0
            for sample_value in fields[9:]:
                values = sample_value.split(":")
                genotype = values[gt_index] if gt_index < len(values) else ""
                statistics["N_GT_total"] += 1
                if "." in genotype:
                    statistics["N_GT_missing"] += 1
                    raise PrivatePanelError(f"private panel line {line_number} has missing GT")
                if "|" not in genotype or "/" in genotype:
                    statistics["N_GT_unphased"] += 1
                    raise PrivatePanelError(f"private panel line {line_number} has unphased GT")
                alleles = genotype.split("|")
                if len(alleles) != 2 or any(not allele.isdigit() for allele in alleles):
                    raise PrivatePanelError(f"private panel line {line_number} has invalid diploid GT")
                statistics["N_GT_phased"] += 1
                for allele in alleles:
                    value = int(allele)
                    if value > alt_count:
                        raise PrivatePanelError(f"private panel line {line_number} has GT outside ALT range")
                    an += 1
                    if value:
                        ac[value - 1] += 1
            info = parse_info(fields[7])
            info["AC"] = ",".join(map(str, ac))
            info["AN"] = str(an)
            info["AF"] = ",".join(f"{count / an:.12g}" for count in ac)
            fields[7] = format_info(info)
            writer.write("\t".join(fields) + "\n")
            for allele_index, alt in enumerate(fields[4].split(","), 1):
                candidate = candidates.get((fields[0], fields[1], fields[3], alt))
                if candidate is not None:
                    native_matches[candidate].append((fields[2], allele_index))
            statistics["native_record_count"] += 1
    if samples is None:
        raise PrivatePanelError("private phased panel has no #CHROM header")
    statistics["sample_count"] = len(samples)
    projection.parent.mkdir(parents=True, exist_ok=True)
    with projection.open("w", encoding="utf-8") as handle:
        handle.write("candidate_id\tnative_record_id\tnative_allele_index\tstatus\n")
        for candidate in sorted(native_matches):
            matches = native_matches[candidate]
            status = "index_addressable" if len(matches) == 1 else (
                "index_unaddressable" if not matches else "ambiguous_projection"
            )
            if len(matches) == 1:
                record, allele_index = matches[0]
                handle.write(f"{candidate}\t{record}\t{allele_index}\t{status}\n")
            else:
                handle.write(f"{candidate}\t.\t.\t{status}\n")
            statistics[status] += 1
    statistics["canonical_total"] = len(candidates)
    gate.parent.mkdir(parents=True, exist_ok=True)
    gate.write_text(json.dumps(dict(statistics), indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return dict(statistics)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-phased-panel", required=True, type=Path)
    parser.add_argument("--canonical-population-source", required=True, type=Path)
    parser.add_argument("--source-graph", required=True, type=Path)
    parser.add_argument("--source-haplotype-manifest", required=True, type=Path)
    parser.add_argument("--provenance", required=True, type=Path)
    parser.add_argument("--canonical-scoring-panel", required=True, type=Path)
    parser.add_argument("--source-cohort-id", required=True)
    parser.add_argument("--exclude-sample", action="append", default=[])
    parser.add_argument("--output-panel", required=True, type=Path)
    parser.add_argument("--output-projection", required=True, type=Path)
    parser.add_argument("--output-gate", required=True, type=Path)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        excluded = set(args.exclude_sample)
        if excluded != FAMILY:
            raise PrivatePanelError("private-panel gate requires the complete HG002 family exclusion")
        for path in (
            args.source_phased_panel, args.canonical_population_source,
            args.source_graph, args.source_haplotype_manifest, args.provenance,
            args.canonical_scoring_panel,
        ):
            if not path.is_file():
                raise PrivatePanelError(f"required private-panel input is missing: {path}")
        provenance = load_provenance(
            args.provenance, canonical_source=args.canonical_population_source,
            source_graph=args.source_graph,
            source_haplotype_manifest=args.source_haplotype_manifest,
            excluded_samples=excluded, cohort_id=args.source_cohort_id,
        )
        statistics = prepare(
            source=args.source_phased_panel, scoring_panel=args.canonical_scoring_panel,
            output=args.output_panel, projection=args.output_projection,
            gate=args.output_gate, excluded_samples=excluded,
        )
        payload = json.loads(args.output_gate.read_text(encoding="utf-8"))
        payload.update({
            "contract": "pgbench_pangenie_private_panel_gate_v1",
            "source_provenance_sha256": sha256(args.provenance),
            "source_graph_hash": sha256(args.source_graph),
            "source_haplotype_manifest_hash": sha256(args.source_haplotype_manifest),
            "family_exclusion_manifest_hash": sha256(args.provenance),
            "pangenie_private_panel_hash": sha256(args.output_panel),
            "canonical_allele_projection_hash": sha256(args.output_projection),
            "private_panel_source_sha256": sha256(args.source_phased_panel),
            "source_cohort_id": args.source_cohort_id,
            "family_exclusion_applied_before_native_panel_generation": provenance["family_exclusion_applied_before_native_panel_generation"],
            "missing_rate": 0.0,
            "unphased_rate": 0.0,
        })
        args.output_gate.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        print(json.dumps(statistics, sort_keys=True))
    except PrivatePanelError as exc:
        print(f"prepare_pangenie_private_panel: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
