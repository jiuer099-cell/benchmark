#!/usr/bin/env python3
"""Freeze and gate PanGenie's native, family-excluded haplotype context.

The canonical scoring panel is intentionally *not* an index input for
PanGenie.  This rule accepts only a panel previously rebuilt from the frozen
HPRC graph/haplotypes after family exclusion, validates its provenance and
phased genotypes, and records the projection back to the fixed scoring
universe.  It never imputes, guesses phase, or removes a sample column
itself.  The sole permitted normalization is ``a/a -> a|a``: a diploid
homozygous slash genotype is phase-equivalent and therefore carries no
unknown haplotype assignment.
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
CONTRACT = "pgbench_pangenie_native_context_v2"
METHOD = "deterministic_family_exclusion_from_official_pangenie_pgin"
MAX_MISSING_HAPLOTYPE_FRACTION = 0.20
PHASE_POLICY = {
    "homozygous_unphased": {
        "phase_ambiguous": False,
        "allowed_to_canonicalize": True,
    },
    "heterozygous_unphased": {
        "phase_ambiguous": True,
        "allowed_to_guess": False,
        "fail_if_unresolved": True,
    },
}


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
    source_graph: Path,
    source_haplotype_manifest: Path,
    source_phased_panel: Path,
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
        "source_graph_sha256",
        "source_haplotype_manifest_sha256",
        "official_source_pgin_sha256",
        "reference_build",
    }
    missing = sorted(required - set(value))
    if missing:
        raise PrivatePanelError("native-context provenance missing: " + ", ".join(missing))
    checks = {
        "source_cohort_id": cohort_id,
        "source_graph_sha256": sha256(source_graph),
        "source_haplotype_manifest_sha256": sha256(source_haplotype_manifest),
        "official_source_pgin_sha256": sha256(source_phased_panel),
        "reference_build": "grch38",
    }
    for key, expected in checks.items():
        if value.get(key) != expected:
            raise PrivatePanelError(f"native-context provenance {key} does not match frozen input")
    return value


def validate_source_identity(
    path: Path, *, frozen_gbz: Path, frozen_sample_manifest: Path, cohort_id: str
) -> dict[str, object]:
    try:
        value = yaml.safe_load(path.read_text(encoding="utf-8"))
    except (OSError, yaml.YAMLError) as exc:
        raise PrivatePanelError(f"cannot read exact-source identity manifest: {exc}") from exc
    if not isinstance(value, dict) or value.get("contract") != "pgbench_hprc_graph_identity_v1":
        raise PrivatePanelError("unsupported exact-source identity manifest")
    required = {
        "status", "source_cohort_id", "current_gbz_sha256",
        "source_sample_manifest_sha256", "source_release", "graph_family",
        "graph_construction_version", "reference_build", "official_manifest_url",
        "official_manifest_sha256", "official_pgin_url", "official_pgin_sha256",
    }
    missing = sorted(required - set(value))
    if missing:
        raise PrivatePanelError("exact-source identity is incomplete: " + ", ".join(missing))
    if value.get("status") != "verified":
        raise PrivatePanelError("exact-source identity must be verified before PanGenie runs")
    checks = {
        "source_cohort_id": cohort_id,
        "current_gbz_sha256": sha256(frozen_gbz),
        "source_sample_manifest_sha256": sha256(frozen_sample_manifest),
        "reference_build": "grch38",
    }
    for key, expected in checks.items():
        if value.get(key) != expected:
            raise PrivatePanelError(f"exact-source identity {key} does not match frozen resource")
    for key in (
        "source_release", "graph_family", "graph_construction_version",
        "official_manifest_url", "official_manifest_sha256", "official_pgin_url",
        "official_pgin_sha256",
    ):
        if not isinstance(value.get(key), str) or not value[key]:
            raise PrivatePanelError(f"exact-source identity {key} must be non-empty")
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


def slash_gt_category(genotype: str) -> str | None:
    """Classify slash GTs before deciding whether their phase is recoverable."""

    if "/" not in genotype:
        return None
    alleles = genotype.split("/")
    if len(alleles) != 2 or any(
        allele != "." and not allele.isdigit() for allele in alleles
    ):
        return "malformed"
    if any(allele == "." for allele in alleles):
        return "partial_missing"
    return "phase_equivalent_homozygous" if alleles[0] == alleles[1] else "phase_ambiguous_heterozygous"


def canonical_diploid_gt(genotype: str, *, line_number: int) -> tuple[list[str], str]:
    """Return a deterministic diploid GT and its phase-audit category.

    A slash is not automatically an unresolved phase error: ``0/0`` and
    ``1/1`` are exactly equivalent to their phased forms.  In contrast,
    heterozygous and partially-missing slash calls do not identify the two
    haplotypes, so this fail-closed gate must reject them.
    """

    if genotype == "./.":
        return [".", "."], "fully_missing"
    if "/" in genotype:
        alleles = genotype.split("/")
        category = slash_gt_category(genotype)
        if category == "malformed":
            raise PrivatePanelError(
                f"private panel line {line_number} has malformed slash GT"
            )
        if category == "partial_missing":
            raise PrivatePanelError(
                f"private panel line {line_number} has partial-missing slash GT"
            )
        if category == "phase_ambiguous_heterozygous":
            raise PrivatePanelError(
                f"private panel line {line_number} has phase-ambiguous heterozygous GT"
            )
        return alleles, "phase_equivalent_homozygous_slash"
    if "|" not in genotype:
        raise PrivatePanelError(f"private panel line {line_number} has invalid diploid GT")
    alleles = genotype.split("|")
    if len(alleles) != 2 or any(
        allele != "." and not allele.isdigit() for allele in alleles
    ):
        raise PrivatePanelError(f"private panel line {line_number} has invalid diploid GT")
    return alleles, "phased"


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
            "N_GT_phase_equivalent_homozygous_slash": 0,
            "N_GT_phase_ambiguous_heterozygous_slash": 0,
            "N_GT_partial_missing_slash": 0,
            "N_GT_malformed_slash": 0,
            "N_GT_phased": 0,
            "native_record_count": 0,
            "index_addressable": 0,
            "index_unaddressable": 0,
            "ambiguous_projection": 0,
            "source_sample_count": 0,
            "retained_sample_count": 0,
            "family_samples_removed": 0,
            "family_samples_not_present": 0,
            "zero_support_records_removed": 0,
            "zero_support_alleles_removed": 0,
            "N_haplotypes_total": 0,
            "N_haplotypes_missing": 0,
            "max_missing_haplotype_fraction": 0.0,
        }
    )
    samples: list[str] | None = None
    retained_indices: list[int] | None = None
    source_sample_count = 0
    output.parent.mkdir(parents=True, exist_ok=True)
    with open_text(source, "rt") as reader, output.open("w", encoding="utf-8") as writer:
        for line_number, line in enumerate(reader, 1):
            if line.startswith("##"):
                writer.write(line)
                continue
            fields = line.rstrip("\n").split("\t")
            if line.startswith("#CHROM"):
                source_samples = fields[9:]
                retained_indices = [
                    index for index, sample in enumerate(source_samples)
                    if sample not in excluded_samples
                ]
                samples = [source_samples[index] for index in retained_indices]
                if not samples:
                    raise PrivatePanelError("private phased panel has no sample columns")
                source_sample_count = len(source_samples)
                statistics["source_sample_count"] = source_sample_count
                statistics["retained_sample_count"] = len(samples)
                statistics["family_samples_removed"] = source_sample_count - len(samples)
                statistics["family_samples_not_present"] = len(excluded_samples - set(source_samples))
                writer.write("\t".join(fields[:9] + samples) + "\n")
                continue
            if not line.strip():
                continue
            if (
                samples is None
                or retained_indices is None
                or len(fields) != 9 + source_sample_count
            ):
                raise PrivatePanelError(f"private panel line {line_number} has invalid sample columns")
            format_fields = fields[8].split(":")
            if "GT" not in format_fields:
                raise PrivatePanelError(f"private panel line {line_number} has no GT")
            gt_index = format_fields.index("GT")
            alt_count = len(fields[4].split(","))
            ac = [0] * alt_count
            an = 0
            retained_genotypes: list[tuple[list[str], list[str]]] = []
            for source_index in retained_indices:
                sample_value = fields[9 + source_index]
                values = sample_value.split(":")
                genotype = values[gt_index] if gt_index < len(values) else ""
                statistics["N_GT_total"] += 1
                try:
                    alleles, phase_category = canonical_diploid_gt(
                        genotype, line_number=line_number
                    )
                except PrivatePanelError as exc:
                    slash_category = slash_gt_category(genotype)
                    if slash_category == "partial_missing":
                        statistics["N_GT_partial_missing_slash"] += 1
                    elif slash_category == "phase_ambiguous_heterozygous":
                        statistics["N_GT_phase_ambiguous_heterozygous_slash"] += 1
                    elif slash_category == "malformed":
                        statistics["N_GT_malformed_slash"] += 1
                    raise exc
                if phase_category == "phase_equivalent_homozygous_slash":
                    statistics["N_GT_phase_equivalent_homozygous_slash"] += 1
                elif phase_category == "phased":
                    statistics["N_GT_phased"] += 1
                for allele in alleles:
                    statistics["N_haplotypes_total"] += 1
                    if allele == ".":
                        statistics["N_haplotypes_missing"] += 1
                        statistics["N_GT_missing"] += 1
                        continue
                    value = int(allele)
                    if value > alt_count:
                        raise PrivatePanelError(f"private panel line {line_number} has GT outside ALT range")
                    an += 1
                    if value:
                        ac[value - 1] += 1
                retained_genotypes.append((values, alleles))
            missing_fraction = sum(
                allele == "." for _values, alleles in retained_genotypes for allele in alleles
            ) / (2 * len(retained_genotypes))
            statistics["max_missing_haplotype_fraction"] = max(
                statistics["max_missing_haplotype_fraction"], missing_fraction
            )
            if missing_fraction > MAX_MISSING_HAPLOTYPE_FRACTION:
                raise PrivatePanelError(
                    f"private panel line {line_number} exceeds official "
                    "prepare-vcf-MC missing-haplotype threshold"
                )
            retained_alt_indexes = [index for index, count in enumerate(ac, start=1) if count]
            if not retained_alt_indexes:
                statistics["zero_support_records_removed"] += 1
                statistics["zero_support_alleles_removed"] += alt_count
                continue
            statistics["zero_support_alleles_removed"] += alt_count - len(retained_alt_indexes)
            old_to_new = {old: new for new, old in enumerate(retained_alt_indexes, start=1)}
            fields[4] = ",".join(
                fields[4].split(",")[index - 1] for index in retained_alt_indexes
            )
            retained_samples: list[str] = []
            for values, alleles in retained_genotypes:
                values[gt_index] = "./." if alleles == [".", "."] else "|".join(
                    "." if allele == "." else (
                        "0" if allele == "0" else str(old_to_new[int(allele)])
                    )
                    for allele in alleles
                )
                retained_samples.append(":".join(values))
            ac = [ac[index - 1] for index in retained_alt_indexes]
            info = parse_info(fields[7])
            info["AC"] = ",".join(map(str, ac))
            info["AN"] = str(an)
            info["AF"] = ",".join(f"{count / an:.12g}" for count in ac)
            fields[7] = format_info(info)
            writer.write("\t".join(fields[:9] + retained_samples) + "\n")
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
    parser.add_argument("--frozen-gbz", required=True, type=Path)
    parser.add_argument("--frozen-sample-manifest", required=True, type=Path)
    parser.add_argument("--source-identity", required=True, type=Path)
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
            args.canonical_scoring_panel, args.frozen_gbz,
            args.frozen_sample_manifest, args.source_identity,
        ):
            if not path.is_file():
                raise PrivatePanelError(f"required private-panel input is missing: {path}")
        provenance = load_provenance(
            args.provenance, source_graph=args.source_graph,
            source_haplotype_manifest=args.source_haplotype_manifest,
            source_phased_panel=args.source_phased_panel,
            cohort_id=args.source_cohort_id,
        )
        identity = validate_source_identity(
            args.source_identity, frozen_gbz=args.frozen_gbz,
            frozen_sample_manifest=args.frozen_sample_manifest,
            cohort_id=args.source_cohort_id,
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
            "source_identity_manifest_hash": sha256(args.source_identity),
            "current_gbz_hash": sha256(args.frozen_gbz),
            "source_sample_manifest_hash": sha256(args.frozen_sample_manifest),
            "source_release": identity["source_release"],
            "graph_family": identity["graph_family"],
            "graph_construction_version": identity["graph_construction_version"],
            "source_graph_hash": sha256(args.source_graph),
            "source_haplotype_manifest_hash": sha256(args.source_haplotype_manifest),
            "family_exclusion_manifest_hash": sha256(args.provenance),
            "pangenie_private_panel_hash": sha256(args.output_panel),
            "canonical_allele_projection_hash": sha256(args.output_projection),
            "private_panel_source_sha256": sha256(args.source_phased_panel),
            "source_cohort_id": args.source_cohort_id,
            "family_exclusion_applied_before_native_panel_generation": True,
            "family_exclusion_transform_version": METHOD,
            "excluded_samples": sorted(excluded),
            "missing_haplotype_policy": "official_prepare_vcf_mc_le_0.20",
            "phase_gate": PHASE_POLICY,
            "phase_equivalent_homozygous_normalization": "a/a_to_a|a_only",
            "max_missing_haplotype_fraction": statistics["max_missing_haplotype_fraction"],
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
