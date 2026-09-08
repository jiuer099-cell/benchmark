#!/usr/bin/env python3
"""Build a truth-blinded challenge panel with frozen tolerant SV matching."""

from __future__ import annotations

import argparse
import bisect
import csv
import gzip
import hashlib
import json
import sys
from collections import defaultdict
from pathlib import Path
from typing import TextIO, cast

from build_pangenome_manifest import format_info, parse_info
from sv_matching import (
    SvRecord,
    SvMatchError,
    genotypes_equal,
    load_evaluator_profile,
    load_vcf,
    one_to_one_match,
    record_from_fields,
)


class ChallengePanelError(ValueError):
    """Raised when a blinded panel cannot be built safely."""


def _open_text(path: Path, mode: str) -> TextIO:
    if path.suffix == ".gz":
        return cast(TextIO, gzip.open(path, mode, encoding="utf-8"))
    return cast(TextIO, path.open(mode, encoding="utf-8"))


AUTOSOMES = {f"chr{number}" for number in range(1, 23)}


def _load_benchmark_regions(
    path: Path,
) -> dict[str, tuple[list[int], list[int]]]:
    raw: dict[str, list[tuple[int, int]]] = defaultdict(list)
    with path.open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                raise ChallengePanelError(
                    f"malformed BED record at {path}:{line_number}"
                )
            try:
                start, end = int(fields[1]), int(fields[2])
            except ValueError as exc:
                raise ChallengePanelError(
                    f"invalid BED coordinates at {path}:{line_number}"
                ) from exc
            if start < 0 or end <= start:
                raise ChallengePanelError(
                    f"invalid BED interval at {path}:{line_number}"
                )
            raw[fields[0]].append((start, end))
    if not raw:
        raise ChallengePanelError(f"benchmark BED contains no regions: {path}")
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


def _fully_contained(
    regions: dict[str, tuple[list[int], list[int]]],
    *,
    contig: str,
    start: int,
    end: int,
) -> bool:
    indexed = regions.get(contig)
    if indexed is None:
        return False
    starts, ends = indexed
    interval_index = bisect.bisect_right(starts, start) - 1
    return interval_index >= 0 and ends[interval_index] >= end


def _candidate_exclusion_reason(
    record: SvRecord,
    profile: dict[str, object],
    regions: dict[str, tuple[list[int], list[int]]] | None,
) -> str | None:
    """Return one exclusive reason why a panel record is outside the universe."""

    universe = profile["universe"]
    if record.chrom not in AUTOSOMES:
        return "non_autosome"
    if "," in record.alt:
        return "multiallelic"
    if record.svtype not in set(universe["allowed_svtypes"]):
        return "svtype"
    if (
        record.alt.startswith("<")
        or record.alt.endswith(">")
        or "[" in record.alt
        or "]" in record.alt
        or (record.svtype == "DEL" and len(record.ref) <= len(record.alt))
        or (record.svtype == "INS" and len(record.alt) <= len(record.ref))
    ):
        return "unresolved_sequence"
    size = abs(record.svlen)
    if not (
        int(universe["minimum_sv_size"])
        <= size
        <= int(universe["maximum_sv_size"])
    ):
        return "size"
    if record.filter_status not in {"PASS", "."}:
        return "filter"
    if regions is not None and not _fully_contained(
        regions,
        contig=record.chrom,
        start=record.pos - 1,
        end=max(record.pos, record.end),
    ):
        return "region"
    return None


def build_challenge_panel(
    *,
    panel_vcf: Path,
    truth_vcf: Path,
    output_vcf: Path,
    hidden_ledger: Path,
    audit_json: Path,
    seed: str,
    evaluator_profile: Path | None = None,
    benchmark_bed: Path | None = None,
    scope_ledger: Path | None = None,
    challenging_vcf: Path | None = None,
) -> dict[str, int | str]:
    profile_path = evaluator_profile or (
        Path(__file__).resolve().parents[2] / "config" / "evaluator_profile.yaml"
    )
    profile = load_evaluator_profile(profile_path)
    truth = load_vcf(truth_vcf, prefix="TRUTH")
    output_vcf.parent.mkdir(parents=True, exist_ok=True)
    hidden_ledger.parent.mkdir(parents=True, exist_ok=True)
    audit_json.parent.mkdir(parents=True, exist_ok=True)

    panel_headers: list[str] = []
    source_entries: list[tuple[list[str], SvRecord, str | None]] = []
    regions = (
        _load_benchmark_regions(benchmark_bed)
        if benchmark_bed is not None
        else None
    )
    exclusions = {
        "excluded_non_autosome_count": 0,
        "excluded_multiallelic_count": 0,
        "excluded_svtype_count": 0,
        "excluded_unresolved_sequence_count": 0,
        "excluded_size_count": 0,
        "excluded_filter_count": 0,
        "excluded_region_count": 0,
    }
    source_candidate_count = 0
    with _open_text(panel_vcf, "rt") as source:
        for line in source:
            if line.startswith("#"):
                panel_headers.append(line)
                continue
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            source_candidate_count += 1
            record = record_from_fields(
                fields,
                prefix="PANEL",
                index=source_candidate_count - 1,
            )
            reason = _candidate_exclusion_reason(record, profile, regions)
            declared_complexity = parse_info(fields[7]).get("GRAPH_COMPLEXITY")
            if reason is None and declared_complexity in {
                "multiallelic",
                "nested_or_overlapping",
                "repeat_ambiguous",
            }:
                reason = str(declared_complexity)
            source_entries.append((fields, record, reason))

    overlapping: set[int] = set()
    active: list[int] = []
    for index in sorted(
        range(len(source_entries)),
        key=lambda item: (
            source_entries[item][1].chrom,
            source_entries[item][1].pos,
            source_entries[item][1].end,
        ),
    ):
        current = source_entries[index][1]
        active = [
            prior
            for prior in active
            if source_entries[prior][1].chrom == current.chrom
            and source_entries[prior][1].end >= current.pos
        ]
        for prior in active:
            overlapping.add(prior)
            overlapping.add(index)
        active.append(index)

    normalized_entries: list[tuple[list[str], SvRecord, str | None]] = []
    for index, (fields, record, reason) in enumerate(source_entries):
        if reason is None and index in overlapping:
            reason = "nested_or_overlapping"
        if reason is not None:
            exclusions.setdefault(f"excluded_{reason}_count", 0)
            exclusions[f"excluded_{reason}_count"] += 1
        normalized_entries.append((fields, record, reason))
    panel_fields = [fields for fields, _, reason in normalized_entries if reason is None]
    panel_records = [record for _, record, reason in normalized_entries if reason is None]
    if not panel_records:
        raise ChallengePanelError(
            "panel contains no candidates in the frozen universe"
        )
    assignments = one_to_one_match(panel_records, truth, profile)

    if scope_ledger is not None:
        scope_ledger.parent.mkdir(parents=True, exist_ok=True)
        with scope_ledger.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=[
                    "source_index",
                    "pangenome_allele_id",
                    "chrom",
                    "pos",
                    "end",
                    "svtype",
                    "svlen",
                    "complexity_class",
                    "main_track_eligible",
                    "exclusion_reason",
                ],
                delimiter="\t",
                lineterminator="\n",
            )
            writer.writeheader()
            for index, (fields, record, reason) in enumerate(normalized_entries, start=1):
                allele_id = parse_info(fields[7]).get("PANGENOME_ALLELE_ID", "")
                complexity = (
                    "nested_or_overlapping"
                    if reason == "nested_or_overlapping"
                    else ("multiallelic" if reason == "multiallelic" else "simple_biallelic")
                )
                writer.writerow(
                    {
                        "source_index": index,
                        "pangenome_allele_id": allele_id,
                        "chrom": record.chrom,
                        "pos": record.pos,
                        "end": record.end,
                        "svtype": record.svtype,
                        "svlen": record.svlen,
                        "complexity_class": complexity,
                        "main_track_eligible": int(reason is None),
                        "exclusion_reason": reason or "",
                    }
                )

    if challenging_vcf is not None:
        challenging_vcf.parent.mkdir(parents=True, exist_ok=True)
        with _open_text(challenging_vcf, "wt") as output:
            for line in panel_headers:
                if line.startswith("##"):
                    if "TRUTH" not in line.upper():
                        output.write(line)
                elif line.startswith("#CHROM"):
                    output.write("\t".join(line.rstrip("\n").split("\t")[:8] + ["FORMAT", "HG002"]) + "\n")
            for fields, _, reason in normalized_entries:
                if reason is None:
                    continue
                info = parse_info(fields[7])
                for key in list(info):
                    if key.upper().startswith("TRUTH"):
                        del info[key]
                output.write("\t".join(fields[:7] + [format_info(info), "GT", "./."]) + "\n")

    counts = {
        "candidate_count": 0,
        "truth_scorable_count": 0,
        "truth_unscorable_count": 0,
        "truth_positive_count": 0,
        "truth_negative_count": 0,
        "truth_event_match_count": len(assignments),
        "genotype_exact_count": 0,
        "source_candidate_count": source_candidate_count,
        **exclusions,
    }
    seen_candidate_ids: set[str] = set()
    saw_column_header = False
    record_index = 0

    with (
        _open_text(output_vcf, "wt") as output,
        hidden_ledger.open("w", encoding="utf-8", newline="") as ledger_handle,
    ):
        ledger = csv.DictWriter(
            ledger_handle,
            fieldnames=[
                "candidate_id",
                "pangenome_allele_id",
                "truth_gt",
                "truth_label",
                "truth_scorable",
                "truth_event_id",
                "match_method",
                "match_score",
                "start_distance",
                "end_distance",
                "size_similarity",
                "sequence_similarity",
                "gt_exact",
            ],
            delimiter="\t",
            lineterminator="\n",
        )
        ledger.writeheader()

        for line in panel_headers:
            if line.startswith("##"):
                if "TRUTH" in line.upper():
                    continue
                output.write(line)
                continue
            if line.startswith("#CHROM"):
                header_fields = line.rstrip("\n").split("\t")
                if len(header_fields) < 8:
                    raise ChallengePanelError("panel VCF has an invalid #CHROM header")
                # PanGenie needs panel haplotypes through pangenome_panel, but
                # those samples must not leak into the shared blinded universe.
                output.write(
                    "\t".join(header_fields[:8] + ["FORMAT", "HG002"]) + "\n"
                )
                saw_column_header = True
                continue
        if not saw_column_header:
            raise ChallengePanelError("panel VCF is missing #CHROM header")

        for original_fields in panel_fields:
            fields = original_fields[:8]
            match = assignments.get(record_index)
            record_index += 1
            info = parse_info(fields[7])
            allele_id = info.get("PANGENOME_ALLELE_ID")
            if not isinstance(allele_id, str) or not allele_id:
                raise ChallengePanelError("panel record is missing PANGENOME_ALLELE_ID")
            # The population allele ID is the public, frozen canonical ID.
            # It is truth-independent and therefore does not need blinding.
            candidate_id = allele_id
            if candidate_id in seen_candidate_ids:
                raise ChallengePanelError(f"candidate ID collision: {candidate_id}")
            seen_candidate_ids.add(candidate_id)

            matched_truth = truth[match.truth_index] if match is not None else None
            truth_gt = matched_truth.gt if matched_truth is not None else "UNSCORABLE"
            positive = bool(
                matched_truth is not None
                and truth_gt not in {"0/0", "0|0", "./.", ".|."}
            )
            truth_scorable = matched_truth is not None and truth_gt not in {
                "./.",
                ".|.",
            }
            truth_label = (
                "positive"
                if truth_scorable and positive
                else ("negative" if truth_scorable else "unscorable")
            )
            match_method = (
                "tolerant_one_to_one" if matched_truth is not None else "no_match"
            )
            gt_exact = (
                genotypes_equal(
                    panel_records[record_index - 1].gt,
                    truth_gt,
                    require_phase=bool(
                        profile["semantics"]["genotype"]["require_phase"]
                    ),
                )
                if matched_truth is not None
                else False
            )
            counts["candidate_count"] += 1
            counts[
                "truth_scorable_count"
                if truth_scorable
                else "truth_unscorable_count"
            ] += 1
            if truth_scorable:
                counts[
                    "truth_positive_count" if positive else "truth_negative_count"
                ] += 1
            counts["genotype_exact_count"] += int(gt_exact)

            for leaked_key in list(info):
                if leaked_key.upper().startswith("TRUTH"):
                    del info[leaked_key]
            fields[2] = candidate_id
            fields[7] = format_info(info)
            output.write("\t".join(fields + ["GT", "./."]) + "\n")
            ledger.writerow(
                {
                    "candidate_id": candidate_id,
                    "pangenome_allele_id": allele_id,
                    "truth_gt": truth_gt,
                    "truth_label": truth_label,
                    "truth_scorable": int(truth_scorable),
                    "truth_event_id": (
                        matched_truth.record_id if matched_truth is not None else ""
                    ),
                    "match_method": match_method,
                    "match_score": f"{match.score:.8f}" if match is not None else "",
                    "start_distance": (
                        str(match.start_distance) if match is not None else ""
                    ),
                    "end_distance": (
                        str(match.end_distance) if match is not None else ""
                    ),
                    "size_similarity": (
                        f"{match.size_similarity:.8f}" if match is not None else ""
                    ),
                    "sequence_similarity": (
                        ""
                        if match is None or match.sequence_similarity is None
                        else f"{match.sequence_similarity:.8f}"
                    ),
                    "gt_exact": int(gt_exact),
                }
            )

    if counts["candidate_count"] == 0:
        raise ChallengePanelError("panel contains no candidates")
    if counts["truth_positive_count"] == 0:
        raise ChallengePanelError("challenge panel has no truth-positive candidate")

    audit: dict[str, int | str] = {
        **counts,
        "matcher_profile": str(profile["profile"]["id"]),
        "evaluator_profile_sha256": str(profile["_sha256"]),
        "seed_sha256": hashlib.sha256(seed.encode("utf-8")).hexdigest(),
        "truth_labels_exposed_to_tool": 0,
        "truth_assignments_one_to_one": 1,
        "candidate_ledger_matches_emitted_vcf": 1,
        "excluded_counts_are_mutually_exclusive": 1,
        "scope_ledger_complete": int(scope_ledger is not None),
        "challenging_track_emitted": int(challenging_vcf is not None),
    }
    temporary = audit_json.with_name(f".{audit_json.name}.tmp")
    temporary.write_text(
        json.dumps(audit, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    temporary.replace(audit_json)
    return audit


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build a blinded synthetic genotype challenge panel."
    )
    parser.add_argument("--panel-vcf", required=True, type=Path)
    parser.add_argument("--truth-vcf", required=True, type=Path)
    parser.add_argument("--benchmark-bed", required=True, type=Path)
    parser.add_argument("--output-vcf", required=True, type=Path)
    parser.add_argument("--hidden-ledger", required=True, type=Path)
    parser.add_argument("--audit-json", required=True, type=Path)
    parser.add_argument("--scope-ledger", required=True, type=Path)
    parser.add_argument("--challenging-vcf", required=True, type=Path)
    parser.add_argument("--seed", required=True)
    parser.add_argument(
        "--evaluator-profile",
        type=Path,
        default=Path("config/evaluator_profile.yaml"),
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        build_challenge_panel(
            panel_vcf=args.panel_vcf,
            truth_vcf=args.truth_vcf,
            benchmark_bed=args.benchmark_bed,
            output_vcf=args.output_vcf,
            hidden_ledger=args.hidden_ledger,
            audit_json=args.audit_json,
            seed=args.seed,
            evaluator_profile=args.evaluator_profile,
            scope_ledger=args.scope_ledger,
            challenging_vcf=args.challenging_vcf,
        )
    except (OSError, ChallengePanelError, SvMatchError) as exc:
        print(f"build_challenge_panel: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
