#!/usr/bin/env python3
"""Build a graph-excluded truth universe for novel pangenome discovery."""

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
from collections.abc import Iterable, Mapping
from dataclasses import dataclass
from pathlib import Path
from typing import Any, TextIO

import yaml  # type: ignore[import-untyped]

from build_pangenome_manifest import parse_info
from sv_matching import (
    SvMatch,
    SvRecord,
    compatibility,
    load_evaluator_profile,
    load_vcf,
    open_text,
    record_shape_in_universe,
    sha256_file,
)


WALK_SEGMENT = re.compile(r"([<>])([^<>]+)")
IUPAC_DNA = re.compile(r"^[ACGTRYSWKMBDHVNacgtryswkmbdhvn]+$")


class NovelTruthError(ValueError):
    """Raised when graph exclusion cannot be proven from frozen assets."""


@dataclass(frozen=True)
class GraphBubble:
    chrom: str
    start: int
    end: int
    reference_allele: int | None
    allele_lengths: tuple[int, ...]
    allele_walks: tuple[str, ...]
    line_number: int


class IndexedReference:
    """Small random-access reader for a registered uncompressed FASTA."""

    def __init__(self, fasta: Path, fai: Path) -> None:
        self.fasta = fasta
        self.entries: dict[str, tuple[int, int, int, int]] = {}
        for raw in fai.read_text(encoding="utf-8").splitlines():
            fields = raw.split("\t")
            if len(fields) >= 5:
                self.entries[fields[0]] = tuple(map(int, fields[1:5]))

    def sequence(self, chrom: str, start: int, end: int) -> str:
        """Return a 0-based half-open reference interval."""

        entry = self.entries.get(chrom)
        if entry is None:
            raise NovelTruthError(f"reference has no contig {chrom}")
        length, offset, line_bases, line_width = entry
        if start < 0 or end < start or end > length:
            raise NovelTruthError(
                f"reference interval is out of bounds: {chrom}:{start}-{end}"
            )
        pieces: list[bytes] = []
        cursor = start
        with self.fasta.open("rb") as handle:
            while cursor < end:
                line_offset = cursor % line_bases
                take = min(end - cursor, line_bases - line_offset)
                byte_offset = (
                    offset
                    + (cursor // line_bases) * line_width
                    + line_offset
                )
                handle.seek(byte_offset)
                pieces.append(handle.read(take))
                cursor += take
        value = b"".join(pieces).decode("ascii").upper()
        if len(value) != end - start or (
            value and not IUPAC_DNA.fullmatch(value)
        ):
            raise NovelTruthError(
                f"reference sequence could not be read at {chrom}:{start}-{end}"
            )
        return value


def _load_profile(path: Path) -> dict[str, Any]:
    try:
        value = yaml.safe_load(path.read_text(encoding="utf-8"))
    except (OSError, yaml.YAMLError) as exc:
        raise NovelTruthError(f"cannot load novel-truth profile {path}") from exc
    if not isinstance(value, dict) or value.get("schema_version") != 2:
        raise NovelTruthError("novel-truth profile schema_version must be 2")
    expected = {
        "source": {
            "format": "minigraph_call_bed_v0.17",
            "reference_sample_column": 0,
            "allele_lengths_info": "ALEN",
            "allele_walks_info": "AWALK",
            "graph_reference_path": "GRCh38",
            "allele_lengths_role": "audit_only_per_sample_call_report",
            "allele_sequence_source": "AWALK_resolved_from_frozen_GFA",
            "sequence_alphabet": "IUPAC_DNA",
            "reference_sequence_source": "frozen_GRCh38_FASTA_interval",
            "reference_sample_genotype_role": "audit_only",
            "graph_link_overlap_contract": "0M",
        },
        "exclusion": {
            "policy": "any_compatible_graph_allele",
            "size_source": (
                "reconstructed_AWALK_length_minus_GRCh38_interval_length"
            ),
            "require_pure_ins_del_after_normalization": True,
            "complex_allele_policy": "exclude_not_decompose",
            "require_deleted_sequence_reference_match": True,
            "fail_on_missing_graph_segment": True,
        },
    }
    for section, contract in expected.items():
        observed = value.get(section)
        if not isinstance(observed, Mapping):
            raise NovelTruthError(f"novel-truth profile is missing {section}")
        for key, expected_value in contract.items():
            if observed.get(key) != expected_value:
                raise NovelTruthError(
                    f"novel-truth profile {section}.{key} is not frozen"
                )
    value["_sha256"] = sha256_file(path)
    return value


def _reference_allele(raw_sample: str) -> int | None:
    gt = raw_sample.split(":", 1)[0]
    alleles = re.split(r"[/|]", gt)
    if alleles and all(allele == "." for allele in alleles):
        return None
    if not alleles or any(not allele.isdigit() for allele in alleles):
        raise NovelTruthError("graph call has no reference-sample genotype")
    unique = {int(allele) for allele in alleles}
    if len(unique) != 1:
        raise NovelTruthError("graph reference sample has a heterozygous call")
    return unique.pop()


def load_graph_bubbles(path: Path) -> list[GraphBubble]:
    bubbles: list[GraphBubble] = []
    with open_text(path) as handle:
        for line_number, raw in enumerate(handle, start=1):
            if not raw.strip() or raw.startswith("#"):
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) < 6:
                raise NovelTruthError(
                    f"malformed minigraph call at {path}:{line_number}"
                )
            try:
                start, end = int(fields[1]), int(fields[2])
            except ValueError as exc:
                raise NovelTruthError("graph call coordinates must be integers") from exc
            info = parse_info(fields[3])
            raw_lengths = info.get("ALEN")
            raw_walks = info.get("AWALK")
            if not isinstance(raw_lengths, str) or not isinstance(raw_walks, str):
                raise NovelTruthError("graph call is missing ALEN or AWALK")
            try:
                lengths = tuple(int(value) for value in raw_lengths.split(","))
            except ValueError as exc:
                raise NovelTruthError("graph call ALEN contains a non-integer") from exc
            walks = tuple(raw_walks.split(","))
            if len(lengths) != len(walks) or len(lengths) < 1:
                raise NovelTruthError("graph call ALEN/AWALK cardinality differs")
            reference_allele = _reference_allele(fields[5])
            if reference_allele is not None and reference_allele >= len(lengths):
                raise NovelTruthError("reference genotype exceeds graph allele count")
            bubbles.append(
                GraphBubble(
                    chrom=fields[0],
                    start=start,
                    end=end,
                    reference_allele=reference_allele,
                    allele_lengths=lengths,
                    allele_walks=walks,
                    line_number=line_number,
                )
            )
    if not bubbles:
        raise NovelTruthError("graph variation call file contains no bubbles")
    return bubbles


def _walk_parts(walk: str) -> tuple[tuple[str, str], ...]:
    if walk == "*":
        return ()
    parts = tuple(WALK_SEGMENT.findall(walk))
    if not parts or "".join(direction + name for direction, name in parts) != walk:
        raise NovelTruthError(f"invalid graph allele walk: {walk}")
    return parts


def load_graph_segments(gfa: Path, required: set[str]) -> dict[str, str]:
    segments: dict[str, str] = {}
    with open_text(gfa) as handle:
        for raw in handle:
            if raw.startswith("L\t"):
                fields = raw.rstrip("\n").split("\t")
                if len(fields) < 6 or fields[5] != "0M":
                    raise NovelTruthError(
                        "graph contains a link outside the frozen 0M overlap contract"
                    )
                continue
            if not raw.startswith("S\t"):
                continue
            fields = raw.rstrip("\n").split("\t", 3)
            if len(fields) >= 3 and fields[1] in required:
                if fields[1] in segments:
                    raise NovelTruthError(
                        f"graph contains duplicate required segment: {fields[1]}"
                    )
                sequence = fields[2].upper()
                if sequence == "*" or not IUPAC_DNA.fullmatch(sequence):
                    raise NovelTruthError(
                        f"required graph segment has no resolved sequence: {fields[1]}"
                    )
                segments[fields[1]] = sequence
    missing = sorted(required - set(segments))
    if missing:
        preview = ", ".join(missing[:5])
        raise NovelTruthError(f"graph is missing required segments: {preview}")
    return segments


def _reverse_complement(sequence: str) -> str:
    return sequence.translate(
        str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN")
    )[::-1]


def _walk_sequence(walk: str, segments: Mapping[str, str]) -> str:
    values = []
    for direction, name in _walk_parts(walk):
        sequence = segments[name]
        values.append(sequence if direction == ">" else _reverse_complement(sequence))
    return "".join(values)


def _eligible_length_delta(
    delta: int,
    evaluator_profile: Mapping[str, Any],
) -> bool:
    """Return whether an exact reconstructed length delta is in scope."""

    universe = evaluator_profile["universe"]
    return (
        delta != 0
        and abs(delta) >= int(universe["minimum_sv_size"])
        and abs(delta) <= int(universe["maximum_sv_size"])
        and ("INS" if delta > 0 else "DEL")
        in set(universe["allowed_svtypes"])
    )


def _normalized_alleles(
    *,
    chrom: str,
    start: int,
    end: int,
    reference_sequence: str,
    alternate_sequence: str,
    reference: IndexedReference,
) -> tuple[int, str, str]:
    if start > 0:
        anchor = reference.sequence(chrom, start - 1, start)
        pos = start
        ref = anchor + reference_sequence
        alt = anchor + alternate_sequence
    else:
        anchor = reference.sequence(chrom, end, end + 1)
        pos = 1
        ref = reference_sequence + anchor
        alt = alternate_sequence + anchor
    suffix = 0
    maximum_suffix = min(len(ref), len(alt)) - 1
    while suffix < maximum_suffix and ref[-1 - suffix] == alt[-1 - suffix]:
        suffix += 1
    remaining_ref = len(ref) - suffix
    remaining_alt = len(alt) - suffix
    prefix = 0
    maximum_prefix = min(remaining_ref, remaining_alt) - 1
    while prefix < maximum_prefix and ref[prefix] == alt[prefix]:
        prefix += 1
    ref_end = len(ref) - suffix if suffix else len(ref)
    alt_end = len(alt) - suffix if suffix else len(alt)
    return pos + prefix, ref[prefix:ref_end], alt[prefix:alt_end]


def graph_alleles(
    *,
    bubbles: list[GraphBubble],
    segments: Mapping[str, str],
    reference: IndexedReference,
    evaluator_profile: dict[str, Any],
) -> tuple[list[SvRecord], dict[str, int]]:
    records: list[SvRecord] = []
    seen: set[tuple[Any, ...]] = set()
    audit = {
        "bubbles": len(bubbles),
        "alleles": sum(len(bubble.allele_walks) for bubble in bubbles),
        "missing_reference_genotype_bubbles": 0,
        "bubbles_with_exact_reference_path": 0,
        "multiple_exact_reference_path_bubbles": 0,
        "declared_reference_length_mismatches": 0,
        "declared_reference_sequence_mismatches": 0,
        "declared_allele_length_mismatches": 0,
        "reference_interval_bubbles": 0,
        "reference_identical_alleles": 0,
        "out_of_universe_length_alleles": 0,
        "in_scope_length_changing_alleles": 0,
        "complex_alleles_excluded": 0,
        "duplicate_pure_sv_alleles": 0,
    }
    for bubble in bubbles:
        reference_sequence = reference.sequence(
            bubble.chrom, bubble.start, bubble.end
        )
        allele_sequences = tuple(
            _walk_sequence(walk, segments) for walk in bubble.allele_walks
        )
        audit["declared_allele_length_mismatches"] += sum(
            len(sequence) != declared
            for sequence, declared in zip(
                allele_sequences, bubble.allele_lengths, strict=True
            )
        )
        exact_reference_paths = tuple(
            index
            for index, sequence in enumerate(allele_sequences)
            if sequence == reference_sequence
        )
        if exact_reference_paths:
            audit["bubbles_with_exact_reference_path"] += 1
        if len(exact_reference_paths) > 1:
            audit["multiple_exact_reference_path_bubbles"] += 1
        if bubble.reference_allele is None:
            audit["missing_reference_genotype_bubbles"] += 1
        else:
            audit["declared_reference_sequence_mismatches"] += (
                allele_sequences[bubble.reference_allele] != reference_sequence
            )
            audit["declared_reference_length_mismatches"] += (
                bubble.allele_lengths[bubble.reference_allele]
                != len(reference_sequence)
            )
        audit["reference_interval_bubbles"] += 1
        for allele_index, alternate_sequence in enumerate(allele_sequences):
            if alternate_sequence == reference_sequence:
                audit["reference_identical_alleles"] += 1
                continue
            delta = len(alternate_sequence) - len(reference_sequence)
            if not _eligible_length_delta(delta, evaluator_profile):
                audit["out_of_universe_length_alleles"] += 1
                continue
            audit["in_scope_length_changing_alleles"] += 1
            svtype = "INS" if delta > 0 else "DEL"
            pos, ref, alt = _normalized_alleles(
                chrom=bubble.chrom,
                start=bubble.start,
                end=bubble.end,
                reference_sequence=reference_sequence,
                alternate_sequence=alternate_sequence,
                reference=reference,
            )
            is_pure_insertion = delta > 0 and len(ref) == 1 and len(alt) > 1
            is_pure_deletion = delta < 0 and len(alt) == 1 and len(ref) > 1
            if not (is_pure_insertion or is_pure_deletion):
                # A length-changing graph allele can also contain substitutions.
                # It is outside this benchmark's pure INS/DEL truth contract.
                audit["complex_alleles_excluded"] += 1
                continue
            if len(alt) - len(ref) != delta:
                raise NovelTruthError(
                    "normalized graph allele length differs from reconstructed AWALK"
                )
            record = SvRecord(
                record_id=(
                    f"GRAPH_{bubble.chrom}_{bubble.start}_{bubble.end}_"
                    f"L{bubble.line_number}_A{allele_index}"
                ),
                chrom=bubble.chrom,
                pos=pos,
                end=pos if svtype == "INS" else pos + len(ref) - 1,
                svtype=svtype,
                svlen=delta,
                ref=ref,
                alt=alt,
                gt="./.",
            )
            if not record_shape_in_universe(record, evaluator_profile):
                continue
            if record.stable_key in seen:
                audit["duplicate_pure_sv_alleles"] += 1
                continue
            seen.add(record.stable_key)
            records.append(record)
    if not records:
        raise NovelTruthError("graph contains no eligible resolved SV alleles")
    audit["eligible_unique_pure_sv_alleles"] = len(records)
    return records, audit


def _graph_index(
    records: list[SvRecord],
) -> dict[tuple[str, str], tuple[list[int], list[int]]]:
    grouped: dict[tuple[str, str], list[tuple[int, int]]] = {}
    for index, record in enumerate(records):
        grouped.setdefault((record.chrom, record.svtype), []).append(
            (record.pos, index)
        )
    result = {}
    for key, rows in grouped.items():
        rows.sort(key=lambda value: (value[0], records[value[1]].stable_key))
        result[key] = ([value[0] for value in rows], [value[1] for value in rows])
    return result


def compatible_graph_matches(
    truth: SvRecord,
    graph_records: list[SvRecord],
    index: Mapping[tuple[str, str], tuple[list[int], list[int]]],
    evaluator_profile: dict[str, Any],
) -> list[tuple[SvRecord, SvMatch]]:
    group = index.get((truth.chrom, truth.svtype))
    if group is None:
        return []
    positions, indices = group
    distance = int(evaluator_profile["sv_match"]["max_breakpoint_distance"])
    lower = bisect.bisect_left(positions, truth.pos - distance)
    upper = bisect.bisect_right(positions, truth.pos + distance)
    matches = []
    for graph_index in indices[lower:upper]:
        graph = graph_records[graph_index]
        match = compatibility(graph, truth, evaluator_profile)
        if match is not None:
            matches.append((graph, match))
    return sorted(matches, key=lambda item: (-item[1].score, item[0].stable_key))


def materialize_plain_novel_truth(
    *,
    truth_vcf: Path,
    graph_records: list[SvRecord],
    evaluator_profile: dict[str, Any],
    output_vcf: Path,
    exclusion_ledger: Path,
) -> dict[str, int]:
    truths = load_vcf(truth_vcf, prefix="truth")
    by_id = {record.record_id: record for record in truths}
    if len(by_id) != len(truths):
        raise NovelTruthError("materialized truth contains duplicate IDs")
    graph_index = _graph_index(graph_records)
    matches_by_truth = {
        truth.record_id: compatible_graph_matches(
            truth, graph_records, graph_index, evaluator_profile
        )
        for truth in truths
    }
    excluded = {
        truth_id for truth_id, matches in matches_by_truth.items() if matches
    }
    counts = {
        "source_truth_records": len(truths),
        "graph_eligible_alleles": len(graph_records),
        "excluded_graph_represented": len(excluded),
        "novel_truth_records": len(truths) - len(excluded),
        "ambiguous_truth_matches": sum(
            len(matches) > 1 for matches in matches_by_truth.values()
        ),
    }
    if counts["novel_truth_records"] <= 0:
        raise NovelTruthError("graph-excluded novel truth universe is empty")

    output_vcf.parent.mkdir(parents=True, exist_ok=True)
    saw_columns = False
    written = 0
    with open_text(truth_vcf) as source, output_vcf.open(
        "w", encoding="utf-8"
    ) as output:
        for raw in source:
            if raw.startswith("##"):
                output.write(raw)
                continue
            if raw.startswith("#CHROM"):
                output.write(
                    "##pgbench_novel_truth=graph_exclusion:"
                    "pgbench_minigraph_novel_truth_v2\n"
                )
                output.write(raw)
                saw_columns = True
                continue
            if not raw.strip():
                continue
            if not saw_columns:
                raise NovelTruthError("truth VCF has no #CHROM header")
            fields = raw.rstrip("\n").split("\t")
            record_id = fields[2]
            if record_id not in by_id:
                raise NovelTruthError(f"truth record ID cannot be resolved: {record_id}")
            if record_id not in excluded:
                output.write(raw)
                written += 1
    if written != counts["novel_truth_records"]:
        raise NovelTruthError("novel truth output count differs from its audit")

    materialized = load_vcf(output_vcf, prefix="novel_truth")
    leakage_count = sum(
        bool(
            compatible_graph_matches(
                truth,
                graph_records,
                graph_index,
                evaluator_profile,
            )
        )
        for truth in materialized
    )
    counts["known_truth_leakage_count"] = leakage_count
    if leakage_count:
        raise NovelTruthError(
            "materialized novel truth retains a compatible graph allele"
        )

    exclusion_ledger.parent.mkdir(parents=True, exist_ok=True)
    with exclusion_ledger.open("w", encoding="utf-8", newline="") as output:
        output.write(
            "truth_event_id\tstatus\tgraph_allele_id\tmatch_score\t"
            "start_distance\tend_distance\tsize_similarity\t"
            "sequence_similarity\tcompatible_graph_alleles\n"
        )
        for truth in truths:
            matches = matches_by_truth[truth.record_id]
            if matches:
                graph, match = matches[0]
                sequence = (
                    "" if match.sequence_similarity is None else match.sequence_similarity
                )
                row: Iterable[object] = (
                    truth.record_id,
                    "graph_represented",
                    graph.record_id,
                    match.score,
                    match.start_distance,
                    match.end_distance,
                    match.size_similarity,
                    sequence,
                    len(matches),
                )
            else:
                row = (truth.record_id, "novel", "", "", "", "", "", "", 0)
            output.write("\t".join(map(str, row)) + "\n")
    return counts


def _verify_graph_lock(
    lock_path: Path,
    *,
    graph_gfa: Path,
    variation_calls: Path,
) -> dict[str, Any]:
    try:
        lock = yaml.safe_load(lock_path.read_text(encoding="utf-8"))
    except (OSError, yaml.YAMLError) as exc:
        raise NovelTruthError("cannot load graph asset lock") from exc
    if not isinstance(lock, dict) or lock.get("schema_version") != 1:
        raise NovelTruthError("graph asset lock schema is invalid")
    if lock.get("profile") != "svarp_minigraph_longread":
        raise NovelTruthError("novel truth requires the SVarp minigraph profile")
    excluded = lock.get("excluded_samples")
    if not isinstance(excluded, list) or not {"HG002", "NA24385"}.issubset(excluded):
        raise NovelTruthError("graph lock does not exclude HG002 aliases")
    assets = lock.get("assets")
    if not isinstance(assets, Mapping) or set(assets) != {
        "gfa",
        "variation_calls",
    }:
        raise NovelTruthError("graph lock has no exact variation-call asset set")
    for name, path in (("gfa", graph_gfa), ("variation_calls", variation_calls)):
        record = assets[name]
        if not isinstance(record, Mapping):
            raise NovelTruthError(f"graph lock {name} record is invalid")
        if (
            record.get("sha256") != sha256_file(path)
            or record.get("size_bytes") != path.stat().st_size
        ):
            raise NovelTruthError(f"graph lock identity mismatch for {name}")
    return lock


def prepare(args: argparse.Namespace) -> dict[str, Any]:
    profile = _load_profile(args.novel_truth_profile)
    evaluator_profile = load_evaluator_profile(args.evaluator_profile)
    if profile.get("matching") != {
        "source": "evaluator_profile.sv_match",
        **evaluator_profile["sv_match"],
    }:
        raise NovelTruthError(
            "novel-truth matcher differs from the frozen evaluator matcher"
        )
    _verify_graph_lock(
        args.graph_asset_lock,
        graph_gfa=args.graph_gfa,
        variation_calls=args.variation_calls,
    )
    bubbles = load_graph_bubbles(args.variation_calls)
    required_segments = {
        name
        for bubble in bubbles
        for walk in bubble.allele_walks
        for _, name in _walk_parts(walk)
    }
    segments = load_graph_segments(args.graph_gfa, required_segments)
    reference = IndexedReference(args.reference, args.reference_index)
    graph_records, graph_reconstruction = graph_alleles(
        bubbles=bubbles,
        segments=segments,
        reference=reference,
        evaluator_profile=evaluator_profile,
    )
    args.output_vcf.parent.mkdir(parents=True, exist_ok=True)
    work = Path(tempfile.mkdtemp(prefix=".novel-truth-", dir=args.output_vcf.parent))
    try:
        plain = work / "novel.truth.vcf"
        staged_ledger = work / "graph-exclusion.tsv"
        counts = materialize_plain_novel_truth(
            truth_vcf=args.truth_vcf,
            graph_records=graph_records,
            evaluator_profile=evaluator_profile,
            output_vcf=plain,
            exclusion_ledger=staged_ledger,
        )
        staged = work / "novel.truth.vcf.gz"
        subprocess.run(
            ["bcftools", "sort", "-Oz", "-o", str(staged), str(plain)],
            check=True,
        )
        subprocess.run(
            ["bcftools", "index", "--tbi", "--force", str(staged)],
            check=True,
        )
        output_index = Path(f"{args.output_vcf}.tbi")
        shutil.copyfile(staged, args.output_vcf)
        shutil.copyfile(Path(f"{staged}.tbi"), output_index)
        shutil.copyfile(staged_ledger, args.exclusion_ledger)
        audit = {
            "schema_version": 2,
            "contract": "pgbench_novel_truth_audit_v2",
            "profile_id": profile["profile"]["id"],
            "counts": counts,
            "graph_reconstruction": graph_reconstruction,
            "truth_vcf_sha256": sha256_file(args.truth_vcf),
            "reference_sha256": sha256_file(args.reference),
            "reference_index_sha256": sha256_file(args.reference_index),
            "graph_gfa_sha256": sha256_file(args.graph_gfa),
            "variation_calls_sha256": sha256_file(args.variation_calls),
            "graph_asset_lock_sha256": sha256_file(args.graph_asset_lock),
            "novel_truth_profile_sha256": profile["_sha256"],
            "evaluator_profile_sha256": evaluator_profile["_sha256"],
            "output_vcf_sha256": sha256_file(args.output_vcf),
            "output_index_sha256": sha256_file(output_index),
            "exclusion_ledger_sha256": sha256_file(args.exclusion_ledger),
            "known_truth_leakage_count": counts["known_truth_leakage_count"],
        }
        args.audit_json.parent.mkdir(parents=True, exist_ok=True)
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
    parser.add_argument("--truth-vcf", required=True, type=Path)
    parser.add_argument("--reference", required=True, type=Path)
    parser.add_argument("--reference-index", required=True, type=Path)
    parser.add_argument("--graph-gfa", required=True, type=Path)
    parser.add_argument("--variation-calls", required=True, type=Path)
    parser.add_argument("--graph-asset-lock", required=True, type=Path)
    parser.add_argument("--novel-truth-profile", required=True, type=Path)
    parser.add_argument("--evaluator-profile", required=True, type=Path)
    parser.add_argument("--output-vcf", required=True, type=Path)
    parser.add_argument("--exclusion-ledger", required=True, type=Path)
    parser.add_argument("--audit-json", required=True, type=Path)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    try:
        prepare(parse_args(argv))
    except (OSError, NovelTruthError, subprocess.CalledProcessError) as exc:
        print(f"prepare_novel_truth: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
