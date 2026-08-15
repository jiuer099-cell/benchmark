#!/usr/bin/env python3
"""Link canonical tool calls to the frozen pangenome allele ledger."""

from __future__ import annotations

import argparse
import bisect
import csv
import re
import sys
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

from build_pangenome_manifest import (
    default_variant_end,
    format_info,
    infer_svtype,
    parse_info,
)

LINK_STATUSES = {
    "in_panel_exact",
    "in_panel_equivalent",
    "out_of_panel",
    "ambiguous",
    "unresolved",
    "invalid_allele_id",
    "wrong_link",
}

LINK_HEADERS = (
    "##INFO=<ID=PANGENOME_LINK_STATUS,Number=1,Type=String,"
    'Description="PGBench pangenome linkage status">',
    "##INFO=<ID=PANGENOME_LINKED_ID,Number=1,Type=String,"
    'Description="Benchmark-inferred stable pangenome allele ID">',
    "##INFO=<ID=PANGENOME_CANDIDATE_COUNT,Number=1,Type=Integer,"
    'Description="Number of compatible pangenome candidates">',
)

INFO_HEADER_PATTERN = re.compile(r"^##INFO=<ID=([^,>]+)")

# Canonical records intentionally retain useful source/panel annotations.  Some
# upstream VCFs omit declarations for those annotations; htslib may tolerate
# that while Truvari correctly rejects the resulting VCF during translation.
# Known fields get their normative VCF types and all other retained fields get
# a conservative String declaration so the linked VCF is self-contained.
KNOWN_INFO_HEADERS = {
    "END": "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position\">",
    "SVTYPE": (
        "##INFO=<ID=SVTYPE,Number=1,Type=String,"
        'Description="Structural variant type">'
    ),
    "SVLEN": (
        "##INFO=<ID=SVLEN,Number=.,Type=Integer,"
        'Description="Structural variant length">'
    ),
    "PANGENOME_ALLELE_ID": (
        "##INFO=<ID=PANGENOME_ALLELE_ID,Number=1,Type=String,"
        'Description="Stable pangenome allele identifier">'
    ),
    "CONFLICT": (
        "##INFO=<ID=CONFLICT,Number=.,Type=String,"
        'Description="Conflicting source annotation retained during normalization">'
    ),
}


def missing_info_headers(
    declared: set[str], used: set[str]
) -> tuple[str, ...]:
    """Return deterministic declarations for INFO tags absent from the header."""

    generated = []
    for field_id in sorted(used - declared):
        generated.append(
            KNOWN_INFO_HEADERS.get(
                field_id,
                f'##INFO=<ID={field_id},Number=.,Type=String,'
                'Description="Retained upstream annotation">',
            )
        )
    return tuple(generated)


class AlleleLinkError(ValueError):
    """Raised when calls or the allele ledger are malformed."""


@dataclass(frozen=True)
class Allele:
    allele_id: str
    chrom: str
    pos: int
    end: int
    svtype: str
    svlen: int
    ref: str
    alt: str
    af: str
    graph_class: str

    @property
    def exact_key(self) -> tuple[str, int, int, str, int, str, str]:
        return (
            self.chrom,
            self.pos,
            self.end,
            self.svtype,
            self.svlen,
            self.ref,
            self.alt,
        )


@dataclass(frozen=True)
class Call:
    canon_id: str
    chrom: str
    pos: int
    end: int
    svtype: str
    svlen: int
    ref: str
    alt: str
    claimed_id: str | None

    @property
    def exact_key(self) -> tuple[str, int, int, str, int, str, str]:
        return (
            self.chrom,
            self.pos,
            self.end,
            self.svtype,
            self.svlen,
            self.ref,
            self.alt,
        )


@dataclass(frozen=True)
class LinkResult:
    status: str
    linked_id: str | None
    candidate_count: int
    method: str
    similarity: float | None


@dataclass(frozen=True)
class AlleleIndex:
    """Indexes that avoid scanning every relevant panel allele for every call."""

    exact: dict[tuple[str, int, int, str, int, str, str], list[Allele]]
    grouped: dict[tuple[str, str], tuple[list[int], list[Allele], float]]


def build_allele_index(alleles: list[Allele]) -> AlleleIndex:
    exact: dict[tuple[str, int, int, str, int, str, str], list[Allele]] = (
        defaultdict(list)
    )
    grouped_values: dict[tuple[str, str], list[Allele]] = defaultdict(list)
    for allele in alleles:
        exact[allele.exact_key].append(allele)
        grouped_values[(allele.chrom, allele.svtype)].append(allele)
    grouped: dict[tuple[str, str], tuple[list[int], list[Allele], float]] = {}
    for key, values in grouped_values.items():
        ordered = sorted(values, key=lambda allele: allele.pos)
        max_distance = max(
            max(100.0, 0.10 * max(1, abs(allele.svlen))) for allele in ordered
        )
        grouped[key] = (
            [allele.pos for allele in ordered],
            ordered,
            max_distance,
        )
    return AlleleIndex(exact=dict(exact), grouped=grouped)


def _candidate_alleles(call: Call, index: AlleleIndex) -> list[Allele]:
    grouped = index.grouped.get((call.chrom, call.svtype))
    if grouped is None:
        return []
    positions, alleles, max_distance = grouped
    left = bisect.bisect_left(positions, call.pos - max_distance)
    right = bisect.bisect_right(positions, call.pos + max_distance)
    return alleles[left:right]


def _allow_large_vcf_alleles() -> None:
    """Raise csv's legacy 128 KiB field cap without assuming C-long width."""

    limit = sys.maxsize
    while True:
        try:
            csv.field_size_limit(limit)
            return
        except OverflowError:
            limit //= 10


def _relevant_call_index(
    calls: list[Call],
) -> dict[tuple[str, str], tuple[list[int], list[Call]]]:
    grouped: dict[tuple[str, str], list[Call]] = defaultdict(list)
    for call in calls:
        grouped[(call.chrom, call.svtype)].append(call)
    return {
        key: (
            [call.pos for call in sorted_calls],
            sorted_calls,
        )
        for key, values in grouped.items()
        for sorted_calls in [sorted(values, key=lambda call: call.pos)]
    }


def _allele_is_relevant(
    allele: Allele,
    *,
    claimed_ids: set[str],
    call_index: dict[tuple[str, str], tuple[list[int], list[Call]]],
) -> bool:
    if allele.allele_id in claimed_ids:
        return True
    indexed = call_index.get((allele.chrom, allele.svtype))
    if indexed is None:
        return False
    positions, calls = indexed
    distance = max(100.0, 0.10 * max(1, abs(allele.svlen)))
    left = bisect.bisect_left(positions, allele.pos - distance)
    right = bisect.bisect_right(positions, allele.pos + distance)
    return any(compatibility(call, allele) is not None for call in calls[left:right])


def load_ledger(
    path: Path,
    *,
    relevant_calls: list[Call] | None = None,
) -> tuple[dict[str, Allele], list[Allele]]:
    _allow_large_vcf_alleles()
    alleles: dict[str, Allele] = {}
    claimed_ids = {
        call.claimed_id for call in (relevant_calls or []) if call.claimed_id
    }
    call_index = _relevant_call_index(relevant_calls or [])
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
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
        }
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise AlleleLinkError("allele ledger is missing required columns")
        for row in reader:
            allele = Allele(
                allele_id=row["PANGENOME_ALLELE_ID"],
                chrom=row["CHROM"],
                pos=int(row["POS"]),
                end=int(row["END"]),
                svtype=row["SVTYPE"],
                svlen=int(row["SVLEN"]),
                ref=row["REF"],
                alt=row["ALT"],
                af=row["AF"],
                graph_class=row["GRAPH_COMPLEXITY_CLASS"],
            )
            if relevant_calls is not None and not _allele_is_relevant(
                allele,
                claimed_ids=claimed_ids,
                call_index=call_index,
            ):
                continue
            if allele.allele_id in alleles:
                raise AlleleLinkError(f"duplicate allele ID {allele.allele_id}")
            alleles[allele.allele_id] = allele
    if not alleles and relevant_calls is None:
        raise AlleleLinkError("allele ledger is empty")
    return alleles, list(alleles.values())


def _int_info(info: dict[str, str | bool], key: str, default: int) -> int:
    value = info.get(key)
    if value is None or value is True:
        return default
    try:
        return int(str(value).split(",", 1)[0])
    except ValueError as exc:
        raise AlleleLinkError(f"{key} must be integer, found {value!r}") from exc


def parse_call(fields: list[str]) -> Call:
    if len(fields) < 8:
        raise AlleleLinkError("canonical VCF record has fewer than 8 fields")
    info = parse_info(fields[7])
    canon_id = str(info.get("CANON_ID", fields[2]))
    svtype = infer_svtype(fields[4], info, ref=fields[3])
    claimed = info.get("PANGENOME_ALLELE_ID")
    pos = int(fields[1])
    return Call(
        canon_id=canon_id,
        chrom=fields[0],
        pos=pos,
        end=_int_info(
            info,
            "END",
            default_variant_end(
                pos=pos,
                ref=fields[3],
                alt=fields[4],
                svtype=svtype,
            ),
        ),
        svtype=svtype,
        svlen=_int_info(info, "SVLEN", len(fields[4]) - len(fields[3])),
        ref=fields[3],
        alt=fields[4],
        claimed_id=str(claimed) if isinstance(claimed, str) else None,
    )


def compatibility(call: Call, allele: Allele) -> float | None:
    if call.chrom != allele.chrom or call.svtype != allele.svtype:
        return None
    truth_length = max(1, abs(allele.svlen))
    query_length = max(1, abs(call.svlen))
    scale = max(100.0, 0.10 * truth_length)
    start_score = max(0.0, 1.0 - abs(call.pos - allele.pos) / scale)
    end_score = max(0.0, 1.0 - abs(call.end - allele.end) / scale)
    size_score = min(query_length, truth_length) / max(query_length, truth_length)
    if start_score <= 0 or end_score <= 0 or size_score < 0.70:
        return None
    sequence_score = 1.0 if call.alt == allele.alt else 0.0
    return (
        0.30 * start_score
        + 0.25 * end_score
        + 0.25 * size_score
        + 0.20 * sequence_score
    )


def link_call(
    call: Call,
    allele_by_id: dict[str, Allele],
    alleles: list[Allele],
    *,
    allele_index: AlleleIndex | None = None,
) -> LinkResult:
    if call.svtype == "BND" and ("[" not in call.alt and "]" not in call.alt):
        return LinkResult("unresolved", None, 0, "unpaired_bnd", None)

    index = allele_index or build_allele_index(alleles)
    exact = index.exact.get(call.exact_key, [])
    compatible = [
        (allele, score)
        for allele in _candidate_alleles(call, index)
        if (score := compatibility(call, allele)) is not None
    ]

    if call.claimed_id:
        claimed = allele_by_id.get(call.claimed_id)
        if claimed is None:
            inferred = exact[0].allele_id if len(exact) == 1 else None
            return LinkResult(
                "invalid_allele_id",
                inferred,
                len(compatible),
                "claimed_id_absent",
                1.0 if inferred else None,
            )
        if call.exact_key == claimed.exact_key:
            return LinkResult(
                "in_panel_exact", claimed.allele_id, 1, "claimed_id_exact", 1.0
            )
        if len(exact) == 1 and exact[0].allele_id != claimed.allele_id:
            return LinkResult(
                "wrong_link",
                exact[0].allele_id,
                len(compatible),
                "claimed_id_conflicts_with_exact",
                1.0,
            )
        claimed_score = compatibility(call, claimed)
        if claimed_score is not None and len(compatible) == 1:
            return LinkResult(
                "in_panel_equivalent",
                claimed.allele_id,
                1,
                "claimed_id_equivalent",
                claimed_score,
            )
        return LinkResult(
            "wrong_link",
            compatible[0][0].allele_id if len(compatible) == 1 else None,
            len(compatible),
            "claimed_id_not_supported",
            compatible[0][1] if len(compatible) == 1 else None,
        )

    if len(exact) == 1:
        return LinkResult(
            "in_panel_exact", exact[0].allele_id, 1, "canonical_exact", 1.0
        )
    if len(exact) > 1 or len(compatible) > 1:
        return LinkResult(
            "ambiguous",
            None,
            max(len(exact), len(compatible)),
            "multiple_candidates",
            None,
        )
    if len(compatible) == 1:
        allele, score = compatible[0]
        return LinkResult(
            "in_panel_equivalent",
            allele.allele_id,
            1,
            "breakpoint_size_equivalent",
            score,
        )
    if call.alt in {"", ".", "<INS>"} and call.svtype == "INS":
        return LinkResult("unresolved", None, 0, "missing_insertion_sequence", None)
    return LinkResult("out_of_panel", None, 0, "no_compatible_panel_allele", None)


def link_vcf(
    canonical_vcf: Path,
    ledger_path: Path,
    output_vcf: Path,
    links_tsv: Path,
) -> int:
    calls: list[Call] = []
    declared_info: set[str] = set()
    used_info: set[str] = set()
    with canonical_vcf.open("r", encoding="utf-8") as source:
        for raw_line in source:
            match = INFO_HEADER_PATTERN.match(raw_line)
            if match:
                declared_info.add(match.group(1))
            if raw_line.strip() and not raw_line.startswith("#"):
                fields = raw_line.rstrip("\n").split("\t")
                calls.append(parse_call(fields))
                used_info.update(parse_info(fields[7]))
    supplemental_headers = missing_info_headers(declared_info, used_info)
    allele_by_id, alleles = load_ledger(
        ledger_path,
        relevant_calls=calls,
    )
    allele_index = build_allele_index(alleles)
    output_vcf.parent.mkdir(parents=True, exist_ok=True)
    links_tsv.parent.mkdir(parents=True, exist_ok=True)
    temporary_vcf = output_vcf.with_name(f".{output_vcf.name}.tmp")
    temporary_tsv = links_tsv.with_name(f".{links_tsv.name}.tmp")
    record_count = 0
    added_headers = False

    with (
        canonical_vcf.open("r", encoding="utf-8") as source,
        temporary_vcf.open("w", encoding="utf-8") as output,
        temporary_tsv.open("w", encoding="utf-8", newline="") as link_handle,
    ):
        writer = csv.DictWriter(
            link_handle,
            fieldnames=[
                "CANON_ID",
                "CLAIMED_PANGENOME_ALLELE_ID",
                "PANGENOME_ALLELE_ID",
                "LINK_STATUS",
                "CANDIDATE_COUNT",
                "MATCH_METHOD",
                "SIMILARITY",
                "PANEL_AF",
                "GRAPH_COMPLEXITY_CLASS",
            ],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for raw_line in source:
            if raw_line.startswith("##"):
                output.write(raw_line)
                continue
            if raw_line.startswith("#CHROM"):
                if not added_headers:
                    for header in supplemental_headers:
                        output.write(header + "\n")
                    for header in LINK_HEADERS:
                        output.write(header + "\n")
                    added_headers = True
                output.write(raw_line)
                continue
            if not raw_line.strip():
                continue
            fields = raw_line.rstrip("\n").split("\t")
            call = parse_call(fields)
            result = link_call(
                call,
                allele_by_id,
                alleles,
                allele_index=allele_index,
            )
            if result.status not in LINK_STATUSES:
                raise AlleleLinkError(f"unexpected link status {result.status}")
            info = parse_info(fields[7])
            info["PANGENOME_LINK_STATUS"] = result.status
            info["PANGENOME_CANDIDATE_COUNT"] = str(result.candidate_count)
            if result.linked_id is not None:
                info["PANGENOME_LINKED_ID"] = result.linked_id
            fields[7] = format_info(info)
            output.write("\t".join(fields) + "\n")
            linked = allele_by_id.get(result.linked_id or "")
            writer.writerow(
                {
                    "CANON_ID": call.canon_id,
                    "CLAIMED_PANGENOME_ALLELE_ID": call.claimed_id or "",
                    "PANGENOME_ALLELE_ID": result.linked_id or "",
                    "LINK_STATUS": result.status,
                    "CANDIDATE_COUNT": result.candidate_count,
                    "MATCH_METHOD": result.method,
                    "SIMILARITY": (
                        "" if result.similarity is None else f"{result.similarity:.12g}"
                    ),
                    "PANEL_AF": linked.af if linked else "",
                    "GRAPH_COMPLEXITY_CLASS": linked.graph_class if linked else "",
                }
            )
            record_count += 1
    temporary_vcf.replace(output_vcf)
    temporary_tsv.replace(links_tsv)
    return record_count


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Link canonical calls to the pangenome allele ledger."
    )
    parser.add_argument("--canonical-vcf", required=True, type=Path)
    parser.add_argument("--allele-ledger", required=True, type=Path)
    parser.add_argument("--output-vcf", required=True, type=Path)
    parser.add_argument("--links-tsv", required=True, type=Path)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        link_vcf(
            args.canonical_vcf,
            args.allele_ledger,
            args.output_vcf,
            args.links_tsv,
        )
    except (OSError, ValueError, AlleleLinkError) as exc:
        print(f"link_pangenome_alleles: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
