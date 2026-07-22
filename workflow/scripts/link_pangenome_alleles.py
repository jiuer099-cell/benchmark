#!/usr/bin/env python3
"""Link canonical tool calls to the frozen pangenome allele ledger."""

from __future__ import annotations

import argparse
import csv
import sys
from dataclasses import dataclass
from pathlib import Path

from build_pangenome_manifest import format_info, parse_info

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


def load_ledger(path: Path) -> tuple[dict[str, Allele], list[Allele]]:
    alleles: dict[str, Allele] = {}
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
            if allele.allele_id in alleles:
                raise AlleleLinkError(f"duplicate allele ID {allele.allele_id}")
            alleles[allele.allele_id] = allele
    if not alleles:
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
    svtype = str(info.get("SVTYPE", "OTHER"))
    claimed = info.get("PANGENOME_ALLELE_ID")
    return Call(
        canon_id=canon_id,
        chrom=fields[0],
        pos=int(fields[1]),
        end=_int_info(info, "END", int(fields[1])),
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
) -> LinkResult:
    if call.svtype == "BND" and ("[" not in call.alt and "]" not in call.alt):
        return LinkResult("unresolved", None, 0, "unpaired_bnd", None)

    exact = [allele for allele in alleles if call.exact_key == allele.exact_key]
    compatible = [
        (allele, score)
        for allele in alleles
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
    allele_by_id, alleles = load_ledger(ledger_path)
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
                    for header in LINK_HEADERS:
                        output.write(header + "\n")
                    added_headers = True
                output.write(raw_line)
                continue
            if not raw_line.strip():
                continue
            fields = raw_line.rstrip("\n").split("\t")
            call = parse_call(fields)
            result = link_call(call, allele_by_id, alleles)
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
