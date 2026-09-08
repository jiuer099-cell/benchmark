#!/usr/bin/env python3
"""Project a tool result onto the frozen canonical candidate universe."""

from __future__ import annotations

import argparse
import csv
import gzip
import io
import json
import sys
from collections import defaultdict
from contextlib import contextmanager
from pathlib import Path
from typing import Iterator, TextIO

from build_pangenome_manifest import format_info, parse_info


class AllSitesError(ValueError):
    """Raised when canonical all-sites projection is ambiguous."""


NATIVE_STATUSES = {
    "addressable",
    "unsupported_representation",
    "adapter_conversion_failure",
    "index_build_failure",
    "linking_failure",
    "ambiguous_mapping",
}


@contextmanager
def open_text(path: Path, mode: str) -> Iterator[TextIO]:
    if path.suffix != ".gz":
        with path.open(mode, encoding="utf-8") as handle:
            yield handle
        return
    if "w" in mode:
        with path.open("wb") as raw:
            with gzip.GzipFile(
                filename="", fileobj=raw, mode="wb", mtime=0
            ) as compressed:
                with io.TextIOWrapper(compressed, encoding="utf-8") as text:
                    yield text
        return
    with gzip.open(path, mode, encoding="utf-8") as handle:
        yield handle


def genotype(fields: list[str]) -> str:
    if len(fields) < 10:
        return "./."
    keys = fields[8].split(":")
    values = fields[9].split(":")
    if "GT" not in keys or keys.index("GT") >= len(values):
        return "./."
    return values[keys.index("GT")]


def normalize_gt(value: str) -> str:
    if value in {".", "./.", ".|."}:
        return "./."
    alleles = value.replace("|", "/").split("/")
    if len(alleles) != 2 or any(not allele.isdigit() for allele in alleles):
        return "./."
    return "/".join(sorted(alleles, key=int))


def load_status_overrides(path: Path | None) -> dict[str, tuple[str, str]]:
    if path is None:
        return {}
    overrides: dict[str, tuple[str, str]] = {}
    with path.open("r", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            candidate_id = str(row.get("candidate_id", ""))
            status = str(row.get("tool_native_status", ""))
            reason = str(row.get("failure_reason", ""))
            if not candidate_id or status not in NATIVE_STATUSES:
                raise AllSitesError("invalid adapter addressability override")
            if candidate_id in overrides:
                raise AllSitesError(f"duplicate addressability override: {candidate_id}")
            overrides[candidate_id] = (status, reason)
    return overrides


def materialize(
    *,
    canonical_panel: Path,
    linked_query: Path,
    output_vcf: Path,
    addressability_tsv: Path,
    audit_json: Path,
    adapter_status: Path | None = None,
) -> dict[str, object]:
    headers: list[str] = []
    candidates: dict[str, list[str]] = {}
    allele_to_candidate: dict[str, str] = {}
    with open_text(canonical_panel, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                headers.append(line)
                continue
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8 or not fields[2] or fields[2] == ".":
                raise AllSitesError("canonical panel record lacks a candidate ID")
            if fields[2] in candidates:
                raise AllSitesError(f"duplicate canonical candidate: {fields[2]}")
            candidates[fields[2]] = fields
            allele_id = parse_info(fields[7]).get("PANGENOME_ALLELE_ID")
            if isinstance(allele_id, str) and allele_id:
                allele_to_candidate[allele_id] = fields[2]
    if not candidates:
        raise AllSitesError("canonical panel contains no candidates")

    mapped: dict[str, list[tuple[str, str, str]]] = defaultdict(list)
    off_panel_records = 0
    with open_text(linked_query, "rt") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                raise AllSitesError("linked query contains a malformed VCF record")
            info = parse_info(fields[7])
            direct = fields[2] if fields[2] in candidates else info.get("ORIG_ID")
            linked = info.get("PANGENOME_LINKED_ID")
            candidate_id = (
                str(direct)
                if isinstance(direct, str) and direct in candidates
                else allele_to_candidate.get(str(linked))
            )
            if candidate_id is None:
                off_panel_records += 1
                continue
            raw_gt = genotype(fields)
            mapped[candidate_id].append((fields[2], raw_gt, normalize_gt(raw_gt)))

    overrides = load_status_overrides(adapter_status)
    unknown_overrides = sorted(set(overrides) - set(candidates))
    if unknown_overrides:
        raise AllSitesError(
            "adapter status references unknown candidates: "
            + ", ".join(unknown_overrides[:5])
        )
    output_vcf.parent.mkdir(parents=True, exist_ok=True)
    addressability_tsv.parent.mkdir(parents=True, exist_ok=True)
    audit_json.parent.mkdir(parents=True, exist_ok=True)
    counts: dict[str, int] = defaultdict(int)
    fieldnames = [
        "candidate_id",
        "canonical_status",
        "tool_native_status",
        "output_status",
        "tool_status",
        "link_status",
        "addressability_status",
        "final_status",
        "raw_gt",
        "canonical_gt",
        "failure_reason",
        "native_variant_id",
    ]
    with (
        open_text(output_vcf, "wt") as output,
        addressability_tsv.open("w", encoding="utf-8", newline="") as status_handle,
    ):
        status_writer = csv.DictWriter(
            status_handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n"
        )
        status_writer.writeheader()
        inserted_header = False
        for line in headers:
            if line.startswith("#CHROM") and not inserted_header:
                output.write(
                    '##INFO=<ID=PGBENCH_OUTPUT_STATUS,Number=1,Type=String,'
                    'Description="Canonical output state">\n'
                )
                inserted_header = True
            output.write(line)

        for candidate_id, original in candidates.items():
            calls = mapped.get(candidate_id, [])
            native_status, failure_reason = overrides.get(
                candidate_id, ("addressable", "")
            )
            native_ids = sorted({native_id for native_id, _, _ in calls})
            raw_gts = {raw_gt for _, raw_gt, _ in calls}
            call_gts = {canonical_gt for _, _, canonical_gt in calls}
            if native_status != "addressable":
                output_status = native_status
                final_status = native_status
                link_status = (
                    "ambiguous_mapping"
                    if native_status == "ambiguous_mapping"
                    else "not_applicable"
                )
                raw_gt = "."
                gt = "./."
            elif len(call_gts) > 1:
                native_status = "ambiguous_mapping"
                output_status = "ambiguous_mapping"
                final_status = "ambiguous_mapping"
                link_status = "ambiguous_mapping"
                failure_reason = "conflicting linked genotypes"
                raw_gt = ",".join(sorted(raw_gts))
                gt = "./."
            elif calls:
                raw_gt = next(iter(raw_gts)) if len(raw_gts) == 1 else ",".join(sorted(raw_gts))
                gt = next(iter(call_gts))
                output_status = "explicit_no_call" if gt == "./." else "called"
                final_status = (
                    "explicit_no_call" if gt == "./." else "addressable_called"
                )
                link_status = "linked"
            else:
                raw_gt = "."
                gt = "./."
                output_status = "missing_output"
                final_status = "missing_output"
                link_status = "no_linked_output"
                failure_reason = failure_reason or "no linked native output record"
            counts[output_status] += 1
            counts[f"native_{native_status}"] += 1
            fields = original[:8]
            info = parse_info(fields[7])
            info["PGBENCH_OUTPUT_STATUS"] = output_status
            fields[7] = format_info(info)
            output.write("\t".join(fields + ["GT", gt]) + "\n")
            status_writer.writerow(
                {
                    "candidate_id": candidate_id,
                    "canonical_status": "main_track",
                    "tool_native_status": native_status,
                    "output_status": output_status,
                    "tool_status": output_status,
                    "link_status": link_status,
                    "addressability_status": native_status,
                    "final_status": final_status,
                    "raw_gt": raw_gt,
                    "canonical_gt": gt,
                    "failure_reason": failure_reason,
                    "native_variant_id": ",".join(native_ids),
                }
            )

    addressable = counts["native_addressable"]
    audit: dict[str, object] = {
        "schema_version": 1,
        "contract": "pgbench_addressability_v1",
        "canonical_candidate_count": len(candidates),
        "addressable_count": addressable,
        "unsupported_representation_count": counts[
            "native_unsupported_representation"
        ],
        "adapter_conversion_failure_count": counts[
            "native_adapter_conversion_failure"
        ],
        "index_build_failure_count": counts["native_index_build_failure"],
        "linking_failure_count": counts["native_linking_failure"],
        "ambiguous_mapping_count": counts["native_ambiguous_mapping"],
        "called_count": counts["called"],
        "explicit_no_call_count": counts["explicit_no_call"],
        "missing_output_count": counts["missing_output"],
        "addressability_rate": addressable / len(candidates),
        "off_panel_native_record_count": off_panel_records,
        "denominator_policy": "canonical_panel_fail_closed",
        "silent_candidate_deletion_count": 0,
    }
    temporary = audit_json.with_name(f".{audit_json.name}.tmp")
    temporary.write_text(
        json.dumps(audit, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    temporary.replace(audit_json)
    return audit


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--canonical-panel", required=True, type=Path)
    parser.add_argument("--linked-query", required=True, type=Path)
    parser.add_argument("--adapter-status", type=Path)
    parser.add_argument("--output-vcf", required=True, type=Path)
    parser.add_argument("--addressability-tsv", required=True, type=Path)
    parser.add_argument("--audit-json", required=True, type=Path)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        materialize(
            canonical_panel=args.canonical_panel,
            linked_query=args.linked_query,
            adapter_status=args.adapter_status,
            output_vcf=args.output_vcf,
            addressability_tsv=args.addressability_tsv,
            audit_json=args.audit_json,
        )
    except (OSError, AllSitesError, csv.Error) as exc:
        print(f"materialize_all_sites: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
