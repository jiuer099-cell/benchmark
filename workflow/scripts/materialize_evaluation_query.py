#!/usr/bin/env python3
"""Derive the non-reference evaluator query from the audited all-sites layer."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import sys
from collections import Counter
from contextlib import contextmanager
from pathlib import Path
from typing import Iterator, TextIO

from build_pangenome_manifest import parse_info
from materialize_all_sites import genotype, normalize_gt


class EvaluationQueryError(ValueError):
    """Raised when all-sites and candidate-status contracts disagree."""


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


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_status(path: Path) -> dict[str, dict[str, str]]:
    rows: dict[str, dict[str, str]] = {}
    required = {
        "candidate_id",
        "tool_status",
        "link_status",
        "addressability_status",
        "raw_gt",
        "canonical_gt",
        "final_status",
    }
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise EvaluationQueryError("candidate-status TSV has an invalid header")
        for row in reader:
            candidate_id = row["candidate_id"]
            if not candidate_id or candidate_id in rows:
                raise EvaluationQueryError(
                    f"empty or duplicate candidate status: {candidate_id!r}"
                )
            rows[candidate_id] = row
    if not rows:
        raise EvaluationQueryError("candidate-status TSV is empty")
    return rows


def load_canonical_panel(path: Path) -> dict[str, tuple[str, str, str, str]]:
    records: dict[str, tuple[str, str, str, str]] = {}
    with open_text(path, "rt") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8 or fields[2] in {"", "."}:
                raise EvaluationQueryError("canonical panel contains a malformed row")
            candidate_id = fields[2]
            if candidate_id in records:
                raise EvaluationQueryError(
                    f"duplicate canonical-panel candidate: {candidate_id}"
                )
            records[candidate_id] = (fields[0], fields[1], fields[3], fields[4])
    if not records:
        raise EvaluationQueryError("canonical panel is empty")
    return records


def materialize_evaluation_query(
    *,
    canonical_panel: Path,
    all_sites_vcf: Path,
    candidate_status_tsv: Path,
    output_vcf: Path,
    audit_json: Path,
) -> dict[str, object]:
    status_rows = load_status(candidate_status_tsv)
    panel_records = load_canonical_panel(canonical_panel)
    if set(panel_records) != set(status_rows):
        raise EvaluationQueryError(
            "candidate-status candidates differ from the frozen canonical panel"
        )
    seen: set[str] = set()
    counters: Counter[str] = Counter()
    output_vcf.parent.mkdir(parents=True, exist_ok=True)
    audit_json.parent.mkdir(parents=True, exist_ok=True)

    with open_text(all_sites_vcf, "rt") as source, open_text(output_vcf, "wt") as target:
        contract_header_written = False
        for line in source:
            if line.startswith("##"):
                target.write(line)
                continue
            if line.startswith("#CHROM"):
                target.write(
                    "##PGBENCH_EVALUATION_QUERY="
                    "non_reference_calls_from_audited_all_sites_v1\n"
                )
                contract_header_written = True
                target.write(line)
                continue
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 10:
                raise EvaluationQueryError("all-sites VCF contains a malformed row")
            candidate_id = fields[2]
            if not candidate_id or candidate_id == "." or candidate_id in seen:
                raise EvaluationQueryError(
                    f"invalid or duplicate all-sites candidate: {candidate_id!r}"
                )
            seen.add(candidate_id)
            row = status_rows.get(candidate_id)
            if row is None:
                raise EvaluationQueryError(
                    f"all-sites candidate lacks status row: {candidate_id}"
                )
            if (fields[0], fields[1], fields[3], fields[4]) != panel_records[candidate_id]:
                raise EvaluationQueryError(
                    f"all-sites allele differs from canonical panel: {candidate_id}"
                )
            canonical_gt = normalize_gt(genotype(fields))
            if canonical_gt != row["canonical_gt"]:
                raise EvaluationQueryError(
                    f"candidate status GT differs from all-sites VCF: {candidate_id}"
                )
            info_status = parse_info(fields[7]).get("PGBENCH_OUTPUT_STATUS")
            if info_status != row["tool_status"]:
                raise EvaluationQueryError(
                    f"candidate status differs from all-sites INFO: {candidate_id}"
                )
            final_status = row["final_status"]
            counters[final_status] += 1
            if final_status != "addressable_called":
                continue
            if canonical_gt == "0/0":
                counters["reference_genotype_excluded"] += 1
                continue
            if canonical_gt not in {"0/1", "1/1"}:
                raise EvaluationQueryError(
                    f"called candidate has unsupported canonical GT: {candidate_id}"
                )
            target.write(line)
            counters["non_reference_query_records"] += 1
        if not contract_header_written:
            raise EvaluationQueryError("all-sites VCF has no #CHROM header")

    missing = sorted(set(status_rows) - seen)
    if missing:
        raise EvaluationQueryError(
            "candidate-status TSV has candidates absent from all-sites VCF: "
            + ", ".join(missing[:5])
        )
    audit: dict[str, object] = {
        "schema_version": 1,
        "contract": "pgbench_evaluation_query_v1",
        "source_contract": "all_sites_plus_candidate_status",
        "canonical_candidate_count": len(status_rows),
        "all_sites_candidate_count": len(seen),
        "non_reference_query_count": counters["non_reference_query_records"],
        "reference_genotype_excluded_count": counters[
            "reference_genotype_excluded"
        ],
        "final_status_counts": {
            key: counters[key]
            for key in sorted(counters)
            if key not in {"non_reference_query_records", "reference_genotype_excluded"}
        },
        "all_sites_sha256": sha256_file(all_sites_vcf),
        "canonical_panel_sha256": sha256_file(canonical_panel),
        "candidate_status_sha256": sha256_file(candidate_status_tsv),
        "output_sha256": sha256_file(output_vcf),
        "unique_candidate_mapping": True,
        "canonical_alleles_verified": True,
        "extra_filtering_applied": False,
        "policy": {
            "0/1": "include",
            "1/1": "include",
            "0/0": "exclude_reference_genotype",
            "./.": "exclude_and_retain_in_candidate_status",
            "non_addressable": "exclude_and_retain_in_candidate_status",
        },
        "status": "valid",
    }
    temporary = audit_json.with_name(f".{audit_json.name}.tmp")
    temporary.write_text(
        json.dumps(audit, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    temporary.replace(audit_json)
    return audit


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--canonical-panel", required=True, type=Path)
    parser.add_argument("--all-sites-vcf", required=True, type=Path)
    parser.add_argument("--candidate-status-tsv", required=True, type=Path)
    parser.add_argument("--output-vcf", required=True, type=Path)
    parser.add_argument("--audit-json", required=True, type=Path)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        materialize_evaluation_query(
            canonical_panel=args.canonical_panel,
            all_sites_vcf=args.all_sites_vcf,
            candidate_status_tsv=args.candidate_status_tsv,
            output_vcf=args.output_vcf,
            audit_json=args.audit_json,
        )
    except (OSError, csv.Error, EvaluationQueryError) as exc:
        print(f"materialize_evaluation_query: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
