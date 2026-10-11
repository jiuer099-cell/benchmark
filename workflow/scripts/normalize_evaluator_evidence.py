#!/usr/bin/env python3
"""Normalize direct evaluator records to PG-F1 per-unit binary evidence."""
from __future__ import annotations

import argparse
import hashlib
import sys
from pathlib import Path

from pgf1 import EVIDENCE_FIELDS, PGF1Error, canonical_alt_dosage, load_units, open_text, read_tsv, sha256_file, write_tsv
from sv_matching import SvMatchError, load_evaluator_profile, load_vcf, query_record_in_universe

VERSION = "pgf1-evaluator-normalizer-v2"


def all_site_gts(path: Path) -> dict[str, str]:
    calls = {}
    with open_text(path, "rt") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) < 10 or fields[2] in {"", "."} or fields[2] in calls:
                raise PGF1Error("CANONICAL_MAPPING_FAILURE")
            keys, values = fields[8].split(":"), fields[9].split(":")
            calls[fields[2]] = values[keys.index("GT")] if "GT" in keys and keys.index("GT") < len(values) else "./."
    return calls


def normalize(*, units: Path, all_sites: Path, raw: Path, evaluator: str,
              release_id: str, run_id: str, sample_id: str, tool_id: str,
              out: Path, query: Path | None = None,
              evaluator_profile: Path | None = None) -> dict[str, int]:
    unit = load_units(units)
    required = {"result_id", "evaluator", "evaluator_accept", "evaluator_resolved",
                "mapping_status", "truth_event_id", "event_eligible"}
    by = {}
    for row in read_tsv(raw, required):
        candidate = row["result_id"]
        if row["evaluator"] != evaluator or not candidate or candidate in by:
            raise PGF1Error("NORMALIZER_RECONCILIATION_FAILURE")
        if row["mapping_status"] in {"partial_complex", "ambiguous", "many_to_many", "one_to_many", "many_to_one"}:
            raise PGF1Error("AMBIGUOUS_COMPLEX_MAPPING")
        if any(row[key] not in {"0", "1"} for key in ("evaluator_accept", "evaluator_resolved", "event_eligible")):
            raise PGF1Error("NORMALIZER_RECONCILIATION_FAILURE")
        if row["event_eligible"] == "0" and (row["evaluator_accept"] != "0" or row["truth_event_id"]):
            raise PGF1Error("NORMALIZER_RECONCILIATION_FAILURE")
        by[candidate] = row
    calls = all_site_gts(all_sites)
    candidates = {row["candidate_id"] for row in unit.values()}
    if set(calls) != candidates or not set(by).issubset(candidates):
        raise PGF1Error("UNIT_SET_MISMATCH")

    # The all-sites table has panel coordinates/filters. Only the actual native
    # evaluator query may establish that a positive call was properly excluded.
    scope = {}
    if (query is None) != (evaluator_profile is None):
        raise PGF1Error("NORMALIZER_RECONCILIATION_FAILURE: query and profile must be supplied together")
    if query is not None:
        profile = load_evaluator_profile(evaluator_profile)
        records = load_vcf(query, prefix="QUERY", allow_empty=True)
        for record in records:
            candidate = record.record_id
            if candidate in scope or candidate not in candidates:
                raise PGF1Error("UNIT_SET_MISMATCH")
            if canonical_alt_dosage(record.gt) != canonical_alt_dosage(calls[candidate]):
                raise PGF1Error("NORMALIZER_RECONCILIATION_FAILURE: query/all-sites dosage mismatch")
            scope[candidate] = query_record_in_universe(record, profile)
        if set(by) != set(scope):
            raise PGF1Error("CONSENSUS_EVIDENCE_INCOMPLETE: query/raw ledger set mismatch")
        for candidate, row in by.items():
            if row.get("scope_eligible") not in {"0", "1"} or (row["scope_eligible"] == "1") != scope[candidate]:
                raise PGF1Error("NORMALIZER_RECONCILIATION_FAILURE: query scope mismatch")
            if row["event_eligible"] == "1" and not scope[candidate]:
                raise PGF1Error("NORMALIZER_RECONCILIATION_FAILURE: excluded query marked as event")

    raw_sha = sha256_file(raw)
    script_sha = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    rows = []
    excluded_positive = set()
    for uid, value in sorted(unit.items()):
        candidate = value["candidate_id"]
        row = by.get(candidate)
        query_positive = canonical_alt_dosage(calls[candidate]) not in {None, 0}
        if row is None:
            if query_positive:
                raise PGF1Error("CONSENSUS_EVIDENCE_INCOMPLETE")
            row = {"truth_event_id": "", "result_id": "", "evaluator_accept": "0", "mapping_status": "core_absence"}
            event = resolved = False
        else:
            event = row["event_eligible"] == "1"
            resolved = row["evaluator_resolved"] == "1"
            if query_positive and not resolved:
                raise PGF1Error("CONSENSUS_EVIDENCE_INCOMPLETE")
            if query_positive and not event:
                # Fail closed without independent query/profile corroboration.
                if scope.get(candidate) is not False or row.get("scope_eligible") != "0" or row["mapping_status"] != "resolved":
                    raise PGF1Error("CONSENSUS_EVIDENCE_INCOMPLETE")
                row = {**row, "mapping_status": "scope_excluded_by_query_policy"}
                excluded_positive.add(candidate)
        relation = "AM" if event and resolved and row["evaluator_accept"] == "1" else "NONE"
        rows.append({"release_id": release_id, "run_id": run_id, "sample_id": sample_id, "tool_id": tool_id,
                     "unit_id": uid, "evaluator": evaluator, "truth_record_id": row.get("truth_event_id", ""),
                     "query_record_id": row.get("result_id", ""), "relation": relation, "binary_vote": int(relation == "AM"),
                     "mapping_cardinality": "one_to_one" if event else "not_applicable",
                     "mapping_status": "resolved" if event and resolved else "not_applicable",
                     "match_group_id": row.get("truth_event_id", ""), "superlocus_id": "",
                     "native_decision": row.get("evaluator_accept", "0"), "native_reason": row.get("mapping_status", ""),
                     "raw_evidence_sha256": raw_sha, "normalizer_version": VERSION, "normalizer_sha256": script_sha})
    write_tsv(out, EVIDENCE_FIELDS, rows)
    return {"units": len(rows), "raw_rows": len(by), "scope_excluded_positive": len(excluded_positive)}


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    for name in ("units", "all-sites", "raw-ledger", "query", "evaluator-profile", "output"):
        parser.add_argument("--" + name, required=True, type=Path)
    parser.add_argument("--evaluator", required=True, choices=("truvari", "aardvark", "vcfdist"))
    for name in ("release-id", "run-id", "sample-id", "tool-id"):
        parser.add_argument("--" + name, required=True)
    args = parser.parse_args(argv)
    try:
        normalize(units=args.units, all_sites=args.all_sites, raw=args.raw_ledger,
                  evaluator=args.evaluator, release_id=args.release_id, run_id=args.run_id,
                  sample_id=args.sample_id, tool_id=args.tool_id, out=args.output,
                  query=args.query, evaluator_profile=args.evaluator_profile)
    except (OSError, PGF1Error, SvMatchError) as error:
        print(f"normalize_evaluator_evidence: {error}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
