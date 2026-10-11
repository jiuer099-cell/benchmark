"""A filtered native query is explicit zero evidence, not a missing vote."""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "workflow/scripts"))
from normalize_evaluator_evidence import VERSION, normalize
from pgf1 import EVIDENCE_FIELDS, UNIT_FIELDS, PGF1Error, read_tsv, write_tsv


def fixture(tmp_path, *, filter_status="lowad", scope="0", event="0", resolved="1",
            accept="0", truth="", query_gt="0/1", all_gt="0/1", svlen=60,
            mapping_status="resolved", missing=False, evaluator="truvari"):
    units, sites, query, raw, out = (tmp_path / name for name in
                                    ("units.tsv", "all.vcf", "query.vcf", "raw.tsv", "evidence.tsv.gz"))
    write_tsv(units, UNIT_FIELDS, [{"unit_id": "C1:A1", "candidate_id": "C1", "allele_id": "C1",
                                  "chrom": "chr1", "pos": "10", "end": "70", "svtype": "DEL",
                                  "canonical_ref": "N", "canonical_alt": "<DEL>",
                                  "truth_gt": "0/1", "scorable_status": "SCORABLE"}])
    header = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG002\n"
    sites.write_text(header + f"chr1\t10\tC1\tN\t<DEL>\t60\tPASS\tEND=70;SVTYPE=DEL;SVLEN=-60\tGT\t{all_gt}\n")
    query.write_text(header + f"chr1\t10\tC1\tN\t<DEL>\t8.6106\t{filter_status}\tEND={10+svlen};SVTYPE=DEL;SVLEN=-{svlen}\tGT:GQ:DP\t{query_gt}:1:8\n")
    fields = ("result_id", "evaluator", "scope_eligible", "event_eligible", "evaluator_resolved",
              "evaluator_accept", "truth_event_id", "mapping_status")
    rows = [] if missing else [{"result_id": "C1", "evaluator": evaluator, "scope_eligible": scope,
                               "event_eligible": event, "evaluator_resolved": resolved,
                               "evaluator_accept": accept, "truth_event_id": truth,
                               "mapping_status": mapping_status}]
    write_tsv(raw, fields, rows)
    return dict(units=units, all_sites=sites, raw=raw, evaluator=evaluator, release_id="R",
                run_id="run", sample_id="HG002", tool_id="tool", out=out,
                query=query, evaluator_profile=ROOT / "config/evaluator_profile.yaml")


@pytest.mark.parametrize("evaluator", ["truvari", "aardvark", "vcfdist"])
def test_filtered_native_positive_is_zero_vote_without_losing_query(tmp_path, evaluator):
    args = fixture(tmp_path, evaluator=evaluator)
    before = args["query"].read_bytes()
    result = normalize(**args)
    assert result["scope_excluded_positive"] == 1
    row, = read_tsv(args["out"], EVIDENCE_FIELDS)
    assert row["binary_vote"] == "0" and row["relation"] == "NONE"
    assert row["query_record_id"] == "C1"
    assert row["native_reason"] == "scope_excluded_by_query_policy"
    assert row["mapping_cardinality"] == row["mapping_status"] == "not_applicable"
    assert row["normalizer_version"] == VERSION
    assert args["query"].read_bytes() == before


def test_shape_policy_exclusion_is_also_independently_checked(tmp_path):
    args = fixture(tmp_path, filter_status="PASS", svlen=49)
    assert normalize(**args)["scope_excluded_positive"] == 1


@pytest.mark.parametrize("changes", [
    {"filter_status": "PASS"}, {"filter_status": "."},
    {"scope": "1"}, {"event": "1"}, {"resolved": "0"},
    {"accept": "1"}, {"truth": "T1"}, {"mapping_status": "unresolved"},
    {"query_gt": "1/1"}, {"missing": True},
])
def test_inconsistent_or_missing_exclusion_evidence_still_fails(tmp_path, changes):
    args = fixture(tmp_path, **changes)
    with pytest.raises(PGF1Error):
        normalize(**args)
    assert not args["out"].exists()


def test_scope_ineligible_flag_alone_cannot_suppress_error(tmp_path):
    args = fixture(tmp_path)
    args.pop("query")
    args.pop("evaluator_profile")
    with pytest.raises(PGF1Error, match="CONSENSUS_EVIDENCE_INCOMPLETE"):
        normalize(**args)


def test_in_scope_positive_without_event_is_still_incomplete(tmp_path):
    with pytest.raises(PGF1Error, match="CONSENSUS_EVIDENCE_INCOMPLETE"):
        normalize(**fixture(tmp_path, filter_status="PASS", scope="1"))


@pytest.mark.parametrize("filter_status", ["PASS", "."])
def test_in_scope_direct_votes_are_unchanged(tmp_path, filter_status):
    args = fixture(tmp_path, filter_status=filter_status, scope="1", event="1", accept="1", truth="T1")
    assert normalize(**args)["scope_excluded_positive"] == 0
    row, = read_tsv(args["out"], EVIDENCE_FIELDS)
    assert row["binary_vote"] == "1" and row["relation"] == "AM"


def test_duplicate_native_query_identity_fails(tmp_path):
    args = fixture(tmp_path)
    args["query"].write_text(args["query"].read_text() * 2)
    with pytest.raises(PGF1Error, match="UNIT_SET_MISMATCH"):
        normalize(**args)


def test_query_and_profile_are_required_together(tmp_path):
    args = fixture(tmp_path)
    args.pop("evaluator_profile")
    with pytest.raises(PGF1Error, match="supplied together"):
        normalize(**args)
