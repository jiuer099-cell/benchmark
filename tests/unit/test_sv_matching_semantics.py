from __future__ import annotations

import csv
import itertools
import sys
from pathlib import Path

import pytest
import yaml

ROOT = Path(__file__).resolve().parents[2]
SCRIPTS = ROOT / "workflow" / "scripts"
sys.path.insert(0, str(SCRIPTS))

from run_formal_evaluator import (  # noqa: E402
    FormalEvaluatorError,
    event_query_indices,
    fully_contained,
    load_regions,
    parse_truvari,
    write_ledger,
)
from sv_matching import (  # noqa: E402
    SvRecord,
    SvMatchError,
    compatibility,
    load_evaluator_profile,
    load_vcf,
    one_to_one_match,
    query_record_in_universe,
)


def record(
    record_id: str,
    *,
    pos: int,
    end: int,
    svtype: str = "DEL",
    svlen: int = -100,
    gt: str = "0/1",
    filter_status: str = "PASS",
) -> SvRecord:
    return SvRecord(
        record_id=record_id,
        chrom="chr1",
        pos=pos,
        end=end,
        svtype=svtype,
        svlen=svlen,
        ref="A",
        alt=f"<{svtype}>",
        gt=gt,
        filter_status=filter_status,
    )


def test_tolerant_matching_is_one_to_one_and_deterministic() -> None:
    profile = load_evaluator_profile(ROOT / "config" / "evaluator_profile.yaml")
    truth = [record("T1", pos=1000, end=1100)]
    queries = [
        record("Q_far", pos=1030, end=1130),
        record("Q_near", pos=1005, end=1105),
    ]
    assignments = one_to_one_match(queries, truth, profile)
    assert set(assignments) == {1}
    assert assignments[1].truth_index == 0
    assert assignments[1].start_distance == 5


def test_matching_rejects_wrong_svtype_even_inside_tolerance() -> None:
    profile = load_evaluator_profile(ROOT / "config" / "evaluator_profile.yaml")
    truth = [record("T1", pos=1000, end=1100)]
    query = [record("Q1", pos=1000, end=1100, svtype="INS", svlen=100)]
    assert one_to_one_match(query, truth, profile) == {}


def test_sequence_similarity_matches_truvari_5_4_edlib_contract() -> None:
    pytest.importorskip("edlib")
    profile = load_evaluator_profile(ROOT / "config" / "evaluator_profile.yaml")
    query = SvRecord(
        "Q",
        "chr1",
        100,
        100,
        "INS",
        50,
        "A",
        "A" + "G" * 50,
        "0/1",
    )
    truth = SvRecord(
        "T",
        "chr1",
        100,
        100,
        "INS",
        50,
        "A",
        "A" + "G" * 47 + "CCC",
        "0/1",
    )
    match = compatibility(query, truth, profile)
    assert match is not None
    assert match.sequence_similarity == pytest.approx(99 / 102)


def test_global_matching_maximizes_cardinality_before_edge_score() -> None:
    """A greedy best-edge choice yields one TP here, while the optimum has two."""

    profile = load_evaluator_profile(ROOT / "config" / "evaluator_profile.yaml")
    truths = [
        record("T1", pos=1000, end=1100),
        record("T2", pos=1400, end=1500),
    ]
    queries = [
        # Q1 strongly prefers T1 but can also match T2.
        record("Q1", pos=1005, end=1105),
        # Q2 can only match T1.
        record("Q2", pos=600, end=700),
    ]

    assignments = one_to_one_match(queries, truths, profile)
    matched_ids = {
        (queries[query_index].record_id, truths[match.truth_index].record_id)
        for query_index, match in assignments.items()
    }
    assert matched_ids == {("Q1", "T2"), ("Q2", "T1")}


def test_global_matching_is_stable_across_input_permutations() -> None:
    profile = load_evaluator_profile(ROOT / "config" / "evaluator_profile.yaml")
    base_truths = [
        record("T1", pos=1000, end=1100),
        record("T2", pos=1400, end=1500),
    ]
    base_queries = [
        record("Q1", pos=1005, end=1105),
        record("Q2", pos=600, end=700),
    ]

    observed: set[tuple[tuple[str, str, float], ...]] = set()
    for queries in itertools.permutations(base_queries):
        for truths in itertools.permutations(base_truths):
            for _ in range(3):
                assignments = one_to_one_match(
                    list(queries),
                    list(truths),
                    profile,
                )
                observed.add(
                    tuple(
                        sorted(
                            (
                                queries[query_index].record_id,
                                truths[match.truth_index].record_id,
                                match.score,
                            )
                            for query_index, match in assignments.items()
                        )
                    )
                )

    assert len(observed) == 1
    assert {
        (query_id, truth_id)
        for query_id, truth_id, _score in next(iter(observed))
    } == {("Q1", "T2"), ("Q2", "T1")}


def test_ledger_separates_detection_genotype_and_no_call(tmp_path: Path) -> None:
    queries = [
        record("Q1", pos=1000, end=1100, gt="0/1"),
        record("Q2", pos=2000, end=2100, gt="./."),
        record("Q3", pos=3000, end=3100, gt="0/0"),
        record("Q4", pos=4000, end=4100, svtype="DUP", gt="0/1"),
    ]
    truths = [record("T1", pos=1000, end=1100, gt="1/1")]
    profile = load_evaluator_profile(ROOT / "config" / "evaluator_profile.yaml")
    assignments = one_to_one_match([queries[0]], truths, profile)
    path = tmp_path / "votes.tsv"
    counts = write_ledger(
        path,
        "truvari",
        queries,
        truths,
        assignments,
        {"Q1": True},
        scope_eligible_ids={"Q1", "Q2", "Q3"},
        event_eligible_ids={"Q1"},
        require_phase=False,
    )
    with path.open(encoding="utf-8", newline="") as handle:
        rows = {
            row["result_id"]: row
            for row in csv.DictReader(handle, delimiter="\t")
        }
    assert rows["Q1"]["detection_correct"] == "1"
    assert rows["Q1"]["genotype_scorable"] == "1"
    assert rows["Q1"]["genotype_correct"] == "0"
    assert rows["Q2"]["event_eligible"] == "0"
    assert rows["Q2"]["no_call"] == "1"
    assert rows["Q3"]["event_eligible"] == "0"
    assert rows["Q4"]["scope_eligible"] == "0"
    assert rows["Q4"]["event_eligible"] == "0"
    assert counts == {
        "rows": 4,
        "scope_eligible": 3,
        "scope_excluded": 1,
        "events": 1,
        "detection_correct": 1,
        "evaluator_resolved": 1,
        "evaluator_unresolved": 0,
        "genotype_scorable": 1,
        "genotype_correct": 0,
        "no_call": 1,
    }


def test_detection_event_semantics_follow_plugin_output_contract() -> None:
    queries = [
        record("VARIANT", pos=1000, end=1100, gt="0/1"),
        record("NO_GT", pos=2000, end=2100, gt="./."),
    ]
    scope = {"VARIANT", "NO_GT"}
    assert event_query_indices(
        queries,
        scope_eligible_ids=scope,
        candidate_output_contract="variant_sites",
    ) == [0, 1]
    assert event_query_indices(
        queries,
        scope_eligible_ids=scope,
        candidate_output_contract="all_sites",
    ) == [0]


def test_variant_sites_contract_rejects_hom_ref_records() -> None:
    queries = [record("HOM_REF", pos=1000, end=1100, gt="0/0")]
    with pytest.raises(FormalEvaluatorError, match="hom-ref"):
        event_query_indices(
            queries,
            scope_eligible_ids={"HOM_REF"},
            candidate_output_contract="variant_sites",
        )


def test_query_scope_is_frozen_by_type_size_filter_and_region(
    tmp_path: Path,
) -> None:
    profile = load_evaluator_profile(ROOT / "config" / "evaluator_profile.yaml")
    assert query_record_in_universe(
        record("PASS", pos=100, end=200),
        profile,
    )
    assert query_record_in_universe(
        record("UNFILTERED", pos=100, end=200, filter_status="."),
        profile,
    )
    assert not query_record_in_universe(
        record("SMALL", pos=100, end=120, svlen=-20),
        profile,
    )
    assert query_record_in_universe(
        record("MAXIMUM", pos=100, end=200, svlen=-10000),
        profile,
    )
    assert not query_record_in_universe(
        record("TOO_LARGE", pos=100, end=200, svlen=-10001),
        profile,
    )
    assert not query_record_in_universe(
        record("OTHER_TYPE", pos=100, end=200, svtype="INV"),
        profile,
    )
    assert not query_record_in_universe(
        record("LOW_QUALITY", pos=100, end=200, filter_status="LowQual"),
        profile,
    )

    bed = tmp_path / "regions.bed"
    bed.write_text("chr1\t99\t200\n", encoding="utf-8")
    regions = load_regions(bed)
    assert fully_contained(regions, record("INSIDE", pos=100, end=200))
    assert not fully_contained(
        regions,
        record("CROSSES_BOUNDARY", pos=100, end=201),
    )


def test_profile_rejects_maximum_below_minimum(tmp_path: Path) -> None:
    source = ROOT / "config" / "evaluator_profile.yaml"
    payload = yaml.safe_load(source.read_text(encoding="utf-8"))
    payload["universe"]["maximum_sv_size"] = 49
    profile = tmp_path / "invalid-maximum.yaml"
    profile.write_text(yaml.safe_dump(payload, sort_keys=False), encoding="utf-8")

    with pytest.raises(SvMatchError, match="maximum_sv_size"):
        load_evaluator_profile(profile)


def test_profile_rejects_evaluator_size_flag_drift(tmp_path: Path) -> None:
    source = ROOT / "config" / "evaluator_profile.yaml"
    payload = yaml.safe_load(source.read_text(encoding="utf-8"))
    payload["universe"]["maximum_sv_size"] = 9999
    profile = tmp_path / "drifted-evaluator-args.yaml"
    profile.write_text(yaml.safe_dump(payload, sort_keys=False), encoding="utf-8")

    with pytest.raises(SvMatchError, match="extra_args must freeze"):
        load_evaluator_profile(profile)


def test_profile_rejects_evaluator_fingerprint_drift(tmp_path: Path) -> None:
    source = ROOT / "config" / "evaluator_profile.yaml"
    payload = yaml.safe_load(source.read_text(encoding="utf-8"))
    payload["evaluators"]["truvari"]["expected_version_sha256"] = "0" * 64
    profile = tmp_path / "drifted-evaluator-fingerprint.yaml"
    profile.write_text(yaml.safe_dump(payload, sort_keys=False), encoding="utf-8")

    with pytest.raises(SvMatchError, match="expected_version_sha256"):
        load_evaluator_profile(profile)


def test_empty_query_vcf_is_allowed_only_when_explicit(tmp_path: Path) -> None:
    path = tmp_path / "empty.vcf"
    path.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
        encoding="utf-8",
    )
    assert load_vcf(path, prefix="query", allow_empty=True) == []


def test_symbolic_deletion_size_falls_back_to_end_span(tmp_path: Path) -> None:
    path = tmp_path / "symbolic.vcf"
    path.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG002\n"
        "chr1\t100\tDEL1\tA\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=175\tGT\t0/1\n",
        encoding="utf-8",
    )
    loaded = load_vcf(path, prefix="query")
    assert loaded[0].svlen == -75
    profile = load_evaluator_profile(ROOT / "config" / "evaluator_profile.yaml")
    assert query_record_in_universe(loaded[0], profile)


def test_resolved_deletion_without_info_uses_reference_span(
    tmp_path: Path,
) -> None:
    path = tmp_path / "resolved.vcf"
    ref = "A" + "C" * 60
    path.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG002\n"
        f"chr1\t100\tDEL1\t{ref}\tA\t.\tPASS\t.\tGT\t0/1\n",
        encoding="utf-8",
    )
    loaded = load_vcf(path, prefix="query")
    assert loaded[0].svtype == "DEL"
    assert loaded[0].end == 160
    assert loaded[0].svlen == -60
    profile = load_evaluator_profile(ROOT / "config" / "evaluator_profile.yaml")
    assert query_record_in_universe(loaded[0], profile)


def test_truvari_parser_accepts_modern_tp_comp_name(tmp_path: Path) -> None:
    query = tmp_path / "query.vcf"
    query.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chr1\t100\tQ1\tA\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=200\n"
        "chr1\t300\tQ2\tA\t<INS>\t.\tPASS\tSVTYPE=INS;END=300\n",
        encoding="utf-8",
    )
    artifacts = tmp_path / "truvari"
    artifacts.mkdir()
    for name, record in (
        (
            "tp-comp.vcf",
            "chr1\t100\tQ1\tA\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=200\n",
        ),
        (
            "fp.vcf",
            "chr1\t300\tQ2\tA\t<INS>\t.\tPASS\tSVTYPE=INS;END=300\n",
        ),
    ):
        (artifacts / name).write_text(
            "##fileformat=VCFv4.2\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            + record,
            encoding="utf-8",
        )
    order, votes = parse_truvari(query, artifacts)
    assert order == ["Q1", "Q2"]
    assert votes == {"Q1": True, "Q2": False}
