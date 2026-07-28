from __future__ import annotations

import csv
import gzip
import sys
from argparse import Namespace
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SCRIPTS = ROOT / "workflow" / "scripts"
sys.path.insert(0, str(SCRIPTS))

from materialize_formal_consensus_metrics import (  # noqa: E402
    ConsensusMetricError,
    materialize,
)
from run_formal_evaluator import parse_aardvark, parse_vcfdist  # noqa: E402


def write_query(path: Path) -> None:
    path.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG002\n"
        "chr1\t10\tCANON_A\tA\tAT\t.\tPASS\tSVTYPE=INS\tGT\t0/1\n"
        "chr1\t20\tCANON_B\tAT\tA\t.\tPASS\tSVTYPE=DEL\tGT\t0/1\n",
        encoding="utf-8",
    )


def write_votes(path: Path, evaluator: str, votes: list[tuple[str, int]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["result_id", "evaluator", "correct"])
        for result_id, correct in votes:
            writer.writerow([result_id, evaluator, correct])


def test_aardvark_and_vcfdist_are_mapped_to_same_query_universe(
    tmp_path: Path,
) -> None:
    query = tmp_path / "query.vcf"
    write_query(query)

    aardvark = tmp_path / "aardvark"
    aardvark.mkdir()
    with gzip.open(aardvark / "query.vcf.gz", "wt", encoding="utf-8") as handle:
        handle.write(
            "##fileformat=VCFv4.2\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG002\n"
            "chr1\t10\tCANON_A\tA\tAT\t.\tPASS\t.\tGT:BD\t0/1:TP\n"
            "chr1\t20\tCANON_B\tAT\tA\t.\tPASS\t.\tGT:BD\t0/1:FP\n"
        )
    order, votes = parse_aardvark(query, aardvark)
    assert order == ["CANON_A", "CANON_B"]
    assert votes == {"CANON_A": True, "CANON_B": False}

    vcfdist = tmp_path / "vcfdist"
    vcfdist.mkdir()
    (vcfdist / "query.tsv").write_text(
        "CONTIG\tPOS\tHAP\tREF\tALT\tCREDIT\n"
        "chr1\t9\t0\tA\tAT\t1.0\n"
        "chr1\t19\t0\tAT\tA\t0.3\n",
        encoding="utf-8",
    )
    order, votes = parse_vcfdist(query, vcfdist, credit_threshold=0.7)
    assert order == ["CANON_A", "CANON_B"]
    assert votes == {"CANON_A": True, "CANON_B": False}


def test_materialize_formal_consensus_counts(tmp_path: Path) -> None:
    paths = {}
    for evaluator, votes in {
        "truvari": [("A", 1), ("B", 1), ("C", 0), ("D", 0)],
        "aardvark": [("A", 1), ("B", 1), ("C", 1), ("D", 0)],
        "vcfdist": [("A", 1), ("B", 0), ("C", 0), ("D", 0)],
    }.items():
        path = tmp_path / f"{evaluator}.tsv"
        write_votes(path, evaluator, votes)
        paths[evaluator] = path
    manifests = []
    for index in range(3):
        path = tmp_path / f"manifest-{index}.json"
        path.write_text(f'{{"id": {index}}}\n', encoding="utf-8")
        manifests.append(path)
    truth = tmp_path / "truth.vcf"
    truth.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chr1\t10\tT1\tA\tAT\t.\tPASS\tSVTYPE=INS\n"
        "chr1\t20\tT2\tAT\tA\t.\tPASS\tSVTYPE=DEL;END=21\n",
        encoding="utf-8",
    )
    regions = tmp_path / "benchmark.bed"
    regions.write_text("chr1\t0\t100\n", encoding="utf-8")
    args = Namespace(
        truvari_ledger=paths["truvari"],
        aardvark_ledger=paths["aardvark"],
        vcfdist_ledger=paths["vcfdist"],
        evaluator_manifests=manifests,
        score_profile=ROOT / "config" / "consensus_scoring.yaml",
        run_id="formal",
        sample_id="HG002",
        tool_id="kanpig",
        official_score_mode="caller_only_shared_alignment",
        primary_truth_profile="giab_hg002_grch38_v5_0q",
        truth_vcf=truth,
        benchmark_bed=regions,
    )
    payload = materialize(args)
    values = {
        record["metric_id"]: record["value"] for record in payload["records"]
    }
    assert values == {
        "benchmark.truth.eligible.count": 2,
        "consensus.all_three_correct.count": 1,
        "consensus.exactly_two_correct.count": 1,
        "consensus.exactly_one_correct.count": 1,
        "consensus.none_correct.count": 1,
    }


def test_materialize_rejects_different_result_universes(tmp_path: Path) -> None:
    paths = {}
    for evaluator, votes in {
        "truvari": [("A", 1)],
        "aardvark": [("A", 1)],
        "vcfdist": [("B", 1)],
    }.items():
        path = tmp_path / f"{evaluator}.tsv"
        write_votes(path, evaluator, votes)
        paths[evaluator] = path
    manifest = tmp_path / "manifest.json"
    manifest.write_text("{}\n", encoding="utf-8")
    args = Namespace(
        truvari_ledger=paths["truvari"],
        aardvark_ledger=paths["aardvark"],
        vcfdist_ledger=paths["vcfdist"],
        evaluator_manifests=[manifest],
        score_profile=ROOT / "config" / "consensus_scoring.yaml",
        run_id="formal",
        sample_id="HG002",
        tool_id="kanpig",
        official_score_mode="caller_only_shared_alignment",
        primary_truth_profile="giab_hg002_grch38_v5_0q",
        truth_vcf=tmp_path / "unused.vcf",
        benchmark_bed=tmp_path / "unused.bed",
    )
    try:
        materialize(args)
    except ConsensusMetricError:
        pass
    else:
        raise AssertionError("different result universes must be rejected")
