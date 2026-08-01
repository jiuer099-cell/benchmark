from __future__ import annotations

import csv
import gzip
import json
import sys
from argparse import Namespace
from pathlib import Path

import pytest
import yaml
from jsonschema import Draft202012Validator

ROOT = Path(__file__).resolve().parents[2]
SCRIPTS = ROOT / "workflow" / "scripts"
sys.path.insert(0, str(SCRIPTS))

from materialize_formal_consensus_metrics import (  # noqa: E402
    candidate_genotype_summary,
    ConsensusMetricError,
    materialize,
    resource_summary,
    semantic_asset_hash,
    stratified_summary,
)
from run_formal_evaluator import (  # noqa: E402
    build_command,
    parse_aardvark,
    parse_vcfdist,
)
from sv_matching import SvRecord, load_evaluator_profile  # noqa: E402


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


def write_semantic_votes(
    path: Path,
    evaluator: str,
    *,
    a_correct: int,
    b_correct: int,
) -> None:
    fields = [
        "result_id",
        "evaluator",
        "scope_eligible",
        "event_eligible",
        "truth_event_id",
        "evaluator_accept",
        "detection_correct",
        "genotype_scorable",
        "genotype_correct",
        "no_call",
        "gt_state",
        "query_svtype",
        "query_svlen",
        "correct",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerow(
            {
                "result_id": "A",
                "evaluator": evaluator,
                "scope_eligible": 1,
                "event_eligible": 1,
                "truth_event_id": "T1",
                "evaluator_accept": a_correct,
                "detection_correct": a_correct,
                "genotype_scorable": 1,
                "genotype_correct": 1,
                "no_call": 0,
                "gt_state": "variant",
                "query_svtype": "INS",
                "query_svlen": 100,
                "correct": a_correct,
            }
        )
        writer.writerow(
            {
                "result_id": "B",
                "evaluator": evaluator,
                "scope_eligible": 1,
                "event_eligible": 1,
                "truth_event_id": "",
                "evaluator_accept": b_correct,
                "detection_correct": b_correct,
                "genotype_scorable": 0,
                "genotype_correct": 0,
                "no_call": 0,
                "gt_state": "variant",
                "query_svtype": "DEL",
                "query_svlen": -60,
                "correct": b_correct,
            }
        )
        writer.writerow(
            {
                "result_id": "C",
                "evaluator": evaluator,
                "scope_eligible": 1,
                "event_eligible": 0,
                "truth_event_id": "",
                "evaluator_accept": 0,
                "detection_correct": 0,
                "genotype_scorable": 0,
                "genotype_correct": 0,
                "no_call": 1,
                "gt_state": "no_call",
                "query_svtype": "INS",
                "query_svlen": 80,
                "correct": 0,
            }
        )


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


def test_evaluator_commands_freeze_the_common_maximum_size(
    tmp_path: Path,
) -> None:
    profile = load_evaluator_profile(ROOT / "config" / "evaluator_profile.yaml")
    common = {
        "query": tmp_path / "query.vcf.gz",
        "truth": tmp_path / "truth.vcf.gz",
        "reference": tmp_path / "reference.fa",
        "regions": tmp_path / "regions.bed",
        "artifacts": tmp_path / "artifacts",
        "threads": 4,
    }
    truvari = build_command(
        "truvari",
        ["truvari"],
        list(profile["evaluators"]["truvari"]["extra_args"]),
        **common,
    )
    vcfdist = build_command(
        "vcfdist",
        ["vcfdist"],
        list(profile["evaluators"]["vcfdist"]["extra_args"]),
        **common,
    )
    aardvark = build_command(
        "aardvark",
        ["aardvark"],
        list(profile["evaluators"]["aardvark"]["extra_args"]),
        **common,
    )
    assert truvari[truvari.index("--sizemax") + 1] == "10000"
    assert vcfdist[vcfdist.index("--largest-variant") + 1] == "10000"
    assert vcfdist[vcfdist.index("--max-supercluster-size") + 1] == "20000"
    assert vcfdist[vcfdist.index("--max-threads") + 1] == "4"
    assert "--sizemax" not in aardvark
    assert "--largest-variant" not in aardvark


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
        "chr1\t20\tT2\tAT\tA\t.\tPASS\tSVTYPE=DEL;END=21\n"
        "chr1\t30\tT3\tA\tAT\t.\tPASS\tSVTYPE=INS\n",
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
        reference=truth,
        pangenome_manifest=manifests[0],
        hidden_truth_ledger=manifests[0],
        graph_asset_lock=None,
    )
    payload = materialize(args)
    values = {
        record["metric_id"]: record["value"] for record in payload["records"]
    }
    assert values == {
        "benchmark.truth.eligible.count": 3,
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
        reference=tmp_path / "unused-reference.fa",
        pangenome_manifest=manifest,
        hidden_truth_ledger=manifest,
        graph_asset_lock=None,
    )
    try:
        materialize(args)
    except ConsensusMetricError:
        pass
    else:
        raise AssertionError("different result universes must be rejected")


def test_semantic_materialization_has_unique_truth_ci_and_strata(
    tmp_path: Path,
) -> None:
    paths: dict[str, Path] = {}
    for evaluator, decisions in {
        "truvari": (1, 0),
        "aardvark": (1, 0),
        "vcfdist": (0, 0),
    }.items():
        paths[evaluator] = tmp_path / f"{evaluator}.tsv"
        write_semantic_votes(
            paths[evaluator],
            evaluator,
            a_correct=decisions[0],
            b_correct=decisions[1],
        )
    manifests = []
    for evaluator in ("truvari", "aardvark", "vcfdist"):
        path = tmp_path / f"{evaluator}.manifest.json"
        path.write_text('{"status":"success"}\n', encoding="utf-8")
        manifests.append(path)
    truth = tmp_path / "truth.vcf"
    truth.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG002\n"
        "chr1\t10\tT1\tA\t<INS>\t.\tPASS\tSVTYPE=INS;END=10;SVLEN=100\tGT\t0/1\n"
        "chr1\t50\tT2\tA\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=110;SVLEN=-60\tGT\t0/1\n",
        encoding="utf-8",
    )
    regions = tmp_path / "benchmark.bed"
    regions.write_text("chr1\t0\t200\n", encoding="utf-8")
    query = tmp_path / "query.vcf"
    query.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG002\n"
        "chr1\t10\tA\tA\t<INS>\t.\tPASS\tSVTYPE=INS;SVLEN=100\tGT\t0/1\n"
        "chr1\t120\tB\tA\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=180;SVLEN=-60\tGT\t0/1\n"
        "chr1\t190\tC\tA\t<INS>\t.\tPASS\tSVTYPE=INS;SVLEN=80\tGT\t./.\n",
        encoding="utf-8",
    )
    resolved_inputs = tmp_path / "resolved-inputs.json"
    resolved_inputs.write_text(
        json.dumps(
            {
                "schema_version": 1,
                "run_id": "formal",
                "tool_id": "kanpig",
                "mode": "caller_only_shared_alignment",
                "alignment_kind": "bam",
                "inputs": [
                    {
                        "name": "shared_alignment",
                        "sha256": "a" * 64,
                        "size_bytes": 100,
                        "path_type": "file",
                    },
                    {
                        "name": "shared_alignment_index",
                        "sha256": "b" * 64,
                        "size_bytes": 10,
                        "path_type": "file",
                    },
                ],
            }
        ),
        encoding="utf-8",
    )
    payload = materialize(
        Namespace(
            truvari_ledger=paths["truvari"],
            aardvark_ledger=paths["aardvark"],
            vcfdist_ledger=paths["vcfdist"],
            evaluator_manifests=manifests,
            evaluator_completions=[],
            evaluator_profile=ROOT / "config" / "evaluator_profile.yaml",
            score_profile=ROOT / "config" / "consensus_scoring.yaml",
            run_id="formal",
            sample_id="HG002",
            tool_id="kanpig",
            official_score_mode="caller_only_shared_alignment",
            primary_truth_profile="giab_hg002_grch38_v5_0q",
            truth_vcf=truth,
            benchmark_bed=regions,
            reference=truth,
            pangenome_manifest=manifests[0],
            hidden_truth_ledger=manifests[0],
            graph_asset_lock=None,
            query_vcf=query,
            resolved_inputs=resolved_inputs,
            sample_technology="pacbio_clr",
            library_id="HG002_CLR",
            source_evidence_id="HG002_GRCh38_BAM",
            coverage_x=30.0,
            read_count=100,
            read_bases=10000,
            downsampling_seed=None,
            resource_benchmark=None,
            cache_policy="isolated_empty_tool_cache",
        )
    )
    analysis = payload["analysis"]
    assert analysis["matching"]["unique_truth_accounting"] is True
    assert analysis["matching"]["credited_truth_events"] == 1
    assert analysis["semantic_summary"]["truvari"]["no_call"] == 1
    assert analysis["comparable_score"] == pytest.approx(100 / 3)
    interval = analysis["comparable_score_confidence_interval"]
    assert interval["replicates"] == 1000
    assert interval["method"] == "paired_genomic_block_bootstrap"
    assert interval["bootstrap_unit"] == "fixed_genomic_block"
    assert interval["block_size_bp"] == 10_000_000
    assert interval["block_count"] == 1
    assert 0 <= interval["lower"] <= interval["upper"] <= 100
    assert {entry["stratum"] for entry in analysis["stratified_summary"]["svtype"]} == {
        "DEL",
        "INS",
    }
    assert analysis["asset_hashes"]["reference"]
    assert analysis["asset_hashes"]["pangenome_manifest"] == semantic_asset_hash(
        manifests[0]
    )
    assert analysis["asset_hashes"]["challenge_hidden_ledger"]
    assert analysis["asset_hashes"]["graph_asset_lock"] is None
    assert analysis["evidence_profile"]["actual_technology"] == "pacbio_clr"
    assert analysis["evidence_profile"]["evidence_kind"] == "shared_alignment"
    assert analysis["evidence_profile"]["quantitative_metadata_complete"] is True
    metrics_schema = yaml.safe_load(
        (ROOT / "workflow" / "schemas" / "metrics.schema.yaml").read_text(
            encoding="utf-8"
        )
    )
    Draft202012Validator(metrics_schema).validate(payload)


def test_resource_repeats_use_median_and_frozen_cold_cache(tmp_path: Path) -> None:
    benchmark = tmp_path / "tool.jsonl"
    benchmark.write_text(
        '{"s":10,"max_rss":100}\n'
        '{"s":30,"max_rss":300}\n'
        '{"s":20,"max_rss":200}\n',
        encoding="utf-8",
    )
    summary = resource_summary(benchmark, "isolated_empty_tool_cache")
    assert summary["cache_policy"] == "isolated_empty_tool_cache"
    assert summary["repeat_count"] == 3
    assert summary["median"]["wall_seconds"] == 20
    assert summary["median"]["max_rss_mb"] == 200


def test_semantic_asset_hash_ignores_paths_but_locks_content(
    tmp_path: Path,
) -> None:
    first = tmp_path / "first.yaml"
    second = tmp_path / "second.yaml"
    first.write_text(
        "generated_at: '2026-01-01T00:00:00Z'\n"
        "panel_vcf:\n"
        "  path: results/run-a/panel.vcf\n"
        f"  sha256: {'a' * 64}\n"
        "graph_assets:\n"
        "  asset_root: /server/a/graph\n"
        f"  gbz: {{path: /server/a/graph/graph.gbz, sha256: {'b' * 64}}}\n",
        encoding="utf-8",
    )
    second.write_text(
        "generated_at: '2026-02-02T00:00:00Z'\n"
        "panel_vcf:\n"
        "  path: results/run-b/panel.vcf\n"
        f"  sha256: {'a' * 64}\n"
        "graph_assets:\n"
        "  asset_root: /server/b/graph\n"
        f"  gbz: {{path: /server/b/graph/graph.gbz, sha256: {'b' * 64}}}\n",
        encoding="utf-8",
    )
    assert semantic_asset_hash(first) == semantic_asset_hash(second)

    second.write_text(
        second.read_text(encoding="utf-8").replace("b" * 64, "c" * 64),
        encoding="utf-8",
    )
    assert semantic_asset_hash(first) != semantic_asset_hash(second)
    pangenome_projection = frozenset(
        {"graph_assets", "graph_build_recipe_sha256"}
    )
    assert semantic_asset_hash(
        first,
        excluded_root_keys=pangenome_projection,
    ) == semantic_asset_hash(
        second,
        excluded_root_keys=pangenome_projection,
    )


def test_matched_query_is_assigned_to_truth_stratum() -> None:
    truth = [
        SvRecord(
            record_id="T1",
            chrom="chr1",
            pos=100,
            end=100,
            svtype="INS",
            svlen=100,
            gt="0/1",
            ref="A",
            alt="<INS>",
        )
    ]
    rows = {
        "Q1": {
            "event_eligible": True,
            "truth_event_id": "T1",
            "query_svtype": "DEL",
            "query_svlen": -10000,
        },
        "Q2": {
            "event_eligible": True,
            "truth_event_id": None,
            "query_svtype": "DEL",
            "query_svlen": -10000,
        },
    }
    summary = stratified_summary(truth, rows, {"Q1": 3, "Q2": 0})
    by_type = {entry["stratum"]: entry for entry in summary["svtype"]}
    by_length = {entry["stratum"]: entry for entry in summary["length_bin"]}
    assert by_type["INS"]["truth_events"] == 1
    assert by_type["INS"]["query_events"] == 1
    assert by_type["INS"]["soft_true_positive_count"] == 1
    assert by_type["DEL"]["truth_events"] == 0
    assert by_type["DEL"]["query_events"] == 1
    assert by_length["100_499"]["truth_events"] == 1
    assert by_length["100_499"]["query_events"] == 1
    assert by_length["ge10000"]["truth_events"] == 0
    assert by_length["ge10000"]["query_events"] == 1


def test_candidate_gt_and_no_call_use_hidden_site_semantics(tmp_path: Path) -> None:
    hidden = tmp_path / "hidden.tsv"
    hidden.write_text(
        "candidate_id\tpangenome_allele_id\ttruth_gt\ttruth_label\t"
        "truth_scorable\ttruth_event_id\n"
        "CAND_A\tPG_A\t0/1\tpositive\t1\tT1\n"
        "CAND_B\tPG_B\t0/0\tnegative\t1\t\n"
        "CAND_C\tPG_C\t1/1\tpositive\t1\tT3\n"
        "CAND_D\tPG_D\t0/0\tnegative\t1\t\n"
        "CAND_E\tPG_E\t./.\tunscorable\t0\t\n",
        encoding="utf-8",
    )
    query = tmp_path / "linked.vcf"
    query.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG002\n"
        "chr1\t10\tCANON_A\tA\t<INS>\t.\tPASS\t"
        "ORIG_ID=CAND_A;SVTYPE=INS;SVLEN=100\tGT\t0/1\n"
        "chr1\t20\tCANON_B\tA\t<DEL>\t.\tPASS\t"
        "ORIG_ID=CAND_B;PANGENOME_LINKED_ID=PG_A;"
        "PANGENOME_LINK_STATUS=in_panel_exact;"
        "SVTYPE=DEL;END=120;SVLEN=-100\tGT\t./.\n"
        "chr1\t30\tCANON_C\tA\t<INS>\t.\tPASS\t"
        "PANGENOME_LINKED_ID=PG_C;"
        "PANGENOME_LINK_STATUS=in_panel_equivalent;"
        "SVTYPE=INS;SVLEN=100\tGT\t1/1\n",
        encoding="utf-8",
    )
    manifest = tmp_path / "tool.yaml"
    manifest.write_text(
        "outputs:\n"
        "  candidate_output_contract: variant_sites\n"
        "  absence_semantics: hom_ref\n",
        encoding="utf-8",
    )
    summary = candidate_genotype_summary(
        query_vcf=query,
        hidden_truth_ledger=hidden,
        tool_manifest=manifest,
        require_phase=False,
    )
    assert summary["contract"] == "hidden_candidate_genotype_v3"
    assert summary["candidate_count"] == 5
    assert summary["truth_scorable"] == 4
    assert summary["truth_positive"] == 2
    assert summary["truth_negative"] == 2
    assert summary["truth_unscorable"] == 1
    assert summary["observed_candidates"] == 3
    assert summary["inferred_absent_candidates"] == 2
    assert summary["inferred_absent_truth_scorable_candidates"] == 1
    assert summary["direct_candidate_records"] == 2
    assert summary["linked_candidate_records"] == 1
    assert summary["candidate_link_conflicts"] == 1
    assert summary["no_call"] == 1
    assert summary["genotype_scorable"] == 4
    assert summary["called_candidates"] == 3
    assert summary["genotype_correct"] == 3
    assert summary["genotype_accuracy"] == 0.75
    assert summary["called_only_genotype_accuracy"] == 1.0
    assert summary["binary_nonref_confusion"] == {
        "true_positive": 2,
        "false_positive": 0,
        "true_negative": 1,
        "false_negative": 0,
        "no_call_truth_positive": 0,
        "no_call_truth_negative": 1,
    }
    assert summary["nonref_precision"] == 1.0
    assert summary["nonref_sensitivity"] == 1.0
    assert summary["specificity"] == 0.5
    assert summary["balanced_accuracy"] == 0.75
    assert summary["nonref_f1"] == 1.0
    assert summary["genotype_macro_f1"] == pytest.approx(
        (2 / 3 + 1.0 + 1.0) / 3
    )
    assert summary["no_call_by_truth_class"] == {
        "heterozygous": 0,
        "hom_alt": 0,
        "hom_ref": 1,
    }
