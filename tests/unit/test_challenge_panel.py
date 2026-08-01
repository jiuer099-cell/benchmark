from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

SCRIPTS = Path(__file__).resolve().parents[2] / "workflow" / "scripts"
sys.path.insert(0, str(SCRIPTS))

from build_challenge_panel import build_challenge_panel  # noqa: E402
from build_pangenome_manifest import assign_stable_alleles  # noqa: E402


def test_challenge_panel_is_blinded_and_has_negative_sites(tmp_path: Path) -> None:
    population = tmp_path / "population.vcf"
    population.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "1\t20\tpositive\tA\t<DEL>\t.\tPASS\t"
        "SVTYPE=DEL;END=30;SVLEN=-10;AF=0.02\n"
        "1\t80\tnegative\tA\tATTT\t.\tPASS\t"
        "SVTYPE=INS;END=80;SVLEN=3;AF=0.10\n",
        encoding="utf-8",
    )
    panel = tmp_path / "panel.vcf"
    allele_ledger = tmp_path / "alleles.tsv"
    assign_stable_alleles(population, panel, allele_ledger, namespace="PGSV")

    truth = tmp_path / "truth.vcf"
    truth.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG002\n"
        "1\t20\ttruth1\tA\t<DEL>\t.\tPASS\t"
        "SVTYPE=DEL;END=30;SVLEN=-10\tGT\t0/1\n",
        encoding="utf-8",
    )
    challenge = tmp_path / "challenge.vcf"
    hidden = tmp_path / "hidden.tsv"
    audit = tmp_path / "audit.json"

    summary = build_challenge_panel(
        panel_vcf=panel,
        truth_vcf=truth,
        output_vcf=challenge,
        hidden_ledger=hidden,
        audit_json=audit,
        seed="fixed-seed",
    )
    assert summary["truth_positive_count"] == 1
    assert summary["truth_negative_count"] == 1
    challenge_text = challenge.read_text(encoding="utf-8")
    assert "positive" not in challenge_text.lower()
    assert "negative" not in challenge_text.lower()
    assert "\tGT\t./.\n" in challenge_text
    assert "0/1" not in challenge_text
    assert "0/0" not in challenge_text

    with hidden.open(encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert {row["truth_gt"] for row in rows} == {"0/0", "0/1"}
    assert all(row["candidate_id"].startswith("CAND_") for row in rows)
    assert json.loads(audit.read_text())["truth_labels_exposed_to_tool"] == 0


def test_candidate_ids_are_seed_deterministic(tmp_path: Path) -> None:
    population = tmp_path / "population.vcf"
    population.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "1\t20\ta\tA\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=30;SVLEN=-10\n"
        "1\t80\tb\tA\t<INS>\t.\tPASS\tSVTYPE=INS;END=80;SVLEN=3\n",
        encoding="utf-8",
    )
    panel = tmp_path / "panel.vcf"
    assign_stable_alleles(population, panel, tmp_path / "alleles.tsv", namespace="PGSV")
    truth = tmp_path / "truth.vcf"
    truth.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG002\n"
        "1\t20\tt\tA\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=30;SVLEN=-10\tGT\t1/1\n",
        encoding="utf-8",
    )

    outputs: list[str] = []
    for prefix in ("first", "second"):
        output = tmp_path / f"{prefix}.vcf"
        build_challenge_panel(
            panel_vcf=panel,
            truth_vcf=truth,
            output_vcf=output,
            hidden_ledger=tmp_path / f"{prefix}.tsv",
            audit_json=tmp_path / f"{prefix}.json",
            seed="same",
        )
        outputs.append(output.read_text())
    assert outputs[0] == outputs[1]


def test_multisample_panel_columns_are_removed_from_blinded_vcf(
    tmp_path: Path,
) -> None:
    panel = tmp_path / "panel.vcf"
    panel.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tPANEL1\tPANEL2\n"
        "chr1\t10\tv1\tA\tAT\t.\tPASS\t"
        "PANGENOME_ALLELE_ID=PGSV_1;SVTYPE=INS;END=10\tGT\t0|1\t1|0\n"
        "chr1\t20\tv2\tAT\tA\t.\tPASS\t"
        "PANGENOME_ALLELE_ID=PGSV_2;SVTYPE=DEL;END=21\tGT\t0|0\t0|1\n",
        encoding="utf-8",
    )
    truth = tmp_path / "truth.vcf"
    truth.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG002\n"
        "chr1\t10\tt1\tA\tAT\t.\tPASS\tSVTYPE=INS;END=10\tGT\t0/1\n",
        encoding="utf-8",
    )
    output = tmp_path / "challenge.vcf"
    build_challenge_panel(
        panel_vcf=panel,
        truth_vcf=truth,
        output_vcf=output,
        hidden_ledger=tmp_path / "hidden.tsv",
        audit_json=tmp_path / "audit.json",
        seed="fixed",
    )

    lines = output.read_text(encoding="utf-8").splitlines()
    header = next(line for line in lines if line.startswith("#CHROM"))
    assert header.split("\t") == [
        "#CHROM",
        "POS",
        "ID",
        "REF",
        "ALT",
        "QUAL",
        "FILTER",
        "INFO",
        "FORMAT",
        "HG002",
    ]
    assert all(
        len(line.split("\t")) == 10
        for line in lines
        if line and not line.startswith("#")
    )
    assert "PANEL1" not in output.read_text(encoding="utf-8")


def test_formal_candidate_truth_only_scores_frozen_bed_universe(
    tmp_path: Path,
) -> None:
    panel = tmp_path / "panel.vcf"
    panel.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chr1\t100\tp1\tA\t<DEL>\t.\tPASS\t"
        "PANGENOME_ALLELE_ID=PGSV_1;SVTYPE=DEL;END=200;SVLEN=-100\n"
        "chr1\t250\tp2\tA\t<INS>\t.\tPASS\t"
        "PANGENOME_ALLELE_ID=PGSV_2;SVTYPE=INS;END=250;SVLEN=80\n"
        "chr1\t500\tp3\tA\t<DEL>\t.\tPASS\t"
        "PANGENOME_ALLELE_ID=PGSV_3;SVTYPE=DEL;END=600;SVLEN=-100\n",
        encoding="utf-8",
    )
    truth = tmp_path / "truth.vcf"
    truth.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG002\n"
        "chr1\t100\tt1\tA\t<DEL>\t.\tPASS\t"
        "SVTYPE=DEL;END=200;SVLEN=-100\tGT\t0/1\n",
        encoding="utf-8",
    )
    regions = tmp_path / "benchmark.bed"
    regions.write_text("chr1\t0\t400\n", encoding="utf-8")
    hidden = tmp_path / "hidden.tsv"

    summary = build_challenge_panel(
        panel_vcf=panel,
        truth_vcf=truth,
        benchmark_bed=regions,
        output_vcf=tmp_path / "challenge.vcf",
        hidden_ledger=hidden,
        audit_json=tmp_path / "audit.json",
        seed="fixed",
    )

    assert summary["candidate_count"] == 3
    assert summary["truth_scorable_count"] == 2
    assert summary["truth_unscorable_count"] == 1
    assert summary["truth_positive_count"] == 1
    assert summary["truth_negative_count"] == 1
    with hidden.open(encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    by_allele = {row["pangenome_allele_id"]: row for row in rows}
    assert by_allele["PGSV_1"]["truth_scorable"] == "1"
    assert by_allele["PGSV_2"]["truth_gt"] == "0/0"
    assert by_allele["PGSV_3"]["truth_scorable"] == "0"
    assert by_allele["PGSV_3"]["truth_label"] == "unscorable"
    assert by_allele["PGSV_3"]["truth_gt"] == "./."
