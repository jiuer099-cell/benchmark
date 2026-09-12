from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import pytest

SCRIPTS = Path(__file__).resolve().parents[2] / "workflow" / "scripts"
sys.path.insert(0, str(SCRIPTS))

from build_challenge_panel import build_challenge_panel  # noqa: E402
from build_pangenome_manifest import assign_stable_alleles  # noqa: E402


def _del(chrom: str, pos: int, record_id: str, size: int, extra: str = "") -> str:
    ref = "A" * (size + 1)
    return (
        f"{chrom}\t{pos}\t{record_id}\t{ref}\tA\t.\tPASS\t"
        f"SVTYPE=DEL;END={pos + size};SVLEN=-{size}{extra}\n"
    )


def _ins(chrom: str, pos: int, record_id: str, size: int, extra: str = "") -> str:
    alt = "A" + "T" * size
    return (
        f"{chrom}\t{pos}\t{record_id}\tA\t{alt}\t.\tPASS\t"
        f"SVTYPE=INS;END={pos};SVLEN={size}{extra}\n"
    )


def test_challenge_panel_is_blinded_and_unmatched_truth_is_unscorable(
    tmp_path: Path,
) -> None:
    population = tmp_path / "population.vcf"
    population.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        + _del("chr1", 20, "positive", 60, ";AF=0.02")
        + _ins("chr1", 180, "unmatched", 60, ";AF=0.10"),
        encoding="utf-8",
    )
    panel = tmp_path / "panel.vcf"
    allele_ledger = tmp_path / "alleles.tsv"
    assign_stable_alleles(population, panel, allele_ledger, namespace="PGSV")

    truth = tmp_path / "truth.vcf"
    truth.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG002\n"
        + _del("chr1", 20, "truth1", 60).rstrip("\n")
        + "\tGT\t0/1\n",
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
    assert summary["truth_negative_count"] == 0
    assert summary["truth_unscorable_count"] == 1
    challenge_text = challenge.read_text(encoding="utf-8")
    assert "positive" not in challenge_text.lower()
    assert "negative" not in challenge_text.lower()
    assert "\tGT\t./.\n" in challenge_text
    assert "0/1" not in challenge_text
    assert "0/0" not in challenge_text

    with hidden.open(encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert {row["truth_gt"] for row in rows} == {"UNSCORABLE", "0/1"}
    assert all(row["candidate_id"].startswith("PGSV_") for row in rows)
    assert json.loads(audit.read_text())["truth_labels_exposed_to_tool"] == 0


def test_candidate_ids_are_seed_deterministic(tmp_path: Path) -> None:
    population = tmp_path / "population.vcf"
    population.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        + _del("chr1", 20, "a", 60)
        + _ins("chr1", 180, "b", 60),
        encoding="utf-8",
    )
    panel = tmp_path / "panel.vcf"
    assign_stable_alleles(population, panel, tmp_path / "alleles.tsv", namespace="PGSV")
    truth = tmp_path / "truth.vcf"
    truth.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG002\n"
        + _del("chr1", 20, "t", 60).rstrip("\n") + "\tGT\t1/1\n",
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


def test_formal_candidate_universe_size_is_enforced(tmp_path: Path) -> None:
    population = tmp_path / "population.vcf"
    population.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        + _del("chr1", 20, "a", 60),
        encoding="utf-8",
    )
    panel = tmp_path / "panel.vcf"
    assign_stable_alleles(population, panel, tmp_path / "alleles.tsv", namespace="PGSV")
    truth = tmp_path / "truth.vcf"
    truth.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG002\n"
        + _del("chr1", 20, "truth", 60).rstrip("\n") + "\tGT\t0/1\n",
        encoding="utf-8",
    )

    with pytest.raises(Exception, match="canonical scoring universe size mismatch"):
        build_challenge_panel(
            panel_vcf=panel,
            truth_vcf=truth,
            output_vcf=tmp_path / "challenge.vcf",
            hidden_ledger=tmp_path / "hidden.tsv",
            audit_json=tmp_path / "audit.json",
            seed="fixed",
            expected_candidate_count=18164,
        )


def test_multisample_panel_columns_are_removed_from_blinded_vcf(
    tmp_path: Path,
) -> None:
    panel = tmp_path / "panel.vcf"
    panel.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tPANEL1\tPANEL2\n"
        + _ins("chr1", 10, "v1", 60, ";PANGENOME_ALLELE_ID=PGSV_1").rstrip("\n")
        + "\tGT\t0|1\t1|0\n"
        + _del("chr1", 120, "v2", 60, ";PANGENOME_ALLELE_ID=PGSV_2").rstrip("\n")
        + "\tGT\t0|0\t0|1\n",
        encoding="utf-8",
    )
    truth = tmp_path / "truth.vcf"
    truth.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG002\n"
        + _ins("chr1", 10, "t1", 60).rstrip("\n") + "\tGT\t0/1\n",
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
        + _del("chr1", 100, "p1", 100, ";PANGENOME_ALLELE_ID=PGSV_1")
        + _ins("chr1", 250, "p2", 80, ";PANGENOME_ALLELE_ID=PGSV_2")
        + _del("chr1", 500, "p3", 100, ";PANGENOME_ALLELE_ID=PGSV_3"),
        encoding="utf-8",
    )
    truth = tmp_path / "truth.vcf"
    truth.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG002\n"
        + _del("chr1", 100, "t1", 100).rstrip("\n") + "\tGT\t0/1\n",
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

    assert summary["candidate_count"] == 2
    assert summary["source_candidate_count"] == 3
    assert summary["truth_scorable_count"] == 1
    assert summary["truth_unscorable_count"] == 1
    assert summary["excluded_region_count"] == 1
    assert summary["truth_positive_count"] == 1
    assert summary["truth_negative_count"] == 0
    with hidden.open(encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    by_allele = {row["pangenome_allele_id"]: row for row in rows}
    assert by_allele["PGSV_1"]["truth_scorable"] == "1"
    assert by_allele["PGSV_2"]["truth_gt"] == "UNSCORABLE"
    assert by_allele["PGSV_2"]["truth_scorable"] == "0"
    assert "PGSV_3" not in by_allele


def test_formal_candidate_universe_exclusions_are_audited(tmp_path: Path) -> None:
    panel = tmp_path / "panel.vcf"
    panel.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        + _del("chr1", 100, "ok1", 100, ";PANGENOME_ALLELE_ID=PGSV_OK1")
        + _ins("chr1", 250, "ok2", 80, ";PANGENOME_ALLELE_ID=PGSV_OK2")
        + "chr1\t300\tmulti\tA\t<DEL>,<INS>\t.\tPASS\t"
        "PANGENOME_ALLELE_ID=PGSV_MULTI;SVTYPE=DEL;END=400;SVLEN=-100\n"
        "chr1\t420\tsymbolic\tA\t<DEL>\t.\tPASS\t"
        "PANGENOME_ALLELE_ID=PGSV_SYMBOLIC;SVTYPE=DEL;END=520;SVLEN=-100\n"
        "chr1\t450\ttype\tA\t<INV>\t.\tPASS\t"
        "PANGENOME_ALLELE_ID=PGSV_TYPE;SVTYPE=INV;END=550;SVLEN=100\n"
        + _del("chr1", 600, "small", 20, ";PANGENOME_ALLELE_ID=PGSV_SMALL")
        + _ins("chr1", 700, "filtered", 80, ";PANGENOME_ALLELE_ID=PGSV_FILTER").replace("\tPASS\t", "\tLowQual\t"),
        encoding="utf-8",
    )
    truth = tmp_path / "truth.vcf"
    truth.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG002\n"
        + _del("chr1", 100, "t1", 100).rstrip("\n") + "\tGT\t0/1\n",
        encoding="utf-8",
    )
    regions = tmp_path / "benchmark.bed"
    regions.write_text("chr1\t0\t1000\n", encoding="utf-8")
    challenge = tmp_path / "challenge.vcf"
    hidden = tmp_path / "hidden.tsv"
    summary = build_challenge_panel(
        panel_vcf=panel,
        truth_vcf=truth,
        benchmark_bed=regions,
        output_vcf=challenge,
        hidden_ledger=hidden,
        audit_json=tmp_path / "audit.json",
        seed="fixed",
    )

    assert summary["source_candidate_count"] == 7
    assert summary["candidate_count"] == 2
    assert summary["excluded_multiallelic_count"] == 1
    assert summary["excluded_svtype_count"] == 1
    assert summary["excluded_unresolved_sequence_count"] == 1
    assert summary["excluded_size_count"] == 1
    assert summary["excluded_filter_count"] == 1
    emitted = [
        line for line in challenge.read_text(encoding="utf-8").splitlines()
        if line and not line.startswith("#")
    ]
    with hidden.open(encoding="utf-8") as handle:
        ledger = list(csv.DictReader(handle, delimiter="\t"))
    assert len(emitted) == len(ledger) == summary["candidate_count"]
    assert all("," not in line.split("\t")[4] for line in emitted)
