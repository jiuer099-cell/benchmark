from __future__ import annotations

import sys
from pathlib import Path

import pytest
import yaml

SCRIPTS = Path(__file__).resolve().parents[2] / "workflow" / "scripts"
sys.path.insert(0, str(SCRIPTS))

from build_pangenome_manifest import (  # noqa: E402
    PangenomeManifestError,
    assign_stable_alleles,
    build_manifest,
)


def _write_fixture(tmp_path: Path) -> tuple[Path, Path]:
    reference = tmp_path / "hs37d5.fa"
    reference.write_text(">1\n" + "A" * 200 + "\n", encoding="utf-8")
    population = tmp_path / "population.vcf"
    population.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "1\t20\tdel1\tA\t<DEL>\t.\tPASS\t"
        "SVTYPE=DEL;END=30;SVLEN=-10;AF=0.02\n"
        "1\t80\tins1\tA\tATTT\t.\tPASS\t"
        "SVTYPE=INS;END=80;SVLEN=3;AF=0.001;"
        "GRAPH_COMPLEXITY=multiallelic\n",
        encoding="utf-8",
    )
    return reference, population


def test_assigns_stable_ids_and_writes_manifest(tmp_path: Path) -> None:
    reference, population = _write_fixture(tmp_path)
    panel = tmp_path / "panel.vcf"
    ledger = tmp_path / "alleles.tsv"

    count = assign_stable_alleles(population, panel, ledger, namespace="PGSV")
    assert count == 2
    panel_text = panel.read_text(encoding="utf-8")
    assert "##INFO=<ID=PANGENOME_ALLELE_ID" in panel_text
    ids = [
        line.split("\t", 1)[0]
        for line in ledger.read_text(encoding="utf-8").splitlines()[1:]
    ]
    assert len(ids) == 2
    assert all(identifier.startswith("PGSV_") for identifier in ids)
    assert len(set(ids)) == 2

    manifest = build_manifest(
        pangenome_id="synthetic_pg_v1",
        backbone_id="synthetic_hs37d5",
        reference=reference,
        population_source_id="synthetic_population",
        population_source=population,
        panel_vcf=panel,
        allele_ledger=ledger,
        namespace="PGSV",
        excluded_truth_samples=[
            "HG002", "NA24385", "HG003", "NA24149", "HG004", "NA24143"
        ],
        graph_build_recipe_sha256=None,
        graph_assets_lock=None,
        generated_at="2026-07-17T00:00:00+00:00",
    )
    assert manifest["panel_vcf"]["record_count"] == 2
    assert manifest["truth_samples_excluded"] == [
        "HG002", "NA24385", "HG003", "NA24149", "HG004", "NA24143"
    ]
    assert len(manifest["target_family_exclusion"]) == 6
    assert len(manifest["panel_vcf"]["sha256"]) == 64

    round_trip = yaml.safe_load(yaml.safe_dump(manifest))
    assert round_trip == manifest


def test_stable_ids_are_reproducible(tmp_path: Path) -> None:
    _, population = _write_fixture(tmp_path)
    first_panel = tmp_path / "first.vcf"
    first_ledger = tmp_path / "first.tsv"
    second_panel = tmp_path / "second.vcf"
    second_ledger = tmp_path / "second.tsv"
    assign_stable_alleles(population, first_panel, first_ledger, namespace="PGSV")
    assign_stable_alleles(population, second_panel, second_ledger, namespace="PGSV")
    assert first_panel.read_text() == second_panel.read_text()
    assert first_ledger.read_text() == second_ledger.read_text()


def test_duplicate_canonical_alleles_are_rejected(tmp_path: Path) -> None:
    population = tmp_path / "duplicate.vcf"
    population.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "1\t20\ta\tA\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=30;SVLEN=-10\n"
        "1\t20\tb\tA\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=30;SVLEN=-10\n",
        encoding="utf-8",
    )
    with pytest.raises(PangenomeManifestError, match="duplicate"):
        assign_stable_alleles(
            population,
            tmp_path / "panel.vcf",
            tmp_path / "ledger.tsv",
            namespace="PGSV",
        )


def test_population_vcf_with_hg002_is_rejected(tmp_path: Path) -> None:
    population = tmp_path / "leaking.vcf"
    population.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG002\n"
        "1\t20\ta\tA\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=30\tGT\t0/1\n",
        encoding="utf-8",
    )
    with pytest.raises(PangenomeManifestError, match="excluded HG002"):
        assign_stable_alleles(
            population,
            tmp_path / "panel.vcf",
            tmp_path / "ledger.tsv",
            namespace="PGSV",
            excluded_samples=["HG002", "NA24385"],
        )
