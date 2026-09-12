from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest
import yaml


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "workflow" / "scripts" / "prepare_pangenie_private_panel.py"
SPEC = importlib.util.spec_from_file_location("pangenie_private", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)

FAMILY = {"HG002", "NA24385", "HG003", "NA24149", "HG004", "NA24143"}


def _context(tmp_path: Path) -> tuple[Path, Path, Path, Path, Path, Path, Path, Path, Path]:
    canonical_source = tmp_path / "canonical.vcf"
    canonical_source.write_text("frozen source\n", encoding="utf-8")
    graph = tmp_path / "source.gfa.gz"
    graph.write_text("graph\n", encoding="utf-8")
    haplotypes = tmp_path / "haplotypes.yaml"
    haplotypes.write_text("samples: [PANEL1]\n", encoding="utf-8")
    frozen_gbz = tmp_path / "frozen.gbz"
    frozen_gbz.write_text("frozen graph\n", encoding="utf-8")
    sample_manifest = tmp_path / "samples.txt"
    sample_manifest.write_text("PANEL1\n", encoding="utf-8")
    private = tmp_path / "private.vcf"
    private.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tPANEL1\n"
        "chr1\t10\tNATIVE1\tA\tAT\t.\tPASS\t.\tGT\t0|1\n",
        encoding="utf-8",
    )
    scoring = tmp_path / "scoring.vcf"
    scoring.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG002\n"
        "chr1\t10\tCAND_A\tA\tAT\t.\tPASS\t.\tGT\t./.\n",
        encoding="utf-8",
    )
    provenance = tmp_path / "provenance.yaml"
    yaml.safe_dump(
        {
            "contract": MODULE.CONTRACT,
            "source_cohort_id": "frozen_hprc",
            "source_graph_sha256": MODULE.sha256(graph),
            "source_haplotype_manifest_sha256": MODULE.sha256(haplotypes),
            "official_source_pgin_sha256": MODULE.sha256(private),
            "reference_build": "grch38",
        },
        provenance.open("w", encoding="utf-8"),
        sort_keys=True,
    )
    identity = tmp_path / "identity.yaml"
    yaml.safe_dump(
        {
            "contract": "pgbench_hprc_graph_identity_v1",
            "status": "verified",
            "source_cohort_id": "frozen_hprc",
            "current_gbz_sha256": MODULE.sha256(frozen_gbz),
            "source_sample_manifest_sha256": MODULE.sha256(sample_manifest),
            "source_release": "test-release",
            "graph_family": "minigraph-cactus",
            "graph_construction_version": "test-v1",
            "reference_build": "grch38",
            "official_manifest_url": "https://example.test/manifest",
            "official_manifest_sha256": "a" * 64,
            "official_pgin_url": "https://example.test/panel.pgin.vcf.gz",
            "official_pgin_sha256": "b" * 64,
        },
        identity.open("w", encoding="utf-8"),
        sort_keys=True,
    )
    return canonical_source, graph, haplotypes, private, scoring, provenance, frozen_gbz, sample_manifest, identity


def test_private_panel_gate_freezes_provenance_and_projection(tmp_path: Path) -> None:
    canonical, graph, haplotypes, private, scoring, provenance, frozen_gbz, sample_manifest, identity = _context(tmp_path)
    MODULE.load_provenance(
        provenance,
        source_graph=graph,
        source_haplotype_manifest=haplotypes,
        source_phased_panel=private,
        cohort_id="frozen_hprc",
    )
    MODULE.validate_source_identity(
        identity,
        frozen_gbz=frozen_gbz,
        frozen_sample_manifest=sample_manifest,
        cohort_id="frozen_hprc",
    )
    result = MODULE.prepare(
        source=private,
        scoring_panel=scoring,
        output=tmp_path / "out.vcf",
        projection=tmp_path / "projection.tsv",
        gate=tmp_path / "gate.json",
        excluded_samples=FAMILY,
    )
    assert result["canonical_total"] == 1
    assert result["index_addressable"] == 1
    assert result["N_GT_missing"] == 0
    assert "AC=1;AN=2;AF=0.5" in (tmp_path / "out.vcf").read_text(encoding="utf-8")
    assert "CAND_A\tNATIVE1\t1\tindex_addressable" in (tmp_path / "projection.tsv").read_text(encoding="utf-8")


def test_private_panel_gate_rejects_phase_ambiguous_heterozygous_gt(tmp_path: Path) -> None:
    _, _, _, private, scoring, _, _, _, _ = _context(tmp_path)
    private.write_text(
        private.read_text(encoding="utf-8").replace("0|1", "0/1"),
        encoding="utf-8",
    )
    with pytest.raises(MODULE.PrivatePanelError, match="phase-ambiguous heterozygous GT"):
        MODULE.prepare(
            source=private,
            scoring_panel=scoring,
            output=tmp_path / "out.vcf",
            projection=tmp_path / "projection.tsv",
            gate=tmp_path / "gate.json",
            excluded_samples=FAMILY,
        )


def test_private_panel_canonicalizes_phase_equivalent_homozygous_slash_gt(
    tmp_path: Path,
) -> None:
    _, _, _, private, scoring, _, _, _, _ = _context(tmp_path)
    private.write_text(
        private.read_text(encoding="utf-8").replace("0|1", "1/1"),
        encoding="utf-8",
    )
    result = MODULE.prepare(
        source=private,
        scoring_panel=scoring,
        output=tmp_path / "out.vcf",
        projection=tmp_path / "projection.tsv",
        gate=tmp_path / "gate.json",
        excluded_samples=FAMILY,
    )
    assert result["N_GT_phase_equivalent_homozygous_slash"] == 1
    assert "\tGT\t1|1\n" in (tmp_path / "out.vcf").read_text(encoding="utf-8")


@pytest.mark.parametrize(
    ("genotype", "category"),
    [
        ("1/1", "phase_equivalent_homozygous"),
        ("0/1", "phase_ambiguous_heterozygous"),
        ("0/.", "partial_missing"),
        ("0/a", "malformed"),
        ("0/1/2", "malformed"),
    ],
)
def test_slash_gt_phase_audit_categories(genotype: str, category: str) -> None:
    assert MODULE.slash_gt_category(genotype) == category


def test_private_panel_derives_family_excluded_retained_haplotypes(tmp_path: Path) -> None:
    _, _, _, private, scoring, _, _, _, _ = _context(tmp_path)
    private.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG002\tPANEL1\n"
        # The first ALT loses all support after HG002 removal; the second ALT
        # stays addressable and its GT must be remapped to the compact ALT set.
        "chr1\t10\tNATIVE1\tA\tAT,ATT\t.\tPASS\t.\tGT\t1|2\t0|2\n"
        # A family-only record is removed from the native index, never from
        # the frozen scoring universe.
        "chr1\t20\tNATIVE2\tC\tCA\t.\tPASS\t.\tGT\t0|1\t0|0\n",
        encoding="utf-8",
    )
    result = MODULE.prepare(
        source=private,
        scoring_panel=scoring,
        output=tmp_path / "out.vcf",
        projection=tmp_path / "projection.tsv",
        gate=tmp_path / "gate.json",
        excluded_samples=FAMILY,
    )
    output = (tmp_path / "out.vcf").read_text(encoding="utf-8")
    assert "\tHG002\t" not in output
    assert "\tAT,ATT\t" not in output
    assert "\tATT\t.\tPASS\tAC=1;AN=2;AF=0.5\tGT\t0|1" in output
    assert "NATIVE2" not in output
    assert result["family_samples_removed"] == 1
    assert result["zero_support_records_removed"] == 1


def test_exact_source_identity_must_be_verified(tmp_path: Path) -> None:
    _, _, _, _, _, _, frozen_gbz, sample_manifest, identity = _context(tmp_path)
    payload = yaml.safe_load(identity.read_text(encoding="utf-8"))
    payload["status"] = "unverified"
    identity.write_text(yaml.safe_dump(payload, sort_keys=True), encoding="utf-8")

    with pytest.raises(MODULE.PrivatePanelError, match="must be verified"):
        MODULE.validate_source_identity(
            identity,
            frozen_gbz=frozen_gbz,
            frozen_sample_manifest=sample_manifest,
            cohort_id="frozen_hprc",
        )
