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


def _context(tmp_path: Path) -> tuple[Path, Path, Path, Path, Path, Path]:
    canonical_source = tmp_path / "canonical.vcf"
    canonical_source.write_text("frozen source\n", encoding="utf-8")
    graph = tmp_path / "source.gfa.gz"
    graph.write_text("graph\n", encoding="utf-8")
    haplotypes = tmp_path / "haplotypes.yaml"
    haplotypes.write_text("samples: [PANEL1]\n", encoding="utf-8")
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
            "canonical_population_source_sha256": MODULE.sha256(canonical_source),
            "source_graph_sha256": MODULE.sha256(graph),
            "source_haplotype_manifest_sha256": MODULE.sha256(haplotypes),
            "reference_build": "grch38",
            "generation_method": MODULE.METHOD,
            "family_exclusion_applied_before_native_panel_generation": True,
            "excluded_samples": sorted(FAMILY),
        },
        provenance.open("w", encoding="utf-8"),
        sort_keys=True,
    )
    return canonical_source, graph, haplotypes, private, scoring, provenance


def test_private_panel_gate_freezes_provenance_and_projection(tmp_path: Path) -> None:
    canonical, graph, haplotypes, private, scoring, provenance = _context(tmp_path)
    MODULE.load_provenance(
        provenance,
        canonical_source=canonical,
        source_graph=graph,
        source_haplotype_manifest=haplotypes,
        excluded_samples=FAMILY,
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


def test_private_panel_gate_rejects_missing_or_unphased_gt(tmp_path: Path) -> None:
    _, _, _, private, scoring, _ = _context(tmp_path)
    private.write_text(
        private.read_text(encoding="utf-8").replace("0|1", "0/1"),
        encoding="utf-8",
    )
    with pytest.raises(MODULE.PrivatePanelError, match="unphased GT"):
        MODULE.prepare(
            source=private,
            scoring_panel=scoring,
            output=tmp_path / "out.vcf",
            projection=tmp_path / "projection.tsv",
            gate=tmp_path / "gate.json",
            excluded_samples=FAMILY,
        )
