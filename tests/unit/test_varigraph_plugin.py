from __future__ import annotations

import gzip
import importlib.util
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[2]
RUNNER = ROOT / "plugins" / "varigraph" / "run.py"
SPEC = importlib.util.spec_from_file_location("varigraph_runner", RUNNER)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def test_manifest_freezes_varigraph_and_paired_reads() -> None:
    manifest = (ROOT / "plugins" / "varigraph" / "tool.yaml").read_text(
        encoding="utf-8"
    )
    environment = (
        ROOT / "plugins" / "varigraph" / "envs" / "environment.yaml"
    ).read_text(encoding="utf-8")
    assert "version: 1.0.8" in manifest
    assert "short_fastq_r1" in manifest
    assert "short_fastq_r2" in manifest
    assert "candidate_output_contract: all_sites" in manifest
    assert "- varigraph=1.0.8" in environment


def _panel(path: Path, sample: str = "PANEL1") -> Path:
    path.write_text(
        "##fileformat=VCFv4.2\n"
        f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}\n"
        "chr1\t10\tv1\tA\tAT\t.\tPASS\t.\tGT\t0|1\n",
        encoding="utf-8",
    )
    return path


def test_panel_header_rejects_benchmark_alias(tmp_path: Path) -> None:
    with pytest.raises(RuntimeError, match="benchmark aliases"):
        MODULE.validate_panel_header(_panel(tmp_path / "panel.vcf", "HG002"))


def test_projection_is_all_sites_and_preserves_no_call(tmp_path: Path) -> None:
    candidates = tmp_path / "candidate.vcf"
    candidates.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chr1\t10\tCAND_1\tA\tAT\t.\tPASS\tPANGENOME_ALLELE_ID=PGSV_1\n"
        "chr1\t20\tCAND_2\tA\tAG\t.\tPASS\tPANGENOME_ALLELE_ID=PGSV_2\n",
        encoding="utf-8",
    )
    generated = tmp_path / "HG002.varigraph.vcf.gz"
    with gzip.open(generated, "wt", encoding="utf-8") as output:
        output.write(
            "##fileformat=VCFv4.2\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG002\n"
            "chr1\t10\tpanel1\tA\tAT\t.\tPASS\t"
            "PANGENOME_ALLELE_ID=PGSV_1\tGT:GQ\t0/1:50\n"
            "chr1\t30\tpanel3\tA\tAC\t.\tPASS\t"
            "PANGENOME_ALLELE_ID=PGSV_3\tGT\t1/1\n"
        )
    destination = tmp_path / "calls.vcf"

    counts = MODULE.project_to_candidates(
        generated, candidates, destination, "HG002"
    )

    records = [
        line.split("\t")
        for line in destination.read_text(encoding="utf-8").splitlines()
        if line and not line.startswith("#")
    ]
    assert counts == (1, 1, 1)
    assert [record[2] for record in records] == ["CAND_1", "CAND_2"]
    assert records[0][9] == "0/1"
    assert records[1][9] == "./."
