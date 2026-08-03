from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[2]
RUNNER = ROOT / "plugins" / "pangenie" / "run.py"
SPEC = importlib.util.spec_from_file_location("pangenie_adapter", RUNNER)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def write_panel(path: Path, *, sample: str = "PANEL1", gt: str = "0|1") -> Path:
    path.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
        f"{sample}\n"
        "chr1\t10\tv1\tA\tAT\t.\tPASS\tEND=10\tGT\t"
        f"{gt}\n",
        encoding="utf-8",
    )
    return path


def test_valid_sequence_resolved_phased_panel(tmp_path: Path) -> None:
    MODULE.validate_pangenie_panel(write_panel(tmp_path / "panel.vcf"))


@pytest.mark.parametrize(
    ("sample", "gt", "alt"),
    [
        ("HG002", "0|1", "AT"),
        ("PANEL1", "0/1", "AT"),
        ("PANEL1", "0|1", "<INS>"),
    ],
)
def test_invalid_pangenie_panel_is_rejected(
    tmp_path: Path,
    sample: str,
    gt: str,
    alt: str,
) -> None:
    panel = write_panel(tmp_path / "panel.vcf", sample=sample, gt=gt)
    panel.write_text(
        panel.read_text(encoding="utf-8").replace("\tAT\t", f"\t{alt}\t"),
        encoding="utf-8",
    )
    with pytest.raises(RuntimeError):
        MODULE.validate_pangenie_panel(panel)


def test_overlapping_records_are_rejected(tmp_path: Path) -> None:
    panel = write_panel(tmp_path / "panel.vcf")
    with panel.open("a", encoding="utf-8") as handle:
        handle.write("chr1\t10\tv2\tA\tAG\t.\tPASS\tEND=10\tGT\t1|0\n")
    with pytest.raises(RuntimeError, match="overlaps"):
        MODULE.validate_pangenie_panel(panel)


def test_unphased_panel_records_are_filtered_without_imputation(
    tmp_path: Path,
) -> None:
    source = write_panel(tmp_path / "source.vcf")
    with source.open("a", encoding="utf-8") as handle:
        handle.write("chr1\t20\tv2\tA\tAG\t.\tPASS\tEND=20\tGT\t./.\n")
    filtered = tmp_path / "filtered.vcf"

    kept, dropped = MODULE.filter_pangenie_panel(source, filtered)

    assert (kept, dropped) == (1, 1)
    assert "\tv1\t" in filtered.read_text(encoding="utf-8")
    assert "\tv2\t" not in filtered.read_text(encoding="utf-8")
    MODULE.validate_pangenie_panel(filtered)


def test_output_ids_are_remapped_to_blinded_candidates(tmp_path: Path) -> None:
    candidate = tmp_path / "candidate.vcf"
    candidate.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG002\n"
        "chr1\t10\tCAND_alpha\tA\tAT\t.\tPASS\t"
        "PANGENOME_ALLELE_ID=PGSV_alpha\tGT\t./.\n",
        encoding="utf-8",
    )
    generated = tmp_path / "generated.vcf"
    generated.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG002\n"
        "chr1\t10\tpanel-id\tA\tAT\t.\tPASS\t"
        "PANGENOME_ALLELE_ID=PGSV_alpha\tGT\t0/1\n",
        encoding="utf-8",
    )
    output = tmp_path / "output.vcf"
    MODULE.remap_to_candidate_ids(generated, candidate, output)
    records = [
        line.rstrip("\n").split("\t")
        for line in output.read_text(encoding="utf-8").splitlines()
        if line and not line.startswith("#")
    ]
    assert records[0][2] == "CAND_alpha"
    assert records[0][9] == "0/1"


def test_filtered_candidate_is_preserved_as_explicit_no_call(tmp_path: Path) -> None:
    candidate = tmp_path / "candidate.vcf"
    candidate.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG002\n"
        "chr1\t10\tCAND_alpha\tA\tAT\t.\tPASS\t"
        "PANGENOME_ALLELE_ID=PGSV_alpha\tGT\t./.\n"
        "chr1\t20\tCAND_beta\tA\tAG\t.\tPASS\t"
        "PANGENOME_ALLELE_ID=PGSV_beta\tGT\t./.\n",
        encoding="utf-8",
    )
    generated = tmp_path / "generated.vcf"
    generated.write_text(
        "##fileformat=VCFv4.2\n"
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG002\n"
        "chr1\t10\tpanel-id\tA\tAT\t.\tPASS\t"
        "PANGENOME_ALLELE_ID=PGSV_alpha\tGT\t0/1\n",
        encoding="utf-8",
    )
    output = tmp_path / "output.vcf"

    MODULE.remap_to_candidate_ids(generated, candidate, output)

    records = [
        line.split("\t")
        for line in output.read_text(encoding="utf-8").splitlines()
        if line and not line.startswith("#")
    ]
    assert [record[2] for record in records] == ["CAND_alpha", "CAND_beta"]
    assert records[0][9] == "0/1"
    assert records[1][9] == "./."


def test_full_panel_records_outside_blinded_universe_are_filtered(
    tmp_path: Path,
) -> None:
    candidate = tmp_path / "candidate.vcf"
    candidate.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG002\n"
        "chr1\t10\tCAND_alpha\tA\tAT\t.\tPASS\t"
        "PANGENOME_ALLELE_ID=PGSV_alpha\tGT\t./.\n",
        encoding="utf-8",
    )
    generated = tmp_path / "generated.vcf"
    generated.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG002\n"
        "chr1\t5\tpanel-extra\tA\tG\t.\tPASS\t"
        "PANGENOME_ALLELE_ID=PGSV_extra\tGT\t0/1\n"
        "chr1\t10\tpanel-id\tA\tAT\t.\tPASS\t"
        "PANGENOME_ALLELE_ID=PGSV_alpha\tGT\t1/1\n",
        encoding="utf-8",
    )
    output = tmp_path / "output.vcf"

    counts = MODULE.remap_to_candidate_ids(generated, candidate, output)

    records = [
        line.split("\t")
        for line in output.read_text(encoding="utf-8").splitlines()
        if line and not line.startswith("#")
    ]
    assert counts == (1, 1, 0)
    assert len(records) == 1
    assert records[0][2] == "CAND_alpha"
    assert records[0][9] == "1/1"
