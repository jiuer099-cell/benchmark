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
