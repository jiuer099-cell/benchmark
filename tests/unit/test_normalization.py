from __future__ import annotations

import sys
from pathlib import Path

SCRIPTS = Path(__file__).resolve().parents[2] / "workflow" / "scripts"
sys.path.insert(0, str(SCRIPTS))

from normalize_sv_vcf import normalize_vcf  # noqa: E402


def test_normalization_sorts_and_preserves_original_representation(
    tmp_path: Path,
) -> None:
    reference = tmp_path / "reference.fa"
    reference.write_text(">1\nAAAA\n>2\nAAAA\n", encoding="utf-8")
    raw = tmp_path / "raw.vcf"
    raw.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ttool_sample\n"
        "2\t9\tb\tA\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=20;SVLEN=-11\tGT\t1/1\n"
        "1\t4\ta\tA\tATTT\t.\tPASS\t.\tGT\t0/1\n",
        encoding="utf-8",
    )
    output = tmp_path / "canonical.vcf"
    assert normalize_vcf(raw, reference, output) == 2
    text = output.read_text(encoding="utf-8")
    assert text.splitlines()[-2].startswith("1\t4\tCANON_")
    assert text.splitlines()[-1].startswith("2\t9\tCANON_")
    assert "ORIG_ID=a" in text
    assert "ORIG_ALT=ATTT" in text
    assert "SVTYPE=INS" in text
    assert text.splitlines()[-3].endswith("\tHG002")
