from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "plugins" / "pangenie" / "mask_partial_missing.py"
SPEC = importlib.util.spec_from_file_location("mask_partial", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def test_masks_partial_missing_without_phase_inference(tmp_path: Path) -> None:
    source = tmp_path / "source.vcf"
    source.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\n"
        "chr20\t10\tv1\tA\tAT\t.\tPASS\t.\tGT:DP\t0/.:5\t1/1:6\t./.:7\n",
        encoding="utf-8",
    )
    output = tmp_path / "output.vcf"
    audit = tmp_path / "audit.json"

    result = MODULE.transform(source, output, audit)

    assert result["partial_missing_slash_masked_to_no_call"] == 1
    assert result["homozygous_slash_canonicalized"] == 1
    assert result["fully_missing_slash_gt"] == 1
    assert "./.:5\t1|1:6\t./.:7" in output.read_text(encoding="utf-8")
    assert json.loads(audit.read_text(encoding="utf-8"))["records"] == 1


@pytest.mark.parametrize("genotype", ["0/1", "1/2"])
def test_rejects_phase_ambiguous_heterozygous_gt(genotype: str) -> None:
    with pytest.raises(MODULE.PhaseGateError, match="phase-ambiguous"):
        MODULE.normalize_gt(genotype, line_number=1)
