from __future__ import annotations

import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
RUNNER = ROOT / "plugins" / "svarp" / "run.py"
SPEC = importlib.util.spec_from_file_location("svarp_runner", RUNNER)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def test_manifest_freezes_long_read_discovery_contract() -> None:
    manifest = (ROOT / "plugins" / "svarp" / "tool.yaml").read_text(
        encoding="utf-8"
    )
    environment = (ROOT / "plugins" / "svarp" / "envs" / "environment.yaml").read_text(
        encoding="utf-8"
    )
    assert "version: 1.2.0" in manifest
    assert "reads: long_fastq" in manifest
    assert "candidate_output_contract: variant_sites" in manifest
    assert "- svarp=1.2.0" in environment
    assert "- wtdbg=2.5" in environment
    assert '"--reads"' not in RUNNER.read_text(encoding="utf-8")


def test_svtig_vcf_adapter_emits_only_large_sequence_resolved_svs(
    tmp_path: Path,
) -> None:
    raw = tmp_path / "paftools.vcf"
    raw.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chr1\t10\t.\tA\tAT\t.\tPASS\t.\n"
        "chr1\t20\t.\t" + "A" * 61 + "\tA\t.\tPASS\t.\n"
        "chr1\t100\t.\tA\t" + "A" + "C" * 60 + "\t.\tPASS\t.\n",
        encoding="utf-8",
    )
    fai = tmp_path / "reference.fa.fai"
    fai.write_text("chr1\t1000\t0\t80\t81\n", encoding="utf-8")
    destination = tmp_path / "calls.vcf"

    count = MODULE.write_discovery_vcf(
        raw, destination, sample="HG002", reference_fai=fai
    )

    records = [
        line.split("\t")
        for line in destination.read_text(encoding="utf-8").splitlines()
        if line and not line.startswith("#")
    ]
    assert count == 2
    assert [record[2] for record in records] == ["SVARP_000000001", "SVARP_000000002"]
    assert "SVTYPE=DEL" in records[0][7]
    assert "SVTYPE=INS" in records[1][7]
    assert records[0][9] == "./."
