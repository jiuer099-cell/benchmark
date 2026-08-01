from __future__ import annotations

import hashlib
import sys
from pathlib import Path

import pytest

SCRIPTS = Path(__file__).resolve().parents[2] / "workflow" / "scripts"
sys.path.insert(0, str(SCRIPTS))

from prepare_sv_truth import (  # noqa: E402
    TruthPreparationError,
    load_filter_contract,
    materialize_plain_truth,
    verify_frozen_inputs,
    verify_universe_contract,
)
from sv_matching import load_evaluator_profile  # noqa: E402


CONTRACT = {
    "pass_only": True,
    "inside_benchmark_bed": True,
    "region_policy": "fully_contained",
    "minimum_sv_size": 50,
    "maximum_sv_size": 10000,
    "allowed_svtypes": ["DEL", "INS"],
    "reject_multiallelic_sv": True,
    "multiallelic_policy": "exclude",
}


def test_materializes_one_frozen_sv_only_universe(tmp_path: Path) -> None:
    source = tmp_path / "joint.vcf"
    source.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG002\n"
        "chr1\t100\t.\tA\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=200;SVLEN=-100\tGT\t0/1\n"
        "chr1\t210\tunflagged\tA\t<DEL>\t.\t.\tSVTYPE=DEL;END=310;SVLEN=-100\tGT\t0/1\n"
        "chr1\t250\tbad_filter\tA\t<DEL>\t.\tq10\tSVTYPE=DEL;END=350;SVLEN=-100\tGT\t0/1\n"
        "chr1\t400\tsmall\tA\t<INS>\t.\tPASS\tSVTYPE=INS;SVLEN=20\tGT\t0/1\n"
        "chr1\t450\ttoo_large\tA\t<INS>\t.\tPASS\t"
        "SVTYPE=INS;SVLEN=10001\tGT\t0/1\n"
        "chr1\t500\tdup\tA\t<DUP>\t.\tPASS\tSVTYPE=DUP;END=600;SVLEN=100\tGT\t0/1\n"
        "chr1\t950\tcrosses\tA\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=1050;SVLEN=-100\tGT\t0/1\n"
        "chr1\t2000\toutside\tA\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=2100;SVLEN=-100\tGT\t0/1\n",
        encoding="utf-8",
    )
    regions = tmp_path / "benchmark.bed"
    regions.write_text("chr1\t0\t1000\n", encoding="utf-8")
    output = tmp_path / "truth.vcf"
    counts = materialize_plain_truth(
        source_vcf=source,
        benchmark_bed=regions,
        output_vcf=output,
        contract=CONTRACT,
    )
    records = [
        line.split("\t")
        for line in output.read_text(encoding="utf-8").splitlines()
        if line and not line.startswith("#")
    ]
    assert len(records) == 2
    assert records[0][2].startswith("TRUTH_")
    assert counts == {
        "source_records": 8,
        "excluded_non_pass": 1,
        "excluded_outside_bed": 2,
        "excluded_below_minimum_size": 1,
        "excluded_above_maximum_size": 1,
        "excluded_svtype": 1,
        "excluded_multiallelic": 0,
        "excluded_duplicate": 0,
        "eligible_records": 2,
    }


def test_excludes_eligible_multiallelic_sv_under_frozen_policy(
    tmp_path: Path,
) -> None:
    source = tmp_path / "multi.vcf"
    source.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chr1\t100\tmulti\tA\t<DEL>,<INS>\t.\tPASS\t"
        "SVTYPE=DEL;END=200;SVLEN=-100,100\n"
        "chr1\t300\tbiallelic\tA\t<DEL>\t.\tPASS\t"
        "SVTYPE=DEL;END=400;SVLEN=-100\n",
        encoding="utf-8",
    )
    regions = tmp_path / "benchmark.bed"
    regions.write_text("chr1\t0\t1000\n", encoding="utf-8")
    output = tmp_path / "truth.vcf"
    counts = materialize_plain_truth(
        source_vcf=source,
        benchmark_bed=regions,
        output_vcf=output,
        contract=CONTRACT,
    )
    assert counts["excluded_multiallelic"] == 1
    assert counts["eligible_records"] == 1
    assert "\tbiallelic\t" in output.read_text(encoding="utf-8")


def test_deduplicates_identical_truth_events(tmp_path: Path) -> None:
    source = tmp_path / "duplicates.vcf"
    record = (
        "chr1\t100\t.\tA\t<DEL>\t.\tPASS\t"
        "SVTYPE=DEL;END=200;SVLEN=-100\n"
    )
    source.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        + record
        + record,
        encoding="utf-8",
    )
    regions = tmp_path / "benchmark.bed"
    regions.write_text("chr1\t0\t1000\n", encoding="utf-8")
    output = tmp_path / "truth.vcf"
    counts = materialize_plain_truth(
        source_vcf=source,
        benchmark_bed=regions,
        output_vcf=output,
        contract=CONTRACT,
    )
    assert counts["eligible_records"] == 1
    assert counts["excluded_duplicate"] == 1


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_formal_catalog_requires_all_predeclared_truth_hashes(
    tmp_path: Path,
) -> None:
    catalog = tmp_path / "truthsets.yaml"
    catalog.write_text(
        "truthsets:\n"
        "  truth:\n"
        "    required_sha256_before_execution: true\n"
        "    files:\n"
        "      vcf: {sha256: null}\n"
        f"      vcf_index: {{sha256: {'1' * 64}}}\n"
        f"      benchmark_bed: {{sha256: {'2' * 64}}}\n"
        "    filters:\n"
        "      pass_only: true\n"
        "      inside_benchmark_bed: true\n"
        "      region_policy: fully_contained\n"
        "      minimum_sv_size: 50\n"
        "      allowed_svtypes: [DEL, INS]\n"
        "      reject_multiallelic_sv: true\n"
        "      multiallelic_policy: exclude\n",
        encoding="utf-8",
    )
    with pytest.raises(TruthPreparationError, match="source_vcf"):
        load_filter_contract(catalog, "truth")


def test_frozen_truth_assets_fail_closed_on_checksum_mismatch(
    tmp_path: Path,
) -> None:
    source = tmp_path / "truth.vcf.gz"
    index = tmp_path / "truth.vcf.gz.tbi"
    regions = tmp_path / "truth.bed"
    for path, content in (
        (source, b"vcf"),
        (index, b"index"),
        (regions, b"chr1\t0\t10\n"),
    ):
        path.write_bytes(content)
    contract = {
        **CONTRACT,
        "required_sha256_before_execution": True,
        "expected_sha256": {
            "source_vcf": _sha(source),
            "source_index": "0" * 64,
            "benchmark_bed": _sha(regions),
        },
    }
    with pytest.raises(TruthPreparationError, match="source_index"):
        verify_frozen_inputs(
            source_vcf=source,
            source_index=index,
            benchmark_bed=regions,
            contract=contract,
        )


def test_truth_catalog_must_match_frozen_evaluator_universe() -> None:
    profile = load_evaluator_profile(
        SCRIPTS.parents[1] / "config" / "evaluator_profile.yaml"
    )
    mismatched = {**CONTRACT, "minimum_sv_size": 30}
    with pytest.raises(TruthPreparationError, match="differ"):
        verify_universe_contract(mismatched, profile)

    mismatched = {**CONTRACT, "maximum_sv_size": 10001}
    with pytest.raises(TruthPreparationError, match="differ"):
        verify_universe_contract(mismatched, profile)
