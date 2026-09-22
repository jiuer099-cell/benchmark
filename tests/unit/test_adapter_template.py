from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[2]
TEMPLATE_PATH = ROOT / "plugins" / "adapter_template" / "run.py"
SPEC = importlib.util.spec_from_file_location("adapter_template_run", TEMPLATE_PATH)
assert SPEC is not None and SPEC.loader is not None
adapter_template = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = adapter_template
SPEC.loader.exec_module(adapter_template)


def _vcf(path: Path, records: list[str]) -> Path:
    path.write_text(
        "##fileformat=VCFv4.2\n"
        '##INFO=<ID=PANGENOME_ALLELE_ID,Number=1,Type=String,Description="x">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG002\n"
        + "".join(records),
        encoding="utf-8",
    )
    return path


def test_template_projection_writes_auditable_native_trace(tmp_path: Path) -> None:
    candidates = _vcf(
        tmp_path / "candidates.vcf",
        ["chr1\t10\tPGSV_1\tA\tT\t.\t.\tPANGENOME_ALLELE_ID=a1\tGT\t./.\n"],
    )
    native = _vcf(
        tmp_path / "native.vcf",
        ["chr1\t10\t.\tA\tT\t.\t.\tPANGENOME_ALLELE_ID=a1\tGT\t0/1\n"],
    )
    output = tmp_path / "raw" / "calls.vcf"

    assert adapter_template.project_all_sites(native, candidates, output, "HG002") == (1, 0)
    assert output.read_text(encoding="utf-8").endswith("GT\t0/1\n")
    trace = output.parent / "tool-work" / "native-to-canonical-projection.tsv"
    assert "PGSV_1\tpangenome_allele_id\tchr1\t10\tA\tT\t0/1" in trace.read_text(encoding="utf-8")


def test_template_projection_traces_native_record_outside_scoring_universe(tmp_path: Path) -> None:
    candidates = _vcf(
        tmp_path / "candidates.vcf",
        ["chr1\t10\tPGSV_1\tA\tT\t.\t.\t.\tGT\t./.\n"],
    )
    native = _vcf(
        tmp_path / "native.vcf",
        ["chr1\t11\t.\tA\tC\t.\t.\t.\tGT\t0/1\n"],
    )
    output = tmp_path / "raw" / "calls.vcf"
    assert adapter_template.project_all_sites(native, candidates, output, "HG002") == (0, 1)
    trace = output.parent / "tool-work" / "native-to-canonical-projection.tsv"
    assert "\t\toutside_canonical_universe\tchr1\t11\tA\tC\t0/1" in trace.read_text(encoding="utf-8")


def test_template_rejects_multiple_primary_read_modes(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv("PGBENCH_INPUT_FASTQ_R1", "/reads/r1.fq.gz")
    monkeypatch.setenv("PGBENCH_INPUT_LONG_READS_FASTQ", "/reads/lr.fq.gz")
    with pytest.raises(RuntimeError, match="more than one primary read evidence"):
        adapter_template.detect_input_mode()


def test_template_exposes_only_resolved_declared_inputs(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    panel = tmp_path / "panel.vcf"
    panel.write_text("x", encoding="utf-8")
    asset = tmp_path / "index"
    asset.mkdir()
    resolved = tmp_path / "resolved-inputs.json"
    resolved.write_text(
        json.dumps(
            {
                "inputs": [
                    {
                        "contract_name": "candidate_panel",
                        "environment_variable": "PGBENCH_CANDIDATE_VCF",
                    },
                    {
                        "contract_name": "adapter_asset.my_index",
                        "environment_variable": "PGBENCH_ADAPTER_ASSET_MY_INDEX",
                    },
                ]
            }
        ),
        encoding="utf-8",
    )
    monkeypatch.setenv("PGBENCH_RESOLVED_INPUTS", str(resolved))
    monkeypatch.setenv("PGBENCH_CANDIDATE_VCF", str(panel))
    monkeypatch.setenv("PGBENCH_ADAPTER_ASSET_MY_INDEX", str(asset))

    assert adapter_template.injected_inputs() == {
        "candidate_panel": panel,
        "adapter_asset.my_index": asset,
    }
