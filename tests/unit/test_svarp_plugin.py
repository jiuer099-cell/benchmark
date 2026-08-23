from __future__ import annotations

import importlib.util
from pathlib import Path

import yaml


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
    parsed = yaml.safe_load(manifest)
    parameters = parsed["parameter_contract"]
    assert parameters["source_release"]["commit"] == (
        "9fc34f848bd6506b4fbe6283e44c55b0a691519d"
    )
    assert parameters["explicit"] == {
        "assembler": "wtdbg2",
        "support": 5,
        "dist_threshold": 100,
        "min_alignment_score": 5000,
        "min_precise_clipping": 0.97,
        "min_map_ratio": 0.90,
        "threads": "PGBENCH_THREADS",
    }
    assert parameters["implicit"]["read_type_option"] == (
        "unsupported_in_v1.2.0"
    )


def test_svtig_paf_adapter_emits_only_large_sequence_resolved_svs(
    tmp_path: Path,
) -> None:
    reverse_insertion = "A" * 20 + "C" * 20 + "G" * 20
    raw = tmp_path / "svtigs.paf"
    raw.write_text(
        "svtig1\t200\t0\t200\t+\tchr1\t1000\t0\t200\t200\t200\t60\t"
        "tp:A:P\tcs:Z::20-" + "A" * 60 + ":120\n"
        "svtig2\t200\t0\t200\t-\tchr1\t1000\t300\t500\t200\t200\t60\t"
        "tp:A:P\tcs:Z::30+" + reverse_insertion + ":110\n"
        "svtig3\t100\t0\t100\t+\tchr1\t1000\t600\t700\t100\t100\t60\t"
        "tp:A:P\tcs:Z::30+" + "T" * 10 + ":60\n"
        "lowmapq\t200\t0\t200\t+\tchr1\t1000\t700\t900\t200\t200\t10\t"
        "tp:A:P\tcs:Z::30+" + "G" * 60 + ":110\n",
        encoding="utf-8",
    )
    reference = tmp_path / "reference.fa"
    reference.write_bytes(b">chr1\n" + b"A" * 1000 + b"\n")
    fai = tmp_path / "reference.fa.fai"
    fai.write_text("chr1\t1000\t6\t1000\t1001\n", encoding="utf-8")
    destination = tmp_path / "calls.vcf"

    count = MODULE.write_discovery_vcf(
        raw,
        destination,
        sample="HG002",
        reference=reference,
        reference_fai=fai,
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
    assert records[1][3:5] == ["A", "A" + reverse_insertion]
    assert records[0][9] == "./."


def test_converter_handles_boundaries_and_rejects_ambiguous_primary_rows(
    tmp_path: Path,
) -> None:
    insertion = "ACGT" * 13
    deletion = "A" * 50
    raw = tmp_path / "svtigs.paf"
    raw.write_text(
        "left_ins\t100\t0\t100\t+\tchr1\t300\t0\t100\t100\t100\t60\t"
        f"tp:A:P\tcs:Z:+{insertion}:48\n"
        "left_del\t100\t0\t100\t+\tchr1\t300\t0\t100\t100\t100\t60\t"
        f"tp:A:P\tcs:Z:-{deletion}:50\n"
        "ambiguous\t100\t0\t100\t+\tchr1\t300\t100\t200\t100\t100\t60\t"
        f"tp:A:P\tcs:Z::20+{insertion}:28\n"
        "ambiguous\t100\t0\t100\t+\tchr1\t300\t150\t250\t100\t100\t60\t"
        f"tp:A:P\tcs:Z::20+{insertion}:28\n",
        encoding="utf-8",
    )
    reference = tmp_path / "reference.fa"
    reference.write_bytes(b">chr1\n" + b"A" * 300 + b"\n")
    fai = tmp_path / "reference.fa.fai"
    fai.write_text("chr1\t300\t6\t300\t301\n", encoding="utf-8")
    destination = tmp_path / "calls.vcf"

    assert (
        MODULE.write_discovery_vcf(
            raw,
            destination,
            sample="HG002",
            reference=reference,
            reference_fai=fai,
        )
        == 2
    )
    records = [
        line.split("\t")
        for line in destination.read_text(encoding="utf-8").splitlines()
        if line and not line.startswith("#")
    ]
    assert records[0][1:5] == ["1", "SVARP_000000001", "A", insertion + "A"]
    assert records[1][1] == "1"
    assert records[1][3] == "A" * 51
    assert records[1][4] == "A"
    assert "END=50" in records[1][7]


def test_runner_forbids_unsealed_cross_attempt_checkpoint_reuse() -> None:
    source = RUNNER.read_text(encoding="utf-8")
    assert "find_reusable_work_dir" not in source
    assert "reusing completed SVarp work checkpoint" not in source
