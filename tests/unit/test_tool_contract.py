from __future__ import annotations

import gzip
import json
import os
import shutil
from pathlib import Path

import pytest
import yaml

from workflow.scripts.pgbench_exec import (
    ResolvedInput,
    ToolContractError,
    ToolExecutionError,
    ToolTimeoutError,
    build_tool_environment,
    execute_tool,
    load_tool_manifest,
    validate_mode_inputs,
)
from workflow.scripts.validate_tool_output import (
    ToolOutputValidationError,
    validate_tool_output,
)

ROOT = Path(__file__).resolve().parents[2]
TOOL_SCHEMA = ROOT / "workflow" / "schemas" / "tool.schema.yaml"
EXAMPLE_PLUGIN = ROOT / "plugins" / "example_genotyper"


def _symlink_or_skip(link: Path, target: Path) -> None:
    try:
        link.symlink_to(target)
    except OSError as exc:
        if os.name == "nt" and getattr(exc, "winerror", None) == 1314:
            pytest.skip("Windows symlink privilege is not enabled")
        raise


def _write_candidate_vcf(path: Path) -> None:
    path.write_text(
        "\n".join(
            [
                "##fileformat=VCFv4.2",
                "##INFO=<ID=PANGENOME_ALLELE_ID,Number=1,Type=String,"
                'Description="Stable pangenome allele ID">',
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
                "1\t100\tCAND_alpha\tN\t<DEL>\t.\tPASS\t"
                "SVTYPE=DEL;END=199;PANGENOME_ALLELE_ID=PGSV0001",
                "1\t300\tCAND_beta\tN\t<INS>\t.\tPASS\t"
                "SVTYPE=INS;END=300;PANGENOME_ALLELE_ID=PGSV0002",
                "",
            ]
        ),
        encoding="utf-8",
    )


def _inputs(tmp_path: Path) -> dict[str, Path]:
    candidate = tmp_path / "challenge.vcf"
    reference = tmp_path / "hs37d5.fa"
    pangenome_manifest = tmp_path / "pangenome-manifest.yaml"
    pangenome_panel = tmp_path / "panel.vcf"
    alignment = tmp_path / "shared.bam"
    _write_candidate_vcf(candidate)
    reference.write_text(">1\nN\n", encoding="utf-8")
    pangenome_manifest.write_text(
        "schema_version: 1\npangenome_id: synthetic\n", encoding="utf-8"
    )
    shutil.copyfile(candidate, pangenome_panel)
    alignment.write_bytes(b"BAM\x01synthetic")
    return {
        "shared_alignment": alignment,
        "reference": reference,
        "pangenome_manifest": pangenome_manifest,
        "pangenome_panel": pangenome_panel,
        "candidate_panel": candidate,
    }


def _execute_example(
    tmp_path: Path,
    *,
    output_name: str,
    supplied_inputs: dict[str, Path] | None = None,
):
    inputs = supplied_inputs if supplied_inputs is not None else _inputs(tmp_path)
    execution_root = tmp_path / output_name
    return execute_tool(
        tool_manifest_path=EXAMPLE_PLUGIN / "tool.yaml",
        schema_path=TOOL_SCHEMA,
        mode="caller_only_shared_alignment",
        run_id=f"run_{output_name}",
        sample_id="HG002",
        supplied_inputs=inputs,
        output_dir=execution_root / "tool",
        resolved_inputs_path=execution_root / "meta" / "resolved_inputs.json",
        attempt_record_path=execution_root / "attempts" / "attempt.json",
        log_path=execution_root / "logs" / "tool.log",
        threads=2,
        memory_mb=1024,
        alignment_kind="bam",
    )


def test_example_tool_yaml_validates_with_runtime_semantics() -> None:
    manifest = load_tool_manifest(EXAMPLE_PLUGIN / "tool.yaml", TOOL_SCHEMA)
    assert manifest["id"] == "example_genotyper"
    assert set(manifest["supported_modes"]) == {
        "caller_only_shared_alignment",
        "end_to_end_from_reads",
    }
    assert manifest["outputs"]["candidate_output_contract"] == "all_sites"


def test_mode_contract_rejects_missing_and_forbidden_inputs(
    tmp_path: Path,
) -> None:
    manifest = load_tool_manifest(EXAMPLE_PLUGIN / "tool.yaml", TOOL_SCHEMA)
    inputs = _inputs(tmp_path)

    missing = dict(inputs)
    missing.pop("candidate_panel")
    with pytest.raises(ToolContractError, match="missing required input"):
        validate_mode_inputs(
            manifest,
            "caller_only_shared_alignment",
            missing,
            alignment_kind="bam",
        )

    forbidden = dict(inputs)
    forbidden_fastq = tmp_path / "forbidden.fastq"
    forbidden_fastq.write_text("@r1\nN\n+\n!\n", encoding="utf-8")
    forbidden["canonical_fastq"] = forbidden_fastq
    with pytest.raises(ToolContractError, match="forbidden input"):
        validate_mode_inputs(
            manifest,
            "caller_only_shared_alignment",
            forbidden,
            alignment_kind="bam",
        )


def test_mode_contract_exposes_only_required_or_optional_inputs(
    tmp_path: Path,
) -> None:
    manifest = load_tool_manifest(EXAMPLE_PLUGIN / "tool.yaml", TOOL_SCHEMA)
    inputs = _inputs(tmp_path)
    graph_assets = tmp_path / "graph"
    graph_assets.mkdir()
    inputs["graph_assets"] = graph_assets
    with pytest.raises(ToolContractError, match="not whitelisted"):
        validate_mode_inputs(
            manifest,
            "caller_only_shared_alignment",
            inputs,
            alignment_kind="bam",
        )


def test_runner_executes_example_and_records_resolved_inputs(
    tmp_path: Path,
) -> None:
    inputs = _inputs(tmp_path)
    result = _execute_example(tmp_path, output_name="first", supplied_inputs=inputs)

    output_vcf = Path(result.output.path)
    assert output_vcf.is_file()
    with gzip.open(output_vcf, "rt", encoding="utf-8") as handle:
        records = [
            line.rstrip("\n").split("\t")
            for line in handle
            if line and not line.startswith("#")
        ]
    assert [record[2] for record in records] == ["CAND_alpha", "CAND_beta"]
    assert all(record[8] == "GT" for record in records)
    assert result.output.candidate_count == 2
    assert result.output.represented_candidate_count == 2

    resolved = json.loads(Path(result.resolved_inputs_path).read_text(encoding="utf-8"))
    assert resolved["mode"] == "caller_only_shared_alignment"
    assert resolved["billable_stages"] == ["genotype", "postprocess"]
    assert {item["name"] for item in resolved["inputs"]} == set(inputs)
    assert all(len(item["sha256"]) == 64 for item in resolved["inputs"])
    assert all(Path(item["path"]).is_absolute() for item in resolved["inputs"])
    assert all(item["read_only"] is True for item in resolved["inputs"])

    attempt = json.loads(Path(result.attempt_record_path).read_text(encoding="utf-8"))
    assert attempt["status"] == "success"
    assert attempt["exit_code"] == 0
    assert attempt["execution_purpose"] == "development_only"
    assert attempt["isolation_status"] == "development_only_contract_enforcement"
    assert attempt["formal_score_eligible"] is False
    assert attempt["output"]["sha256"] == result.output.sha256
    assert attempt["output"]["full_compression_validation_status"] == (
        "infrastructure_required"
    )
    assert Path(result.log_path).is_file()
    assert Path(result.attempt_archive_path).is_file()
    assert Path(result.log_archive_path).is_file()
    assert Path(result.attempt_archive_path) != Path(result.attempt_record_path)
    assert Path(result.log_archive_path) != Path(result.log_path)


def test_example_genotyper_output_is_byte_deterministic(tmp_path: Path) -> None:
    inputs = _inputs(tmp_path)
    first = _execute_example(tmp_path, output_name="first", supplied_inputs=inputs)
    second = _execute_example(tmp_path, output_name="second", supplied_inputs=inputs)
    assert Path(first.output.path).read_bytes() == Path(second.output.path).read_bytes()


def test_parent_pgbench_environment_is_scrubbed(tmp_path: Path) -> None:
    manifest = load_tool_manifest(EXAMPLE_PLUGIN / "tool.yaml", TOOL_SCHEMA)
    resolved = (
        ResolvedInput(
            name="candidate_panel",
            contract_name="candidate_panel",
            path="/inputs/challenge.vcf",
            sha256="a" * 64,
            size_bytes=10,
            mtime_ns=1,
            path_type="file",
            read_only=True,
            environment_variable="PGBENCH_CANDIDATE_VCF",
        ),
    )
    environment = build_tool_environment(
        base_environment={
            "PATH": "/usr/bin",
            "PGBENCH_TRUTH_VCF": "/secret/truth.vcf",
            "PGBENCH_INPUT_FASTQ": "/stale/reads.fastq",
        },
        manifest=manifest,
        run_id="run_env",
        sample_id="HG002",
        resolved_inputs=resolved,
        resolved_inputs_path=Path("/meta/resolved.json"),
        output_dir=Path("/output"),
        output_vcf=Path("/output/raw/calls.vcf.gz"),
        threads=1,
        memory_mb=512,
        alignment_kind=None,
        attempt_work_dir=tmp_path / "attempt",
    )
    assert environment["PATH"] == "/usr/bin"
    assert "PGBENCH_TRUTH_VCF" not in environment
    assert "PGBENCH_INPUT_FASTQ" not in environment
    assert environment["PGBENCH_CANDIDATE_VCF"] == "/inputs/challenge.vcf"
    assert {key for key in environment if key.startswith("PGBENCH_")} <= {
        "PGBENCH_RUN_ID",
        "PGBENCH_SAMPLE_ID",
        "PGBENCH_RESOLVED_INPUTS",
        "PGBENCH_OUTPUT_DIR",
        "PGBENCH_OUTPUT_VCF",
        "PGBENCH_THREADS",
        "PGBENCH_MEMORY_MB",
        "PGBENCH_ALLELE_NAMESPACE",
        "PGBENCH_CANDIDATE_VCF",
    }


def test_preexisting_vcf_is_rejected_before_current_execution(
    tmp_path: Path,
) -> None:
    plugin = tmp_path / "plugin"
    shutil.copytree(EXAMPLE_PLUGIN, plugin)
    no_output_runner = plugin / "no_output.py"
    no_output_runner.write_text("print('runner executed')\n", encoding="utf-8")
    manifest = yaml.safe_load((plugin / "tool.yaml").read_text(encoding="utf-8"))
    manifest["execution"]["runner"] = "no_output.py"
    (plugin / "tool.yaml").write_text(
        yaml.safe_dump(manifest, sort_keys=False), encoding="utf-8"
    )

    inputs = _inputs(tmp_path)
    output_dir = tmp_path / "execution" / "tool"
    stale_output = output_dir / "raw" / "calls.vcf.gz"
    stale_output.parent.mkdir(parents=True)
    stale_output.write_bytes(b"precomputed result")
    attempt_path = tmp_path / "execution" / "attempt.json"

    with pytest.raises(ToolContractError, match="already exists"):
        execute_tool(
            tool_manifest_path=plugin / "tool.yaml",
            schema_path=TOOL_SCHEMA,
            mode="caller_only_shared_alignment",
            run_id="run_stale",
            sample_id="HG002",
            supplied_inputs=inputs,
            output_dir=output_dir,
            resolved_inputs_path=tmp_path / "execution" / "resolved.json",
            attempt_record_path=attempt_path,
            log_path=tmp_path / "execution" / "tool.log",
            threads=1,
            memory_mb=512,
            alignment_kind="bam",
        )

    assert stale_output.read_bytes() == b"precomputed result"
    assert not attempt_path.exists()


def test_verified_previous_vcf_and_index_are_archived_on_rerun(
    tmp_path: Path,
) -> None:
    inputs = _inputs(tmp_path)
    first = _execute_example(
        tmp_path,
        output_name="rerun",
        supplied_inputs=inputs,
    )
    final_output = Path(first.output.path)
    previous_bytes = final_output.read_bytes()
    previous_index = Path(f"{final_output}.tbi")
    previous_index.write_bytes(b"synthetic index")
    previous_stat = final_output.stat()
    os.utime(
        final_output,
        ns=(previous_stat.st_atime_ns, previous_stat.st_mtime_ns + 1_000_000),
    )

    second = _execute_example(
        tmp_path,
        output_name="rerun",
        supplied_inputs=inputs,
    )

    archive = final_output.parents[1] / "meta" / "previous_outputs" / second.attempt_id
    archived_vcf = archive / "raw" / "calls.vcf.gz"
    archived_index = archive / "raw" / "calls.vcf.gz.tbi"
    assert archived_vcf.read_bytes() == previous_bytes
    assert archived_index.read_bytes() == b"synthetic index"
    assert final_output.is_file()
    assert final_output.read_bytes() == previous_bytes
    assert not previous_index.exists()

    attempt = json.loads(Path(second.attempt_record_path).read_text(encoding="utf-8"))
    archived = attempt["previous_output_archive"]
    assert archived["path"] == str(archive)
    assert archived["source_attempt_id"] == first.attempt_id
    assert archived["output_sha256"] == first.output.sha256
    assert archived["files"] == [
        "raw/calls.vcf.gz",
        "raw/calls.vcf.gz.tbi",
    ]


def test_verified_previous_output_survives_cross_platform_repo_move(
    tmp_path: Path,
) -> None:
    inputs = _inputs(tmp_path)
    first = _execute_example(
        tmp_path,
        output_name="relocated",
        supplied_inputs=inputs,
    )
    attempt_path = Path(first.attempt_record_path)
    previous = json.loads(attempt_path.read_text(encoding="utf-8"))
    old_root = "/legacy/benchmark/results/HG002/tool"
    old_output = f"{old_root}/raw/calls.vcf.gz"
    previous["output_dir"] = old_root
    previous["expected_output_vcf"] = old_output
    previous["output"]["path"] = old_output
    attempt_path.write_text(
        json.dumps(previous, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )

    second = _execute_example(
        tmp_path,
        output_name="relocated",
        supplied_inputs=inputs,
    )

    assert Path(second.output.path).is_file()
    assert second.output.sha256 == first.output.sha256


def test_failed_rerun_keeps_previous_vcf_only_in_recoverable_archive(
    tmp_path: Path,
) -> None:
    plugin = tmp_path / "plugin"
    shutil.copytree(EXAMPLE_PLUGIN, plugin)
    inputs = _inputs(tmp_path)
    execution = tmp_path / "execution"
    output_dir = execution / "tool"
    attempt_path = execution / "attempt.json"
    common = {
        "tool_manifest_path": plugin / "tool.yaml",
        "schema_path": TOOL_SCHEMA,
        "mode": "caller_only_shared_alignment",
        "sample_id": "HG002",
        "supplied_inputs": inputs,
        "output_dir": output_dir,
        "resolved_inputs_path": execution / "resolved.json",
        "attempt_record_path": attempt_path,
        "log_path": execution / "tool.log",
        "threads": 1,
        "memory_mb": 512,
        "alignment_kind": "bam",
    }
    first = execute_tool(run_id="first", **common)
    final_output = Path(first.output.path)
    previous_bytes = final_output.read_bytes()

    no_output_runner = plugin / "no_output.py"
    no_output_runner.write_text("print('runner executed')\n", encoding="utf-8")
    manifest = yaml.safe_load((plugin / "tool.yaml").read_text(encoding="utf-8"))
    manifest["execution"]["runner"] = "no_output.py"
    (plugin / "tool.yaml").write_text(
        yaml.safe_dump(manifest, sort_keys=False), encoding="utf-8"
    )

    with pytest.raises(ToolExecutionError, match="output validation failed"):
        execute_tool(run_id="failed_rerun", **common)

    failed_attempt = json.loads(attempt_path.read_text(encoding="utf-8"))
    archive = Path(failed_attempt["previous_output_archive"]["path"])
    assert failed_attempt["status"] == "failed"
    assert failed_attempt["output"] is None
    assert not final_output.exists()
    assert (archive / "raw" / "calls.vcf.gz").read_bytes() == previous_bytes


@pytest.mark.parametrize("existing_kind", ["symlink", "directory"])
def test_preexisting_non_regular_vcf_is_rejected_without_archival(
    tmp_path: Path,
    existing_kind: str,
) -> None:
    output_dir = tmp_path / "execution" / "tool"
    final_output = output_dir / "raw" / "calls.vcf.gz"
    final_output.parent.mkdir(parents=True)
    if existing_kind == "symlink":
        outside = tmp_path / "outside.vcf.gz"
        outside.write_bytes(b"outside")
        _symlink_or_skip(final_output, outside)
    else:
        final_output.mkdir()

    with pytest.raises(ToolContractError, match="symlink|regular file"):
        execute_tool(
            tool_manifest_path=EXAMPLE_PLUGIN / "tool.yaml",
            schema_path=TOOL_SCHEMA,
            mode="caller_only_shared_alignment",
            run_id="unsafe_existing_output",
            sample_id="HG002",
            supplied_inputs=_inputs(tmp_path),
            output_dir=output_dir,
            resolved_inputs_path=tmp_path / "execution" / "resolved.json",
            attempt_record_path=tmp_path / "execution" / "attempt.json",
            log_path=tmp_path / "execution" / "tool.log",
            threads=1,
            memory_mb=512,
            alignment_kind="bam",
        )

    if existing_kind == "symlink":
        assert final_output.is_symlink()
    else:
        assert final_output.is_dir()


def test_previous_output_index_symlink_is_rejected_before_vcf_moves(
    tmp_path: Path,
) -> None:
    inputs = _inputs(tmp_path)
    first = _execute_example(
        tmp_path,
        output_name="unsafe_index",
        supplied_inputs=inputs,
    )
    final_output = Path(first.output.path)
    previous_bytes = final_output.read_bytes()
    outside = tmp_path / "outside.tbi"
    outside.write_bytes(b"outside index")
    index = Path(f"{final_output}.tbi")
    _symlink_or_skip(index, outside)

    with pytest.raises(ToolContractError, match="must not be a symlink"):
        _execute_example(
            tmp_path,
            output_name="unsafe_index",
            supplied_inputs=inputs,
        )

    assert final_output.read_bytes() == previous_bytes
    assert index.is_symlink()
    previous_attempt = json.loads(
        Path(first.attempt_record_path).read_text(encoding="utf-8")
    )
    assert previous_attempt["attempt_id"] == first.attempt_id


def test_exit_zero_without_new_vcf_records_failed_attempt(tmp_path: Path) -> None:
    plugin = tmp_path / "plugin"
    shutil.copytree(EXAMPLE_PLUGIN, plugin)
    no_output_runner = plugin / "no_output.py"
    no_output_runner.write_text("print('runner executed')\n", encoding="utf-8")
    manifest = yaml.safe_load((plugin / "tool.yaml").read_text(encoding="utf-8"))
    manifest["execution"]["runner"] = "no_output.py"
    (plugin / "tool.yaml").write_text(
        yaml.safe_dump(manifest, sort_keys=False), encoding="utf-8"
    )
    attempt_path = tmp_path / "execution" / "attempt.json"
    log_path = tmp_path / "execution" / "tool.log"
    with pytest.raises(ToolExecutionError, match="output validation failed"):
        execute_tool(
            tool_manifest_path=plugin / "tool.yaml",
            schema_path=TOOL_SCHEMA,
            mode="caller_only_shared_alignment",
            run_id="run_missing_output",
            sample_id="HG002",
            supplied_inputs=_inputs(tmp_path),
            output_dir=tmp_path / "execution" / "tool",
            resolved_inputs_path=tmp_path / "execution" / "resolved.json",
            attempt_record_path=attempt_path,
            log_path=log_path,
            threads=1,
            memory_mb=512,
            alignment_kind="bam",
        )
    attempt = json.loads(attempt_path.read_text(encoding="utf-8"))
    assert attempt["status"] == "failed"
    assert attempt["exit_code"] == 0
    assert "runner executed" in log_path.read_text(encoding="utf-8")


def test_manifest_rejects_output_path_escape(tmp_path: Path) -> None:
    plugin = tmp_path / "plugin"
    shutil.copytree(EXAMPLE_PLUGIN, plugin)
    manifest = yaml.safe_load((plugin / "tool.yaml").read_text(encoding="utf-8"))
    manifest["outputs"]["vcf"] = "../submitted.vcf"
    (plugin / "tool.yaml").write_text(
        yaml.safe_dump(manifest, sort_keys=False), encoding="utf-8"
    )
    with pytest.raises(ToolContractError, match="relative path"):
        load_tool_manifest(plugin / "tool.yaml", TOOL_SCHEMA)


def test_manifest_rejects_runner_outside_plugin_root(tmp_path: Path) -> None:
    plugin = tmp_path / "plugin"
    shutil.copytree(EXAMPLE_PLUGIN, plugin)
    outside = tmp_path / "outside.py"
    outside.write_text("raise SystemExit(0)\n", encoding="utf-8")
    manifest = yaml.safe_load((plugin / "tool.yaml").read_text(encoding="utf-8"))
    manifest["execution"]["runner"] = "../outside.py"
    (plugin / "tool.yaml").write_text(
        yaml.safe_dump(manifest, sort_keys=False), encoding="utf-8"
    )
    with pytest.raises(ToolContractError, match="plugin root"):
        load_tool_manifest(plugin / "tool.yaml", TOOL_SCHEMA)


def test_formal_execution_refuses_unsandboxed_external_tool(
    tmp_path: Path,
) -> None:
    with pytest.raises(ToolContractError, match="refuses sandbox_backend=none"):
        execute_tool(
            tool_manifest_path=EXAMPLE_PLUGIN / "tool.yaml",
            schema_path=TOOL_SCHEMA,
            mode="caller_only_shared_alignment",
            run_id="formal_without_sandbox",
            sample_id="HG002",
            supplied_inputs=_inputs(tmp_path),
            output_dir=tmp_path / "formal" / "tool",
            resolved_inputs_path=tmp_path / "formal" / "resolved.json",
            attempt_record_path=tmp_path / "formal" / "attempt.json",
            log_path=tmp_path / "formal" / "tool.log",
            threads=1,
            memory_mb=512,
            alignment_kind="bam",
            execution_purpose="formal",
        )


def test_timeout_terminates_runner_and_records_failure(tmp_path: Path) -> None:
    plugin = tmp_path / "plugin"
    shutil.copytree(EXAMPLE_PLUGIN, plugin)
    slow_runner = plugin / "slow.py"
    slow_runner.write_text(
        "import subprocess, sys, time\n"
        "subprocess.Popen([sys.executable, '-c', 'import time; time.sleep(60)'])\n"
        "time.sleep(60)\n",
        encoding="utf-8",
    )
    manifest = yaml.safe_load((plugin / "tool.yaml").read_text(encoding="utf-8"))
    manifest["execution"]["runner"] = "slow.py"
    (plugin / "tool.yaml").write_text(
        yaml.safe_dump(manifest, sort_keys=False), encoding="utf-8"
    )
    attempt_path = tmp_path / "timeout" / "attempt.json"
    with pytest.raises(ToolTimeoutError, match="exceeded timeout"):
        execute_tool(
            tool_manifest_path=plugin / "tool.yaml",
            schema_path=TOOL_SCHEMA,
            mode="caller_only_shared_alignment",
            run_id="run_timeout",
            sample_id="HG002",
            supplied_inputs=_inputs(tmp_path),
            output_dir=tmp_path / "timeout" / "tool",
            resolved_inputs_path=tmp_path / "timeout" / "resolved.json",
            attempt_record_path=attempt_path,
            log_path=tmp_path / "timeout" / "tool.log",
            threads=1,
            memory_mb=512,
            alignment_kind="bam",
            timeout_seconds=1,
        )
    attempt = json.loads(attempt_path.read_text(encoding="utf-8"))
    assert attempt["status"] == "failed"
    assert "ToolTimeoutError" in attempt["error"]
    assert attempt["timed_out"] is True
    if os.name == "nt":
        assert attempt["exit_code"] is not None
        assert attempt["termination_signal"] is None
    else:
        assert attempt["exit_code"] in {-15, -9}
        assert attempt["termination_signal"] in {15, 9}


def test_output_validator_rejects_file_that_predates_attempt(
    tmp_path: Path,
) -> None:
    output_dir = tmp_path / "output"
    output_dir.mkdir()
    candidate = tmp_path / "candidate.vcf"
    _write_candidate_vcf(candidate)
    output = output_dir / "calls.vcf"
    output.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG002\n"
        "1\t100\tCAND_alpha\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL\tGT\t0/1\n"
        "1\t300\tCAND_beta\tN\t<INS>\t.\tPASS\tSVTYPE=INS\tGT\t0/0\n",
        encoding="utf-8",
    )
    started_after_file = output.stat().st_mtime_ns + 1
    with pytest.raises(ToolOutputValidationError, match="predates"):
        validate_tool_output(
            output_vcf=output,
            output_dir=output_dir,
            started_at_ns=started_after_file,
            sample_id="HG002",
            candidate_output_contract="all_sites",
            candidate_vcf=candidate,
        )


def test_output_validator_enforces_all_sites_candidate_coverage(
    tmp_path: Path,
) -> None:
    output_dir = tmp_path / "output"
    output_dir.mkdir()
    candidate = tmp_path / "candidate.vcf"
    _write_candidate_vcf(candidate)
    output = output_dir / "calls.vcf"
    output.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG002\n"
        "1\t100\tCAND_alpha\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL\tGT\t0/1\n",
        encoding="utf-8",
    )
    with pytest.raises(ToolOutputValidationError, match="missing 1 candidate"):
        validate_tool_output(
            output_vcf=output,
            output_dir=output_dir,
            started_at_ns=0,
            sample_id="HG002",
            candidate_output_contract="all_sites",
            candidate_vcf=candidate,
        )


def test_output_validator_rejects_hardlinked_result(tmp_path: Path) -> None:
    output_dir = tmp_path / "output"
    output_dir.mkdir()
    candidate = tmp_path / "candidate.vcf"
    _write_candidate_vcf(candidate)
    source = tmp_path / "precomputed.vcf"
    source.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG002\n"
        "1\t100\tCAND_alpha\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL\tGT\t0/1\n"
        "1\t300\tCAND_beta\tN\t<INS>\t.\tPASS\tSVTYPE=INS\tGT\t0/0\n",
        encoding="utf-8",
    )
    output = output_dir / "calls.vcf"
    os.link(source, output)
    with pytest.raises(ToolOutputValidationError, match="hardlink"):
        validate_tool_output(
            output_vcf=output,
            output_dir=output_dir,
            started_at_ns=0,
            sample_id="HG002",
            candidate_output_contract="all_sites",
            candidate_vcf=candidate,
        )
