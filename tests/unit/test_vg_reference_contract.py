from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest

from workflow.scripts.pgbench_exec import (
    ResolvedInput,
    ToolContractError,
    build_tool_environment,
)


ROOT = Path(__file__).resolve().parents[2]
RUNNER = ROOT / "plugins" / "vg" / "run.py"
SPEC = importlib.util.spec_from_file_location("vg_adapter", RUNNER)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def _resolved_graph(tmp_path: Path, *, write_lock: bool = True) -> ResolvedInput:
    graph_dir = tmp_path / "graph"
    graph_dir.mkdir()
    if write_lock:
        (graph_dir / "graph-assets.lock.yaml").write_text(
            "schema_version: 1\nreference_path: GRCh38\n",
            encoding="utf-8",
        )
    return ResolvedInput(
        name="graph_assets",
        contract_name="graph_assets",
        path=str(graph_dir.resolve()),
        sha256="a" * 64,
        size_bytes=1,
        mtime_ns=1,
        path_type="directory",
        read_only=True,
        environment_variable="PGBENCH_GRAPH_DIR",
    )


def _environment(tmp_path: Path, graph: ResolvedInput) -> dict[str, str]:
    return build_tool_environment(
        base_environment={
            "PATH": "/usr/bin",
            "PGBENCH_GRAPH_REFERENCE_PATH": "stale-reference",
        },
        manifest={"capabilities": {"allele_namespace": ["PGSV"]}},
        run_id="vg_reference_contract",
        sample_id="HG002",
        resolved_inputs=(graph,),
        resolved_inputs_path=tmp_path / "resolved.json",
        output_dir=tmp_path / "output",
        output_vcf=tmp_path / "output" / "raw" / "calls.vcf",
        threads=4,
        memory_mb=4096,
        alignment_kind=None,
        attempt_work_dir=tmp_path / "attempt",
    )


def test_executor_exposes_reference_from_graph_lock(tmp_path: Path) -> None:
    graph = _resolved_graph(tmp_path)
    environment = _environment(tmp_path, graph)

    assert environment["PGBENCH_GRAPH_DIR"] == graph.path
    assert environment["PGBENCH_GRAPH_REFERENCE_PATH"] == "GRCh38"


def test_executor_fails_closed_without_graph_reference_lock(tmp_path: Path) -> None:
    graph = _resolved_graph(tmp_path, write_lock=False)

    with pytest.raises(ToolContractError, match="frozen reference selector"):
        _environment(tmp_path, graph)


def _set_runner_environment(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    *,
    reference_path: str | None,
) -> None:
    graph_dir = tmp_path / "runner-graph"
    graph_dir.mkdir()
    for name in (
        "graph.gbz",
        "graph.shortread.withzip.min",
        "graph.shortread.zipcodes",
        "graph.dist",
    ):
        (graph_dir / name).write_bytes(name.encode())
    reads_r1 = tmp_path / "reads.R1.fastq"
    reads_r2 = tmp_path / "reads.R2.fastq"
    reads_r1.write_text("@read/1\nAC\n+\n!!\n", encoding="utf-8")
    reads_r2.write_text("@read/2\nGT\n+\n!!\n", encoding="utf-8")
    candidates = tmp_path / "candidates.vcf"
    candidates.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chr1\t10\tCAND_1\tA\tAT\t.\tPASS\tPANGENOME_ALLELE_ID=PGSV_1\n",
        encoding="utf-8",
    )
    output_dir = tmp_path / "runner-output"
    monkeypatch.setenv("PGBENCH_INPUT_FASTQ_R1", str(reads_r1))
    monkeypatch.setenv("PGBENCH_INPUT_FASTQ_R2", str(reads_r2))
    monkeypatch.setenv("PGBENCH_GRAPH_DIR", str(graph_dir))
    monkeypatch.setenv("PGBENCH_OUTPUT_DIR", str(output_dir))
    monkeypatch.setenv(
        "PGBENCH_OUTPUT_VCF",
        str(output_dir / "raw" / "calls.vcf"),
    )
    monkeypatch.setenv("PGBENCH_SAMPLE_ID", "HG002")
    monkeypatch.setenv("PGBENCH_THREADS", "4")
    monkeypatch.setenv("PGBENCH_CANDIDATE_VCF", str(candidates))
    if reference_path is None:
        monkeypatch.delenv("PGBENCH_GRAPH_REFERENCE_PATH", raising=False)
    else:
        monkeypatch.setenv("PGBENCH_GRAPH_REFERENCE_PATH", reference_path)


def test_vg_call_is_restricted_to_frozen_reference_sample(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    _set_runner_environment(monkeypatch, tmp_path, reference_path="GRCh38")
    commands: list[list[str]] = []

    def fake_run(command: list[str], *, stdout: Path | None = None) -> None:
        commands.append(command)
        if stdout is not None:
            stdout.parent.mkdir(parents=True, exist_ok=True)
            stdout.write_bytes(b"synthetic\n")

    monkeypatch.setattr(MODULE, "run", fake_run)
    monkeypatch.setattr(MODULE, "project_all_sites", lambda *args: (1, 0, 0))
    assert MODULE.main() == 0

    call = next(command for command in commands if command[:2] == ["vg", "call"])
    selector_index = call.index("--ref-sample")
    assert call[selector_index + 1] == "GRCh38"

    giraffe = next(
        command for command in commands if command[:2] == ["vg", "giraffe"]
    )
    assert giraffe.count("-f") == 2
    assert "graph.shortread.withzip.min" in giraffe[giraffe.index("-m") + 1]
    assert "graph.shortread.zipcodes" in giraffe[giraffe.index("-z") + 1]


def test_vg_runner_fails_closed_without_reference_selector(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    _set_runner_environment(monkeypatch, tmp_path, reference_path=None)

    with pytest.raises(
        RuntimeError,
        match="PGBENCH_GRAPH_REFERENCE_PATH",
    ):
        MODULE.main()
