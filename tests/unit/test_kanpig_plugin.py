from __future__ import annotations

import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
RUNNER = ROOT / "plugins" / "kanpig" / "run.py"
SPEC = importlib.util.spec_from_file_location("kanpig_adapter", RUNNER)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def test_runner_passes_frozen_thread_allocation(
    tmp_path: Path,
    monkeypatch,
) -> None:
    candidate = tmp_path / "candidate.vcf"
    alignment = tmp_path / "reads.bam"
    reference = tmp_path / "reference.fasta"
    output = tmp_path / "output" / "calls.vcf"
    for path in (candidate, alignment, reference):
        path.write_text("fixture\n", encoding="utf-8")
    monkeypatch.setenv("PGBENCH_CANDIDATE_VCF", str(candidate))
    monkeypatch.setenv("PGBENCH_SHARED_ALIGNMENT", str(alignment))
    monkeypatch.setenv("PGBENCH_REFERENCE_FASTA", str(reference))
    monkeypatch.setenv("PGBENCH_OUTPUT_VCF", str(output))
    monkeypatch.setenv("PGBENCH_THREADS", "16")
    commands: list[list[str]] = []

    class Completed:
        returncode = 0

    def fake_run(command: list[str], *, check: bool):
        assert check is False
        commands.append(command)
        return Completed()

    monkeypatch.setattr(MODULE.subprocess, "run", fake_run)
    assert MODULE.main() == 0
    assert commands == [[
        "kanpig", "gt",
        "--input", str(candidate),
        "--reads", str(alignment),
        "--reference", str(reference),
        "--out", str(output),
        "--threads", "16",
    ]]
