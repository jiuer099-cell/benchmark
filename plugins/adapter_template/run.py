#!/usr/bin/env python3
"""External adapter scaffold. Implement run_tool(), then pass reviewed gates.

Core owns public-input validation, scoring, admission and final sealing.
This stdlib-only runner never imports native Python extensions into Core's
interpreter. Use native_python() to invoke helpers in the tool environment.
"""
from __future__ import annotations

import dataclasses
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from adapter_projection import digest, open_text, project_all_sites  # noqa: E402


def required(name: str) -> str:
    value = os.environ.get(name)
    if not value:
        raise RuntimeError("required environment variable is absent: " + name)
    return value


def require_readonly_input(path: Path, label: str) -> None:
    """Injection check only: do not repeat Core read-count/hash validation."""
    if not path.exists() or (path.is_file() and path.stat().st_size == 0):
        raise RuntimeError(f"Core injected missing/empty validated {label}: {path}")


def detect_input_mode() -> str:
    modes = {
        "shared_bam": bool(os.environ.get("PGBENCH_SHARED_ALIGNMENT_BAM")),
        "short_fastq": bool(os.environ.get("PGBENCH_INPUT_FASTQ_R1")),
        "long_fastq": bool(os.environ.get("PGBENCH_INPUT_LONG_READS_FASTQ")),
    }
    selected = [name for name, present in modes.items() if present]
    if len(selected) != 1:
        raise RuntimeError("no primary read mode or more than one primary read evidence mode")
    return selected[0]


def injected_inputs() -> dict[str, Path]:
    payload = json.loads(Path(required("PGBENCH_RESOLVED_INPUTS")).read_text(encoding="utf-8"))
    entries = payload.get("inputs") if isinstance(payload, dict) else None
    if not isinstance(entries, list):
        raise RuntimeError("Core resolved-input manifest has no input list")
    inputs = {}
    for entry in entries:
        if not isinstance(entry, dict):
            raise RuntimeError("malformed Core resolved input")
        name, variable = entry.get("contract_name"), entry.get("environment_variable")
        if not isinstance(name, str) or not isinstance(variable, str) or name in inputs:
            raise RuntimeError("missing/duplicate Core input binding")
        path = Path(required(variable))
        if entry.get("path") is not None and path != Path(entry["path"]):
            raise RuntimeError("injected path differs from the Core resolved-input binding")
        require_readonly_input(path, name)
        inputs[name] = path
    return inputs


@dataclasses.dataclass(frozen=True)
class AdapterContext:
    input_mode: str
    sample: str
    threads: int
    reference: Path
    candidates: Path
    output_dir: Path
    output_vcf: Path
    work_dir: Path
    inputs: dict[str, Path] = dataclasses.field(default_factory=dict)
    parameters: dict = dataclasses.field(default_factory=dict)
    bam: Path | None = None
    bai: Path | None = None
    fastq_r1: Path | None = None
    fastq_r2: Path | None = None
    long_fastq: Path | None = None
    optional_inputs: dict[str, Path] = dataclasses.field(default_factory=dict)


def native_python(tool_executable: str) -> Path:
    """Use the interpreter beside an executable in the frozen native env."""
    executable = shutil.which(tool_executable)
    if not executable:
        raise RuntimeError("native executable missing from declared environment: " + tool_executable)
    interpreter = Path(executable).parent / "python"
    if not interpreter.is_file():
        raise RuntimeError("native helper requires a declared Python beside " + executable)
    return interpreter


def run_native(context: AdapterContext, argv: list[str]) -> None:
    """Record real invocations; errors are fatal, never successful no-calls.

    Build argv from the pre-reviewed parameter contract. Generic code cannot
    identify the biological meaning of arbitrary tool flags: conformance and
    human review must check official defaults/recommendations and exceptions.
    """
    if not argv or not all(isinstance(value, str) and value for value in argv):
        raise ValueError("native command must be a nonempty argv list")
    log = context.work_dir / "native-commands.jsonl"
    with log.open("a", encoding="utf-8") as stream:
        stream.write(json.dumps({"argv": argv, "state": "started"}) + "\n")
    result = subprocess.run(argv, cwd=context.work_dir, check=False)
    with log.open("a", encoding="utf-8") as stream:
        stream.write(json.dumps({"argv": argv, "state": "finished", "returncode": result.returncode}) + "\n")
    result.check_returncode()


def run_tool(context: AdapterContext) -> Path:
    """Implement the official frozen native chain using run_native().

    Return a nonempty, newly owned biallelic sequence-resolved native VCF inside
    work_dir. If necessary, use an official ALT-specific converter and retain
    BOTH original input and converted output here. Never copy GT=1/2 into a
    biallelic record, infer reference calls for silence, tune to truth/score,
    or add non-default genotype filters/rejection switches. Declare and review
    compatibility/resource/conversion exceptions before running.

    Tool-specific chunking must preserve its algorithm and complete evidence;
    do not assume an arbitrary per-contig split is valid for every method.
    Clean disposable temporary indexes before returning, but retain all native
    outputs and conversion/mapping evidence needed to reproduce the calls.
    """
    raise NotImplementedError("implement and review run_tool() for the external tool")


def main() -> int:
    mode = detect_input_mode()
    sample, threads = required("PGBENCH_SAMPLE_ID"), int(required("PGBENCH_THREADS"))
    if threads < 1:
        raise ValueError("thread budget must be positive")
    output_dir = Path(required("PGBENCH_OUTPUT_DIR")).resolve()
    output = Path(required("PGBENCH_OUTPUT_VCF")).resolve()
    output.relative_to(output_dir)
    if output.suffix != ".vcf":
        raise ValueError("template emits plain .vcf; do not declare .vcf.gz")
    inputs = injected_inputs()
    purpose = required("PGBENCH_EXECUTION_PURPOSE")
    if purpose not in {"formal", "development_only"}:
        raise ValueError("invalid Core execution purpose")
    parameters = json.loads(required("PGBENCH_TOOL_PARAMETERS"))
    policy = json.loads(required("PGBENCH_FORMAL_PRIMARY_POLICY"))
    tool_id = required("PGBENCH_TOOL_ID")
    report_path = None
    report = None
    if purpose == "formal":
        report_path = Path(required("PGBENCH_CONFORMANCE_REPORT"))
        report = json.loads(report_path.read_text(encoding="utf-8"))
        if (report.get("contract") != "pgbench-adapter-conformance-v2"
                or report.get("production_eligibility") != "PASS" or report.get("tool_id") != tool_id
                or report.get("parameter_contract") != parameters or report.get("formal_primary_policy") != policy):
            raise ValueError("a Core-verified reviewed conformance receipt is required")
    work = output_dir / "tool-work"
    work.mkdir(parents=True, exist_ok=True)
    context = AdapterContext(mode, sample, threads, inputs["reference"], inputs["candidate_panel"],
        output_dir, output, work, inputs=inputs, parameters=parameters)
    bindings = {
        "shared_bam": {"bam": "shared_shortread_alignment", "bai": "shared_shortread_alignment_index"},
        "short_fastq": {"fastq_r1": "short_fastq_r1", "fastq_r2": "short_fastq_r2"},
        "long_fastq": {"long_fastq": "long_reads_fastq"},
    }
    context = dataclasses.replace(context, **{name: inputs[token] for name, token in bindings[mode].items()})
    context.optional_inputs.update({name: inputs[name] for name in ("reference_index", "graph_assets", "tool_index") if name in inputs})
    native_path = run_tool(context)
    if native_path.is_symlink():
        raise ValueError("native output must not be a symlink")
    native = native_path.resolve(strict=True)
    native.relative_to(work.resolve())
    if not native.is_file() or native.stat().st_size == 0 or native.stat().st_nlink != 1:
        raise ValueError("native output must be a newly owned nonempty file under tool-work")
    commands = work / "native-commands.jsonl"
    entries = [json.loads(line) for line in commands.read_text().splitlines() if line.strip()]
    starts = [row["argv"] for row in entries if row.get("state") == "started"]
    finishes = [row for row in entries if row.get("state") == "finished"]
    if not starts or [row.get("argv") for row in finishes] != starts or any(row.get("returncode") != 0 for row in finishes):
        raise ValueError("native invocations are incomplete or failed")
    matched, no_call = project_all_sites(native, context.candidates, output, sample,
        reference=context.reference, reference_index=inputs["reference_index"])
    # tool-work remains retained in this attempt; Core and the final seal bind it.
    frozen = work / "frozen-native-policy.json"
    frozen.write_text(json.dumps({"tool_id": tool_id, "parameters": parameters,
        "formal_primary_policy": policy, "execution_purpose": purpose}, sort_keys=True) + "\n", encoding="utf-8")
    artifacts = {"native_commands": commands, "frozen_parameters": frozen}
    if report_path is not None:
        artifacts["conformance_receipt"] = report_path
        plan = report["gate_plan"]
        plan_path = Path(plan["path"])
        if digest(plan_path) != plan["sha256"]:
            raise ValueError("reviewed gate plan changed during native execution")
        artifacts["conformance_gate_plan"] = plan_path
    for gate, row in (report["gates"] if report else {}).items():
        for index, entry in enumerate([row["log"], *row["evidence"]]):
            path = Path(entry["path"])
            if digest(path) != entry["sha256"]:
                raise ValueError("conformance evidence changed during native execution")
            artifacts[f"conformance/{gate}/{index}"] = path
    for path in work.rglob("*"):
        if path.is_symlink():
            raise ValueError("retained native artifacts must not be symlinks")
        if path.is_file():
            artifacts["tool_work/" + path.relative_to(work).as_posix()] = path
    audit = {"contract": "pgbench-native-output-audit-v1", "truth_used": False,
        "tool_id": tool_id, "output_sha256": digest(output),
        "execution_policy": ("pre_frozen_official_defaults_or_recommendations_no_extra_genotype_filtering_one_formal_pgf1"
                             if purpose == "formal" else "development_only_unsealed_conformance_probe"),
        "execution_purpose": purpose, "formal_primary_policy": policy, "parameter_contract": context.parameters,
        "artifacts": {name: {"path": str(path.resolve()), "sha256": digest(path)} for name, path in artifacts.items()}}
    output.with_suffix(".native-audit.json").write_text(json.dumps(audit, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"adapter projection: mapped={matched} no_call={no_call}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
