from __future__ import annotations

import copy
import json
import shutil
import sys
from pathlib import Path

import jsonschema
import pytest
import yaml

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "workflow/scripts"))
from validate_adapter_conformance import (
    AdapterConformanceError, GATES, PRIMARY_POLICY, SCHEMA,
    validate, verify_receipt, validate_primary_configuration,
)


def _kit(tmp_path):
    root = tmp_path / "test-tool"
    shutil.copytree(ROOT / "plugins/adapter_template", root, ignore=shutil.ignore_patterns("__pycache__"))
    manifest = yaml.safe_load((root / "tool.template.yaml").read_text())
    manifest["id"] = "fixture-tool"
    manifest["execution"]["trust_level"] = "trusted_reviewed"
    manifest["execution"]["sandbox_backend"] = "bwrap"
    manifest["parameter_contract"]["source_release"] = {
        "tag": "fixture-v1", "commit": "1" * 40, "source_url": "https://example.invalid/fixture",
        "archive_sha256": "2" * 64,
    }
    (root / "tool.template.yaml").unlink()
    path = root / "tool.yaml"
    path.write_text(yaml.safe_dump(manifest))
    conformance = root / "conformance"
    conformance.mkdir()
    plan = {"contract": "pgbench-adapter-gate-plan-v2", "gates": {}}
    for gate in GATES:
        proof = "conformance/" + gate + ".json"
        # Synthetic gate-runner fixture only; production tools need real assertions.
        command = [sys.executable, "-c", f"from pathlib import Path; Path({proof!r}).write_text('fixture evidence'); print('fixture executed')"]
        plan["gates"][gate] = {"command": command, "evidence": [proof]}
    checks = conformance / "plan.json"
    checks.write_text(json.dumps(plan))
    return root, path, checks, conformance / "report.json"


def test_template_schema_requires_parameters_audits_and_declared_reference_index(tmp_path):
    _, manifest_path, _, _ = _kit(tmp_path)
    manifest = yaml.safe_load(manifest_path.read_text())
    validator = jsonschema.Draft202012Validator(yaml.safe_load(SCHEMA.read_text()))
    validator.validate(manifest)
    for field in ("parameter_contract", "formal_primary_policy"):
        broken = copy.deepcopy(manifest)
        broken.pop(field)
        assert list(validator.iter_errors(broken))
    for flag in ("native_projection_audit", "native_output_audit"):
        broken = copy.deepcopy(manifest)
        broken["outputs"][flag] = False
        assert list(validator.iter_errors(broken))
    broken = copy.deepcopy(manifest)
    broken["supported_modes"]["end_to_end_from_reads"]["required_inputs"].remove("reference_index")
    assert list(validator.iter_errors(broken))


def test_executed_gate_receipt_binds_commands_sources_and_evidence(tmp_path):
    _, manifest, checks, out = _kit(tmp_path)
    result = validate(checks, out, tool_manifest=manifest, execute=True)
    assert result["production_eligibility"] == "PASS"
    assert result["formal_primary_policy"] == PRIMARY_POLICY
    assert set(result["gates"]) == set(GATES)
    assert verify_receipt(out, manifest) == result


def test_handwritten_pass_labels_cannot_enable_production(tmp_path):
    _, manifest, checks, out = _kit(tmp_path)
    checks.write_text(json.dumps({name: "PASS" for name in GATES}))
    with pytest.raises(AdapterConformanceError, match="--execute"):
        validate(checks, out, tool_manifest=manifest)
    with pytest.raises(AdapterConformanceError, match="execution plan"):
        validate(checks, out, tool_manifest=manifest, execute=True)
    assert not out.exists()


@pytest.mark.parametrize("change", ["source", "gate_helper", "plan", "evidence", "parameters", "policy"])
def test_receipt_rejects_any_changed_source_evidence_or_configuration(tmp_path, change):
    root, manifest, checks, out = _kit(tmp_path)
    validate(checks, out, tool_manifest=manifest, execute=True)
    if change == "source":
        with (root / "run.py").open("a") as stream:
            stream.write("\n# changed implementation\n")
    elif change == "gate_helper":
        (root / "conformance/test_probe.py").write_text("print('changed gate implementation')\n")
    elif change == "plan":
        checks.write_text("{}")
    elif change == "evidence":
        (root / "conformance/real_input.json").write_text("changed")
    else:
        value = yaml.safe_load(manifest.read_text())
        if change == "parameters":
            value["parameter_contract"]["explicit"]["threads"] = 32
        else:
            value["formal_primary_policy"]["extra_genotype_filtering_or_rejection"] = "allowed"
        manifest.write_text(yaml.safe_dump(value))
    with pytest.raises((ValueError, jsonschema.ValidationError)):
        verify_receipt(out, manifest)


def test_failed_gate_never_publishes_eligibility(tmp_path):
    _, manifest, checks, out = _kit(tmp_path)
    plan = json.loads(checks.read_text())
    plan["gates"]["environment"]["command"] = [sys.executable, "-c", "raise SystemExit(7)"]
    checks.write_text(json.dumps(plan))
    with pytest.raises(AdapterConformanceError, match="failed executed gate"):
        validate(checks, out, tool_manifest=manifest, execute=True)
    assert not out.exists()


def test_missing_gate_evidence_never_publishes_eligibility(tmp_path):
    _, manifest, checks, out = _kit(tmp_path)
    plan = json.loads(checks.read_text())
    plan["gates"]["environment"]["evidence"] = ["conformance/not-created.json"]
    checks.write_text(json.dumps(plan))
    with pytest.raises(AdapterConformanceError, match="evidence"):
        validate(checks, out, tool_manifest=manifest, execute=True)
    assert not out.exists()


def test_gate_cannot_reuse_stale_evidence_with_a_noop_command(tmp_path):
    import os
    root, manifest, checks, out = _kit(tmp_path)
    stale = root / "conformance/environment.json"
    stale.write_text("old successful probe")
    os.utime(stale, (100, 100))
    plan = json.loads(checks.read_text())
    plan["gates"]["environment"]["command"] = [sys.executable, "-c", "print('no probe')"]
    checks.write_text(json.dumps(plan))
    with pytest.raises(AdapterConformanceError, match="fresh evidence"):
        validate(checks, out, tool_manifest=manifest, execute=True)
    assert not out.exists()


def test_gate_plan_cannot_change_while_tests_execute(tmp_path):
    root, manifest, checks, out = _kit(tmp_path)
    plan = json.loads(checks.read_text())
    plan["gates"]["schema"]["command"] = [sys.executable, "-c",
        "from pathlib import Path; Path('conformance/schema.json').write_text('proof'); "
        "Path('conformance/plan.json').write_text('{}')"]
    checks.write_text(json.dumps(plan))
    with pytest.raises(AdapterConformanceError, match="plan changed while running"):
        validate(checks, out, tool_manifest=manifest, execute=True)
    assert not out.exists()
    assert not (root / "conformance/environment.json").exists()


@pytest.mark.parametrize("field", ["command", "evidence", "timeout_seconds", "log"])
def test_receipt_rows_must_match_the_frozen_execution_plan(tmp_path, field):
    _, manifest, checks, out = _kit(tmp_path)
    validate(checks, out, tool_manifest=manifest, execute=True)
    receipt = json.loads(out.read_text())
    row = receipt["gates"]["schema"]
    if field == "command":
        row[field] = [sys.executable, "-c", "print('different')"]
    elif field == "timeout_seconds":
        row[field] = 1
    elif field == "evidence":
        row[field] = receipt["gates"]["environment"]["evidence"]
    else:
        row[field] = receipt["gates"]["environment"]["log"]
    out.write_text(json.dumps(receipt))
    with pytest.raises(AdapterConformanceError, match="differs from frozen plan"):
        verify_receipt(out, manifest)


@pytest.mark.parametrize("proof", ["conformance/plan.json", "conformance/report.json", "conformance/schema.log", "conformance/environment.json"])
def test_gate_proofs_cannot_alias_plan_receipt_logs_or_other_proofs(tmp_path, proof):
    _, manifest, checks, out = _kit(tmp_path)
    plan = json.loads(checks.read_text())
    plan["gates"]["schema"]["evidence"] = [proof]
    checks.write_text(json.dumps(plan))
    with pytest.raises(AdapterConformanceError, match="aliases"):
        validate(checks, out, tool_manifest=manifest, execute=True)
    assert not out.exists()


@pytest.mark.parametrize("filename", ["test-data.json", "test-helper.js"])
def test_conformance_fixtures_and_non_python_helpers_are_frozen(tmp_path, filename):
    root, manifest, checks, out = _kit(tmp_path)
    helper = root / "conformance" / filename
    helper.write_text("original fixture")
    validate(checks, out, tool_manifest=manifest, execute=True)
    helper.write_text("changed fixture")
    with pytest.raises(AdapterConformanceError, match="stale/unbound"):
        verify_receipt(out, manifest)


def test_no_pass_receipt_is_left_if_a_later_gate_changes_previous_evidence(tmp_path):
    _, manifest, checks, out = _kit(tmp_path)
    plan = json.loads(checks.read_text())
    plan["gates"]["environment"]["command"] = [sys.executable, "-c",
        "from pathlib import Path; Path('conformance/environment.json').write_text('proof'); "
        "Path('conformance/schema.json').write_text('tampered')"]
    checks.write_text(json.dumps(plan))
    with pytest.raises(AdapterConformanceError, match="changed conformance evidence"):
        validate(checks, out, tool_manifest=manifest, execute=True)
    assert not out.exists()


@pytest.mark.skipif(sys.platform == "win32", reason="symlink fixture requires Linux")
@pytest.mark.parametrize("target", ["plan", "proof", "parent"])
def test_conformance_rejects_symlink_plans_proofs_and_parent_directories(tmp_path, target):
    root, manifest, checks, out = _kit(tmp_path)
    if target == "plan":
        original = checks.with_name("original.json")
        checks.rename(original)
        checks.symlink_to(original)
    elif target == "proof":
        actual = root / "conformance/original.json"
        actual.write_text("fixture")
        (root / "conformance/schema.json").symlink_to(actual)
    else:
        actual = root / "original-conformance"
        (root / "conformance").rename(actual)
        (root / "conformance").symlink_to(actual, target_is_directory=True)
    with pytest.raises(AdapterConformanceError, match="symlink|escapes"):
        validate(checks, out, tool_manifest=manifest, execute=True)
    assert not out.exists()


def test_primary_policy_rejects_validation_tuning_and_release_placeholders(tmp_path):
    _, manifest, _, _ = _kit(tmp_path)
    value = yaml.safe_load(manifest.read_text())
    value["information_contract"]["tuning_mode"] = "validation_tuned"
    value["information_contract"]["parameter_source"] = "independent_validation_sample"
    with pytest.raises(AdapterConformanceError, match="official"):
        validate_primary_configuration(value)
    value = yaml.safe_load(manifest.read_text())
    value["parameter_contract"]["source_release"]["commit"] = "0" * 40
    with pytest.raises(AdapterConformanceError, match="placeholder"):
        validate_primary_configuration(value)


def test_unreviewed_extra_exception_is_forbidden(tmp_path):
    _, manifest, _, _ = _kit(tmp_path)
    value = yaml.safe_load(manifest.read_text())
    value["parameter_contract"]["adapter_exceptions"] = [
        {"category": "extra_genotype_filter", "rationale": "higher score", "parameters": {"min_gq": 20}}
    ]
    with pytest.raises(AdapterConformanceError, match="exceptions"):
        validate_primary_configuration(value)


def test_production_executor_checks_receipt_before_resolving_large_inputs(tmp_path, monkeypatch):
    import pgbench_exec
    _, manifest_path, _, _ = _kit(tmp_path)
    manifest = yaml.safe_load(manifest_path.read_text())
    monkeypatch.setattr(pgbench_exec, "load_tool_manifest", lambda *args: manifest)
    monkeypatch.setattr(pgbench_exec, "validate_mode_inputs", lambda *args: None)
    monkeypatch.setattr(pgbench_exec, "resolve_inputs", lambda **kwargs: pytest.fail("must reject before input hashing"))
    with pytest.raises(pgbench_exec.ToolContractError, match="PRODUCTION_NOT_ELIGIBLE"):
        pgbench_exec.execute_tool(tool_manifest_path=manifest_path, schema_path=SCHEMA,
            mode="end_to_end_from_reads", run_id="test", sample_id="HG002", supplied_inputs={},
            output_dir=tmp_path / "output", resolved_inputs_path=tmp_path / "resolved.json",
            attempt_record_path=tmp_path / "attempt.json", log_path=tmp_path / "tool.log",
            threads=1, memory_mb=1024, execution_purpose="formal")


def _native_main(tmp_path, monkeypatch, purpose):
    from test_adapter_template import adapter_template
    _, manifest, checks, receipt = _kit(tmp_path)
    report = validate(checks, receipt, tool_manifest=manifest, execute=True)
    data = tmp_path / "data"
    data.mkdir()
    inputs = {}
    variables = {
        "short_fastq_r1": "PGBENCH_INPUT_FASTQ_R1", "short_fastq_r2": "PGBENCH_INPUT_FASTQ_R2",
        "reference": "PGBENCH_REFERENCE_FASTA", "reference_index": "PGBENCH_REFERENCE_INDEX",
        "candidate_panel": "PGBENCH_CANDIDATE_VCF",
    }
    for token in variables:
        path = data / token
        path.write_text("fixture input\n")
        inputs[token] = path
    inputs["reference"].write_text(">chr1\n" + "A" * 100 + "\n")
    inputs["reference_index"].write_text("chr1\t100\t6\t100\t101\n")
    vcf = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG002\nchr1\t10\tPGSV_1\tA\tAT\t20\tLowQual\tUK=0\tGT:GQ\t0/1:8\n"
    inputs["candidate_panel"].write_text(vcf.replace("0/1:8", "./.:."))
    resolved = data / "resolved.json"
    resolved.write_text(json.dumps({"inputs": [{"contract_name": token, "environment_variable": variable,
        "path": str(inputs[token])} for token, variable in variables.items()]}))
    output_dir = tmp_path / "attempt-output"
    output_dir.mkdir()
    output = output_dir / "raw/calls.vcf"
    monkeypatch.delenv("PGBENCH_SHARED_ALIGNMENT_BAM", raising=False)
    monkeypatch.delenv("PGBENCH_INPUT_LONG_READS_FASTQ", raising=False)
    environment = {
        "PGBENCH_SAMPLE_ID": "HG002", "PGBENCH_THREADS": "1", "PGBENCH_OUTPUT_DIR": str(output_dir),
        "PGBENCH_OUTPUT_VCF": str(output), "PGBENCH_RESOLVED_INPUTS": str(resolved),
        "PGBENCH_EXECUTION_PURPOSE": purpose, "PGBENCH_TOOL_ID": report["tool_id"],
        "PGBENCH_TOOL_PARAMETERS": json.dumps(report["parameter_contract"]),
        "PGBENCH_FORMAL_PRIMARY_POLICY": json.dumps(report["formal_primary_policy"]),
    }
    if purpose == "formal":
        environment["PGBENCH_CONFORMANCE_REPORT"] = str(receipt)
    else:
        monkeypatch.delenv("PGBENCH_CONFORMANCE_REPORT", raising=False)
    environment.update({variable: str(inputs[token]) for token, variable in variables.items()})
    for name, value in environment.items():
        monkeypatch.setenv(name, value)
    def native_tool(context):
        native = context.work_dir / "native.vcf"
        adapter_template.run_native(context, [sys.executable, "-c",
            f"from pathlib import Path; Path({str(native)!r}).write_text({vcf!r})"])
        return native
    monkeypatch.setattr(adapter_template, "run_tool", native_tool)
    assert adapter_template.main() == 0
    return output, inputs, report


def test_template_real_runner_emits_complete_formal_native_and_projection_audits(tmp_path, monkeypatch):
    output, _, report = _native_main(tmp_path, monkeypatch, "formal")
    native = json.loads(output.with_suffix(".native-audit.json").read_text())
    projection = json.loads(output.with_suffix(".projection.json").read_text())
    assert native["execution_purpose"] == "formal"
    assert native["formal_primary_policy"] == PRIMARY_POLICY
    assert native["parameter_contract"] == report["parameter_contract"]
    assert projection["reference_equivalence_checked"]
    for entry in native["artifacts"].values():
        from validate_adapter_conformance import sha
        assert sha(Path(entry["path"])) == entry["sha256"]
    assert any(name.startswith("conformance/fixed_mini_panel/") for name in native["artifacts"])


def test_development_probe_runs_without_a_circular_receipt_requirement(tmp_path, monkeypatch):
    output, _, _ = _native_main(tmp_path, monkeypatch, "development_only")
    audit = json.loads(output.with_suffix(".native-audit.json").read_text())
    assert audit["execution_purpose"] == "development_only"
    assert "conformance_receipt" not in audit["artifacts"]
    assert audit["execution_policy"] == "development_only_unsealed_conformance_probe"


def test_template_audits_survive_core_publication_and_enter_latest_seal(tmp_path, monkeypatch):
    from test_finalize_pgf1_release import _fixture, _seal, _sha, _json
    output, inputs, _ = _native_main(tmp_path, monkeypatch, "formal")
    # Model Core publication: move the trace and update only its declared path.
    final_dir = tmp_path / "published/raw"
    final_dir.mkdir(parents=True)
    final_output = final_dir / "calls.vcf"
    shutil.copyfile(output, final_output)
    trace = final_output.with_suffix(".native-trace.jsonl.gz")
    shutil.move(output.with_suffix(".native-trace.jsonl.gz"), trace)
    projection = json.loads(output.with_suffix(".projection.json").read_text())
    projection["trace_path"] = str(trace)
    projected_audit = final_output.with_suffix(".projection.json")
    _json(projected_audit, projection)
    native_audit = final_output.with_suffix(".native-audit.json")
    shutil.move(output.with_suffix(".native-audit.json"), native_audit)
    fixture_dir = tmp_path / "seal-fixture"
    fixture_dir.mkdir()
    paths = _fixture(fixture_dir)
    score = json.loads(paths["score"].read_text())
    score["tuple"]["tool"] = "fixture-tool"
    _json(paths["score"], score)
    _json(paths["admission"], {"status": "PASS", "score_sha256": _sha(paths["score"])})
    _json(paths["invariants"], {"canonical_unit_set_sha256": _sha(paths["units"]),
        "canonical_panel_sha256": _sha(inputs["candidate_panel"]), "reference_sha256": _sha(inputs["reference"])})
    package = _seal(paths, require_projection=True, projection_audit=projected_audit, projection_vcf=final_output,
        require_native_output_audit=True, native_output_audit=native_audit, native_output_vcf=final_output)
    assert package["sealed"]
    lineage = package["artifact_hashes"]["native_output_lineage"]
    assert "conformance_receipt" in lineage and "tool_work/native.vcf" in lineage
    assert package["artifact_hashes"]["native_trace"]["path"] == str(trace)


@pytest.mark.skipif(shutil.which("bwrap") is None, reason="real sandbox requires Linux bwrap")
def test_reviewed_template_runs_in_real_sandbox_and_publishes_both_audits(tmp_path, monkeypatch):
    import pgbench_exec
    from validate_adapter_conformance import sha
    root, manifest, checks, receipt = _kit(tmp_path)
    runner = root / "run.py"
    vcf = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG002\nchr1\t10\tPGSV_1\tA\tAT\t20\tLowQual\tUK=0\tGT:GQ\t0/1:8\n"
    replacement = ("native = context.work_dir / 'native.vcf'\n"
        "    native_data = " + repr(vcf) + "\n"
        "    run_native(context, [sys.executable, '-c', "
        "f\"from pathlib import Path; Path({str(native)!r}).write_text({native_data!r})\"])\n"
        "    return native")
    runner.write_text(runner.read_text().replace(
        'raise NotImplementedError("implement and review run_tool() for the external tool")', replacement))
    validate(checks, receipt, tool_manifest=manifest, execute=True)
    data = tmp_path / "public-inputs"
    data.mkdir()
    supplied = {}
    for token in yaml.safe_load(manifest.read_text())["supported_modes"]["end_to_end_from_reads"]["required_inputs"]:
        path = data / token
        path.write_text("fixture opaque public asset\n")
        supplied[token] = path
    supplied["reference"].write_text(">chr1\n" + "A" * 100 + "\n")
    supplied["reference_index"].write_text("chr1\t100\t6\t100\t101\n")
    supplied["candidate_panel"].write_text(vcf.replace("0/1:8", "./.:."))
    supplied["short_fastq_r1"].write_text("@fixture/1\nAAAA\n+\nIIII\n")
    supplied["short_fastq_r2"].write_text("@fixture/2\nTTTT\n+\nIIII\n")
    output = tmp_path / "formal-output"
    output.mkdir()
    monkeypatch.setattr(pgbench_exec, "_input_validation_cache_root", lambda *_: tmp_path / "input-validation-cache")
    result = pgbench_exec.execute_tool(tool_manifest_path=manifest, schema_path=SCHEMA,
        mode="end_to_end_from_reads", run_id="isolated-fixture", sample_id="HG002", supplied_inputs=supplied,
        output_dir=output, resolved_inputs_path=output / "meta/resolved.json",
        attempt_record_path=output / "meta/attempt.json", log_path=tmp_path / "fixture.log",
        threads=1, memory_mb=1024, execution_purpose="formal", timeout_seconds=30)
    raw = Path(result.output.path)
    assert result.output.record_count == 1 and result.output.sample_id == "HG002"
    projection = json.loads(raw.with_suffix(".projection.json").read_text())
    native = json.loads(raw.with_suffix(".native-audit.json").read_text())
    assert sha(Path(projection["trace_path"])) == projection["trace_sha256"]
    assert native["output_sha256"] == sha(raw)
    assert "\t20\tLowQual\t" in raw.read_text()
    for entry in native["artifacts"].values():
        assert sha(Path(entry["path"])) == entry["sha256"]
