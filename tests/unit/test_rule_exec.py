from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Any

import pytest

SCRIPTS = Path(__file__).resolve().parents[2] / "workflow" / "scripts"
sys.path.insert(0, str(SCRIPTS))

from pgbench_provenance import (  # noqa: E402
    sha256_bytes,
    sha256_path,
    validate_manifest,
)
from pgbench_rule_exec import main as rule_exec_main  # noqa: E402

SHA_A = sha256_bytes(b"a")
SHA_B = sha256_bytes(b"b")


def _write_json(path: Path, payload: dict[str, Any]) -> None:
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def _prepare_workspace(root: Path) -> dict[str, Path]:
    paths = {
        "rule_source": root / "test_rule.smk",
        "producer": root / "produce.py",
        "config": root / "config.json",
        "score_profile": root / "score-profile.yaml",
        "run_context": root / "run-context.json",
        "input": root / "input.txt",
        "upstream": root / "upstream.json",
        "output": root / "output.txt",
        "manifest": root / "manifest.json",
    }
    paths["rule_source"].write_text(
        'rule test_rule:\n    shell: "python produce.py"\n',
        encoding="utf-8",
    )
    paths["producer"].write_text(
        "from pathlib import Path\n"
        "import sys\n"
        "source, destination = map(Path, sys.argv[1:])\n"
        "destination.write_text(source.read_text(encoding='utf-8').upper(), "
        "encoding='utf-8')\n",
        encoding="utf-8",
    )
    _write_json(paths["config"], {"benchmark": "synthetic"})
    paths["score_profile"].write_text(
        "profile_id: pgbench_test_v1\n",
        encoding="utf-8",
    )
    _write_json(
        paths["run_context"],
        {
            "git_head": "1" * 40,
            "git_dirty": False,
            "git_diff_sha256": SHA_A,
            "snakemake_version": "9.23.1",
            "execution_profile": "pytest",
            "hardware_fingerprint_sha256": SHA_B,
            "random_seed": 42,
        },
    )
    paths["input"].write_text("synthetic input\n", encoding="utf-8")
    _write_json(
        paths["upstream"],
        {"manifest_id": sha256_bytes(b"upstream-manifest")},
    )
    return paths


def _base_argv(root: Path, paths: dict[str, Path]) -> list[str]:
    return [
        "--rule-name",
        "test_rule",
        "--job-key",
        "HG002",
        "--run-id",
        "synthetic-run",
        "--module-or-tool-id",
        "pgbench-core",
        "--snakefile-path",
        paths["rule_source"].name,
        "--rule-source-path",
        paths["rule_source"].name,
        "--script-or-wrapper-path",
        paths["producer"].name,
        "--config-snapshot",
        paths["config"].name,
        "--score-profile",
        paths["score_profile"].name,
        "--run-context",
        paths["run_context"].name,
        "--truth-profile",
        "synthetic_truth",
        "--snakemake-version",
        "9.23.1",
        "--execution-profile",
        "pytest",
        "--random-seed",
        "42",
        "--threads",
        "1",
        "--resource",
        "mem_mb=128",
        "--param",
        "profile=pgbench_test_v1",
        "--wildcard",
        "sample=HG002",
        "--input",
        paths["input"].name,
        "--output",
        paths["output"].name,
        "--upstream-manifest",
        paths["upstream"].name,
        "--manifest-output",
        paths["manifest"].name,
        "--repo-root",
        str(root),
    ]


def test_successful_argv_execution_writes_valid_hashed_manifest(
    tmp_path: Path,
) -> None:
    paths = _prepare_workspace(tmp_path)
    command = [
        sys.executable,
        paths["producer"].name,
        paths["input"].name,
        paths["output"].name,
    ]

    exit_code = rule_exec_main([*_base_argv(tmp_path, paths), "--", *command])

    assert exit_code == 0
    assert paths["output"].read_text(encoding="utf-8") == "SYNTHETIC INPUT\n"
    manifest = json.loads(paths["manifest"].read_text(encoding="utf-8"))
    validate_manifest(
        manifest,
        verify_paths=True,
        base_dir=tmp_path,
        require_success=True,
    )
    assert manifest["command"] == command
    assert manifest["input_paths"] == [
        paths["input"].name,
        paths["upstream"].name,
    ]
    assert manifest["input_sha256"] == {
        paths["input"].name: sha256_path(paths["input"]),
        paths["upstream"].name: sha256_path(paths["upstream"]),
    }
    assert manifest["output_sha256"] == {
        paths["output"].name: sha256_path(paths["output"])
    }
    assert manifest["upstream_manifest_ids"] == [sha256_bytes(b"upstream-manifest")]


def test_failed_command_does_not_write_success_manifest(tmp_path: Path) -> None:
    paths = _prepare_workspace(tmp_path)
    failing_command = [sys.executable, "-c", "raise SystemExit(7)"]

    exit_code = rule_exec_main([*_base_argv(tmp_path, paths), "--", *failing_command])

    assert exit_code == 7
    assert not paths["manifest"].exists()


def test_failed_rerun_archives_old_manifest_and_leaves_no_canonical_manifest(
    tmp_path: Path,
) -> None:
    paths = _prepare_workspace(tmp_path)
    successful_command = [
        sys.executable,
        paths["producer"].name,
        paths["input"].name,
        paths["output"].name,
    ]
    assert (
        rule_exec_main([*_base_argv(tmp_path, paths), "--", *successful_command]) == 0
    )
    previous_bytes = paths["manifest"].read_bytes()

    failing_command = [sys.executable, "-c", "raise SystemExit(7)"]
    assert rule_exec_main([*_base_argv(tmp_path, paths), "--", *failing_command]) == 7

    assert not paths["manifest"].exists()
    archived = list((paths["manifest"].parent / ".previous").glob("*/manifest.json"))
    assert len(archived) == 1
    assert archived[0].read_bytes() == previous_bytes


def test_successful_rerun_records_recoverable_previous_manifest(
    tmp_path: Path,
) -> None:
    paths = _prepare_workspace(tmp_path)
    command = [
        sys.executable,
        paths["producer"].name,
        paths["input"].name,
        paths["output"].name,
    ]
    assert rule_exec_main([*_base_argv(tmp_path, paths), "--", *command]) == 0
    previous_bytes = paths["manifest"].read_bytes()
    previous = json.loads(previous_bytes)

    assert rule_exec_main([*_base_argv(tmp_path, paths), "--", *command]) == 0

    current = json.loads(paths["manifest"].read_text(encoding="utf-8"))
    archived = Path(current["previous_manifest_archive"])
    assert archived.read_bytes() == previous_bytes
    assert archived.parent.name == current["attempt_id"]
    assert current["attempt_id"] != previous["attempt_id"]


def test_declared_pangenome_manifest_output_is_hashed_after_execution(
    tmp_path: Path,
) -> None:
    paths = _prepare_workspace(tmp_path)
    paths["output"].write_text("stale pangenome manifest\n", encoding="utf-8")
    command = [
        sys.executable,
        paths["producer"].name,
        paths["input"].name,
        paths["output"].name,
    ]

    exit_code = rule_exec_main(
        [
            *_base_argv(tmp_path, paths),
            "--pangenome-manifest",
            paths["output"].name,
            "--",
            *command,
        ]
    )

    assert exit_code == 0
    manifest = json.loads(paths["manifest"].read_text(encoding="utf-8"))
    assert manifest["pangenome_manifest_sha256"] == sha256_path(paths["output"])


def test_input_change_during_execution_fails_without_manifest(
    tmp_path: Path,
) -> None:
    paths = _prepare_workspace(tmp_path)
    command = [
        sys.executable,
        "-c",
        (
            "from pathlib import Path; import sys; "
            "source, destination = map(Path, sys.argv[1:]); "
            "destination.write_text('output\\n'); "
            "source.write_text('mutated input\\n')"
        ),
        paths["input"].name,
        paths["output"].name,
    ]

    exit_code = rule_exec_main([*_base_argv(tmp_path, paths), "--", *command])

    assert exit_code == 3
    assert not paths["manifest"].exists()


def test_code_change_during_execution_fails_without_manifest(
    tmp_path: Path,
) -> None:
    paths = _prepare_workspace(tmp_path)
    paths["producer"].write_text(
        "from pathlib import Path\n"
        "import sys\n"
        "source, destination = map(Path, sys.argv[1:])\n"
        "destination.write_text(source.read_text(encoding='utf-8'), "
        "encoding='utf-8')\n"
        "Path(__file__).write_text('# changed during execution\\n', "
        "encoding='utf-8')\n",
        encoding="utf-8",
    )
    command = [
        sys.executable,
        paths["producer"].name,
        paths["input"].name,
        paths["output"].name,
    ]

    exit_code = rule_exec_main([*_base_argv(tmp_path, paths), "--", *command])

    assert exit_code == 3
    assert not paths["manifest"].exists()


@pytest.mark.parametrize(
    "guarded_path",
    ["config", "score_profile", "run_context", "upstream"],
)
def test_special_provenance_input_change_fails_without_manifest(
    tmp_path: Path,
    guarded_path: str,
) -> None:
    paths = _prepare_workspace(tmp_path)
    command = [
        sys.executable,
        "-c",
        (
            "from pathlib import Path; import sys; "
            "target, destination = map(Path, sys.argv[1:]); "
            "destination.write_text('output\\n'); "
            "target.write_text('changed\\n')"
        ),
        paths[guarded_path].name,
        paths["output"].name,
    ]

    exit_code = rule_exec_main([*_base_argv(tmp_path, paths), "--", *command])

    assert exit_code == 3
    assert not paths["manifest"].exists()


@pytest.mark.parametrize("existing_kind", ["symlink", "directory"])
def test_non_regular_existing_manifest_is_rejected_before_command(
    tmp_path: Path,
    existing_kind: str,
) -> None:
    paths = _prepare_workspace(tmp_path)
    if existing_kind == "symlink":
        outside = tmp_path / "outside-manifest.json"
        outside.write_text("{}\n", encoding="utf-8")
        paths["manifest"].symlink_to(outside)
    else:
        paths["manifest"].mkdir()
    command = [
        sys.executable,
        paths["producer"].name,
        paths["input"].name,
        paths["output"].name,
    ]

    exit_code = rule_exec_main([*_base_argv(tmp_path, paths), "--", *command])

    assert exit_code == 3
    assert not paths["output"].exists()


def test_nullable_environment_and_pangenome_fields_have_reasons(
    tmp_path: Path,
) -> None:
    paths = _prepare_workspace(tmp_path)
    command = [
        sys.executable,
        paths["producer"].name,
        paths["input"].name,
        paths["output"].name,
    ]
    assert rule_exec_main([*_base_argv(tmp_path, paths), "--", *command]) == 0

    manifest = json.loads(paths["manifest"].read_text(encoding="utf-8"))
    reasons = manifest["not_applicable_reason"]
    nullable_fields = (
        "pangenome_manifest_sha256",
        "conda_lock_sha256",
        "container_uri",
        "container_digest",
    )
    for field in nullable_fields:
        assert manifest[field] is None
        assert isinstance(reasons[field], str)
        assert reasons[field].strip()
