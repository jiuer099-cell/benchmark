from __future__ import annotations

import subprocess
import sys
from pathlib import Path

SCRIPTS = Path(__file__).resolve().parents[2] / "workflow" / "scripts"
sys.path.insert(0, str(SCRIPTS))

from snapshot_run_context import capture_run_context  # noqa: E402


def test_run_context_records_dirty_state_and_profile_hash(tmp_path: Path) -> None:
    root = tmp_path / "repository"
    root.mkdir()
    profile = root / "me_f1_scoring.yaml"
    profile.write_text("schema_version: 1\n", encoding="utf-8")
    subprocess.run(["git", "init", "-q"], cwd=root, check=True)
    subprocess.run(
        ["git", "config", "user.email", "pgbench@example.invalid"],
        cwd=root,
        check=True,
    )
    subprocess.run(
        ["git", "config", "user.name", "PGBench Test"],
        cwd=root,
        check=True,
    )
    subprocess.run(["git", "add", profile.name], cwd=root, check=True)
    subprocess.run(
        ["git", "commit", "-q", "-m", "fixture"], cwd=root, check=True
    )
    profile.write_text("schema_version: 1\nchanged: true\n", encoding="utf-8")
    context = capture_run_context(
        repo_root=root,
        score_profile=profile,
        random_seed=20260717,
        snakemake_version="9.23.1",
        execution_profile="local",
    )
    assert len(context["git_head"]) == 40
    assert context["git_dirty"] is True
    assert len(context["repo_state_sha256"]) == 64
    assert len(context["score_profile_sha256"]) == 64
    assert context["random_seed"] == 20260717
    assert context["snakemake_version"] == "9.23.1"
