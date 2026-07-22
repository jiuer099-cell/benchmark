from __future__ import annotations

import sys
from pathlib import Path

SCRIPTS = Path(__file__).resolve().parents[2] / "workflow" / "scripts"
sys.path.insert(0, str(SCRIPTS))

from snapshot_run_context import capture_run_context  # noqa: E402


def test_run_context_records_dirty_state_and_profile_hash() -> None:
    root = Path(__file__).resolve().parents[2]
    profile = root / "config" / "consensus_scoring.yaml"
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
