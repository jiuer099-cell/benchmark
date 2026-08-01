from __future__ import annotations

import hashlib
import subprocess
import sys
from pathlib import Path

SCRIPTS = Path(__file__).resolve().parents[2] / "workflow" / "scripts"
sys.path.insert(0, str(SCRIPTS))

from obsidian_sync import (  # noqa: E402
    BEGIN_MARKER,
    END_MARKER,
    inspect_sync,
    resolve_paths,
    synchronize_design,
)


def _git(root: Path, *args: str) -> None:
    subprocess.run(["git", *args], cwd=root, check=True)


def _workspace(tmp_path: Path) -> tuple[Path, Path, Path]:
    repo = tmp_path / "repo"
    vault = tmp_path / "vault"
    note = vault / "Writing" / "design.md"
    repo.mkdir()
    note.parent.mkdir(parents=True)
    (repo / "docs").mkdir()
    (repo / "config").mkdir()
    (repo / "workflow" / "schemas").mkdir(parents=True)
    (repo / "plugins" / "example").mkdir(parents=True)
    (repo / ".claude" / "project-memory").mkdir(parents=True)

    (repo / "docs" / "design.md").write_text(
        "# Confirmed design\n\n"
        "Exact generated body.\n\n"
        "```text\n"
        f"{BEGIN_MARKER}\n"
        "example\n"
        f"{END_MARKER}\n"
        "```\n",
        encoding="utf-8",
    )
    (repo / "config" / "score.yaml").write_text(
        "profile:\n  id: pgbench_v1\n",
        encoding="utf-8",
    )
    (repo / "config" / "config.yaml").write_text(
        "execution:\n"
        "  random_seed: 42\n"
        "pangenome:\n"
        "  id: synthetic_pg\n"
        "score:\n"
        "  profile: pgbench_v1\n"
        "catalogs:\n"
        "  score_weights: config/score.yaml\n"
        "obsidian:\n"
        "  vault_root: null\n"
        "  vault_root_env: PGBENCH_TEST_VAULT\n"
        "  project_relpath: Research/hg002-grch38-pangenome-sv-benchmark\n"
        "  canonical_note: Writing/design.md\n",
        encoding="utf-8",
    )
    (repo / "workflow" / "rule-registry.yaml").write_text(
        "schema_version: 1\nrules: []\n",
        encoding="utf-8",
    )
    (repo / "workflow" / "schemas" / "pangenome-manifest.schema.yaml").write_text(
        "type: object\n",
        encoding="utf-8",
    )
    (repo / "plugins" / "example" / "tool.yaml").write_text(
        "id: example\nversion: 1\n",
        encoding="utf-8",
    )
    binding = (
        repo / ".claude" / "project-memory" / "hg002-grch38-pangenome-sv-benchmark.md"
    )
    binding.write_text(
        "---\n"
        "project_slug: hg002-grch38-pangenome-sv-benchmark\n"
        f"vault_root: {vault}\n"
        "canonical_note: Research/hg002-grch38-pangenome-sv-benchmark/"
        "Writing/design.md\n"
        "---\n\nBinding only.\n",
        encoding="utf-8",
    )
    note.write_text(
        "---\n"
        "type: writing\n"
        "status: active\n"
        "---\n\n"
        "Human annotation.\n\n"
        f"{BEGIN_MARKER}\nold\n{END_MARKER}\n",
        encoding="utf-8",
    )

    _git(repo, "init", "-q")
    _git(repo, "config", "user.name", "PGBench Test")
    _git(repo, "config", "user.email", "pgbench@example.invalid")
    _git(repo, "add", ".")
    _git(repo, "commit", "-qm", "fixture")
    return repo, vault, binding


def _paths(repo: Path, vault: Path, binding: Path):
    return resolve_paths(
        repo_root=repo,
        design=Path("docs/design.md"),
        config=Path("config/config.yaml"),
        rule_registry=Path("workflow/rule-registry.yaml"),
        tool_root=Path("plugins"),
        binding=binding.relative_to(repo),
        vault_root=vault,
        pangenome_schema=Path("workflow/schemas/pangenome-manifest.schema.yaml"),
    )


def test_clean_repo_sync_writes_exact_design_and_current_hashes(
    tmp_path: Path,
) -> None:
    repo, vault, binding = _workspace(tmp_path)
    paths = _paths(repo, vault, binding)
    result = synchronize_design(
        paths,
        code_state="partial",
        snakemake_version="9.23.1",
    )
    assert result["status"] == "current"
    assert result["previous_canonical_note_sha256"]
    assert result["canonical_note_sha256"]
    text = paths.canonical_note.read_text(encoding="utf-8")
    assert "Human annotation." in text
    assert "# Confirmed design\n\nExact generated body." in text
    assert text.count(BEGIN_MARKER) == 2
    assert text.count(END_MARKER) == 2
    assert "code_state: partial" in text
    expected_design_sha = hashlib.sha256(
        (repo / "docs" / "design.md").read_bytes()
    ).hexdigest()
    assert f"design_sha256: {expected_design_sha}" in text


def test_uncommitted_design_change_is_reported_stale(tmp_path: Path) -> None:
    repo, vault, binding = _workspace(tmp_path)
    paths = _paths(repo, vault, binding)
    synchronize_design(
        paths,
        code_state="partial",
        snakemake_version="9.23.1",
    )
    (repo / "docs" / "design.md").write_text(
        "# Changed without sync\n", encoding="utf-8"
    )
    result = inspect_sync(
        paths,
        code_state="partial",
        snakemake_version="9.23.1",
    )
    assert result["status"] == "stale"
    fields = {mismatch["field"] for mismatch in result["mismatches"]}
    assert "design_sha256" in fields
    assert "generated_design" in fields
    assert "git_dirty" in fields


def test_vault_root_uses_configured_project_relpath(tmp_path: Path) -> None:
    repo, project_root, binding = _workspace(tmp_path)
    vault_root = tmp_path / "memory"
    nested_project = vault_root / "Research" / "hg002-grch38-pangenome-sv-benchmark"
    nested_project.parent.mkdir(parents=True)
    project_root.rename(nested_project)
    paths = resolve_paths(
        repo_root=repo,
        design=Path("docs/design.md"),
        config=Path("config/config.yaml"),
        rule_registry=Path("workflow/rule-registry.yaml"),
        tool_root=Path("plugins"),
        binding=binding.relative_to(repo),
        vault_root=vault_root,
        pangenome_schema=Path("workflow/schemas/pangenome-manifest.schema.yaml"),
    )
    assert paths.vault_root == nested_project
    assert paths.canonical_note == nested_project / "Writing" / "design.md"
