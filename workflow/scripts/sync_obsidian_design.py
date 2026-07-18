#!/usr/bin/env python3
"""Synchronize the repository design into the bound Obsidian canonical note."""

from __future__ import annotations

import argparse
import sys
from collections.abc import Sequence
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path

from obsidian_sync import (
    CODE_STATES,
    ObsidianSyncError,
    resolve_paths,
    synchronize_design,
    write_check,
)


def _snakemake_version() -> str:
    try:
        return version("snakemake")
    except PackageNotFoundError:
        return "not-installed"


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path, default=Path.cwd())
    parser.add_argument(
        "--design",
        type=Path,
        default=Path(
            "docs/superpowers/specs/2026-07-16-hg002-grch37-sv-benchmark-design.md"
        ),
    )
    parser.add_argument(
        "--config", type=Path, default=Path("config/config.example.yaml")
    )
    parser.add_argument(
        "--rule-registry",
        type=Path,
        default=Path("workflow/rule-registry.yaml"),
    )
    parser.add_argument("--tool-root", type=Path, default=Path("plugins"))
    parser.add_argument(
        "--binding",
        type=Path,
        default=Path(".claude/project-memory/hg002-grch37-pangenome-sv-benchmark.md"),
    )
    parser.add_argument("--vault-root", type=Path)
    parser.add_argument("--canonical-note", type=Path)
    parser.add_argument(
        "--pangenome-schema",
        type=Path,
        default=Path("workflow/schemas/pangenome-manifest.schema.yaml"),
    )
    parser.add_argument("--code-state", choices=sorted(CODE_STATES), default="partial")
    parser.add_argument("--output-check", type=Path)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        paths = resolve_paths(
            repo_root=args.repo_root,
            design=args.design,
            config=args.config,
            rule_registry=args.rule_registry,
            tool_root=args.tool_root,
            binding=args.binding,
            vault_root=args.vault_root,
            canonical_note=args.canonical_note,
            pangenome_schema=args.pangenome_schema,
        )
        result = synchronize_design(
            paths,
            code_state=args.code_state,
            snakemake_version=_snakemake_version(),
        )
        if args.output_check is not None:
            write_check(args.output_check, result)
        print(f"sync_obsidian_design: {result['status']} {result['canonical_note']}")
    except (OSError, ObsidianSyncError, ValueError) as exc:
        print(f"sync_obsidian_design: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
