"""Build deterministic JSON and TSV rule lineage artifacts."""

from __future__ import annotations

import argparse
from collections.abc import Sequence
from pathlib import Path

try:
    from pgbench_provenance import (
        ProvenanceError,
        atomic_write_json,
        atomic_write_text,
        build_lineage,
        discover_manifest_paths,
        lineage_tsv,
        load_manifests,
    )
except ModuleNotFoundError:  # pragma: no cover - package-style invocation
    from .pgbench_provenance import (
        ProvenanceError,
        atomic_write_json,
        atomic_write_text,
        build_lineage,
        discover_manifest_paths,
        lineage_tsv,
        load_manifests,
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", action="append", default=[], type=Path)
    parser.add_argument("--manifest-dir", action="append", default=[], type=Path)
    parser.add_argument("--target-manifest-id", action="append", default=[])
    parser.add_argument("--output-json", required=True, type=Path)
    parser.add_argument("--output-tsv", required=True, type=Path)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        paths = discover_manifest_paths(args.manifest, args.manifest_dir)
        if not paths:
            raise ProvenanceError("at least one manifest is required")
        manifests = load_manifests(paths)
        lineage = build_lineage(
            manifests,
            target_manifest_ids=args.target_manifest_id or None,
        )
        atomic_write_json(args.output_json, lineage)
        atomic_write_text(args.output_tsv, lineage_tsv(lineage))
    except ProvenanceError as exc:
        raise SystemExit(f"build_rule_lineage: {exc}") from exc
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
