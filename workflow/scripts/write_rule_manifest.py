"""CLI for preparing and atomically writing one formal rule manifest."""

from __future__ import annotations

import argparse
from collections.abc import Sequence
from pathlib import Path

try:
    from pgbench_provenance import (
        ProvenanceError,
        load_json,
        write_rule_manifest,
    )
except ModuleNotFoundError:  # pragma: no cover - package-style invocation
    from .pgbench_provenance import (
        ProvenanceError,
        load_json,
        write_rule_manifest,
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--payload", required=True, type=Path, help="Input JSON payload"
    )
    parser.add_argument("--output", required=True, type=Path, help="Manifest JSON path")
    parser.add_argument(
        "--base-dir",
        type=Path,
        default=Path.cwd(),
        help="Base directory for relative declared paths",
    )
    parser.add_argument(
        "--no-path-hashes",
        action="store_true",
        help="Use path metadata already present in the payload",
    )
    parser.add_argument(
        "--allow-non-success",
        action="store_true",
        help="Allow a non-success manifest (intended for validation fixtures only)",
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        payload = load_json(args.payload)
        if not isinstance(payload, dict):
            raise ProvenanceError("manifest payload must be a JSON object")
        write_rule_manifest(
            payload,
            args.output,
            base_dir=args.base_dir,
            calculate_path_metadata=not args.no_path_hashes,
            require_success=not args.allow_non_success,
        )
    except ProvenanceError as exc:
        raise SystemExit(f"write_rule_manifest: {exc}") from exc
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
