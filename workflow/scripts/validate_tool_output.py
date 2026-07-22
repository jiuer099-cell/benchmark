#!/usr/bin/env python3
"""Validate an external tool VCF under the PGBench output contract."""

from __future__ import annotations

import argparse
import gzip
import re
from collections.abc import Sequence
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import TextIO

try:
    from pgbench_provenance import atomic_write_json, fingerprint_path
except ModuleNotFoundError:  # pragma: no cover - package-style invocation
    from .pgbench_provenance import atomic_write_json, fingerprint_path


_GT_RE = re.compile(r"^(?:\.|\d+)(?:(?:/|\|)(?:\.|\d+))*$")


class ToolOutputValidationError(ValueError):
    """Raised when a tool output violates the declared VCF contract."""


@dataclass(frozen=True)
class ToolOutputValidation:
    """Auditable result of validating one tool-produced VCF."""

    path: str
    sha256: str
    size_bytes: int
    mtime_ns: int
    record_count: int
    candidate_count: int | None
    represented_candidate_count: int | None
    sample_id: str
    candidate_output_contract: str
    compression_kind: str
    full_compression_validation_status: str

    def to_dict(self) -> dict[str, str | int | None]:
        """Return a JSON-serializable representation."""

        return asdict(self)


def _open_text(path: Path) -> TextIO:
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("r", encoding="utf-8")


def _compression_status(path: Path) -> tuple[str, str]:
    if path.suffix != ".gz":
        return "plain_vcf", "not_applicable"
    try:
        with path.open("rb") as handle:
            header = handle.read(18)
    except OSError as exc:
        raise ToolOutputValidationError(
            f"cannot inspect compressed output {path}: {exc}"
        ) from exc
    if len(header) < 10 or header[:2] != b"\x1f\x8b":
        raise ToolOutputValidationError(
            "output ending in .gz does not have a gzip header"
        )
    is_bgzf_header = (
        len(header) >= 18
        and header[3] & 0x04
        and header[12:14] == b"BC"
        and header[14:16] == b"\x02\x00"
    )
    return (
        "bgzf_header" if is_bgzf_header else "gzip",
        "infrastructure_required",
    )


def _ensure_within_output_dir(path: Path, output_dir: Path) -> Path:
    if path.is_symlink():
        raise ToolOutputValidationError("tool output must not be a symlink")
    try:
        root = output_dir.resolve(strict=True)
    except OSError as exc:
        raise ToolOutputValidationError(
            f"output directory cannot be resolved: {output_dir}"
        ) from exc
    try:
        resolved = path.resolve(strict=True)
    except OSError as exc:
        raise ToolOutputValidationError(f"tool output does not exist: {path}") from exc
    try:
        resolved.relative_to(root)
    except ValueError as exc:
        raise ToolOutputValidationError(
            f"tool output escapes the assigned output directory: {path}"
        ) from exc
    if not resolved.is_file():
        raise ToolOutputValidationError(f"tool output is not a regular file: {path}")
    if resolved.stat().st_nlink != 1:
        raise ToolOutputValidationError(
            "tool output must be a newly owned file, not a hardlink"
        )
    return resolved


def _candidate_ids(path: Path) -> set[str]:
    candidate_ids: set[str] = set()
    saw_header = False
    try:
        with _open_text(path) as handle:
            for line in handle:
                if line.startswith("#CHROM"):
                    saw_header = True
                    continue
                if not line.strip() or line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 8:
                    raise ToolOutputValidationError(
                        "candidate VCF contains a record with fewer than 8 columns"
                    )
                candidate_id = fields[2]
                if not candidate_id or candidate_id == ".":
                    raise ToolOutputValidationError(
                        "candidate VCF records must have blinded candidate IDs"
                    )
                if candidate_id in candidate_ids:
                    raise ToolOutputValidationError(
                        f"duplicate candidate ID in candidate VCF: {candidate_id}"
                    )
                candidate_ids.add(candidate_id)
    except (OSError, EOFError, gzip.BadGzipFile) as exc:
        raise ToolOutputValidationError(
            f"cannot read candidate VCF {path}: {exc}"
        ) from exc
    if not saw_header:
        raise ToolOutputValidationError("candidate VCF is missing a #CHROM header")
    if not candidate_ids:
        raise ToolOutputValidationError("candidate VCF contains no records")
    return candidate_ids


def _validate_vcf_records(
    path: Path,
    *,
    sample_id: str,
    required_candidate_ids: set[str] | None,
) -> tuple[int, set[str]]:
    saw_fileformat = False
    saw_header = False
    sample_index: int | None = None
    record_ids: set[str] = set()
    record_count = 0

    try:
        with _open_text(path) as handle:
            for line in handle:
                if line.startswith("##fileformat=VCF"):
                    saw_fileformat = True
                    continue
                if line.startswith("##"):
                    continue
                if line.startswith("#CHROM"):
                    fields = line.rstrip("\n").split("\t")
                    if len(fields) < 10:
                        raise ToolOutputValidationError(
                            "tool VCF must contain FORMAT and a sample column"
                        )
                    samples = fields[9:]
                    if samples != [sample_id]:
                        raise ToolOutputValidationError(
                            "tool VCF must contain exactly one sample named "
                            f"{sample_id!r}"
                        )
                    sample_index = 9
                    saw_header = True
                    continue
                if not line.strip() or line.startswith("#"):
                    continue
                if not saw_header or sample_index is None:
                    raise ToolOutputValidationError(
                        "tool VCF record occurs before the #CHROM header"
                    )
                fields = line.rstrip("\n").split("\t")
                if len(fields) <= sample_index or len(fields) < 10:
                    raise ToolOutputValidationError(
                        "tool VCF record does not contain the declared sample"
                    )
                record_id = fields[2]
                if required_candidate_ids is not None and (
                    not record_id or record_id == "."
                ):
                    raise ToolOutputValidationError(
                        "all-sites tool VCF records must retain candidate IDs"
                    )
                if required_candidate_ids is not None and record_id in record_ids:
                    raise ToolOutputValidationError(
                        f"duplicate record ID in tool VCF: {record_id}"
                    )
                format_keys = fields[8].split(":")
                if "GT" not in format_keys:
                    raise ToolOutputValidationError(
                        f"tool VCF record {record_id} has no GT FORMAT field"
                    )
                gt_index = format_keys.index("GT")
                sample_values = fields[sample_index].split(":")
                if gt_index >= len(sample_values):
                    raise ToolOutputValidationError(
                        f"tool VCF record {record_id} has no GT sample value"
                    )
                genotype = sample_values[gt_index]
                if not _GT_RE.fullmatch(genotype):
                    raise ToolOutputValidationError(
                        f"tool VCF record {record_id} has invalid GT {genotype!r}"
                    )
                record_ids.add(record_id)
                record_count += 1
    except (OSError, EOFError, UnicodeError, gzip.BadGzipFile) as exc:
        raise ToolOutputValidationError(f"cannot read tool VCF {path}: {exc}") from exc

    if not saw_fileformat:
        raise ToolOutputValidationError("tool VCF is missing ##fileformat=VCF")
    if not saw_header:
        raise ToolOutputValidationError("tool VCF is missing a #CHROM header")
    if required_candidate_ids is not None:
        missing = sorted(required_candidate_ids - record_ids)
        if missing:
            preview = ", ".join(missing[:5])
            raise ToolOutputValidationError(
                f"all-sites output is missing {len(missing)} candidate(s): {preview}"
            )
        extra = sorted(record_ids - required_candidate_ids)
        if extra:
            preview = ", ".join(extra[:5])
            raise ToolOutputValidationError(
                f"all-sites output contains {len(extra)} off-panel record(s): {preview}"
            )
    return record_count, record_ids


def validate_tool_output(
    *,
    output_vcf: Path,
    output_dir: Path,
    started_at_ns: int,
    sample_id: str,
    candidate_output_contract: str,
    candidate_vcf: Path | None = None,
) -> ToolOutputValidation:
    """Validate path provenance, creation time, VCF syntax, and all-sites coverage."""

    if candidate_output_contract not in {"all_sites", "variant_sites"}:
        raise ToolOutputValidationError(
            f"unsupported candidate output contract: {candidate_output_contract}"
        )
    resolved = _ensure_within_output_dir(output_vcf, output_dir)
    stat_result = resolved.stat()
    if stat_result.st_size <= 0:
        raise ToolOutputValidationError("tool output VCF is empty")
    if stat_result.st_mtime_ns < started_at_ns:
        raise ToolOutputValidationError(
            "tool output VCF predates the current execution attempt"
        )
    compression_kind, compression_validation = _compression_status(resolved)

    required_candidate_ids: set[str] | None = None
    if candidate_output_contract == "all_sites":
        if candidate_vcf is None:
            raise ToolOutputValidationError(
                "all-sites validation requires the blinded candidate VCF"
            )
        required_candidate_ids = _candidate_ids(candidate_vcf)

    record_count, represented_ids = _validate_vcf_records(
        resolved,
        sample_id=sample_id,
        required_candidate_ids=required_candidate_ids,
    )
    fingerprint = fingerprint_path(resolved)
    return ToolOutputValidation(
        path=str(resolved),
        sha256=fingerprint.sha256,
        size_bytes=fingerprint.size,
        mtime_ns=fingerprint.mtime_ns,
        record_count=record_count,
        candidate_count=(
            len(required_candidate_ids) if required_candidate_ids is not None else None
        ),
        represented_candidate_count=(
            len(required_candidate_ids & represented_ids)
            if required_candidate_ids is not None
            else None
        ),
        sample_id=sample_id,
        candidate_output_contract=candidate_output_contract,
        compression_kind=compression_kind,
        full_compression_validation_status=compression_validation,
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-vcf", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--started-at-ns", required=True, type=int)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument(
        "--candidate-output-contract",
        required=True,
        choices=("all_sites", "variant_sites"),
    )
    parser.add_argument("--candidate-vcf", type=Path)
    parser.add_argument("--validation-json", type=Path)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        result = validate_tool_output(
            output_vcf=args.output_vcf,
            output_dir=args.output_dir,
            started_at_ns=args.started_at_ns,
            sample_id=args.sample_id,
            candidate_output_contract=args.candidate_output_contract,
            candidate_vcf=args.candidate_vcf,
        )
    except ToolOutputValidationError as exc:
        raise SystemExit(f"validate_tool_output: {exc}") from exc
    if args.validation_json is not None:
        atomic_write_json(args.validation_json, result.to_dict())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
