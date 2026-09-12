#!/usr/bin/env python3
"""Conservatively prepare a Minigraph-Cactus VCF for PanGenie P0.

Only an unresolved *partial* missing slash GT (for example ``0/.``) is
changed, and it is changed to the less informative ``./.``.  This does not
infer genotype or phase.  Fully observed slash heterozygotes remain a hard
failure; homozygous slash GTs are phase-equivalent and canonicalized.  A JSON
audit makes every transformation reviewable before official preprocessing.
"""

from __future__ import annotations

import argparse
import gzip
import json
import sys
from collections import Counter
from pathlib import Path
from typing import TextIO


class PhaseGateError(ValueError):
    """A source genotype cannot be safely represented for PanGenie."""


def open_text(path: Path, mode: str) -> TextIO:
    if path.suffix == ".gz":
        return gzip.open(path, mode, encoding="utf-8")  # type: ignore[return-value]
    return path.open(mode, encoding="utf-8")


def normalize_gt(genotype: str, *, line_number: int) -> tuple[str, str]:
    """Return an allowed GT and the deterministic action taken.

    Masking partial-missing calls is deliberately one-way: it removes an
    uncertain haplotype assignment instead of manufacturing one.  The caller
    will account for the missing allele under the official MC <=20% filter.
    """

    if "/" not in genotype:
        return genotype, "unchanged"
    alleles = genotype.split("/")
    if len(alleles) != 2 or any(
        allele != "." and not allele.isdigit() for allele in alleles
    ):
        raise PhaseGateError(f"line {line_number}: malformed slash GT {genotype!r}")
    if alleles == [".", "."]:
        return "./.", "fully_missing_unchanged"
    if "." in alleles:
        return "./.", "partial_missing_masked_to_no_call"
    if alleles[0] == alleles[1]:
        return f"{alleles[0]}|{alleles[1]}", "homozygous_canonicalized"
    raise PhaseGateError(
        f"line {line_number}: phase-ambiguous heterozygous GT {genotype!r}"
    )


def transform(source: Path, destination: Path, audit: Path) -> dict[str, int]:
    statistics: Counter[str] = Counter(
        {
            "records": 0,
            "sample_genotypes": 0,
            "fully_missing_slash_gt": 0,
            "partial_missing_slash_masked_to_no_call": 0,
            "homozygous_slash_canonicalized": 0,
            "phase_ambiguous_heterozygous_slash_gt": 0,
            "malformed_slash_gt": 0,
        }
    )
    samples: list[str] | None = None
    destination.parent.mkdir(parents=True, exist_ok=True)
    with open_text(source, "rt") as reader, open_text(destination, "wt") as writer:
        for line_number, line in enumerate(reader, start=1):
            if line.startswith("##") or not line.strip():
                writer.write(line)
                continue
            fields = line.rstrip("\n").split("\t")
            if line.startswith("#CHROM"):
                samples = fields[9:]
                writer.write(line)
                continue
            if samples is None or len(fields) != 9 + len(samples):
                raise PhaseGateError(f"line {line_number}: invalid VCF sample columns")
            format_fields = fields[8].split(":")
            if "GT" not in format_fields:
                raise PhaseGateError(f"line {line_number}: missing GT field")
            gt_index = format_fields.index("GT")
            for index, sample_value in enumerate(fields[9:], start=9):
                values = sample_value.split(":")
                if gt_index >= len(values):
                    raise PhaseGateError(f"line {line_number}: missing GT value")
                original = values[gt_index]
                statistics["sample_genotypes"] += 1
                try:
                    normalized, action = normalize_gt(original, line_number=line_number)
                except PhaseGateError as exc:
                    if "/" in original:
                        if original.count("/") != 1 or any(
                            part != "." and not part.isdigit()
                            for part in original.split("/")
                        ):
                            statistics["malformed_slash_gt"] += 1
                        else:
                            statistics["phase_ambiguous_heterozygous_slash_gt"] += 1
                    raise exc
                if action == "fully_missing_unchanged":
                    statistics["fully_missing_slash_gt"] += 1
                elif action == "partial_missing_masked_to_no_call":
                    statistics["partial_missing_slash_masked_to_no_call"] += 1
                elif action == "homozygous_canonicalized":
                    statistics["homozygous_slash_canonicalized"] += 1
                values[gt_index] = normalized
                fields[index] = ":".join(values)
            writer.write("\t".join(fields) + "\n")
            statistics["records"] += 1
    audit.parent.mkdir(parents=True, exist_ok=True)
    audit.write_text(json.dumps(dict(statistics), indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return dict(statistics)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--audit", required=True, type=Path)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        print(json.dumps(transform(args.input, args.output, args.audit), sort_keys=True))
    except PhaseGateError as exc:
        print(f"mask_pangenie_partial_missing: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
