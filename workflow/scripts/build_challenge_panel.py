#!/usr/bin/env python3
"""Build a blinded genotyping challenge panel.

Phase 1 uses an exact canonical matcher for synthetic fixtures. The matcher
profile is recorded so these artifacts cannot be mistaken for the formal
multi-representation HG002 panel.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import sys
from pathlib import Path
from typing import TextIO, cast

from build_pangenome_manifest import format_info, infer_svtype, parse_info


class ChallengePanelError(ValueError):
    """Raised when a blinded panel cannot be built safely."""


def _open_text(path: Path, mode: str) -> TextIO:
    if path.suffix == ".gz":
        return cast(TextIO, gzip.open(path, mode, encoding="utf-8"))
    return cast(TextIO, path.open(mode, encoding="utf-8"))


def _integer_info(info: dict[str, str | bool], key: str, default: int) -> int:
    value = info.get(key)
    if value is None or value is True:
        return default
    try:
        return int(str(value).split(",", 1)[0])
    except ValueError as exc:
        raise ChallengePanelError(f"{key} must be integer, found {value!r}") from exc


def canonical_key(fields: list[str]) -> tuple[str, int, str, str, str, int]:
    if len(fields) < 8:
        raise ChallengePanelError("VCF record must contain at least 8 columns")
    chrom, pos_raw, _, ref, alt = fields[:5]
    try:
        pos = int(pos_raw)
    except ValueError as exc:
        raise ChallengePanelError(f"invalid POS {pos_raw!r}") from exc
    info = parse_info(fields[7])
    svtype = infer_svtype(alt, info)
    end = _integer_info(info, "END", pos)
    return chrom, pos, ref, alt, svtype, end


def _truth_gt(fields: list[str]) -> str:
    if len(fields) < 10:
        return "0/1"
    format_keys = fields[8].split(":")
    sample_values = fields[9].split(":")
    if "GT" not in format_keys:
        return "0/1"
    index = format_keys.index("GT")
    if index >= len(sample_values):
        return "./."
    return sample_values[index]


def load_truth(path: Path) -> dict[tuple[str, int, str, str, str, int], str]:
    truth: dict[tuple[str, int, str, str, str, int], str] = {}
    with _open_text(path, "rt") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            key = canonical_key(fields)
            if key in truth:
                raise ChallengePanelError(f"duplicate truth key: {key!r}")
            truth[key] = _truth_gt(fields)
    if not truth:
        raise ChallengePanelError("truth VCF contains no records")
    return truth


def _candidate_id(pangenome_allele_id: str, seed: str) -> str:
    digest = hashlib.sha256(f"{seed}\0{pangenome_allele_id}".encode()).hexdigest()
    return f"CAND_{digest[:16]}"


def build_challenge_panel(
    *,
    panel_vcf: Path,
    truth_vcf: Path,
    output_vcf: Path,
    hidden_ledger: Path,
    audit_json: Path,
    seed: str,
) -> dict[str, int | str]:
    truth = load_truth(truth_vcf)
    output_vcf.parent.mkdir(parents=True, exist_ok=True)
    hidden_ledger.parent.mkdir(parents=True, exist_ok=True)
    audit_json.parent.mkdir(parents=True, exist_ok=True)

    counts = {
        "candidate_count": 0,
        "truth_positive_count": 0,
        "truth_negative_count": 0,
    }
    seen_candidate_ids: set[str] = set()
    saw_column_header = False

    with (
        _open_text(panel_vcf, "rt") as source,
        _open_text(output_vcf, "wt") as output,
        hidden_ledger.open("w", encoding="utf-8", newline="") as ledger_handle,
    ):
        ledger = csv.DictWriter(
            ledger_handle,
            fieldnames=[
                "candidate_id",
                "pangenome_allele_id",
                "truth_gt",
                "truth_label",
                "match_method",
            ],
            delimiter="\t",
            lineterminator="\n",
        )
        ledger.writeheader()

        for line in source:
            if line.startswith("##"):
                if "TRUTH" in line.upper():
                    continue
                output.write(line)
                continue
            if line.startswith("#CHROM"):
                header_fields = line.rstrip("\n").split("\t")
                if len(header_fields) < 8:
                    raise ChallengePanelError("panel VCF has an invalid #CHROM header")
                # PanGenie needs panel haplotypes through pangenome_panel, but
                # those samples must not leak into the shared blinded universe.
                output.write(
                    "\t".join(header_fields[:8] + ["FORMAT", "HG002"]) + "\n"
                )
                saw_column_header = True
                continue
            if not line.strip():
                continue
            if not saw_column_header:
                raise ChallengePanelError("panel VCF is missing #CHROM header")

            fields = line.rstrip("\n").split("\t")
            fields = fields[:8]
            info = parse_info(fields[7])
            allele_id = info.get("PANGENOME_ALLELE_ID")
            if not isinstance(allele_id, str) or not allele_id:
                raise ChallengePanelError("panel record is missing PANGENOME_ALLELE_ID")
            candidate_id = _candidate_id(allele_id, seed)
            if candidate_id in seen_candidate_ids:
                raise ChallengePanelError(f"candidate ID collision: {candidate_id}")
            seen_candidate_ids.add(candidate_id)

            key = canonical_key(fields)
            truth_gt = truth.get(key, "0/0")
            positive = truth_gt not in {"0/0", "0|0", "./.", ".|."}
            truth_label = "positive" if positive else "negative"
            match_method = "synthetic_exact" if key in truth else "no_exact_match"
            counts["candidate_count"] += 1
            counts["truth_positive_count" if positive else "truth_negative_count"] += 1

            for leaked_key in list(info):
                if leaked_key.upper().startswith("TRUTH"):
                    del info[leaked_key]
            fields[2] = candidate_id
            fields[7] = format_info(info)
            output.write("\t".join(fields + ["GT", "./."]) + "\n")
            ledger.writerow(
                {
                    "candidate_id": candidate_id,
                    "pangenome_allele_id": allele_id,
                    "truth_gt": truth_gt,
                    "truth_label": truth_label,
                    "match_method": match_method,
                }
            )

    if counts["candidate_count"] == 0:
        raise ChallengePanelError("panel contains no candidates")
    if counts["truth_positive_count"] == 0:
        raise ChallengePanelError("challenge panel has no truth-positive candidate")
    if counts["truth_negative_count"] == 0:
        raise ChallengePanelError("challenge panel has no truth-negative candidate")

    audit: dict[str, int | str] = {
        **counts,
        "matcher_profile": "synthetic_exact_v1",
        "seed_sha256": hashlib.sha256(seed.encode("utf-8")).hexdigest(),
        "truth_labels_exposed_to_tool": 0,
    }
    temporary = audit_json.with_name(f".{audit_json.name}.tmp")
    temporary.write_text(
        json.dumps(audit, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    temporary.replace(audit_json)
    return audit


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build a blinded synthetic genotype challenge panel."
    )
    parser.add_argument("--panel-vcf", required=True, type=Path)
    parser.add_argument("--truth-vcf", required=True, type=Path)
    parser.add_argument("--output-vcf", required=True, type=Path)
    parser.add_argument("--hidden-ledger", required=True, type=Path)
    parser.add_argument("--audit-json", required=True, type=Path)
    parser.add_argument("--seed", required=True)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        build_challenge_panel(
            panel_vcf=args.panel_vcf,
            truth_vcf=args.truth_vcf,
            output_vcf=args.output_vcf,
            hidden_ledger=args.hidden_ledger,
            audit_json=args.audit_json,
            seed=args.seed,
        )
    except (OSError, ChallengePanelError) as exc:
        print(f"build_challenge_panel: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
