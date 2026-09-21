#!/usr/bin/env python3
"""Build a deterministic, non-production PG-F1 replay panel from frozen assets."""
from __future__ import annotations

import argparse
import csv
from collections import Counter
from pathlib import Path

from pgf1 import PGF1Error, canonical_alt_dosage, json_dump, open_text, sha256_file


def parse_info(value: str) -> dict[str, str]:
    return dict(item.split("=", 1) for item in value.split(";") if "=" in item)


def record_end(fields: list[str]) -> int:
    info = parse_info(fields[7])
    return max(int(fields[1]) + len(fields[3]) - 1, int(info.get("END", fields[1])))


def genotype_class(fields: list[str]) -> str:
    keys = fields[8].split(":") if len(fields) > 8 else []
    values = fields[9].split(":") if len(fields) > 9 else []
    gt = values[keys.index("GT")] if "GT" in keys and keys.index("GT") < len(values) else "."
    dosage = canonical_alt_dosage(gt)
    if dosage is None:
        return "no_call"
    return "positive" if dosage > 0 else "reference"


def read_vcf(path: Path) -> tuple[list[str], list[list[str]]]:
    headers: list[str] = []
    records: list[list[str]] = []
    with open_text(path, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                headers.append(line)
            elif line.strip():
                records.append(line.rstrip("\r\n").split("\t"))
    return headers, records


def write_vcf(path: Path, headers: list[str], records: list[list[str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open_text(path, "wt") as handle:
        handle.writelines(headers)
        for fields in records:
            handle.write("\t".join(fields) + "\n")


def compact_selection(records: list[dict[str, object]], targets: dict[str, int]) -> list[dict[str, object]]:
    have: Counter[str] = Counter()
    left = 0
    best: tuple[int, int, int] | None = None
    for right, row in enumerate(records):
        have[str(row["class"])] += 1
        while all(have[key] >= value for key, value in targets.items()):
            span = int(records[right]["end"]) - int(records[left]["pos"])
            candidate = (span, left, right)
            if best is None or candidate < best:
                best = candidate
            have[str(records[left]["class"])] -= 1
            left += 1
    if best is None:
        raise PGF1Error("small replay selection cannot satisfy requested genotype strata")
    selected: list[dict[str, object]] = []
    remaining = dict(targets)
    for row in records[best[1] : best[2] + 1]:
        category = str(row["class"])
        if remaining[category] > 0:
            selected.append(row)
            remaining[category] -= 1
    if any(remaining.values()):
        raise PGF1Error("small replay selection is incomplete")
    return selected


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    for name in ("canonical-panel", "hidden-truth-ledger", "linked-vcf", "all-sites-vcf", "truth-vcf"):
        parser.add_argument(f"--{name}", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--chrom", default="chr1")
    parser.add_argument("--positive-units", type=int, default=20)
    parser.add_argument("--reference-units", type=int, default=5)
    parser.add_argument("--no-call-units", type=int, default=5)
    parser.add_argument("--padding-bp", type=int, default=500)
    args = parser.parse_args(argv)
    try:
        args.output_dir.mkdir(parents=True, exist_ok=True)
        with args.hidden_truth_ledger.open(encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            hidden_fields = reader.fieldnames
            hidden = {row["candidate_id"]: row for row in reader}
        if not hidden_fields:
            raise PGF1Error("invalid hidden-truth ledger")

        linked_headers, linked_records = read_vcf(args.linked_vcf)
        eligible: list[dict[str, object]] = []
        for fields in linked_records:
            info = parse_info(fields[7])
            candidate_id = info.get("PANGENOME_LINKED_ID", "")
            if (
                fields[0] != args.chrom
                or info.get("PANGENOME_CANDIDATE_COUNT") != "1"
                or hidden.get(candidate_id, {}).get("truth_scorable") != "1"
            ):
                continue
            eligible.append(
                {
                    "candidate_id": candidate_id,
                    "class": genotype_class(fields),
                    "pos": int(fields[1]),
                    "end": record_end(fields),
                    "fields": fields,
                }
            )
        eligible.sort(key=lambda row: (int(row["pos"]), str(row["candidate_id"])))
        targets = {
            "positive": args.positive_units,
            "reference": args.reference_units,
            "no_call": args.no_call_units,
        }
        selected = compact_selection(eligible, targets)
        selected_ids = {str(row["candidate_id"]) for row in selected}
        region_start = max(0, min(int(row["pos"]) - 1 for row in selected) - args.padding_bp)
        region_end = max(int(row["end"]) for row in selected) + args.padding_bp

        panel_headers, panel_records = read_vcf(args.canonical_panel)
        panel_subset = [fields for fields in panel_records if fields[2] in selected_ids]
        all_headers, all_records = read_vcf(args.all_sites_vcf)
        all_subset = [fields for fields in all_records if fields[2] in selected_ids]
        linked_subset = [row["fields"] for row in selected]
        if not (len(panel_subset) == len(all_subset) == len(linked_subset) == len(selected_ids)):
            raise PGF1Error("UNIT_SET_MISMATCH: replay asset subsets do not reconcile")

        truth_headers, truth_records = read_vcf(args.truth_vcf)
        truth_subset = [
            fields
            for fields in truth_records
            if fields[0] == args.chrom
            and int(fields[1]) - 1 >= region_start
            and record_end(fields) <= region_end
        ]
        if not truth_subset:
            raise PGF1Error("small replay truth subset is empty")

        outputs = {
            "canonical_panel": args.output_dir / "canonical-panel.vcf",
            "hidden_truth": args.output_dir / "hidden-truth.tsv",
            "linked": args.output_dir / "linked.vcf",
            "all_sites": args.output_dir / "all-sites.vcf.gz",
            "truth": args.output_dir / "truth.vcf.gz",
            "bed": args.output_dir / "benchmark.bed",
        }
        write_vcf(outputs["canonical_panel"], panel_headers, panel_subset)
        write_vcf(outputs["linked"], linked_headers, linked_subset)
        write_vcf(outputs["all_sites"], all_headers, all_subset)
        write_vcf(outputs["truth"], truth_headers, truth_subset)
        with outputs["hidden_truth"].open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=hidden_fields, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            for candidate_id in sorted(selected_ids):
                writer.writerow(hidden[candidate_id])
        outputs["bed"].write_text(f"{args.chrom}\t{region_start}\t{region_end}\n", encoding="utf-8")

        counts = Counter(str(row["class"]) for row in selected)
        manifest = {
            "artifact_type": "pgf1_small_real_replay_panel",
            "production_eligible": False,
            "leaderboard_admissible": False,
            "selection_policy": "smallest_genomic_window_satisfying_fixed_gt_strata_v1",
            "chrom": args.chrom,
            "region_start_0based": region_start,
            "region_end_0based_exclusive": region_end,
            "selected_unit_count": len(selected_ids),
            "selected_genotype_counts": dict(sorted(counts.items())),
            "selected_candidate_ids": sorted(selected_ids),
            "truth_record_count": len(truth_subset),
            "source_sha256": {
                "canonical_panel": sha256_file(args.canonical_panel),
                "hidden_truth_ledger": sha256_file(args.hidden_truth_ledger),
                "linked_vcf": sha256_file(args.linked_vcf),
                "all_sites_vcf": sha256_file(args.all_sites_vcf),
                "truth_vcf": sha256_file(args.truth_vcf),
            },
            "output_sha256": {name: sha256_file(path) for name, path in outputs.items()},
        }
        json_dump(args.output_dir / "selection-manifest.json", manifest)
    except (OSError, ValueError, PGF1Error) as error:
        print(f"prepare_pgf1_small_replay: {error}", file=__import__("sys").stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
