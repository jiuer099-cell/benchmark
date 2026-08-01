#!/usr/bin/env python3
"""Run one frozen formal evaluator and emit an auditable semantic ledger."""

from __future__ import annotations

import argparse
import bisect
import csv
import gzip
import hashlib
import json
import os
import shutil
import subprocess
import sys
import tempfile
from collections import defaultdict
from pathlib import Path
from typing import Any, TextIO

import yaml  # type: ignore[import-untyped]

from sv_matching import (
    SvMatchError,
    SvRecord,
    genotypes_equal,
    gt_state,
    load_evaluator_profile,
    load_vcf,
    one_to_one_match,
    query_record_in_universe,
    truth_record_in_universe,
)


class FormalEvaluatorError(RuntimeError):
    """Raised when an evaluator cannot produce a complete semantic ledger."""


LEDGER_FIELDS = [
    "result_id",
    "evaluator",
    "scope_eligible",
    "event_eligible",
    "truth_event_id",
    "evaluator_accept",
    "detection_correct",
    "genotype_scorable",
    "genotype_correct",
    "no_call",
    "gt_state",
    "query_svtype",
    "query_svlen",
    "match_score",
    "start_distance",
    "end_distance",
    "size_similarity",
    "sequence_similarity",
    "correct",
]


def open_text(path: Path) -> TextIO:
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("r", encoding="utf-8")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def vcf_records(path: Path) -> list[list[str]]:
    records: list[list[str]] = []
    with open_text(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                raise FormalEvaluatorError(f"malformed VCF record in {path}")
            records.append(fields)
    return records


def record_key(fields: list[str]) -> tuple[str, int, str, str]:
    return fields[0], int(fields[1]), fields[3], fields[4]


def query_universe(path: Path) -> tuple[list[str], dict[str, tuple[str, int, str, str]]]:
    order: list[str] = []
    keys: dict[str, tuple[str, int, str, str]] = {}
    for fields in vcf_records(path):
        result_id = fields[2]
        if not result_id or result_id == ".":
            raise FormalEvaluatorError("canonical query VCF contains an empty ID")
        if result_id in keys:
            raise FormalEvaluatorError(f"duplicate canonical query ID: {result_id}")
        order.append(result_id)
        keys[result_id] = record_key(fields)
    return order, keys


def _resolve_votes(
    query: Path,
    decisions_by_id: dict[str, bool],
    decisions_by_key: dict[tuple[str, int, str, str], list[bool]],
) -> tuple[list[str], dict[str, bool]]:
    order, keys = query_universe(query)
    votes: dict[str, bool] = {}
    for result_id in order:
        if result_id in decisions_by_id:
            votes[result_id] = decisions_by_id[result_id]
            continue
        key_votes = decisions_by_key.get(keys[result_id], [])
        if not key_votes:
            raise FormalEvaluatorError(
                f"evaluator did not classify query result {result_id}"
            )
        if len(set(key_votes)) != 1:
            raise FormalEvaluatorError(
                f"evaluator produced conflicting classifications for {result_id}"
            )
        votes[result_id] = key_votes[0]
    return order, votes


def parse_truvari(query: Path, artifacts: Path) -> tuple[list[str], dict[str, bool]]:
    by_id: dict[str, bool] = {}
    by_key: dict[tuple[str, int, str, str], list[bool]] = defaultdict(list)
    artifact_candidates = (
        (("tp-comp.vcf.gz", "tp-comp.vcf", "tp-call.vcf.gz", "tp-call.vcf"), True),
        (("fp.vcf.gz", "fp.vcf"), False),
    )
    for names, correct in artifact_candidates:
        path = next(
            (artifacts / name for name in names if (artifacts / name).is_file()),
            None,
        )
        if path is None:
            searched = ", ".join(str(artifacts / name) for name in names)
            raise FormalEvaluatorError(
                f"Truvari output is absent; searched: {searched}"
            )
        for fields in vcf_records(path):
            if fields[2] and fields[2] != ".":
                by_id[fields[2]] = correct
            by_key[record_key(fields)].append(correct)
    return _resolve_votes(query, by_id, by_key)


def _format_value(fields: list[str], name: str) -> str:
    if len(fields) < 10:
        raise FormalEvaluatorError("evaluator VCF has no sample column")
    keys = fields[8].split(":")
    values = fields[9].split(":")
    if name not in keys:
        raise FormalEvaluatorError(f"evaluator VCF is missing FORMAT/{name}")
    index = keys.index(name)
    if index >= len(values):
        raise FormalEvaluatorError(f"evaluator VCF has no FORMAT/{name} value")
    return values[index]


def parse_aardvark(query: Path, artifacts: Path) -> tuple[list[str], dict[str, bool]]:
    path = artifacts / "query.vcf.gz"
    if not path.is_file():
        raise FormalEvaluatorError(f"Aardvark output is absent: {path}")
    by_id: dict[str, bool] = {}
    by_key: dict[tuple[str, int, str, str], list[bool]] = defaultdict(list)
    for fields in vcf_records(path):
        decision = _format_value(fields, "BD")
        if decision not in {"TP", "FP"}:
            raise FormalEvaluatorError(f"unexpected Aardvark BD={decision!r}")
        correct = decision == "TP"
        if fields[2] and fields[2] != ".":
            by_id[fields[2]] = correct
        by_key[record_key(fields)].append(correct)
    return _resolve_votes(query, by_id, by_key)


def parse_vcfdist(
    query: Path,
    artifacts: Path,
    *,
    credit_threshold: float,
) -> tuple[list[str], dict[str, bool]]:
    path = artifacts / "query.tsv"
    if not path.is_file():
        raise FormalEvaluatorError(f"vcfdist output is absent: {path}")
    credits: dict[tuple[str, int, str, str], list[float]] = defaultdict(list)
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"CONTIG", "POS", "REF", "ALT", "CREDIT"}
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise FormalEvaluatorError("vcfdist query.tsv has an unsupported header")
        for row in reader:
            key = (
                row["CONTIG"],
                int(row["POS"]) + 1,
                row["REF"],
                row["ALT"],
            )
            credits[key].append(float(row["CREDIT"]))
    decisions = {
        key: [sum(values) / len(values) >= credit_threshold]
        for key, values in credits.items()
    }
    return _resolve_votes(query, {}, decisions)


def load_regions(path: Path) -> dict[str, tuple[list[int], list[int]]]:
    raw: dict[str, list[tuple[int, int]]] = defaultdict(list)
    with path.open("r", encoding="utf-8") as handle:
        for number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                raise FormalEvaluatorError(f"malformed BED record at {path}:{number}")
            start, end = int(fields[1]), int(fields[2])
            if start < 0 or end <= start:
                raise FormalEvaluatorError(f"invalid BED interval at {path}:{number}")
            raw[fields[0]].append((start, end))
    merged: dict[str, tuple[list[int], list[int]]] = {}
    for chrom, intervals in raw.items():
        compact: list[list[int]] = []
        for start, end in sorted(intervals):
            if compact and start <= compact[-1][1]:
                compact[-1][1] = max(compact[-1][1], end)
            else:
                compact.append([start, end])
        merged[chrom] = (
            [interval[0] for interval in compact],
            [interval[1] for interval in compact],
        )
    if not merged:
        raise FormalEvaluatorError(f"benchmark BED contains no intervals: {path}")
    return merged


def fully_contained(
    regions: dict[str, tuple[list[int], list[int]]], record: SvRecord
) -> bool:
    indexed = regions.get(record.chrom)
    if indexed is None:
        return False
    starts, ends = indexed
    start = record.pos - 1
    end = max(record.pos, record.end)
    index = bisect.bisect_right(starts, start) - 1
    return index >= 0 and ends[index] >= end


def write_variant_query(source: Path, destination: Path, allowed_ids: set[str]) -> None:
    with open_text(source) as reader, destination.open("w", encoding="utf-8") as writer:
        for line in reader:
            if line.startswith("#"):
                writer.write(line)
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 3 and fields[2] in allowed_ids:
                writer.write(line)


def command_prefix(config: dict[str, Any], evaluator: str) -> list[str]:
    value = config["evaluation"]["commands"][evaluator]
    if not isinstance(value, list) or not all(
        isinstance(item, str) and item for item in value
    ):
        raise FormalEvaluatorError(f"invalid command prefix for {evaluator}")
    command = list(value)
    # `mamba run -n ...` may contend on shared wrapper/lock state when the
    # three formal evaluators launch concurrently. Resolve the conventional
    # named-environment form to its immutable executable whenever available.
    if len(command) >= 5 and command[:3] == ["mamba", "run", "-n"]:
        environment_name = command[3]
        executable = command[4]
        roots = []
        configured_root = os.environ.get("MAMBA_ROOT_PREFIX")
        if configured_root:
            roots.append(Path(configured_root))
        roots.append(Path.home() / ".local" / "share" / "mamba")
        for root in roots:
            candidate = root / "envs" / environment_name / "bin" / executable
            if candidate.is_file():
                return [str(candidate.resolve()), *command[5:]]
    return command


def build_command(
    evaluator: str,
    prefix: list[str],
    extra_args: list[str],
    *,
    query: Path,
    truth: Path,
    reference: Path,
    regions: Path,
    artifacts: Path,
    threads: int,
) -> list[str]:
    if evaluator == "truvari":
        return [
            *prefix,
            "bench",
            "-b",
            str(truth),
            "-c",
            str(query),
            "-f",
            str(reference),
            "--includebed",
            str(regions),
            "-o",
            str(artifacts),
            *extra_args,
        ]
    if evaluator == "aardvark":
        return [
            *prefix,
            "compare",
            "--reference",
            str(reference),
            "--truth-vcf",
            str(truth),
            "--query-vcf",
            str(query),
            "--regions",
            str(regions),
            "--output-dir",
            str(artifacts),
            "--threads",
            str(threads),
            *extra_args,
        ]
    if evaluator == "vcfdist":
        return [
            *prefix,
            str(query),
            str(truth),
            str(reference),
            "-b",
            str(regions),
            "-p",
            str(artifacts) + os.sep,
            "--max-threads",
            str(threads),
            *extra_args,
        ]
    raise FormalEvaluatorError(f"unsupported evaluator: {evaluator}")


def evaluator_version(
    prefix: list[str], version_args: list[str]
) -> dict[str, Any]:
    command = [*prefix, *version_args]
    completed = subprocess.run(
        command,
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
    )
    text = completed.stdout.strip()
    if not text:
        raise FormalEvaluatorError("evaluator version command returned no text")
    return {
        "command": command,
        "output": text,
        "sha256": hashlib.sha256(text.encode("utf-8")).hexdigest(),
    }


def event_query_indices(
    queries: list[SvRecord],
    *,
    scope_eligible_ids: set[str],
    candidate_output_contract: str,
) -> list[int]:
    """Select submitted detection events without conflating GT and detection."""

    if candidate_output_contract not in {"all_sites", "variant_sites"}:
        raise FormalEvaluatorError(
            f"invalid candidate_output_contract: {candidate_output_contract!r}"
        )
    selected: list[int] = []
    for index, record in enumerate(queries):
        if record.record_id not in scope_eligible_ids:
            continue
        state = gt_state(record.gt)
        if candidate_output_contract == "variant_sites":
            if state == "hom_ref":
                raise FormalEvaluatorError(
                    "variant_sites output contains an in-scope hom-ref record: "
                    f"{record.record_id}"
                )
            # Discovery VCFs commonly omit GT. Presence in a variant-sites VCF
            # is the detection assertion; no-call remains a separate GT
            # diagnostic and must not erase the detection.
            selected.append(index)
        elif state == "variant":
            selected.append(index)
    return selected


def write_ledger(
    path: Path,
    evaluator: str,
    queries: list[SvRecord],
    truths: list[SvRecord],
    assignments: dict[int, Any],
    evaluator_votes: dict[str, bool],
    *,
    scope_eligible_ids: set[str],
    event_eligible_ids: set[str],
    require_phase: bool,
) -> dict[str, int]:
    counts = {
        "rows": 0,
        "scope_eligible": 0,
        "scope_excluded": 0,
        "events": 0,
        "detection_correct": 0,
        "genotype_scorable": 0,
        "genotype_correct": 0,
        "no_call": 0,
    }
    seen_truth: set[str] = set()
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=LEDGER_FIELDS,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for index, query in enumerate(queries):
            state = gt_state(query.gt)
            scope_eligible = query.record_id in scope_eligible_ids
            event_eligible = query.record_id in event_eligible_ids
            if event_eligible and not scope_eligible:
                raise FormalEvaluatorError(
                    f"out-of-scope query was marked as a detection: {query.record_id}"
                )
            match = assignments.get(index) if event_eligible else None
            truth = truths[match.truth_index] if match is not None else None
            truth_id = truth.record_id if truth is not None else ""
            if truth_id:
                if truth_id in seen_truth:
                    raise FormalEvaluatorError(
                        f"truth event was assigned more than once: {truth_id}"
                    )
                seen_truth.add(truth_id)
            accepted = evaluator_votes.get(query.record_id, False) if event_eligible else False
            detection_correct = bool(accepted and truth is not None)
            genotype_scorable = truth is not None and state != "no_call"
            genotype_correct = bool(
                genotype_scorable
                and genotypes_equal(
                    query.gt,
                    truth.gt,
                    require_phase=require_phase,
                )
            )
            no_call = scope_eligible and state == "no_call"
            writer.writerow(
                {
                    "result_id": query.record_id,
                    "evaluator": evaluator,
                    "scope_eligible": int(scope_eligible),
                    "event_eligible": int(event_eligible),
                    "truth_event_id": truth_id,
                    "evaluator_accept": int(accepted),
                    "detection_correct": int(detection_correct),
                    "genotype_scorable": int(genotype_scorable),
                    "genotype_correct": int(genotype_correct),
                    "no_call": int(no_call),
                    "gt_state": state,
                    "query_svtype": query.svtype,
                    "query_svlen": query.svlen,
                    "match_score": "" if match is None else f"{match.score:.12g}",
                    "start_distance": "" if match is None else match.start_distance,
                    "end_distance": "" if match is None else match.end_distance,
                    "size_similarity": (
                        "" if match is None else f"{match.size_similarity:.12g}"
                    ),
                    "sequence_similarity": (
                        ""
                        if match is None or match.sequence_similarity is None
                        else f"{match.sequence_similarity:.12g}"
                    ),
                    # Compatibility alias: the primary vote remains detection only.
                    "correct": int(detection_correct),
                }
            )
            counts["rows"] += 1
            counts["scope_eligible"] += int(scope_eligible)
            counts["scope_excluded"] += int(not scope_eligible)
            counts["events"] += int(event_eligible)
            counts["detection_correct"] += int(detection_correct)
            counts["genotype_scorable"] += int(genotype_scorable)
            counts["genotype_correct"] += int(genotype_correct)
            counts["no_call"] += int(no_call)
    return counts


def run_evaluator(args: argparse.Namespace) -> None:
    config = yaml.safe_load(args.config.read_text(encoding="utf-8"))
    if not isinstance(config, dict):
        raise FormalEvaluatorError("configuration must be a YAML mapping")
    profile = load_evaluator_profile(args.evaluator_profile)
    tool_manifest = yaml.safe_load(args.tool_manifest.read_text(encoding="utf-8"))
    if not isinstance(tool_manifest, dict) or not isinstance(
        tool_manifest.get("outputs"), dict
    ):
        raise FormalEvaluatorError("tool manifest has no outputs contract")
    candidate_output_contract = tool_manifest["outputs"].get(
        "candidate_output_contract"
    )
    evaluator_profile = profile["evaluators"][args.evaluator]
    prefix = command_prefix(config, args.evaluator)
    version = evaluator_version(prefix, list(evaluator_profile["version_args"]))
    regions = load_regions(args.regions)
    all_queries = load_vcf(args.query, prefix="query", allow_empty=True)
    queries = [
        record for record in all_queries if fully_contained(regions, record)
    ]
    all_truths = load_vcf(args.truth, prefix="truth")
    outside_truth = [
        record.record_id
        for record in all_truths
        if not fully_contained(regions, record)
    ]
    if outside_truth:
        raise FormalEvaluatorError(
            "materialized truth contains events not fully contained in the "
            f"benchmark BED: {', '.join(outside_truth[:5])}"
        )
    out_of_scope_truth = [
        record.record_id
        for record in all_truths
        if not truth_record_in_universe(record, profile)
    ]
    if out_of_scope_truth:
        raise FormalEvaluatorError(
            "materialized truth violates the frozen SV universe: "
            f"{', '.join(out_of_scope_truth[:5])}"
        )
    truths = all_truths
    if not truths:
        raise FormalEvaluatorError("truth VCF has no events in the benchmark regions")
    if len({record.record_id for record in queries}) != len(queries):
        raise FormalEvaluatorError("query result IDs must be unique")
    scope_eligible_ids = {
        record.record_id
        for record in queries
        if query_record_in_universe(record, profile)
    }
    detection_indices = event_query_indices(
        queries,
        scope_eligible_ids=scope_eligible_ids,
        candidate_output_contract=str(candidate_output_contract),
    )
    event_queries = [queries[index] for index in detection_indices]
    event_assignments = one_to_one_match(event_queries, truths, profile)
    assignments = {
        detection_indices[event_index]: match
        for event_index, match in event_assignments.items()
    }
    event_ids = {queries[index].record_id for index in detection_indices}

    args.output_dir.parent.mkdir(parents=True, exist_ok=True)
    work = Path(
        tempfile.mkdtemp(prefix=f".{args.output_dir.name}.", dir=args.output_dir.parent)
    )
    command: list[str] = []
    try:
        evaluator_votes: dict[str, bool] = {}
        if event_ids:
            plain = work / "input" / "query.events.vcf"
            plain.parent.mkdir(parents=True)
            write_variant_query(args.query, plain, event_ids)
            prepared = work / "input" / "query.events.vcf.gz"
            bcftools = command_prefix(config, "bcftools")
            tabix = command_prefix(config, "tabix")
            subprocess.run(
                [*bcftools, "view", "-Oz", "-o", str(prepared), str(plain)],
                check=True,
            )
            subprocess.run([*tabix, "-f", "-p", "vcf", str(prepared)], check=True)
            artifacts = work / "artifacts"
            command = build_command(
                args.evaluator,
                prefix,
                list(evaluator_profile.get("extra_args", [])),
                query=prepared,
                truth=args.truth,
                reference=args.reference,
                regions=args.regions,
                artifacts=artifacts,
                threads=args.threads,
            )
            with (work / "evaluator.log").open("w", encoding="utf-8") as log:
                subprocess.run(
                    command,
                    check=True,
                    stdout=log,
                    stderr=subprocess.STDOUT,
                    text=True,
                )
            if args.evaluator == "truvari":
                _, evaluator_votes = parse_truvari(prepared, artifacts)
            elif args.evaluator == "aardvark":
                _, evaluator_votes = parse_aardvark(prepared, artifacts)
            else:
                _, evaluator_votes = parse_vcfdist(
                    prepared,
                    artifacts,
                    credit_threshold=float(evaluator_profile["credit_threshold"]),
                )
        else:
            (work / "evaluator.log").write_text(
                "No non-reference query events in benchmark regions; evaluator skipped.\n",
                encoding="utf-8",
            )

        ledger = work / "votes.tsv"
        counts = write_ledger(
            ledger,
            args.evaluator,
            queries,
            truths,
            assignments,
            evaluator_votes,
            scope_eligible_ids=scope_eligible_ids,
            event_eligible_ids=event_ids,
            require_phase=bool(
                profile["semantics"]["genotype"].get("require_phase", False)
            ),
        )
        complete = {
            "schema_version": 3,
            "evaluator": args.evaluator,
            "command": command,
            "version": version,
            "evaluator_profile_id": profile["profile"]["id"],
            "evaluator_profile_sha256": profile["_sha256"],
            "evaluator_profile_path": str(args.evaluator_profile),
            "query_sha256": sha256_file(args.query),
            "truth_sha256": sha256_file(args.truth),
            "benchmark_bed_sha256": sha256_file(args.regions),
            "reference_sha256": sha256_file(args.reference),
            "vote_ledger_sha256": sha256_file(ledger),
            "matching": {
                "algorithm": (
                    "deterministic_max_cardinality_max_weight_one_to_one_v1"
                ),
                "truth_assignments": len(assignments),
                "truth_ids_unique": True,
            },
            "semantics": {
                "primary_vote": "detection",
                "candidate_output_contract": candidate_output_contract,
                "genotype_separate": True,
                "no_call_separate": True,
            },
            "universe": {
                **profile["universe"],
                "query_records_outside_regions": len(all_queries) - len(queries),
            },
            **counts,
        }
        (work / ".complete.json").write_text(
            json.dumps(complete, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        backup = args.output_dir.with_name(f".{args.output_dir.name}.previous")
        if backup.exists():
            shutil.rmtree(backup)
        if args.output_dir.exists():
            args.output_dir.replace(backup)
        try:
            work.replace(args.output_dir)
        except BaseException:
            if backup.exists() and not args.output_dir.exists():
                backup.replace(args.output_dir)
            raise
        shutil.rmtree(backup, ignore_errors=True)
    except BaseException:
        shutil.rmtree(work, ignore_errors=True)
        raise


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--evaluator",
        required=True,
        choices=["truvari", "aardvark", "vcfdist"],
    )
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--evaluator-profile", required=True, type=Path)
    parser.add_argument("--tool-manifest", required=True, type=Path)
    parser.add_argument("--query", required=True, type=Path)
    parser.add_argument("--truth", required=True, type=Path)
    parser.add_argument("--reference", required=True, type=Path)
    parser.add_argument("--regions", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--threads", type=int, default=4)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    try:
        run_evaluator(parse_args(argv))
    except (
        OSError,
        ValueError,
        KeyError,
        TypeError,
        subprocess.CalledProcessError,
        FormalEvaluatorError,
        SvMatchError,
    ) as exc:
        print(f"run_formal_evaluator: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
