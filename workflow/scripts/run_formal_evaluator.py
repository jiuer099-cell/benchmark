#!/usr/bin/env python3
"""Run one formal evaluator and emit one binary vote per query result."""

from __future__ import annotations

import argparse
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
from typing import TextIO

import yaml  # type: ignore[import-untyped]


class FormalEvaluatorError(RuntimeError):
    """Raised when an evaluator cannot produce a complete query vote ledger."""


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
    if not order:
        raise FormalEvaluatorError("canonical query VCF contains no results")
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
    for name, correct in (("tp-call.vcf.gz", True), ("fp.vcf.gz", False)):
        path = artifacts / name
        if not path.is_file():
            raise FormalEvaluatorError(f"Truvari output is absent: {path}")
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


def write_votes(
    path: Path,
    evaluator: str,
    order: list[str],
    votes: dict[str, bool],
) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["result_id", "evaluator", "correct"],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for result_id in order:
            writer.writerow(
                {
                    "result_id": result_id,
                    "evaluator": evaluator,
                    "correct": int(votes[result_id]),
                }
            )


def command_prefix(config: dict, evaluator: str) -> list[str]:
    value = config["evaluation"]["commands"][evaluator]
    if not isinstance(value, list) or not all(
        isinstance(item, str) and item for item in value
    ):
        raise FormalEvaluatorError(f"invalid command prefix for {evaluator}")
    return list(value)


def build_command(
    evaluator: str,
    prefix: list[str],
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
        ]
    raise FormalEvaluatorError(f"unsupported evaluator: {evaluator}")


def run_evaluator(args: argparse.Namespace) -> None:
    config = yaml.safe_load(args.config.read_text(encoding="utf-8"))
    if not isinstance(config, dict):
        raise FormalEvaluatorError("configuration must be a YAML mapping")
    prefix = command_prefix(config, args.evaluator)
    args.output_dir.parent.mkdir(parents=True, exist_ok=True)
    work = Path(
        tempfile.mkdtemp(
            prefix=f".{args.output_dir.name}.",
            dir=args.output_dir.parent,
        )
    )
    try:
        prepared = work / "input" / "query.vcf.gz"
        prepared.parent.mkdir(parents=True)
        bcftools = command_prefix(config, "bcftools")
        tabix = command_prefix(config, "tabix")
        subprocess.run(
            [*bcftools, "view", "-Oz", "-o", str(prepared), str(args.query)],
            check=True,
        )
        subprocess.run([*tabix, "-f", "-p", "vcf", str(prepared)], check=True)

        artifacts = work / "artifacts"
        command = build_command(
            args.evaluator,
            prefix,
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

        threshold = float(config["evaluation"]["vcfdist_credit_threshold"])
        if args.evaluator == "truvari":
            order, votes = parse_truvari(args.query, artifacts)
        elif args.evaluator == "aardvark":
            order, votes = parse_aardvark(args.query, artifacts)
        else:
            order, votes = parse_vcfdist(
                args.query,
                artifacts,
                credit_threshold=threshold,
            )
        ledger = work / "votes.tsv"
        write_votes(ledger, args.evaluator, order, votes)
        complete = {
            "schema_version": 1,
            "evaluator": args.evaluator,
            "command": command,
            "query_sha256": sha256_file(args.query),
            "truth_sha256": sha256_file(args.truth),
            "vote_ledger_sha256": sha256_file(ledger),
            "result_count": len(order),
            "correct_count": sum(votes.values()),
        }
        (work / ".complete.json").write_text(
            json.dumps(complete, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        # Snakemake may pre-create the parent of declared output files, which
        # is the evaluator directory itself. Commit the completed work tree
        # with a rollback directory so a failed rerun never exposes a partial
        # vote ledger.
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
    except (OSError, ValueError, KeyError, subprocess.CalledProcessError, FormalEvaluatorError) as exc:
        print(f"run_formal_evaluator: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
