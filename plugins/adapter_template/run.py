#!/usr/bin/env python3
"""PGBench generic adapter template — copy this directory, fill in YOUR tool.

The benchmark core knows nothing about your tool.  It only:
  1. injects benchmark-supplied inputs as PGBENCH_* environment variables
     (driven by the ``required_inputs`` you declare in ``tool.yaml``),
  2. executes this ``run.py`` inside the declared sandbox/environment, and
  3. scores the all-sites VCF you produce at ``PGBENCH_OUTPUT_VCF``.

Your only mandatory code change is ``run_tool()``.  Everything else — input
preflight, the all-sites output contract, no-call semantics — is provided so
you cannot accidentally violate the benchmark contract.

INPUT MODES (auto-detected from which environment variables the core injected):
  - shared_bam : the benchmark's frozen shared alignment (BAM+BAI), e.g.
                 Paragraph / GraphTyper2 style short-read graph genotypers.
  - short_fastq: paired Illumina FASTQ (R1/R2), e.g. vg / BayesTyper style.
  - long_fastq : one long-read FASTQ (PacBio CLR / HiFi, or the frozen ONT
                 R9.4.1 Guppy5 SUP pass dataset), e.g. long-read callers.

HARD RULES (violations invalidate the run):
  - Output must cover EVERY candidate in PGBENCH_CANDIDATE_VCF (all_sites).
  - A candidate your tool did not genotype is ``./.`` — never ``0/0``.
  - Never read truth VCFs, assemblies, relatives, or undeclared callsets.
  - If your tool builds genome-wide graphs in memory, CHUNK your input (see
    plugins/paragraph/README.md for the 389 GB OOM lesson) and keep peak RSS
    inside the memory you declare in the run config.
  - Never subset or re-filter the frozen ONT pass reads (read selection belongs
    to the benchmark: official_qscore_pass_only).
"""

from __future__ import annotations

import dataclasses
import json
import os
import subprocess
from pathlib import Path
from typing import Any, TextIO


# ---------------------------------------------------------------------------
# Environment plumbing (the stable core -> adapter transport contract)
# ---------------------------------------------------------------------------
def required(name: str) -> str:
    value = os.environ.get(name)
    if not value:
        raise RuntimeError(f"required environment variable is absent: {name}")
    return value


def optional(name: str) -> str | None:
    return os.environ.get(name) or None


def open_text(path: Path, mode: str) -> TextIO:
    import gzip
    return gzip.open(path, mode, encoding="utf-8") if path.name.endswith(".gz") else path.open(mode, encoding="utf-8")


# ---------------------------------------------------------------------------
# Input boundary helpers
# ---------------------------------------------------------------------------
def require_readonly_input(path: Path, label: str) -> None:
    """Check injection only; Core owns content validation and its cache.

    Do not add SHA256, gzip passes, read counting, `samtools quickcheck`, or
    generic read-group validation here. Those public-input semantics are
    performed once by Core before this external adapter starts.
    """
    if not path.exists() or (path.is_file() and path.stat().st_size == 0):
        raise RuntimeError(f"Core injected missing/empty validated {label}: {path}")


# ---------------------------------------------------------------------------
# All-sites projection — the output contract, provided for you
# ---------------------------------------------------------------------------
def _allele_id(info: str) -> str | None:
    for item in info.split(";"):
        if item.startswith("PANGENOME_ALLELE_ID="):
            return item.split("=", 1)[1]
    return None


def _genotype(fields: list[str]) -> str:
    if len(fields) < 10:
        return "./."
    names = fields[8].split(":")
    values = fields[9].split(":")
    if "GT" not in names or names.index("GT") >= len(values):
        return "./."
    value = values[names.index("GT")]
    return value if value not in {"", "."} else "./."


def project_all_sites(native: Path, candidates: Path, destination: Path, sample: str) -> tuple[int, int]:
    """Project your tool's native VCF onto every frozen candidate.

    Matching precedence: blinded candidate ID -> PANGENOME_ALLELE_ID ->
    (CHROM, POS, REF, ALT).  Unmatched candidates are emitted as ``./.``.
    Returns (matched, no_call) for your run log.
    """

    order: list[str] = []
    records: dict[str, list[str]] = {}
    by_allele: dict[str, str] = {}
    by_key: dict[tuple[str, str, str, str], str] = {}
    meta: list[str] = []
    with open_text(candidates, "rt") as handle:
        for line in handle:
            if line.startswith("##"):
                meta.append(line)
                continue
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8 or not fields[2]:
                raise RuntimeError("candidate panel contains a malformed record")
            candidate_id = fields[2]
            if candidate_id in records:
                raise RuntimeError(f"candidate panel duplicates {candidate_id}")
            order.append(candidate_id)
            records[candidate_id] = fields
            key = (fields[0], fields[1], fields[3], fields[4])
            if key in by_key:
                raise RuntimeError(f"ambiguous duplicate canonical allele key: {key}")
            by_key[key] = candidate_id
            allele_id = _allele_id(fields[7])
            if allele_id:
                if allele_id in by_allele:
                    raise RuntimeError(f"ambiguous duplicate pangenome allele ID: {allele_id}")
                by_allele[allele_id] = candidate_id
    if not order:
        raise RuntimeError("candidate panel contains no records")
    calls: dict[str, str] = {}
    trace_rows: list[tuple[int, str, str, str, str, str, str, str]] = []
    with open_text(native, "rt") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                raise RuntimeError("tool output contains a malformed record")
            candidate_id = fields[2] if fields[2] in records else None
            match_strategy = "candidate_id" if candidate_id is not None else None
            if candidate_id is None:
                allele_id = _allele_id(fields[7])
                candidate_id = by_allele.get(allele_id) if allele_id else None
                if candidate_id is not None:
                    match_strategy = "pangenome_allele_id"
            if candidate_id is None:
                candidate_id = by_key.get((fields[0], fields[1], fields[3], fields[4]))
                if candidate_id is not None:
                    match_strategy = "exact_variant_key"
            if candidate_id is None or match_strategy is None:
                # A tool-native panel/graph may contain calls outside the fixed
                # canonical scoring universe. Retain that fact in the trace,
                # but never manufacture a candidate link or score it.
                trace_rows.append(
                    (
                        line_number,
                        "",
                        "outside_canonical_universe",
                        fields[0],
                        fields[1],
                        fields[3],
                        fields[4],
                        _genotype(fields),
                    )
                )
                continue
            if candidate_id in calls:
                raise RuntimeError(f"tool output duplicates {candidate_id}")
            genotype = _genotype(fields)
            calls[candidate_id] = genotype
            trace_rows.append(
                (
                    line_number,
                    candidate_id,
                    match_strategy,
                    fields[0],
                    fields[1],
                    fields[3],
                    fields[4],
                    genotype,
                )
            )
    destination.parent.mkdir(parents=True, exist_ok=True)
    with open_text(destination, "wt") as output:
        output.write("##fileformat=VCFv4.2\n")
        for line in meta:
            if not line.startswith("##fileformat") and not line.startswith("##FORMAT=<ID=GT,"):
                output.write(line)
        output.write("##source=PGBench-adapter-template-all-sites\n")
        output.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        output.write(f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}\n")
        for candidate_id in order:
            output.write(
                "\t".join(records[candidate_id][:8] + ["GT", calls.get(candidate_id, "./.")]) + "\n"
            )
    trace_path = destination.parent / "tool-work" / "native-to-canonical-projection.tsv"
    trace_path.parent.mkdir(parents=True, exist_ok=True)
    with trace_path.open("w", encoding="utf-8", newline="\n") as trace:
        trace.write(
            "native_record_line\tcandidate_id\tmatch_strategy\tnative_chrom\t"
            "native_pos\tnative_ref\tnative_alt\tnative_gt\n"
        )
        for row in trace_rows:
            trace.write("\t".join(str(value) for value in row) + "\n")
    return len(calls), len(order) - len(calls)


# ---------------------------------------------------------------------------
# The part YOU implement
# ---------------------------------------------------------------------------
@dataclasses.dataclass(frozen=True)
class AdapterContext:
    input_mode: str                 # "shared_bam" | "short_fastq" | "long_fastq"
    sample: str
    threads: int
    reference: Path
    candidates: Path
    output_dir: Path
    output_vcf: Path
    work_dir: Path
    # Every declared Core input is keyed by its contract token, e.g.
    # ``context.inputs[\"pangenome_manifest\"]`` or
    # ``context.inputs[\"adapter_asset.my_index\"]``.
    inputs: dict[str, Path] = dataclasses.field(default_factory=dict)
    # populated per input mode:
    bam: Path | None = None
    bai: Path | None = None
    fastq_r1: Path | None = None
    fastq_r2: Path | None = None
    long_fastq: Path | None = None
    optional_inputs: dict[str, Path] = dataclasses.field(default_factory=dict)


def run_tool(context: AdapterContext) -> Path:
    """Invoke YOUR tool and return the path of its native VCF.

    TODO: replace the body of this function.

    Guidance:
      - Keep every artifact inside ``context.work_dir``.
      - Pass ``context.threads`` to your tool; do not hardcode thread counts.
      - If your tool builds genome-wide graphs in memory, chunk your input
        (one contig per chunk) and merge before projecting — see
        plugins/paragraph/run.py for a worked, memory-safe example.
      - If your tool needs a tool-specific index/asset, declare it in
        tool.yaml as ``adapter_asset.<name>`` so the core manifests it.
    """

    raise NotImplementedError(
        "adapter_template is a scaffold: implement run_tool() for your tool"
    )


# ---------------------------------------------------------------------------
# Wiring — you normally do not need to change anything below
# ---------------------------------------------------------------------------
def detect_input_mode() -> str:
    modes = {
        "shared_bam": bool(os.environ.get("PGBENCH_SHARED_ALIGNMENT_BAM")),
        "short_fastq": bool(os.environ.get("PGBENCH_INPUT_FASTQ_R1")),
        "long_fastq": bool(os.environ.get("PGBENCH_INPUT_LONG_READS_FASTQ")),
    }
    selected = [name for name, present in modes.items() if present]
    if len(selected) == 1:
        return selected[0]
    if not selected:
        raise RuntimeError(
            "no benchmark read input was injected; check supported_modes."
            "required_inputs in tool.yaml"
        )
    raise RuntimeError(
        "Core injected more than one primary read evidence mode: " + ", ".join(selected)
    )


def injected_inputs() -> dict[str, Path]:
    """Expose only the Core-resolved input set, including adapter assets."""

    manifest_path = Path(required("PGBENCH_RESOLVED_INPUTS"))
    try:
        payload: Any = json.loads(manifest_path.read_text(encoding="utf-8"))
        entries = payload["inputs"]
    except (OSError, ValueError, KeyError, TypeError) as exc:
        raise RuntimeError(f"cannot read Core resolved-input manifest: {manifest_path}") from exc
    if not isinstance(entries, list):
        raise RuntimeError("Core resolved-input manifest has no input list")
    inputs: dict[str, Path] = {}
    for entry in entries:
        if not isinstance(entry, dict):
            raise RuntimeError("Core resolved-input manifest contains a malformed input entry")
        name = entry.get("contract_name")
        variable = entry.get("environment_variable")
        if not isinstance(name, str) or not isinstance(variable, str):
            raise RuntimeError("Core resolved-input manifest lacks a transport binding")
        path = Path(required(variable))
        require_readonly_input(path, name)
        if name in inputs:
            raise RuntimeError(f"Core resolved-input manifest duplicates {name}")
        inputs[name] = path
    return inputs


def main() -> int:
    input_mode = detect_input_mode()
    sample = required("PGBENCH_SAMPLE_ID")
    threads = int(required("PGBENCH_THREADS"))
    output_dir = Path(required("PGBENCH_OUTPUT_DIR")).resolve()
    output_vcf = Path(required("PGBENCH_OUTPUT_VCF")).resolve()
    output_vcf.relative_to(output_dir)  # refuse to escape the run directory
    work_dir = output_dir / "tool-work"
    work_dir.mkdir(parents=True, exist_ok=True)
    context = AdapterContext(
        input_mode=input_mode,
        sample=sample,
        threads=threads,
        reference=Path(required("PGBENCH_REFERENCE_FASTA")),
        candidates=Path(required("PGBENCH_CANDIDATE_VCF")),
        output_dir=output_dir,
        output_vcf=output_vcf,
        work_dir=work_dir,
        inputs=injected_inputs(),
    )
    if input_mode == "shared_bam":
        context = dataclasses.replace(
            context,
            bam=Path(required("PGBENCH_SHARED_ALIGNMENT_BAM")),
            bai=Path(required("PGBENCH_SHARED_ALIGNMENT_BAI")),
        )
        require_readonly_input(context.bam, "shared BAM")  # type: ignore[arg-type]
        require_readonly_input(context.bai, "shared BAM index")  # type: ignore[arg-type]
    elif input_mode == "short_fastq":
        context = dataclasses.replace(
            context,
            fastq_r1=Path(required("PGBENCH_INPUT_FASTQ_R1")),
            fastq_r2=Path(required("PGBENCH_INPUT_FASTQ_R2")),
        )
        require_readonly_input(context.fastq_r1, "R1 FASTQ")  # type: ignore[arg-type]
        require_readonly_input(context.fastq_r2, "R2 FASTQ")  # type: ignore[arg-type]
    elif input_mode == "long_fastq":
        context = dataclasses.replace(
            context,
            long_fastq=Path(required("PGBENCH_INPUT_LONG_READS_FASTQ")),
        )
        require_readonly_input(context.long_fastq, "long-read FASTQ")  # type: ignore[arg-type]
    optional_variables = {
        "reference_index": "PGBENCH_REFERENCE_INDEX",
        "graph_assets": "PGBENCH_GRAPH_DIR",
        "tool_index": "PGBENCH_INDEX_DIR",
    }
    for name, variable in optional_variables.items():
        value = optional(variable)
        if value:
            context.optional_inputs[name] = Path(value)

    native = run_tool(context)
    matched, no_call = project_all_sites(native, context.candidates, context.output_vcf, sample)
    print(f"adapter-template projection: matched={matched} no_call={no_call}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
