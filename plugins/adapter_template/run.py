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
import os
import subprocess
from pathlib import Path
from typing import TextIO


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
# Input preflight helpers — reuse them, they encode benchmark policy
# ---------------------------------------------------------------------------
def preflight_shared_bam(bam: Path, bai: Path, sample: str) -> None:
    """Validate the frozen shared BAM and its read group (SM must match)."""
    if not bai.is_file() or bai.stat().st_size == 0:
        raise RuntimeError("shared alignment requires a non-empty BAI index")
    subprocess.run(["samtools", "quickcheck", "-v", str(bam)], check=True)
    header = subprocess.check_output(["samtools", "view", "-H", str(bam)], text=True)
    read_groups: dict[str, str] = {}
    for line in header.splitlines():
        if line.startswith("@RG\t"):
            tags = dict(
                part.split(":", 1) for part in line.split("\t")[1:] if ":" in part
            )
            if tags.get("ID"):
                read_groups[tags["ID"]] = tags.get("SM", "")
    if not read_groups:
        raise RuntimeError("shared BAM header declares no @RG read group")
    mismatched = {rg: sm for rg, sm in read_groups.items() if sm != sample}
    if mismatched:
        raise RuntimeError(f"@RG SM mismatch for sample {sample}: {mismatched}")


def preflight_fastq(path: Path) -> None:
    """Existence + gzip integrity for a (possibly compressed) FASTQ input."""
    if not path.is_file() or path.stat().st_size == 0:
        raise RuntimeError(f"input FASTQ is missing or empty: {path}")
    if path.name.endswith(".gz"):
        import gzip
        with gzip.open(path, "rb") as handle:
            while handle.read(1 << 20):
                pass  # full pass; raises gzip.BadGzipFile on truncation/corruption


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
            by_key[(fields[0], fields[1], fields[3], fields[4])] = candidate_id
            allele_id = _allele_id(fields[7])
            if allele_id:
                by_allele[allele_id] = candidate_id
    if not order:
        raise RuntimeError("candidate panel contains no records")
    calls: dict[str, str] = {}
    with open_text(native, "rt") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                raise RuntimeError("tool output contains a malformed record")
            candidate_id = fields[2] if fields[2] in records else None
            if candidate_id is None:
                allele_id = _allele_id(fields[7])
                candidate_id = by_allele.get(allele_id) if allele_id else None
            if candidate_id is None:
                candidate_id = by_key.get((fields[0], fields[1], fields[3], fields[4]))
            if candidate_id is not None:
                if candidate_id in calls:
                    raise RuntimeError(f"tool output duplicates {candidate_id}")
                calls[candidate_id] = _genotype(fields)
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
    if os.environ.get("PGBENCH_SHARED_ALIGNMENT_BAM"):
        return "shared_bam"
    if os.environ.get("PGBENCH_INPUT_FASTQ_R1"):
        return "short_fastq"
    if os.environ.get("PGBENCH_INPUT_LONG_READS_FASTQ"):
        return "long_fastq"
    raise RuntimeError(
        "no benchmark read input was injected; check supported_modes."
        "required_inputs in tool.yaml"
    )


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
    )
    if input_mode == "shared_bam":
        context = dataclasses.replace(
            context,
            bam=Path(required("PGBENCH_SHARED_ALIGNMENT_BAM")),
            bai=Path(required("PGBENCH_SHARED_ALIGNMENT_BAI")),
        )
        preflight_shared_bam(context.bam, context.bai, sample)  # type: ignore[arg-type]
    elif input_mode == "short_fastq":
        context = dataclasses.replace(
            context,
            fastq_r1=Path(required("PGBENCH_INPUT_FASTQ_R1")),
            fastq_r2=Path(required("PGBENCH_INPUT_FASTQ_R2")),
        )
        preflight_fastq(context.fastq_r1)  # type: ignore[arg-type]
        preflight_fastq(context.fastq_r2)  # type: ignore[arg-type]
    elif input_mode == "long_fastq":
        context = dataclasses.replace(
            context,
            long_fastq=Path(required("PGBENCH_INPUT_LONG_READS_FASTQ")),
        )
        preflight_fastq(context.long_fastq)  # type: ignore[arg-type]
    for name in ("reference_index", "graph_assets", "tool_index"):
        value = optional(f"PGBENCH_{name.upper()}")
        if value:
            context.optional_inputs[name] = Path(value)

    native = run_tool(context)
    matched, no_call = project_all_sites(native, context.candidates, context.output_vcf, sample)
    print(f"adapter-template projection: matched={matched} no_call={no_call}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
