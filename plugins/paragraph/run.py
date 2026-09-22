#!/usr/bin/env python3
"""PGBench adapter for Paragraph from a frozen shared BAM through all-sites VCF.

Memory contract (revised after the 2026-09-21 OOM incident): ``multigrmpy.py``
builds one pangenome graph per invocation and its peak RSS scales with the
number of graph sites it receives.  Feeding the whole 18,164-candidate genome
panel in a single invocation peaked at ~389 GB and was killed by the kernel.
This adapter now chunks the candidate panel (one contig per chunk, optionally
split further by a candidate-count budget) and genotypes chunks with bounded
concurrency, so peak memory scales with the largest chunk instead of the whole
genome.  Chunking never changes results: chunks are disjoint by contig and the
final projection deterministically covers every canonical candidate.
"""

from __future__ import annotations

import dataclasses
import json
import os
import subprocess
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Mapping, Sequence, TextIO

DEFAULT_MAX_CANDIDATES_PER_CHUNK = 1000
DEFAULT_CHUNK_PARALLELISM = 2


def required(name: str) -> str:
    value = os.environ.get(name)
    if not value:
        raise RuntimeError(f"required environment variable is absent: {name}")
    return value


def open_text(path: Path, mode: str) -> TextIO:
    import gzip
    return gzip.open(path, mode, encoding="utf-8") if path.name.endswith(".gz") else path.open(mode, encoding="utf-8")


def read_length_from_bam(path: Path) -> int:
    """Read a representative sequenced-read length without buffering the BAM.

    The previous implementation captured the full ``samtools view`` stdout in
    memory before scanning it; this version streams one line and stops.
    """

    process = subprocess.Popen(
        ["samtools", "view", str(path)],
        stdout=subprocess.PIPE,
        stderr=subprocess.DEVNULL,
        text=True,
    )
    assert process.stdout is not None
    try:
        for line in process.stdout:
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 10 and fields[9] not in {"", "*"}:
                return len(fields[9])
        raise RuntimeError("frozen shared BAM contains no sequenced reads")
    finally:
        process.kill()
        process.wait()


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


@dataclasses.dataclass(frozen=True)
class CandidatePanel:
    """Deterministic in-memory index of the frozen canonical candidate panel."""

    meta: tuple[str, ...]
    header: str
    order: tuple[str, ...]
    records: dict[str, list[str]]
    by_allele: dict[str, str]
    by_key: dict[tuple[str, str, str, str], str]
    contigs: tuple[str, ...]
    contig_records: dict[str, tuple[str, ...]]


@dataclasses.dataclass(frozen=True)
class FastaIndexEntry:
    """One random-access entry from a standard FASTA ``.fai`` index."""

    length: int
    offset: int
    bases_per_line: int
    bytes_per_line: int


def load_fasta_index(reference: Path) -> dict[str, FastaIndexEntry]:
    """Read the adjacent FAI without loading the multi-gigabase reference."""

    entries: dict[str, FastaIndexEntry] = {}
    index = Path(f"{reference}.fai")
    with index.open("rt", encoding="utf-8") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 5:
                raise RuntimeError(f"malformed FASTA index record: {line!r}")
            name = fields[0]
            if not name or name in entries:
                raise RuntimeError(f"duplicate or empty FASTA index contig: {name!r}")
            entries[name] = FastaIndexEntry(
                length=int(fields[1]),
                offset=int(fields[2]),
                bases_per_line=int(fields[3]),
                bytes_per_line=int(fields[4]),
            )
    if not entries:
        raise RuntimeError(f"FASTA index is empty: {index}")
    return entries


def reference_base(
    reference: Path,
    index: Mapping[str, FastaIndexEntry],
    contig: str,
    position: int,
) -> str:
    """Return one 1-based reference base using the frozen FAI geometry."""

    entry = index.get(contig)
    if entry is None or not 1 <= position <= entry.length:
        raise RuntimeError(f"reference position is unavailable: {contig}:{position}")
    zero_based = position - 1
    offset = (
        entry.offset
        + (zero_based // entry.bases_per_line) * entry.bytes_per_line
        + zero_based % entry.bases_per_line
    )
    with reference.open("rb") as handle:
        handle.seek(offset)
        base = handle.read(1).decode("ascii").upper()
    if base not in {"A", "C", "G", "T", "N"}:
        raise RuntimeError(f"invalid reference base at {contig}:{position}: {base!r}")
    return base


def paragraph_compatible_records(
    panel: CandidatePanel,
    reference: Path,
) -> tuple[dict[str, list[str]], int]:
    """Add a shared left anchor only where Paragraph rejects the VCF syntax.

    Paragraph 2.3 requires REF and ALT to begin with the same padding base.
    Some valid sequence-resolved canonical SV records do not use that encoding.
    Prefixing the preceding reference base is an equivalent VCF representation:
    it preserves the biological allele, stable candidate ID, END and SVLEN.
    The original canonical records remain untouched and are used for output
    projection and all downstream scoring.
    """

    index = load_fasta_index(reference)
    prepared: dict[str, list[str]] = {}
    changed = 0
    for candidate_id in panel.order:
        fields = list(panel.records[candidate_id])
        ref, alt = fields[3], fields[4]
        if (
            not ref
            or not alt
            or alt.startswith("<")
            or ref[0] == alt[0]
        ):
            prepared[candidate_id] = fields
            continue
        position = int(fields[1])
        if position <= 1:
            raise RuntimeError(
                f"Paragraph cannot left-anchor {candidate_id} at {fields[0]}:{position}"
            )
        anchor = reference_base(reference, index, fields[0], position - 1)
        fields[1] = str(position - 1)
        fields[3] = anchor + ref
        fields[4] = anchor + alt
        prepared[candidate_id] = fields
        changed += 1
    return prepared, changed


def load_candidate_panel(path: Path) -> CandidatePanel:
    order: list[str] = []
    records: dict[str, list[str]] = {}
    by_allele: dict[str, str] = {}
    by_key: dict[tuple[str, str, str, str], str] = {}
    meta: list[str] = []
    header = ""
    with open_text(path, "rt") as handle:
        for line in handle:
            if line.startswith("##"):
                meta.append(line)
                continue
            if line.startswith("#CHROM"):
                header = line.rstrip("\n")
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
    if not order or not header:
        raise RuntimeError("candidate panel contains no records")
    contig_records: dict[str, list[str]] = {}
    for candidate_id in order:
        contig_records.setdefault(records[candidate_id][0], []).append(candidate_id)
    return CandidatePanel(
        meta=tuple(meta),
        header=header,
        order=tuple(order),
        records=records,
        by_allele=by_allele,
        by_key=by_key,
        contigs=tuple(contig_records),
        contig_records={contig: tuple(ids) for contig, ids in contig_records.items()},
    )


def write_chunks(
    panel: CandidatePanel,
    chunk_dir: Path,
    max_per_chunk: int,
    contig_filter: frozenset[str] | None = None,
    records: Mapping[str, Sequence[str]] | None = None,
) -> list[Path]:
    """Write disjoint VCF chunks; each chunk is one contig (or a part of it)."""

    if contig_filter is not None:
        missing = sorted(contig_filter - set(panel.contigs))
        if missing:
            raise RuntimeError(f"requested contigs absent from candidate panel: {missing}")
    selected = panel.contigs if contig_filter is None else [
        contig for contig in panel.contigs if contig in contig_filter
    ]
    chunk_dir.mkdir(parents=True, exist_ok=True)
    chunk_records = panel.records if records is None else records
    chunks: list[Path] = []
    index = 0
    for contig in selected:
        ids = panel.contig_records[contig]
        for start in range(0, len(ids), max_per_chunk):
            part_ids = ids[start : start + max_per_chunk]
            index += 1
            path = chunk_dir / f"chunk-{index:03d}-{contig}-part{start // max_per_chunk + 1:03d}.vcf"
            with path.open("wt", encoding="utf-8") as out:
                out.writelines(panel.meta)
                out.write(panel.header + "\n")
                for candidate_id in part_ids:
                    out.write("\t".join(chunk_records[candidate_id]) + "\n")
            chunks.append(path)
    if not chunks:
        raise RuntimeError("chunking produced zero chunks")
    return chunks


def run_chunk(
    chunk: Path,
    manifest: Path,
    reference: Path,
    work: Path,
    threads: int,
    max_depth: int,
) -> Path | None:
    """Genotype one chunk with multigrmpy; return its native genotypes VCF."""

    out = work / f"native-{chunk.stem}"
    subprocess.run(
        [
            "multigrmpy.py", "-i", str(chunk), "-m", str(manifest),
            "-r", str(reference), "-o", str(out), "-t", str(threads),
            "-M", str(max_depth),
        ],
        check=True,
    )
    genotypes = out / "genotypes.vcf.gz"
    return genotypes if genotypes.is_file() else None


def project_all_sites(
    native_paths: Sequence[Path],
    panel: CandidatePanel,
    destination: Path,
    sample: str,
    pilot_contigs: Sequence[str] | None = None,
) -> tuple[int, int]:
    """Project native Paragraph genotypes (possibly chunked) onto every candidate."""

    calls: dict[str, str] = {}
    for native in native_paths:
        with open_text(native, "rt") as handle:
            for line in handle:
                if not line.strip() or line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 8:
                    raise RuntimeError("Paragraph output contains a malformed record")
                candidate_id = fields[2] if fields[2] in panel.records else None
                if candidate_id is None:
                    allele_id = _allele_id(fields[7])
                    candidate_id = panel.by_allele.get(allele_id) if allele_id else None
                if candidate_id is None:
                    candidate_id = panel.by_key.get((fields[0], fields[1], fields[3], fields[4]))
                if candidate_id is not None:
                    if candidate_id in calls:
                        raise RuntimeError(f"Paragraph output duplicates {candidate_id}")
                    calls[candidate_id] = _genotype(fields)
    destination.parent.mkdir(parents=True, exist_ok=True)
    with open_text(destination, "wt") as output:
        output.write("##fileformat=VCFv4.2\n")
        for line in panel.meta:
            if not line.startswith("##fileformat") and not line.startswith("##FORMAT=<ID=GT,"):
                output.write(line)
        output.write("##source=PGBench-Paragraph-all-sites-adapter\n")
        if pilot_contigs:
            output.write(f"##PGBENCH_Paragraph_PilotContigs={','.join(pilot_contigs)}\n")
        output.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        output.write(f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}\n")
        for candidate_id in panel.order:
            output.write(
                "\t".join(panel.records[candidate_id][:8] + ["GT", calls.get(candidate_id, "./.")]) + "\n"
            )
    return len(calls), len(panel.order) - len(calls)


def _positive_int_env(name: str, default: int) -> int:
    raw = os.environ.get(name)
    if not raw:
        return default
    value = int(raw)
    if value < 1:
        raise RuntimeError(f"{name} must be a positive integer, got {raw}")
    return value


def main() -> int:
    bam = Path(required("PGBENCH_SHARED_ALIGNMENT_BAM"))
    bai = Path(required("PGBENCH_SHARED_ALIGNMENT_BAI"))
    reference = Path(required("PGBENCH_REFERENCE_FASTA"))
    candidates = Path(required("PGBENCH_CANDIDATE_VCF"))
    sample = required("PGBENCH_SAMPLE_ID")
    threads = _positive_int_env("PGBENCH_THREADS", 1)
    output_dir = Path(required("PGBENCH_OUTPUT_DIR")).resolve()
    output_vcf = Path(required("PGBENCH_OUTPUT_VCF")).resolve()
    output_vcf.relative_to(output_dir)
    work = output_dir / "paragraph-work"
    work.mkdir(parents=True, exist_ok=True)

    contigs_raw = os.environ.get("PGBENCH_PARAGRAPH_CONTIGS")
    contig_filter = (
        frozenset(part.strip() for part in contigs_raw.split(",") if part.strip())
        if contigs_raw
        else None
    )
    max_per_chunk = _positive_int_env(
        "PGBENCH_PARAGRAPH_MAX_CANDIDATES_PER_CHUNK", DEFAULT_MAX_CANDIDATES_PER_CHUNK
    )
    parallelism = min(
        _positive_int_env("PGBENCH_PARAGRAPH_CHUNK_PARALLELISM", DEFAULT_CHUNK_PARALLELISM),
        threads,
    )

    if not bai.is_file() or bai.stat().st_size == 0:
        raise RuntimeError("Paragraph requires the frozen shared BAM index")
    subprocess.run(["samtools", "quickcheck", "-v", str(bam)], check=True)
    coverage = subprocess.check_output(
        ["samtools", "coverage", str(bam)], text=True
    ).splitlines()
    depths = [float(line.split("\t")[6]) for line in coverage if line and not line.startswith("#")]
    depth = sum(depths) / len(depths) if depths else 1.0
    max_depth = max(20, round(depth * 20))
    manifest = work / "sample.tsv"
    manifest.write_text(
        "id\tpath\tread length\tdepth\n"
        f"{sample}\t{bam}\t{read_length_from_bam(bam)}\t{depth:.6f}\n",
        encoding="utf-8",
    )

    panel = load_candidate_panel(candidates)
    prepared_records, left_anchored = paragraph_compatible_records(panel, reference)
    chunks = write_chunks(
        panel,
        work / "chunks",
        max_per_chunk,
        contig_filter,
        records=prepared_records,
    )
    threads_per_chunk = max(1, threads // parallelism)
    print(
        f"Paragraph chunked genotyping: candidates={len(panel.order)} "
        f"chunks={len(chunks)} max_per_chunk={max_per_chunk} "
        f"parallelism={parallelism} threads_per_chunk={threads_per_chunk} "
        f"pilot_contigs={','.join(sorted(contig_filter)) if contig_filter else 'none'} "
        f"left_anchored_records={left_anchored}"
    )
    with ThreadPoolExecutor(max_workers=parallelism) as pool:
        futures = [
            pool.submit(run_chunk, chunk, manifest, reference, work, threads_per_chunk, max_depth)
            for chunk in chunks
        ]
        native = [future.result() for future in futures]
    empty = [chunk for chunk, path in zip(chunks, native) if path is None]
    native_paths = [path for path in native if path is not None]
    if empty:
        print(f"WARNING: {len(empty)} chunk(s) produced no genotypes.vcf.gz: {[c.name for c in empty]}")

    pilot_contigs = sorted(contig_filter) if contig_filter else None
    matched, no_call = project_all_sites(native_paths, panel, output_vcf, sample, pilot_contigs)
    if contig_filter:
        marker = output_dir / "paragraph-pilot.json"
        marker.write_text(
            json.dumps(
                {
                    "pilot": True,
                    "genotyped_contigs": pilot_contigs,
                    "chunks": len(chunks),
                    "matched": matched,
                    "no_call": no_call,
                },
                sort_keys=True,
            )
            + "\n",
            encoding="utf-8",
        )
        print(
            "WARNING: PGBENCH_PARAGRAPH_CONTIGS is set, so this is a PILOT run. "
            "The output covers only the pilot contigs; it must not be used for a formal score."
        )
    print(f"Paragraph candidate projection: matched={matched} no_call={no_call}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
