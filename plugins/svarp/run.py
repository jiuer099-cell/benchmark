#!/usr/bin/env python3
"""PGBench external adapter for SVarp long-read pangenome SV discovery.

SVarp's native result is a collection of assembled SV allele sequences
("svtigs").  The adapter makes the conversion to the benchmark's VCF exchange
format explicit and reproducible: svtigs are aligned against the registered
linear reference and large indels are extracted directly from minimap2's
reference-oriented ``cs`` strings.  ``paftools call`` is intentionally not
used: it expects a set of genome-scale, non-overlapping assembly contigs and
therefore discards SVarp's local alternate-allele contigs.
"""

from __future__ import annotations

import os
import re
import shutil
import subprocess
from pathlib import Path
from typing import BinaryIO, Iterable


MINIMUM_SV_SIZE = 50
MINIMUM_MAPQ = 20
CS_OPERATION = re.compile(
    r"(?::(?P<match>\d+)|=(?P<equal>[A-Za-z]+)|\*(?P<sub>[A-Za-z]{2})|"
    r"\+(?P<insertion>[A-Za-z]+)|-(?P<deletion>[A-Za-z]+)|"
    r"~[A-Za-z]{2}(?P<skip>\d+)[A-Za-z]{2})"
)


def required(name: str) -> str:
    value = os.environ.get(name)
    if not value:
        raise RuntimeError(f"required environment variable is absent: {name}")
    return value


def run(command: list[str], *, cwd: Path, stdout: BinaryIO | None = None) -> None:
    """Run one immutable tool command and show it in the attempt log."""

    print("+ " + " ".join(command), flush=True)
    subprocess.run(command, cwd=cwd, stdout=stdout, check=True)


def convert_fastq_to_bgzip_fasta(source: Path, destination: Path, *, cwd: Path) -> None:
    """Convert a registered FASTQ to bgzip FASTA without materializing plain text."""

    print(f"+ seqtk seq -A {source} | bgzip -c > {destination}", flush=True)
    destination.parent.mkdir(parents=True, exist_ok=True)
    with destination.open("wb") as output:
        seqtk = subprocess.Popen(
            ["seqtk", "seq", "-A", str(source)],
            cwd=cwd,
            stdout=subprocess.PIPE,
        )
        assert seqtk.stdout is not None
        bgzip = subprocess.Popen(
            ["bgzip", "-c"], cwd=cwd, stdin=seqtk.stdout, stdout=output
        )
        seqtk.stdout.close()
        bgzip_code = bgzip.wait()
        seqtk_code = seqtk.wait()
    if seqtk_code != 0 or bgzip_code != 0:
        destination.unlink(missing_ok=True)
        raise subprocess.CalledProcessError(
            seqtk_code or bgzip_code,
            ["seqtk", "seq", "-A", str(source), "|", "bgzip", "-c"],
        )
    if not destination.is_file() or destination.stat().st_size == 0:
        raise RuntimeError("FASTQ-to-FASTA conversion produced an empty file")


def find_svtigs(output_dir: Path, sample: str) -> Path | None:
    """Return SVarp's supported unphased or phased svtig output, if present."""

    candidates = (
        output_dir / f"{sample}_svtigs.fa",
        output_dir / f"{sample}_svtigs_H1.fa",
        output_dir / f"{sample}_svtigs_H2.fa",
        output_dir / f"{sample}_svtigs_untagged.fa",
    )
    existing = [path for path in candidates if path.is_file() and path.stat().st_size]
    if not existing:
        return None
    if len(existing) == 1:
        return existing[0]
    merged = output_dir / f"{sample}_svtigs.merged.fa"
    with merged.open("wb") as output:
        for path in existing:
            with path.open("rb") as source:
                shutil.copyfileobj(source, output)
    return merged


def _header_from_fai(reference_fai: Path) -> Iterable[str]:
    for raw_line in reference_fai.read_text(encoding="utf-8").splitlines():
        fields = raw_line.split("\t")
        if len(fields) >= 2:
            yield f"##contig=<ID={fields[0]},length={fields[1]}>\n"


class IndexedFasta:
    """Minimal random-access FASTA reader using the registered .fai index."""

    def __init__(self, fasta: Path, fai: Path) -> None:
        self.fasta = fasta
        self.entries: dict[str, tuple[int, int, int]] = {}
        for raw_line in fai.read_text(encoding="utf-8").splitlines():
            fields = raw_line.split("\t")
            if len(fields) >= 5:
                self.entries[fields[0]] = (
                    int(fields[2]),
                    int(fields[3]),
                    int(fields[4]),
                )

    def base(self, contig: str, position: int) -> str | None:
        """Return one 1-based reference base, or None outside a contig."""

        entry = self.entries.get(contig)
        if entry is None or position < 1:
            return None
        offset, line_bases, line_width = entry
        zero_based = position - 1
        byte_offset = offset + (zero_based // line_bases) * line_width + (
            zero_based % line_bases
        )
        with self.fasta.open("rb") as handle:
            handle.seek(byte_offset)
            value = handle.read(1).decode("ascii").upper()
        return value if value in {"A", "C", "G", "T", "N"} else None


def _tag_value(tags: list[str], prefix: str) -> str | None:
    for tag in tags:
        if tag.startswith(prefix):
            return tag[len(prefix) :]
    return None


def _paf_events(
    paf_line: str,
    reference: IndexedFasta,
) -> Iterable[tuple[str, int, str, str, str, int]]:
    """Yield sequence-resolved INS/DEL calls from one primary minimap2 PAF row."""

    fields = paf_line.rstrip("\n").split("\t")
    if len(fields) < 12:
        return
    try:
        mapq = int(fields[11])
        target_position = int(fields[7])  # 0-based reference coordinate.
    except ValueError:
        return
    if mapq < MINIMUM_MAPQ or "tp:A:P" not in fields[12:]:
        return
    cs = _tag_value(fields[12:], "cs:Z:")
    if cs is None:
        return
    contig = fields[5]
    cursor = 0
    for match in CS_OPERATION.finditer(cs):
        if match.start() != cursor:
            return
        cursor = match.end()
        kind = match.lastgroup
        value = match.group(kind) if kind else ""
        if kind in {"match", "equal"}:
            target_position += int(value) if kind == "match" else len(value)
        elif kind == "sub":
            target_position += 1
        elif kind == "skip":
            target_position += int(value)
        elif kind == "insertion":
            if len(value) < MINIMUM_SV_SIZE or target_position == 0:
                continue
            # cs is in reference alignment orientation, including reverse PAF rows.
            anchor = reference.base(contig, target_position)
            if anchor is not None:
                yield (
                    contig,
                    target_position,
                    anchor,
                    anchor + value.upper(),
                    "INS",
                    len(value),
                )
        elif kind == "deletion":
            if len(value) >= MINIMUM_SV_SIZE and target_position > 0:
                anchor = reference.base(contig, target_position)
                if anchor is not None:
                    yield (
                        contig,
                        target_position,
                        anchor + value.upper(),
                        anchor,
                        "DEL",
                        -len(value),
                    )
            target_position += len(value)
    if cursor != len(cs):
        return


def write_discovery_vcf(
    paf: Path,
    destination: Path,
    *,
    sample: str,
    reference: Path,
    reference_fai: Path,
) -> int:
    """Write non-redundant, sequence-resolved SVs directly from svtig PAF rows."""

    destination.parent.mkdir(parents=True, exist_ok=True)
    emitted = 0
    with destination.open("w", encoding="utf-8") as output:
        output.write("##fileformat=VCFv4.2\n")
        output.write("##source=PGBench-SVarp-1.2.0+minimap2-cs-v1\n")
        output.write(
            '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">\n'
        )
        output.write(
            '##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant">\n'
        )
        output.write(
            '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">\n'
        )
        output.write(
            '##INFO=<ID=SVARP_SOURCE,Number=1,Type=String,Description="SVarp evidence source">\n'
        )
        output.write(
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype; SVarp discovery output does not genotype">\n'
        )
        for line in _header_from_fai(reference_fai):
            output.write(line)
        output.write(
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
            + sample
            + "\n"
        )
        indexed_reference = IndexedFasta(reference, reference_fai)
        seen: set[tuple[str, int, str, str]] = set()
        with paf.open("r", encoding="utf-8") as source:
            for raw_line in source:
                for chrom, pos, ref, alt, svtype, svlen in _paf_events(
                    raw_line, indexed_reference
                ):
                    key = (chrom, pos, ref, alt)
                    if key in seen:
                        continue
                    seen.add(key)
                    end = pos if svtype == "INS" else pos + len(ref) - 1
                    record = [
                        chrom,
                        str(pos),
                        f"SVARP_{emitted + 1:09d}",
                        ref,
                        alt,
                        ".",
                        "PASS",
                        f"SVTYPE={svtype};END={end};SVLEN={svlen};SVARP_SOURCE=svtig",
                        "GT",
                        "./.",
                    ]
                    output.write("\t".join(record) + "\n")
                    emitted += 1
    return emitted


def main() -> int:
    fastq = Path(required("PGBENCH_INPUT_FASTQ"))
    reference = Path(required("PGBENCH_REFERENCE_FASTA"))
    reference_fai = Path(required("PGBENCH_REFERENCE_INDEX"))
    graph_dir = Path(required("PGBENCH_GRAPH_DIR"))
    output_dir = Path(required("PGBENCH_OUTPUT_DIR"))
    output_vcf = Path(required("PGBENCH_OUTPUT_VCF"))
    sample = required("PGBENCH_SAMPLE_ID")
    threads = required("PGBENCH_THREADS")

    graph = graph_dir / "graph.gfa.gz"
    if not graph.is_file() or graph.stat().st_size == 0:
        raise RuntimeError(f"SVarp requires a frozen rGFA at {graph}")
    if not reference_fai.is_file():
        raise RuntimeError(f"SVarp adapter requires reference index {reference_fai}")

    work = output_dir / "work"
    work.mkdir(parents=True, exist_ok=True)
    reads_fasta = work / f"{sample}.clr.fasta.gz"
    reads_fai = Path(f"{reads_fasta}.fai")
    if not reads_fasta.is_file() or not reads_fai.is_file():
        reads_fasta.unlink(missing_ok=True)
        reads_fai.unlink(missing_ok=True)
        convert_fastq_to_bgzip_fasta(fastq, reads_fasta, cwd=work)
        run(["samtools", "faidx", str(reads_fasta)], cwd=work)

    gaf = work / f"{sample}.svarp.gaf"
    if not gaf.is_file() or gaf.stat().st_size == 0:
        with gaf.open("wb") as output:
            run(
                [
                    "minigraph",
                    "-cx",
                    "lr",
                    "--vc",
                    "-t",
                    threads,
                    str(graph),
                    str(reads_fasta),
                ],
                cwd=work,
                stdout=output,
            )

    svarp_dir = work / "svarp"
    svarp_dir.mkdir(exist_ok=True)
    svtigs = find_svtigs(svarp_dir, sample)
    if svtigs is None:
        run(
            [
                "svarp",
                "--gaf",
                str(gaf),
                "--graph",
                str(graph),
                "--fasta",
                str(reads_fasta),
                "--sample",
                sample,
                "--out",
                str(svarp_dir),
                "--threads",
                threads,
            ],
            cwd=work,
        )
        svtigs = find_svtigs(svarp_dir, sample)

    if svtigs is None:
        paf = work / f"{sample}.svtigs.paf"
        paf.write_text("", encoding="utf-8")
    else:
        paf = work / f"{sample}.svtigs.paf"
        with paf.open("wb") as output:
            run(
                [
                    "minimap2",
                    "-cx",
                    "asm20",
                    "--cs",
                    "-t",
                    threads,
                    str(reference),
                    str(svtigs),
                ],
                cwd=work,
                stdout=output,
            )

    count = write_discovery_vcf(
        paf,
        output_vcf,
        sample=sample,
        reference=reference,
        reference_fai=reference_fai,
    )
    print(
        f"SVarp discovery VCF records={count}; genotype is intentionally no-call",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
