#!/usr/bin/env python3
"""PGBench external adapter for SVarp long-read pangenome SV discovery.

SVarp's native result is a collection of assembled SV allele sequences
("svtigs").  The adapter makes the conversion to the benchmark's VCF exchange
format explicit and reproducible: svtigs are aligned against the registered
linear reference and `paftools.js call` provides sequence-resolved INS/DEL
records.  These are discovery records, not invented genotypes.
"""

from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path
from typing import BinaryIO, Iterable


MINIMUM_SV_SIZE = 50


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


def _sv_record(fields: list[str], ordinal: int) -> list[str] | None:
    """Turn one paftools VCF record into a benchmark discovery record."""

    if len(fields) < 8 or "," in fields[4]:
        return None
    chrom, position, _, ref, alt, _, filter_value, _ = fields[:8]
    if ref in {"", "."} or alt in {"", "."} or alt.startswith("<"):
        return None
    try:
        pos = int(position)
    except ValueError:
        return None
    difference = len(alt) - len(ref)
    if abs(difference) < MINIMUM_SV_SIZE:
        return None
    if difference > 0:
        svtype, end, svlen = "INS", pos, difference
    else:
        svtype, end, svlen = "DEL", pos + len(ref) - 1, difference
    return [
        chrom,
        str(pos),
        f"SVARP_{ordinal:09d}",
        ref,
        alt,
        ".",
        filter_value if filter_value in {"PASS", "."} else "PASS",
        f"SVTYPE={svtype};END={end};SVLEN={svlen};SVARP_SOURCE=svtig",
        "GT",
        "./.",
    ]


def write_discovery_vcf(
    raw_vcf: Path,
    destination: Path,
    *,
    sample: str,
    reference_fai: Path,
) -> int:
    """Filter paftools VCF to sequence-resolved structural discovery records."""

    destination.parent.mkdir(parents=True, exist_ok=True)
    emitted = 0
    with destination.open("w", encoding="utf-8") as output:
        output.write("##fileformat=VCFv4.2\n")
        output.write("##source=PGBench-SVarp-1.2.0+minimap2-2.31-paftools\n")
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
        for raw_line in raw_vcf.read_text(encoding="utf-8").splitlines():
            if not raw_line or raw_line.startswith("#"):
                continue
            record = _sv_record(raw_line.split("\t"), emitted + 1)
            if record is None:
                continue
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

    raw_vcf = work / f"{sample}.svtigs.paftools.vcf"
    if svtigs is None:
        raw_vcf.write_text("##fileformat=VCFv4.2\n", encoding="utf-8")
    else:
        paf = work / f"{sample}.svtigs.paf"
        sorted_paf = work / f"{sample}.svtigs.sorted.paf"
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
        with sorted_paf.open("wb") as output:
            run(["sort", "-k6,6", "-k8,8n", str(paf)], cwd=work, stdout=output)
        with raw_vcf.open("wb") as output:
            run(
                ["paftools.js", "call", "-f", str(reference), str(sorted_paf)],
                cwd=work,
                stdout=output,
            )

    count = write_discovery_vcf(
        raw_vcf,
        output_vcf,
        sample=sample,
        reference_fai=reference_fai,
    )
    print(
        f"SVarp discovery VCF records={count}; genotype is intentionally no-call",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
