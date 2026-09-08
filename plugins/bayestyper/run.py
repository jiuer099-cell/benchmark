#!/usr/bin/env python3
"""BayesTyper 1.5 adapter for the fixed short-read candidate panel."""

from __future__ import annotations

import gzip
import os
import subprocess
from pathlib import Path
from typing import TextIO


def required(name: str) -> str:
    value = os.environ.get(name)
    if not value:
        raise RuntimeError(f"required environment variable is absent: {name}")
    return value


def run(command: list[str], *, cwd: Path | None = None) -> None:
    print("+ " + " ".join(command), flush=True)
    subprocess.run(command, check=True, cwd=cwd)


def open_text(path: Path, mode: str) -> TextIO:
    if path.name.endswith(".gz"):
        return gzip.open(path, mode, encoding="utf-8")
    return path.open(mode, encoding="utf-8")


def allele_id(info: str) -> str | None:
    for item in info.split(";"):
        if item.startswith("PANGENOME_ALLELE_ID="):
            return item.split("=", 1)[1]
    return None


def genotype(fields: list[str]) -> str:
    if len(fields) < 10:
        return "./."
    keys, values = fields[8].split(":"), fields[9].split(":")
    if "GT" not in keys or keys.index("GT") >= len(values):
        return "./."
    gt = values[keys.index("GT")]
    return gt if gt not in {"", "."} else "./."


def project_all_sites(native: Path, panel: Path, output: Path, sample: str) -> None:
    order: list[str] = []
    records: dict[str, list[str]] = {}
    by_key: dict[tuple[str, str, str, str], str] = {}
    by_allele: dict[str, str] = {}
    meta: list[str] = []
    with open_text(panel, "rt") as handle:
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
                raise RuntimeError(f"duplicate candidate ID: {candidate_id}")
            order.append(candidate_id)
            records[candidate_id] = fields
            by_key[(fields[0], fields[1], fields[3], fields[4])] = candidate_id
            source_id = allele_id(fields[7])
            if source_id:
                by_allele[source_id] = candidate_id
    calls: dict[str, str] = {}
    with open_text(native, "rt") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                raise RuntimeError("BayesTyper output contains a malformed record")
            candidate_id = fields[2] if fields[2] in records else None
            if candidate_id is None:
                source_id = allele_id(fields[7])
                candidate_id = by_allele.get(source_id) if source_id else None
            if candidate_id is None:
                candidate_id = by_key.get((fields[0], fields[1], fields[3], fields[4]))
            if candidate_id:
                if candidate_id in calls:
                    raise RuntimeError(f"BayesTyper output duplicates {candidate_id}")
                calls[candidate_id] = genotype(fields)
    output.parent.mkdir(parents=True, exist_ok=True)
    with open_text(output, "wt") as handle:
        handle.write("##fileformat=VCFv4.2\n")
        for line in meta:
            if not line.startswith("##fileformat") and not line.startswith("##FORMAT=<ID=GT,"):
                handle.write(line)
        handle.write("##source=PGBench-BayesTyper-all-sites-adapter\n")
        handle.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        handle.write(f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}\n")
        for candidate_id in order:
            handle.write("\t".join(records[candidate_id][:8] + ["GT", calls.get(candidate_id, "./.")]) + "\n")
    print(f"BayesTyper candidate projection: matched={len(calls)} no_call={len(order) - len(calls)}")


def main() -> int:
    r1 = Path(required("PGBENCH_INPUT_FASTQ_R1"))
    r2 = Path(required("PGBENCH_INPUT_FASTQ_R2"))
    panel = Path(required("PGBENCH_CANDIDATE_VCF"))
    index = Path(required("PGBENCH_INDEX_DIR"))
    sample = required("PGBENCH_SAMPLE_ID")
    threads = required("PGBENCH_THREADS")
    output_dir = Path(required("PGBENCH_OUTPUT_DIR")).resolve()
    output_vcf = Path(required("PGBENCH_OUTPUT_VCF")).resolve()
    output_vcf.relative_to(output_dir)
    canon, decoy = index / "canon.fa", index / "decoy.fa"
    if not canon.is_file() or not decoy.is_file():
        raise RuntimeError("BayesTyper tool_index must contain canon.fa and decoy.fa")
    work = (output_dir / "bayestyper-work").resolve()
    kmers, temporary = work / "kmers", work / "kmc-tmp"
    kmers.mkdir(parents=True, exist_ok=True)
    temporary.mkdir(parents=True, exist_ok=True)
    read_list = work / "reads.list"
    read_list.write_text(f"{r1}\n{r2}\n", encoding="utf-8")
    prefix = kmers / sample
    run(["kmc", "-k55", "-ci1", "-fq", f"@{read_list}", str(prefix), str(temporary)])
    run(["bayesTyperTools", "makeBloom", "-k", str(prefix), "-p", threads])
    sample_file = work / "samples.tsv"
    sample_file.write_text(f"{sample}\tmale\t{prefix}\n", encoding="utf-8")
    candidates = work / "candidates.vcf.gz"
    with open_text(panel, "rt") as source, gzip.GzipFile(filename="", mode="wb", fileobj=candidates.open("wb"), mtime=0) as raw:
        for line in source:
            raw.write(line.encode("utf-8"))
    run(["bayesTyper", "cluster", "-v", str(candidates), "-s", str(sample_file), "-g", str(canon), "-d", str(decoy), "-p", threads], cwd=work)
    units = sorted(work.glob("bayestyper_unit_*/variant_clusters.bin"), key=lambda path: int(path.parent.name.rsplit("_", 1)[1]))
    if not units:
        raise RuntimeError("BayesTyper cluster produced no variant units")
    vcfs: list[Path] = []
    for unit in units:
        out_prefix = unit.parent / "calls"
        run(["bayesTyper", "genotype", "-v", str(unit), "-c", str(work / "bayestyper_cluster_data"), "-s", str(sample_file), "-g", str(canon), "-d", str(decoy), "-p", threads, "-z", "-o", str(out_prefix)])
        vcfs.append(Path(f"{out_prefix}.vcf.gz"))
    native = work / "native.vcf.gz"
    run(["bcftools", "concat", "-O", "z", "-o", str(native), *map(str, vcfs)])
    project_all_sites(native, panel, output_vcf, sample)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
