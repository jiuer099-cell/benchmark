#!/usr/bin/env python3
"""PGBench external adapter for GraphTyper2 SV genotyping."""

from __future__ import annotations

import gzip
import os
import shutil
import subprocess
import threading
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import TextIO


def required(name: str) -> str:
    value = os.environ.get(name)
    if not value:
        raise RuntimeError(f"required environment variable is absent: {name}")
    return value


def open_text(path: Path, mode: str) -> TextIO:
    if path.suffix.casefold() == ".gz":
        return gzip.open(path, mode, encoding="utf-8")
    return path.open(mode, encoding="utf-8")


def run(command: list[str], *, stdout: Path | None = None) -> None:
    print("+ " + " ".join(command), flush=True)
    if stdout is None:
        subprocess.run(command, check=True)
        return
    with stdout.open("wb") as handle:
        subprocess.run(command, check=True, stdout=handle)


def candidate_keys(path: Path) -> tuple[list[str], dict[tuple[str, str, str, str], str]]:
    order: list[str] = []
    keys: dict[tuple[str, str, str, str], str] = {}
    with open_text(path, "rt") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8 or not fields[2].startswith("CAND_"):
                raise RuntimeError(f"invalid blinded candidate at line {line_number}")
            key = (fields[0], fields[1], fields[3], fields[4])
            if key in keys:
                raise RuntimeError(f"duplicate blinded candidate coordinates: {key}")
            keys[key] = fields[2]
            order.append(fields[2])
    if not order:
        raise RuntimeError("blinded candidate panel contains no records")
    return order, keys


def stage_candidate(source: Path, destination: Path) -> None:
    run(["bcftools", "sort", "-Oz", "-o", str(destination), str(source)])
    run(["tabix", "-f", "-p", "vcf", str(destination)])


def write_regions(candidate: Path, destination: Path, chunk_size: int = 10_000_000) -> None:
    regions: set[tuple[str, int, int]] = set()
    with open_text(candidate, "rt") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.split("\t", 3)
            pos = int(fields[1])
            start = ((pos - 1) // chunk_size) * chunk_size + 1
            regions.add((fields[0], start, start + chunk_size - 1))
    with destination.open("w", encoding="utf-8") as handle:
        for chrom, start, end in sorted(regions, key=lambda item: (item[0], item[1])):
            handle.write(f"{chrom}:{start}-{end}\n")


def write_region_shards(
    regions: Path,
    destination: Path,
    shard_count: int,
) -> list[Path]:
    """Split ordered regions deterministically across balanced worker shards."""

    ordered = [
        line.strip()
        for line in regions.read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]
    if not ordered:
        raise RuntimeError("GraphTyper2 region list is empty")
    workers = min(max(1, shard_count), len(ordered))
    destination.mkdir(parents=True, exist_ok=False)
    buckets: list[list[str]] = [[] for _ in range(workers)]
    for index, region in enumerate(ordered):
        buckets[index % workers].append(region)

    shards: list[Path] = []
    for index, bucket in enumerate(buckets):
        shard = destination / f"regions.{index:03d}.txt"
        shard.write_text("\n".join(bucket) + "\n", encoding="utf-8")
        shards.append(shard)
    return shards


def run_graphtyper_shards(
    reference: Path,
    candidate: Path,
    bam: Path,
    regions: Path,
    output: Path,
    threads: int,
) -> None:
    """Run one single-BAM GraphTyper process per deterministic region shard.

    GraphTyper's internal SAM-reader parallelism is bounded by the number of
    input BAMs.  A one-sample benchmark therefore uses process-level regional
    parallelism while keeping the aggregate worker count within PGBENCH_THREADS.
    """

    workers = max(1, threads)
    shard_root = output.parent / "region-shards"
    if shard_root.exists():
        shutil.rmtree(shard_root)
    shards = write_region_shards(regions, shard_root, workers)
    output.mkdir(parents=True, exist_ok=False)

    active: dict[int, subprocess.Popen[bytes]] = {}
    active_lock = threading.Lock()

    def execute(index_and_shard: tuple[int, Path]) -> None:
        index, shard = index_and_shard
        shard_output = output / f"shard_{index:03d}"
        command = [
            "graphtyper",
            "genotype_sv",
            str(reference),
            str(candidate),
            "--sam=" + str(bam),
            "--region_file=" + str(shard),
            "--threads=1",
            "--output=" + str(shard_output),
            "--force_no_copy_reference",
        ]
        print("+ " + " ".join(command), flush=True)
        with active_lock:
            process = subprocess.Popen(command)
            active[index] = process
        returncode = process.wait()
        with active_lock:
            active.pop(index, None)
            if returncode != 0:
                for sibling in active.values():
                    sibling.terminate()
        if returncode != 0:
            raise subprocess.CalledProcessError(returncode, command)

    failures: list[BaseException] = []
    with ThreadPoolExecutor(max_workers=len(shards)) as executor:
        futures = [
            executor.submit(execute, item)
            for item in enumerate(shards)
        ]
        for future in as_completed(futures):
            try:
                future.result()
            except BaseException as exc:  # propagate after all children settle
                failures.append(exc)
                for pending in futures:
                    pending.cancel()
    if failures:
        raise RuntimeError(
            f"{len(failures)} GraphTyper2 region shard(s) failed"
        ) from failures[0]


def align_reads(
    reference: Path,
    read1: Path,
    read2: Path,
    bam: Path,
    sample: str,
    threads: str,
) -> None:
    index_prefix = bam.parent / "grch38.bwa-mem2"
    index_sentinel = Path(f"{index_prefix}.bwt.2bit.64")
    if not index_sentinel.is_file():
        run(["bwa-mem2", "index", "-p", str(index_prefix), str(reference)])
    read_group = f"@RG\\tID:{sample}.illumina\\tSM:{sample}\\tPL:ILLUMINA"
    mapper = subprocess.Popen(
        [
            "bwa-mem2", "mem", "-t", threads, "-R", read_group,
            str(index_prefix), str(read1), str(read2),
        ],
        stdout=subprocess.PIPE,
    )
    assert mapper.stdout is not None
    sorter = subprocess.run(
        ["samtools", "sort", "-@", threads, "-o", str(bam), "-"],
        stdin=mapper.stdout,
        check=False,
    )
    mapper.stdout.close()
    mapper_returncode = mapper.wait()
    if mapper_returncode != 0 or sorter.returncode != 0:
        raise subprocess.CalledProcessError(
            mapper_returncode or sorter.returncode,
            ["bwa-mem2", "mem", "|", "samtools", "sort"],
        )
    run(["samtools", "index", "-@", threads, str(bam)])
    run(["samtools", "quickcheck", "-v", str(bam)])


def generated_vcfs(root: Path) -> list[Path]:
    result = sorted(
        path for path in root.rglob("*.vcf.gz")
        if path.is_file() and path.stat().st_size > 0
    )
    if not result:
        result = sorted(
            path for path in root.rglob("*.vcf")
            if path.is_file() and path.stat().st_size > 0
        )
    if not result:
        raise RuntimeError("GraphTyper2 produced no VCF files")
    return result


def project_calls(candidate: Path, generated: list[Path], output: Path, sample: str) -> None:
    order, key_to_id = candidate_keys(candidate)
    headers: list[str] = []
    by_candidate: dict[str, list[list[str]]] = defaultdict(list)
    for source in generated:
        with open_text(source, "rt") as handle:
            for line in handle:
                if line.startswith("##"):
                    if line not in headers:
                        headers.append(line)
                    continue
                if not line.strip() or line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 10:
                    continue
                candidate_id = key_to_id.get((fields[0], fields[1], fields[3], fields[4]))
                if candidate_id is not None:
                    by_candidate[candidate_id].append(fields)

    candidate_records: dict[str, list[str]] = {}
    with open_text(candidate, "rt") as handle:
        for line in handle:
            if line.strip() and not line.startswith("#"):
                fields = line.rstrip("\n").split("\t")
                candidate_records[fields[2]] = fields

    output.parent.mkdir(parents=True, exist_ok=True)
    matched = conflicts = no_calls = 0
    with output.open("w", encoding="utf-8") as handle:
        handle.write("##fileformat=VCFv4.2\n")
        for header in headers:
            if not header.startswith("##fileformat") and not header.startswith(
                "##FORMAT=<ID=GT,"
            ):
                handle.write(header if header.endswith("\n") else header + "\n")
        handle.write("##source=PGBench-GraphTyper2-2.7.7-adapter\n")
        handle.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        handle.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + sample + "\n")
        for candidate_id in order:
            records = by_candidate.get(candidate_id, [])
            unique = {tuple(record) for record in records}
            if len(unique) == 1:
                fields = list(next(iter(unique)))
                fields[2] = candidate_id
                handle.write("\t".join(fields[:10]) + "\n")
                matched += 1
                continue
            base = candidate_records[candidate_id]
            handle.write("\t".join(base[:8] + ["GT", "./."]) + "\n")
            no_calls += 1
            if len(unique) > 1:
                conflicts += 1
    print(
        f"GraphTyper2 candidate projection: matched={matched} "
        f"no_call={no_calls} conflicts={conflicts}",
        flush=True,
    )


def main() -> int:
    read1 = Path(required("PGBENCH_INPUT_FASTQ_R1"))
    read2 = Path(required("PGBENCH_INPUT_FASTQ_R2"))
    reference = Path(required("PGBENCH_REFERENCE_FASTA"))
    reference_index = Path(required("PGBENCH_REFERENCE_INDEX"))
    expected_reference_index = Path(f"{reference}.fai")
    if reference_index.resolve() != expected_reference_index.resolve():
        raise RuntimeError(
            "GraphTyper2 requires reference_index to be the .fai adjacent to "
            f"the reference FASTA: expected {expected_reference_index}, got "
            f"{reference_index}"
        )
    candidate = Path(required("PGBENCH_CANDIDATE_VCF"))
    output = Path(required("PGBENCH_OUTPUT_VCF"))
    output_dir = Path(required("PGBENCH_OUTPUT_DIR"))
    sample = required("PGBENCH_SAMPLE_ID")
    threads = required("PGBENCH_THREADS")

    work = output_dir / "work"
    work.mkdir(parents=True, exist_ok=True)
    bam = work / f"{sample}.illumina.sorted.bam"
    if not bam.is_file() or not Path(f"{bam}.bai").is_file():
        align_reads(reference, read1, read2, bam, sample, threads)
    else:
        run(["samtools", "quickcheck", "-v", str(bam)])

    staged_candidate = work / "candidate.vcf.gz"
    stage_candidate(candidate, staged_candidate)
    regions = work / "regions.txt"
    write_regions(candidate, regions)
    graphtyper_output = work / "graphtyper"
    if graphtyper_output.exists():
        shutil.rmtree(graphtyper_output)
    run_graphtyper_shards(
        reference,
        staged_candidate,
        bam,
        regions,
        graphtyper_output,
        int(threads),
    )
    project_calls(candidate, generated_vcfs(graphtyper_output), output, sample)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
