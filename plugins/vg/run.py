#!/usr/bin/env python3
"""PGBench paired-short-read adapter for vg Giraffe + pack + call."""

from __future__ import annotations

import gzip
import hashlib
import json
import os
import shutil
import subprocess
import time
from contextlib import ExitStack
from pathlib import Path
from threading import Event, Thread
from typing import TextIO


SHARD_COUNT = 8
FRAGMENT_MEAN = "596.974"
FRAGMENT_STDEV = "153.29"
CHECKPOINT_VERSION = 1
GIRAFFE_RECIPE = ["vg", "giraffe", "-Z", "{gbz}", "-d", "{distance}", "-m", "{minimizer}", "-z", "{zipcodes}", "-f", "{r1}", "-f", "{r2}", "-t", "{threads}", "-p", "--fragment-mean", "{mean}", "--fragment-stdev", "{stdev}"]
# Immutable recipe of the completed 20261009 checkpoints; never broaden this.
LEGACY_GIRAFFE_RECIPE = ["vg", "giraffe", "-Z", "{gbz}", "-d", "{distance}", "-m", "{minimizer}", "-z", "{zipcodes}", "-f", "{r1}", "-f", "{r2}", "-t", "{threads}", "-p", "--fragment-mean", "{mean}", "--fragment-stdev", "{stdev}"]


def required(name: str) -> str:
    value = os.environ.get(name)
    if not value:
        raise RuntimeError(f"required environment variable is absent: {name}")
    return value


def run(command: list[str], *, stdout: Path | None = None) -> None:
    if stdout is None:
        subprocess.run(command, check=True)
        return
    temporary = stdout.with_name(stdout.name + ".tmp")
    wait_for_memory()
    with temporary.open("wb") as handle:
        process = subprocess.Popen(command, stdout=handle)
        while process.poll() is None:
            try:
                process.wait(timeout=10)
            except subprocess.TimeoutExpired:
                if available_memory_gib() < 96:
                    process.terminate()
                    try:
                        process.wait(timeout=30)
                    except subprocess.TimeoutExpired:
                        process.kill()
                        process.wait()
                    raise RuntimeError("available memory fell below 96 GiB during " + " ".join(command[:2]))
        if process.returncode:
            raise RuntimeError(f"{command[:2]} failed: {process.returncode}")
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temporary, stdout)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def identity(value: object) -> str:
    encoded = json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(encoded).hexdigest()


def write_json_atomic(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", encoding="utf-8") as handle:
        json.dump(value, handle, indent=2, sort_keys=True)
        handle.write("\n")
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temporary, path)


def read_fastq_record(handle: object) -> bytes | None:
    header = handle.readline()
    if not header:
        return None
    sequence = handle.readline()
    separator = handle.readline()
    quality = handle.readline()
    if not sequence or not separator.startswith(b"+") or not quality or not header.startswith(b"@"):
        raise RuntimeError("FASTQ input ends with a malformed record")
    if len(sequence.rstrip(b"\r\n")) != len(quality.rstrip(b"\r\n")):
        raise RuntimeError("FASTQ sequence and quality lengths differ")
    return header + sequence + separator + quality


def pair_name(record: bytes) -> bytes:
    header = record.split(b"\n", 1)[0].split(None, 1)[0]
    return header.removesuffix(b"/1").removesuffix(b"/2")


def shard_paths(root: Path, index: int) -> tuple[Path, Path]:
    return root / f"part-{index:03d}.R1.fastq", root / f"part-{index:03d}.R2.fastq"


def valid_shards(root: Path, split_key: str) -> list[tuple[Path, Path]] | None:
    manifest = root / "manifest.json"
    if not manifest.is_file() or manifest.is_symlink():
        return None
    try:
        data = json.loads(manifest.read_text(encoding="utf-8"))
        if data["split_key"] != split_key or len(data["parts"]) != SHARD_COUNT:
            return None
        pairs = []
        for index, item in enumerate(data["parts"]):
            r1, r2 = shard_paths(root, index)
            if item["index"] != index or item["pairs"] < 1:
                return None
            for path, label in ((r1, "r1"), (r2, "r2")):
                if not path.is_file() or path.is_symlink():
                    return None
                if path.stat().st_size != item[f"{label}_size"]:
                    return None
                if sha256_file(path) != item[f"{label}_sha256"]:
                    return None
            pairs.append((r1, r2))
        return pairs
    except (OSError, KeyError, TypeError, ValueError, json.JSONDecodeError):
        return None


def split_pairs(r1: Path, r2: Path, root: Path, split_key: str) -> list[tuple[Path, Path]]:
    root.mkdir(parents=True, exist_ok=True)
    counts = [0] * SHARD_COUNT
    pairs = [shard_paths(root, index) for index in range(SHARD_COUNT)]
    print(f"vg split: writing {SHARD_COUNT} paired FASTQ shards", flush=True)
    started = time.monotonic()
    with ExitStack() as stack:
        source1 = stack.enter_context(gzip.open(r1, "rb"))
        source2 = stack.enter_context(gzip.open(r2, "rb"))
        outputs = [
            (
                stack.enter_context(a.open("wb", buffering=1024 * 1024)),
                stack.enter_context(b.open("wb", buffering=1024 * 1024)),
            )
            for a, b in pairs
        ]
        pair_number = 0
        while True:
            first = read_fastq_record(source1)
            second = read_fastq_record(source2)
            if first is None and second is None:
                break
            if first is None or second is None or pair_name(first) != pair_name(second):
                raise RuntimeError(f"paired FASTQs disagree at pair {pair_number + 1}")
            index = pair_number % SHARD_COUNT
            outputs[index][0].write(first)
            outputs[index][1].write(second)
            counts[index] += 1
            pair_number += 1
            if pair_number % 10_000_000 == 0:
                print(f"vg split: {pair_number} pairs in {time.monotonic() - started:.0f}s", flush=True)
        for first, second in outputs:
            for handle in (first, second):
                handle.flush()
                os.fsync(handle.fileno())
    if any(count == 0 for count in counts):
        raise RuntimeError("a FASTQ shard is empty")
    parts = []
    for index, (first, second) in enumerate(pairs):
        parts.append({
            "index": index,
            "pairs": counts[index],
            "r1_size": first.stat().st_size,
            "r1_sha256": sha256_file(first),
            "r2_size": second.stat().st_size,
            "r2_sha256": sha256_file(second),
        })
    write_json_atomic(root / "manifest.json", {"split_key": split_key, "parts": parts})
    print(f"vg split: {pair_number} pairs complete in {time.monotonic() - started:.0f}s", flush=True)
    return pairs


def prior_work_roots(output_dir: Path) -> list[Path]:
    history = output_dir.parent.parent
    roots = [output_dir / "work"]
    for attempt in sorted(history.iterdir(), key=lambda path: path.stat().st_mtime, reverse=True):
        if attempt == output_dir.parent or attempt.is_symlink() or not attempt.is_dir():
            continue
        work = attempt / "output" / "work"
        if work.is_dir() and not work.is_symlink():
            roots.append(work)
    return roots


def cached_gam(roots: list[Path], index: int, mapping_key: str) -> Path | None:
    for root in roots:
        directory = root / "checkpoints"
        path = directory / f"part-{index:03d}.gam"
        manifest = directory / f"part-{index:03d}.json"
        if not manifest.is_file() or manifest.is_symlink() or not path.is_file() or path.is_symlink():
            continue
        try:
            data = json.loads(manifest.read_text(encoding="utf-8"))
            if data["mapping_key"] != mapping_key or data["index"] != index:
                continue
            if path.stat().st_size != data["size"] or sha256_file(path) != data["sha256"]:
                continue
            print(f"vg checkpoint: reused completed GAM shard {index + 1}/{SHARD_COUNT}", flush=True)
            return path
        except (OSError, KeyError, TypeError, ValueError, json.JSONDecodeError):
            continue
    return None


def available_memory_gib() -> float:
    with Path("/proc/meminfo").open(encoding="ascii") as handle:
        for line in handle:
            if line.startswith("MemAvailable:"):
                return int(line.split()[1]) / (1024 * 1024)
    raise RuntimeError("cannot read available system memory")


def wait_for_memory() -> None:
    while available_memory_gib() < 240:
        print("vg map: waiting for 240 GiB available memory before another shard", flush=True)
        time.sleep(60)


def map_shard(
    index: int,
    reads: tuple[Path, Path],
    *,
    gbz: Path,
    distance: Path,
    minimizer: Path,
    zipcodes: Path,
    work: Path,
    mapping_threads: int,
    mapping_key: str,
) -> Path:
    wait_for_memory()
    directory = work / "checkpoints"
    directory.mkdir(parents=True, exist_ok=True)
    gam = directory / f"part-{index:03d}.gam"
    temporary = directory / f"part-{index:03d}.gam.tmp"
    print(f"vg map: shard {index + 1}/{SHARD_COUNT} started", flush=True)
    started = time.monotonic()
    with temporary.open("wb") as output:
        giraffe = subprocess.Popen(
            [token.format(gbz=gbz, distance=distance, minimizer=minimizer,
                          zipcodes=zipcodes, r1=reads[0], r2=reads[1],
                          threads=mapping_threads, mean=FRAGMENT_MEAN,
                          stdev=FRAGMENT_STDEV) for token in GIRAFFE_RECIPE],
            stdout=output,
        )
        while giraffe.poll() is None:
            try:
                giraffe.wait(timeout=10)
            except subprocess.TimeoutExpired:
                if available_memory_gib() < 96:
                    print("vg map: aborting current shard to avoid global OOM", flush=True)
                    giraffe.terminate()
                    try:
                        giraffe.wait(timeout=30)
                    except subprocess.TimeoutExpired:
                        giraffe.kill()
                        giraffe.wait()
                    raise RuntimeError("available memory fell below 96 GiB")
        if giraffe.returncode != 0:
            raise RuntimeError(f"vg shard {index + 1} failed: giraffe={giraffe.returncode}")
        output.flush()
        os.fsync(output.fileno())
    if not temporary.is_file() or temporary.stat().st_size == 0:
        raise RuntimeError(f"vg shard {index + 1} produced no GAM")
    os.replace(temporary, gam)
    write_json_atomic(directory / f"part-{index:03d}.json", {
        "mapping_key": mapping_key,
        "index": index,
        "size": gam.stat().st_size,
        "sha256": sha256_file(gam),
    })
    print(f"vg map: shard {index + 1}/{SHARD_COUNT} complete in {time.monotonic() - started:.0f}s", flush=True)
    return gam


def cached_final_pack(roots: list[Path], pack_key: str) -> Path | None:
    for root in roots:
        pack = root / "final-pack" / "HG002.pack"
        manifest = root / "final-pack" / "HG002.pack.json"
        if not pack.is_file() or pack.is_symlink() or not manifest.is_file() or manifest.is_symlink():
            continue
        try:
            data = json.loads(manifest.read_text(encoding="utf-8"))
            if data["pack_key"] == pack_key and data["size"] == pack.stat().st_size:
                if data["sha256"] == sha256_file(pack):
                    print("vg checkpoint: reused completed whole-sample pack", flush=True)
                    return pack
        except (OSError, KeyError, TypeError, ValueError, json.JSONDecodeError):
            continue
    return None


def gam_fingerprint(path: Path) -> dict:
    # cached_gam has already verified this digest; map_shard writes it after exit 0.
    # Avoid rereading 256 GB multiple times solely to build the pack cache key.
    manifest = path.with_suffix(".json")
    if manifest.is_file():
        data = json.loads(manifest.read_text())
        if data.get("size") == path.stat().st_size and len(data.get("sha256", "")) == 64:
            return {"sha256": data["sha256"], "size": data["size"]}
    return {"sha256": sha256_file(path), "size": path.stat().st_size}


def pack_gams(gams: list[Path], gbz: Path, work: Path, pack_key: str) -> Path:
    wait_for_memory()
    directory = work / "final-pack"
    directory.mkdir(parents=True, exist_ok=True)
    temporary = directory / "HG002.pack.tmp"
    pack = directory / "HG002.pack"
    print("vg pack: processing all completed GAM shards once", flush=True)
    process = subprocess.Popen(
        ["vg", "pack", "-x", str(gbz), "-g", "-", "-o", str(temporary), "-t", "4", "-Q", "5"],
        stdin=subprocess.PIPE,
    )
    assert process.stdin is not None
    stop = Event()
    memory_abort = Event()

    def guard_memory() -> None:
        while not stop.wait(10):
            if process.poll() is not None:
                return
            if available_memory_gib() < 96:
                memory_abort.set()
                process.terminate()
                return

    watcher = Thread(target=guard_memory, daemon=True)
    watcher.start()
    try:
        for gam in gams:
            print(f"vg pack: reading {gam.name}", flush=True)
            with gam.open("rb") as source:
                shutil.copyfileobj(source, process.stdin, length=8 * 1024 * 1024)
        process.stdin.close()
        status = process.wait()
    except BaseException:
        if process.poll() is None:
            process.terminate()
        process.wait()
        raise
    finally:
        stop.set()
        watcher.join()
    if memory_abort.is_set():
        raise RuntimeError("vg pack stopped because available memory fell below 96 GiB")
    if status != 0 or not temporary.is_file() or temporary.stat().st_size == 0:
        raise RuntimeError(f"vg pack failed with status {status}")
    os.replace(temporary, pack)
    write_json_atomic(directory / "HG002.pack.json", {
        "pack_key": pack_key,
        "source_gams": [gam_fingerprint(path) for path in gams],
        "min_quality": 5,
        "size": pack.stat().st_size,
        "sha256": sha256_file(pack),
    })
    print("vg pack: whole-sample pack complete", flush=True)
    return pack


def open_text(path: Path, mode: str) -> TextIO:
    return gzip.open(path, mode, encoding="utf-8") if path.name.endswith(".gz") else path.open(mode, encoding="utf-8")


def graph_asset(graph_dir: Path, canonical_name: str, pattern: str) -> Path:
    """Resolve a frozen graph asset without making an index or guessing."""

    canonical = graph_dir / canonical_name
    if canonical.is_file():
        return canonical
    matches = sorted(path for path in graph_dir.glob(pattern) if path.is_file())
    if len(matches) != 1:
        raise RuntimeError(
            f"required vg graph asset {canonical_name} is absent or ambiguous "
            f"under {graph_dir} (pattern {pattern!r}, matches={len(matches)})"
        )
    return matches[0]


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


def project_all_sites(native: Path, candidates: Path, destination: Path, sample: str) -> tuple[int, int, int]:
    # Import by file path as tests load this adapter outside its package.
    import importlib.util
    spec = importlib.util.spec_from_file_location("vg_native_projection", Path(__file__).with_name("projection.py"))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.project(
        native, candidates, destination, sample,
        reference=Path(required("PGBENCH_REFERENCE_FASTA")),
        native_policy="official_pack_q5_all_reference_snarls_v1",
    )

def main() -> int:
    reads_r1 = Path(required("PGBENCH_INPUT_FASTQ_R1"))
    reads_r2 = Path(required("PGBENCH_INPUT_FASTQ_R2"))
    graph_dir = Path(required("PGBENCH_GRAPH_DIR"))
    output_dir = Path(required("PGBENCH_OUTPUT_DIR"))
    output_vcf = Path(required("PGBENCH_OUTPUT_VCF"))
    sample = required("PGBENCH_SAMPLE_ID")
    candidates = Path(required("PGBENCH_CANDIDATE_VCF"))
    threads = int(required("PGBENCH_THREADS"))
    if threads < 8:
        raise RuntimeError("vg sharded mapping requires at least eight assigned threads")
    reference_path = required("PGBENCH_GRAPH_REFERENCE_PATH")

    gbz = graph_asset(graph_dir, "graph.gbz", "*.gbz")
    minimizer = graph_asset(graph_dir, "graph.shortread.withzip.min", "*.shortread.withzip.min")
    zipcodes = graph_asset(graph_dir, "graph.shortread.zipcodes", "*.shortread.zipcodes")
    distance = graph_asset(graph_dir, "graph.dist", "*.dist")

    work = output_dir / "work"
    work.mkdir(parents=True, exist_ok=True)
    output_vcf.parent.mkdir(parents=True, exist_ok=True)
    frozen_snarls = graph_dir / "graph.snarls"
    if not frozen_snarls.is_file():
        snarl_candidates = sorted(path for path in graph_dir.glob("*.snarls") if path.is_file())
        frozen_snarls = snarl_candidates[0] if len(snarl_candidates) == 1 else None
    snarls = frozen_snarls or (work / "graph.snarls")
    native_vcf = work / "native.calls.vcf"

    resolved = Path(required("PGBENCH_RESOLVED_INPUTS"))
    with resolved.open(encoding="utf-8") as handle:
        input_manifest = json.load(handle)
    input_hashes = {item["name"]: item["sha256"] for item in input_manifest["inputs"]}
    split_key = identity({
        "format": CHECKPOINT_VERSION,
        "r1": input_hashes["short_fastq_r1"],
        "r2": input_hashes["short_fastq_r2"],
        "shards": SHARD_COUNT,
        "strategy": "paired-round-robin",
    })
    vg_version = subprocess.check_output(["vg", "version"], text=True).splitlines()[0]
    mapping_identity = {
        "format": CHECKPOINT_VERSION,
        "split_key": split_key,
        "graph": input_hashes["graph_assets"],
        "vg_version": vg_version,
        "mapping_threads": threads,
        "fragment_mean": FRAGMENT_MEAN,
        "fragment_stdev": FRAGMENT_STDEV,
        "giraffe_preset": "default",
    }
    legacy_mapping_key = identity(mapping_identity)
    mapping_key = identity({**mapping_identity, "argv_recipe": GIRAFFE_RECIPE})
    pack_key = identity({
        "format": 2,
        "pack_command_recipe": ["vg", "pack", "-x", "{gbz}", "-g", "-", "-o", "{pack}", "-t", "4", "-Q", "5"],
        "mapping_key": mapping_key,
        "graph": input_hashes["graph_assets"],
        "pack_threads": 4,
    })
    # Recipe changes cannot silently reuse the pre-Q5 pack. Mapping is unchanged.
    roots = prior_work_roots(output_dir)
    gams = [cached_gam(roots, index, mapping_key) for index in range(SHARD_COUNT)]
    if GIRAFFE_RECIPE == LEGACY_GIRAFFE_RECIPE:
        gams = [path or cached_gam(roots, index, legacy_mapping_key)
                for index, path in enumerate(gams)]
    pack = None
    if all(gams):
        source_gams = [gam_fingerprint(path) for path in gams]
        pack_key = identity({"recipe_identity": pack_key, "source_gams": source_gams})
        pack = cached_final_pack(roots, pack_key)
    if pack is None:
        if any(path is None for path in gams):
            shards = None
            for root in roots:
                shards = valid_shards(root / "shards", split_key)
                if shards is not None:
                    print(f"vg split: reused verified shards from {root}", flush=True)
                    break
            if shards is None:
                shards = split_pairs(reads_r1, reads_r2, work / "shards", split_key)
            for index, existing in enumerate(gams):
                if existing is None:
                    gams[index] = map_shard(
                        index, shards[index], gbz=gbz, distance=distance,
                        minimizer=minimizer, zipcodes=zipcodes, work=work,
                        mapping_threads=threads, mapping_key=mapping_key,
                    )
        if any(path is None for path in gams):
            raise RuntimeError("not every vg shard has a verified GAM")
        if not all(gams):
            raise RuntimeError("pack requires all source GAMs")
        # Compute the same source-bound identity whether shards were cached or new.
        recipe_key = identity({"format": 2, "pack_command_recipe": ["vg", "pack", "-x", "{gbz}", "-g", "-", "-o", "{pack}", "-t", "4", "-Q", "5"], "mapping_key": mapping_key, "graph": input_hashes["graph_assets"], "pack_threads": 4})
        pack_key = identity({"recipe_identity": recipe_key, "source_gams": [gam_fingerprint(path) for path in gams]})
        pack = pack_gams(gams, gbz, work, pack_key)
    if frozen_snarls is None:
        run(["vg", "snarls", str(gbz), "-t", str(threads)], stdout=snarls)
    run(
        [
            "vg",
            "call",
            str(gbz),
            "-r",
            str(snarls),
            "-k",
            str(pack),
            "-s",
            sample,
            "-z",
            "--ref-sample",
            reference_path,
            "-t", str(threads),
            "-a",
        ],
        stdout=native_vcf,
    )
    matched, no_call, outside = project_all_sites(native_vcf, candidates, output_vcf, sample)
    print(f"vg candidate projection: matched={matched} no_call={no_call} outside={outside}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
