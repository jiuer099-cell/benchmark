"""Re-genotype retained Paragraph breakpoint counts using upstream grmpy.

No truth is read. Reconstructed counts are usable only after an old-depth
round trip reproduces the saved native GT, filters, likelihoods and counts.
Any non-equivalent site falls back to native BAM genotyping, never zero-fill.
"""
from __future__ import annotations

import copy
import gzip
import json
import math
import shutil
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import subprocess
import sys


def restore_alignment(graph: dict, native: dict, sample: str) -> dict:
    if graph.get("ID") != native.get("graphinfo", {}).get("ID"):
        raise ValueError("graph/native ID mismatch")
    stats = native["samples"][sample]
    counts = {}
    breakpoints = stats.get("breakpoints")
    if not isinstance(breakpoints, dict) or not breakpoints:
        raise ValueError("no complete saved breakpoint counts")
    for bp in breakpoints.values():
        edges = bp.get("counts", {}).get("edges")
        if not isinstance(edges, dict) or not edges:
            raise ValueError("incomplete saved breakpoint counts")
        for edge, count in edges.items():
            if not isinstance(count, int) or isinstance(count, bool) or count < 0:
                raise ValueError("invalid saved edge count")
            if edge in counts and counts[edge] != count:
                raise ValueError("conflicting saved edge counts")
            counts[edge] = count
    restored = copy.deepcopy(graph)
    restored["read_counts_by_edge"] = counts
    restored["alignment_statistics"] = {k: v for k, v in stats.items() if k not in {"gt", "breakpoints"}}
    restored["fragment_statistics"] = {}
    return restored


def equivalent(a, b) -> bool:
    """Exact discrete comparison, tolerance only for printed floating GLs."""
    if isinstance(a, dict) and isinstance(b, dict):
        return a.keys() == b.keys() and all(equivalent(a[k], b[k]) for k in a)
    if isinstance(a, list) and isinstance(b, list):
        return len(a) == len(b) and all(equivalent(x, y) for x, y in zip(a, b))
    if isinstance(a, int) and isinstance(b, int):
        return a == b
    if isinstance(a, (int, float)) and isinstance(b, (int, float)):
        return math.isclose(a, b, rel_tol=1e-7, abs_tol=1e-8)
    return a == b


def assert_roundtrip(saved: dict, replayed: dict, sample: str) -> None:
    if saved["graphinfo"]["ID"] != replayed["graphinfo"]["ID"]:
        raise ValueError("roundtrip graph ID mismatch")
    for key in ("gt", "breakpoints"):
        if not equivalent(saved["samples"][sample].get(key), replayed["samples"][sample].get(key)):
            raise ValueError("old-depth native roundtrip mismatch: " + key)


def _grmpy(alignment: Path, graph: Path, manifest: Path, reference: Path, output: Path) -> dict:
    subprocess.run(["grmpy", "-m", str(manifest), "-g", str(graph), "-r", str(reference),
                    "-o", str(output), "-t", "1", "--progress=0", "--log-level=error"],
                   check=True, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
    value = json.loads(output.read_text())
    if isinstance(value, list):
        if len(value) != 1:
            raise ValueError("unexpected single-site grmpy output")
        value = value[0]
    return value


def run(bundle: dict, *, panel_ids: set[str], work: Path, bam: Path, bai: Path,
        reference: Path, sample: str, depth: float, read_length: int, workers: int) -> tuple[list[Path], dict]:
    if bundle.get("contract") != "paragraph-native-readcount-replay-v1" or bundle.get("sample") != sample:
        raise RuntimeError("invalid Paragraph count replay bundle")
    chunks = bundle["chunks"]
    seen = set()
    events = []
    for ci, chunk in enumerate(chunks):
        for ei, event in enumerate(chunk["events"]):
            cid = event["vcf_record"].split("\t")[2]
            if cid not in panel_ids or cid in seen:
                raise RuntimeError("Paragraph replay candidate set mismatch/duplicate")
            seen.add(cid)
            events.append((ci, ei, event, chunk["vcf_header"]))
    if seen != panel_ids:
        raise RuntimeError("Paragraph replay does not cover the complete frozen panel")
    max_reads = int(bundle["legacy_max_reads_per_event"])
    legacy_depth = float(bundle["legacy_depth"])
    legacy_length = int(bundle["legacy_read_length"])
    if max_reads <= 0 or legacy_depth <= 0 or legacy_length != read_length:
        raise RuntimeError("invalid legacy count acquisition parameters")
    site_root = work / "readcount-replay"
    site_root.mkdir(parents=True, exist_ok=True)

    def process(item):
        ci, ei, event, header = item
        site = site_root / f"{ci:03d}-{ei:04d}"
        site.mkdir(exist_ok=True)
        graph = event["graph"]
        graph_path = site / "graph.json"
        graph_path.write_text(json.dumps(graph))
        alignment = site / "alignment.json"
        old_manifest = site / "sample-old.tsv"
        new_manifest = site / "sample-new.tsv"
        # BAM is kept as the authentic primary evidence; pre-aligned data
        # prevents read extraction. Its index is explicitly registered.
        columns = "id\tpath\tindex_path\tparagraph\tread length\tdepth\n"
        def manifest(path, d, length):
            path.write_text(columns + f"{sample}\t{bam}\t{bai}\t{alignment}\t{length}\t{d:.6f}\n")
        try:
            restored = restore_alignment(graph, event["native"], sample)
            alignment.write_text(json.dumps(restored))
            manifest(old_manifest, legacy_depth, legacy_length)
            old = _grmpy(alignment, graph_path, old_manifest, reference, site / "roundtrip.json")
            assert_roundtrip(event["native"], old, sample)
            manifest(new_manifest, depth, read_length)
            new = _grmpy(alignment, graph_path, new_manifest, reference, site / "corrected.json")
            return ci, ei, new, "reused_counts"
        except (ValueError, KeyError) as error:
            print(f"Paragraph count replay fallback {graph.get('ID')}: {error}", flush=True)
        # A failed grmpy command is not masked as a successful no-call.
        # Only incomplete/non-equivalent source statistics permit fallback.
        vcf = site / "input.vcf"
        vcf.write_text(header + event["vcf_record"].rstrip("\n") + "\n")
        fresh_manifest = site / "sample-bam.tsv"
        fresh_manifest.write_text("id\tpath\tindex_path\tread length\tdepth\n" +
                                  f"{sample}\t{bam}\t{bai}\t{read_length}\t{depth:.6f}\n")
        out = site / "fresh"
        subprocess.run(["multigrmpy.py", "-i", str(vcf), "-m", str(fresh_manifest),
                        "-r", str(reference), "-o", str(out), "-t", "1", "-M", str(max_reads)],
                       check=True, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
        with gzip.open(out / "genotypes.json.gz", "rt") as h:
            new = json.load(h)
        if len(new) != 1:
            raise RuntimeError("fallback produced unexpected native genotype set")
        return ci, ei, new[0], "fresh_alignment"

    completed = {}
    counts = {"reused_counts": 0, "fresh_alignment": 0}
    with ThreadPoolExecutor(max_workers=max(1, workers)) as pool:
        for ci, ei, result, policy in pool.map(process, events):
            completed[ci, ei] = result
            counts[policy] += 1
            total = sum(counts.values())
            if total % 100 == 0 or total == len(events):
                print(f"Paragraph native re-genotyping {total}/{len(events)}; {counts}", flush=True)
    # Use upstream's own VCF writer, not a reimplementation of GT conversion.
    native_runner = shutil.which("multigrmpy.py")
    if native_runner is None:
        raise RuntimeError("native Paragraph environment is unavailable")
    native_bin = Path(native_runner).resolve().parent
    native_python = native_bin / "python"
    if not native_python.is_file():
        raise RuntimeError("native Paragraph Python interpreter is unavailable")
    result_paths = []
    for ci, chunk in enumerate(chunks):
        out = work / f"replayed-chunk-{ci:03d}"
        out.mkdir(exist_ok=True)
        native_json = out / "genotypes.json.gz"
        with gzip.open(native_json, "wt") as h:
            json.dump([completed[ci, ei] for ei in range(len(chunk["events"]))], h)
        source_vcf = out / "variants.vcf"
        source_vcf.write_text(chunk["vcf_header"] + "".join(e["vcf_record"].rstrip("\n") + "\n" for e in chunk["events"]))
        result_vcf = out / "genotypes.vcf.gz"
        # Core may use Python 3.12 while native pysam is built for 3.11.
        # Execute the released writer in its own declared native environment.
        subprocess.run([str(native_python), "-B", "-c",
            "import sys; sys.path.insert(0,sys.argv[1]); import findgrm; "
            "from grm.vcfgraph import vcfupdate; "
            "vcfupdate.update_vcf_from_grmpy(sys.argv[2],vcfupdate.read_grmpy(sys.argv[3]),sys.argv[4],[sys.argv[5]])",
            str(native_bin), str(source_vcf), str(native_json), str(result_vcf), sample], check=True)
        result_paths.append(result_vcf)
    return result_paths, {**counts, "old_depth_roundtrip_required": True,
                          "legacy_depth": legacy_depth, "new_depth": depth,
                          "read_length": read_length, "count_acquisition_max_reads": max_reads}
