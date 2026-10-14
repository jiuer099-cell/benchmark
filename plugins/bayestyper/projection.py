"""Truth-blind, representation-preserving projection of native BayesTyper alleles."""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import gzip
import hashlib
import json
from pathlib import Path
from urllib.parse import quote

CONTRACT = "pgbench-native-projection-v2"


def digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def text(path, mode="rt"):
    return gzip.open(path, mode, encoding="utf-8") if str(path).endswith(".gz") else Path(path).open(mode, encoding="utf-8")


def info(value):
    return dict(item.split("=", 1) for item in value.split(";") if "=" in item)


def trimmed(chrom, pos, ref, alt):
    pos = int(pos)
    ref, alt = ref.upper(), alt.upper()
    while len(ref) > 1 and len(alt) > 1 and ref[-1] == alt[-1]:
        ref, alt = ref[:-1], alt[:-1]
    while len(ref) > 1 and len(alt) > 1 and ref[0] == alt[0]:
        pos += 1
        ref, alt = ref[1:], alt[1:]
    return chrom, pos, ref, alt


class Reference:
    def __init__(self, path):
        self.stream = Path(path).open("rb")
        with Path(str(path) + ".fai").open() as handle:
            self.index = {f[0]: tuple(map(int, f[1:5])) for f in (line.rstrip().split("\t") for line in handle)}

    def close(self):
        self.stream.close()

    def fetch(self, chrom, start, end):
        length, offset, bases, width = self.index[chrom]
        if not 0 <= start <= end <= length:
            raise ValueError("reference interval is out of bounds")
        if start == end:
            return ""
        first = offset + start // bases * width + start % bases
        last = offset + (end - 1) // bases * width + (end - 1) % bases + 1
        self.stream.seek(first)
        result = self.stream.read(last - first).replace(b"\n", b"").replace(b"\r", b"").decode("ascii").upper()
        if len(result) != end - start:
            raise ValueError("reference FASTA/index mismatch")
        return result

    def equivalent(self, a, b):
        # Do not infer equivalence from graph IDs alone, size, or proximity.
        if a[0] != b[0]:
            return False
        left = min(int(a[1]), int(b[1])) - 1
        right = max(int(a[1]) - 1 + len(a[3]), int(b[1]) - 1 + len(b[3]))
        if right - left > 1_000_000:
            return False
        seq = self.fetch(a[0], left, right)
        haplotypes = []
        for f in (a, b):
            offset = int(f[1]) - 1 - left
            if seq[offset:offset + len(f[3])] != f[3].upper():
                return False
            if not set(f[3].upper() + f[4].upper()) <= set("ACGT"):
                return False
            haplotypes.append(seq[:offset] + f[4].upper() + seq[offset + len(f[3]):])
        return haplotypes[0] == haplotypes[1]


def dosage_gt(gt, allele, alt_count, alts=None):
    tokens = gt.replace("|", "/").split("/")
    if len(tokens) != 2 or any(not x.isdigit() or int(x) > alt_count for x in tokens):
        return "./."
    # BayesTyper writes a '*' allele for a nested variant not covered by
    # a haplotype. It cannot supply a complete diploid sequence genotype.
    if alts is not None and any(int(x) and not set(alts[int(x) - 1].upper()) <= set("ACGT") for x in tokens):
        return "./."
    dosage = sum(int(x) == allele for x in tokens)
    return ("0/0", "0/1", "1/1")[dosage]


def project(native, candidates, destination, sample, *, reference=None, audit_path=None, trace_path=None, native_policy="not_declared"):
    native, candidates, destination = map(Path, (native, candidates, destination))
    audit_path = Path(audit_path or destination.with_suffix(".projection.json"))
    trace_path = Path(trace_path or destination.with_suffix(".native-trace.jsonl.gz"))
    panel = {}; headers = []; keys = defaultdict(set); paths = defaultdict(set); claimed = defaultdict(set)
    with text(candidates) as stream:
        for line in stream:
            if line.startswith("##"):
                headers.append(line); continue
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\r\n").split("\t")
            if len(f) < 8 or f[2] in panel or f[2] in {"", "."} or "," in f[4]:
                raise ValueError("projection requires a unique-ID biallelic panel")
            panel[f[2]] = f
            keys[trimmed(f[0], f[1], f[3], f[4])].add(f[2])
            data = info(f[7]); at = data.get("AT", "").split(",")
            if len(at) == 2 and all(at): paths[(f[0], at[0], at[1])].add(f[2])
            if data.get("PANGENOME_ALLELE_ID"): claimed[data["PANGENOME_ALLELE_ID"]].add(f[2])
    if not panel: raise ValueError("empty candidate panel")
    ref = Reference(reference) if reference else None
    selected = {}; ambiguous = set(); counts = Counter(); native_headers = []; saw_header = False
    for p in (destination, audit_path, trace_path): p.parent.mkdir(parents=True, exist_ok=True)
    try:
        with text(native) as stream, gzip.open(trace_path, "wt", encoding="utf-8") as trace:
            for line in stream:
                if line.startswith("##"):
                    native_headers.append(line); continue
                if line.startswith("#CHROM"):
                    if saw_header:
                        raise ValueError("duplicate native VCF sample header")
                    if line.rstrip().split("\t")[9:] != [sample]:
                        raise ValueError("native VCF must contain the requested single sample")
                    saw_header = True
                    continue
                if line.startswith("#") or not line.strip(): continue
                if not saw_header:
                    raise ValueError("native VCF record occurs before the sample header")
                f = line.rstrip("\r\n").split("\t")
                if len(f) != 10: raise ValueError("malformed native VCF record")
                counts["native_records"] += 1
                data = info(f[7]); alts = f[4].split(","); at = data.get("AT", "").split(",")
                values = dict(zip(f[8].split(":"), f[9].split(":")))
                record_matched = False
                for index, alt in enumerate(alts, 1):
                    n = list(f); n[4] = alt
                    hits = set(keys.get(trimmed(f[0], f[1], f[3], alt), ()))
                    method = "normalized_sequence"
                    path_hits = paths.get((f[0], at[0], at[index]), ()) if len(at) == len(alts) + 1 else ()
                    hits.update(path_hits)
                    hits.update(cid for cid in f[2].split(";") if cid in panel)
                    hits.update(claimed.get(data.get("PANGENOME_ALLELE_ID"), ()))
                    # Validate every proposed mapping, including claimed IDs.
                    verified = {cid for cid in hits if (ref.equivalent(n, panel[cid]) if ref else trimmed(*[n[i] for i in (0,1,3,4)]) == trimmed(*[panel[cid][i] for i in (0,1,3,4)]))}
                    if len(verified) > 1:
                        ambiguous.update(verified); counts["ambiguous_native_alleles"] += 1; continue
                    if not verified: continue
                    cid = next(iter(verified)); record_matched = True
                    method = "graph_path_and_reference" if cid in path_hits else method
                    gt = dosage_gt(values.get("GT", "./."), index, len(alts), alts)
                    trace_record = {"candidate_id": cid, "native_record": f, "native_alt_index": index, "projected_gt": gt, "method": method}
                    trace.write(json.dumps(trace_record, separators=(",", ":")) + "\n")
                    if cid in selected:
                        # Never pick a record based on quality, truth, or scan order.
                        old = selected[cid]
                        if old[1] != index or old[0] != f:
                            ambiguous.add(cid)
                        counts["duplicate_native_mappings"] += 1
                    else:
                        selected[cid] = (f, index, gt, method)
                        counts[method] += 1
                if not record_matched: counts["unmatched_native_records"] += 1
        if not saw_header:
            raise ValueError("native VCF lacks a sample header")
        for cid in ambiguous: selected.pop(cid, None)
        output_headers = {}
        for line in headers + native_headers:
            if line.startswith(("##fileformat", "##source=", "##FORMAT=<ID=GT,")): continue
            if line.startswith("##FORMAT=<") and not line.startswith(("##FORMAT=<ID=GQ,", "##FORMAT=<ID=DP,")):
                continue
            key = line.split(",", 1)[0] if line.startswith(("##INFO=<", "##FORMAT=<", "##FILTER=<", "##contig=<")) else line
            output_headers[key] = line
        temporary = destination.with_name(".tmp." + destination.name)
        with text(temporary, "wt") as out:
            out.write("##fileformat=VCFv4.2\n")
            out.writelines(output_headers.values())
            out.write("##source=PGBench-BayesTyper-native-projection-v2\n")
            out.write("##PGBENCH_NATIVE_PROJECTION=reference_haplotype_verified_v2\n")
            for key, typ, description in (("SVTYPE", "String", "Allele-specific structural variant type"), ("SVLEN", "Integer", "ALT minus REF sequence length"), ("END", "Integer", "Last REF position")):
                if f"##INFO=<ID={key}" not in output_headers:
                    out.write(f'##INFO=<ID={key},Number=1,Type={typ},Description="{description}">\n')
            for key, description in {"PGBENCH_NATIVE_STATUS": "Native projection state, not inferred from output completeness", "PGBENCH_PROJECTION_METHOD": "Truth-blind verified mapping method", "PGBENCH_NATIVE_RECORD_ID": "Percent-encoded original native record ID", "PGBENCH_NATIVE_ALT_INDEX": "One-based source ALT index", "PGBENCH_NATIVE_TRACE_SHA256": "Digest of full native-record trace"}.items():
                out.write(f'##INFO=<ID={key},Number=1,Type=String,Description="{description}">\n')
            out.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Unphased ALT-specific diploid dosage">\n')
            out.write(f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}\n")
            trace_sha = digest(trace_path)
            for cid, candidate in panel.items():
                if cid not in selected:
                    f = candidate[:8].copy(); f[5:7] = [".", "."]
                    status = "ambiguous_mapping" if cid in ambiguous else "not_observed"
                    f[7] += f";PGBENCH_NATIVE_STATUS={status}"
                    out.write("\t".join(f + ["GT", "./."]) + "\n"); continue
                native_fields, index, gt, method = selected[cid]
                f = native_fields[:8].copy(); f[2] = cid; f[4] = native_fields[4].split(",")[index - 1]
                # Rebuild allele-specific INFO; untouched original fields live in trace.
                allele_id = info(candidate[7]).get("PANGENOME_ALLELE_ID", cid)
                length = len(f[4]) - len(f[3]); svtype = "INS" if length > 0 else "DEL"
                f[7] = f"PANGENOME_ALLELE_ID={allele_id};SVTYPE={svtype};SVLEN={length};END={int(f[1]) + len(f[3]) - 1};PGBENCH_NATIVE_STATUS=addressable;PGBENCH_PROJECTION_METHOD={method};PGBENCH_NATIVE_RECORD_ID={quote(native_fields[2], safe='')};PGBENCH_NATIVE_ALT_INDEX={index};PGBENCH_NATIVE_TRACE_SHA256={trace_sha}"
                values = dict(zip(native_fields[8].split(":"), native_fields[9].split(":")))
                fmt = ["GT"] + [key for key in ("GQ", "DP") if key in values]
                sample_value = ":".join([gt] + [values[key] for key in fmt[1:]])
                out.write("\t".join(f + [":".join(fmt), sample_value]) + "\n")
        temporary.replace(destination)
        audit = {"contract": CONTRACT, "truth_used": False, "reference_equivalence_checked": ref is not None, "caller_alleles_preserved": True, "quality_and_filter_preserved": True, "candidate_count": len(panel), "matched_candidates": len(selected), "ambiguous_candidates": len(ambiguous), "not_observed_candidates": len(panel) - len(selected) - len(ambiguous), "native_vcf": str(native), "native_sha256": digest(native), "candidate_sha256": digest(candidates), "reference_path": str(reference) if reference else None, "reference_fai_sha256": digest(str(reference) + ".fai") if reference else None, "output_sha256": digest(destination), "trace_path": str(trace_path), "trace_sha256": trace_sha, "counts": dict(counts), "biological_performance_validated": False}
        audit["native_execution_policy"] = native_policy
        # The release contract binds FASTA bytes, not merely its path or FAI.
        audit["reference_sha256"] = digest(reference) if reference else None
        audit["native_vcf"] = str(native.resolve())
        audit["reference_hash_owner"] = "benchmark_core"
        audit_path.write_text(json.dumps(audit, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        return len(selected), len(panel) - len(selected), counts["unmatched_native_records"]
    finally:
        if ref: ref.close()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("native", "candidates", "output", "reference", "audit", "trace"):
        parser.add_argument("--" + name, required=True, type=Path)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--native-policy", default="not_declared")
    a = parser.parse_args()
    print(project(a.native, a.candidates, a.output, a.sample, reference=a.reference, audit_path=a.audit, trace_path=a.trace, native_policy=a.native_policy), flush=True)


if __name__ == "__main__": main()
