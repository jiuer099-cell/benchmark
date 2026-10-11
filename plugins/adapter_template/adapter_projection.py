"""Truth-blind, bounded-memory projection for the external adapter kit.

Only identity, exact normalized sequence, or reference-verified haplotype
equivalence is accepted. This is not a truth-aware matching/evaluation engine.
Unsupported native representations fail the entire run, never become no-calls.
"""
from __future__ import annotations

import gzip
import hashlib
import json
import re
from pathlib import Path
from typing import TextIO
from urllib.parse import quote


def open_text(path: Path, mode: str) -> TextIO:
    if "w" in mode and path.suffix == ".gz":
        raise ValueError("formal adapter output must be plain .vcf; Core owns BGZF validation")
    return gzip.open(path, mode, encoding="utf-8") if path.suffix == ".gz" else path.open(mode, encoding="utf-8")


def digest(path: Path) -> str:
    result = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            result.update(block)
    return result.hexdigest()


def info(value: str) -> dict[str, str]:
    return dict(item.split("=", 1) for item in value.split(";") if "=" in item)


def key(fields: list[str]) -> tuple[str, int, str, str]:
    pos, ref, alt = int(fields[1]), fields[3].upper(), fields[4].upper()
    while len(ref) > 1 and len(alt) > 1 and ref[-1] == alt[-1]:
        ref, alt = ref[:-1], alt[:-1]
    while len(ref) > 1 and len(alt) > 1 and ref[0] == alt[0]:
        pos += 1
        ref, alt = ref[1:], alt[1:]
    return fields[0], pos, ref, alt


def validate_allele(fields: list[str]) -> None:
    if len(fields) < 8 or int(fields[1]) < 1 or not fields[3] or not fields[4]:
        raise ValueError("malformed native/candidate VCF record")
    if "," in fields[4]:
        raise ValueError("multiallelic native VCF requires a declared official ALT-specific conversion before projection")
    if not set((fields[3] + fields[4]).upper()) <= set("ACGTN"):
        raise ValueError("symbolic/breakend alleles require a declared sequence-resolving conversion")


def genotype(fields: list[str], sample_index: int) -> str:
    names, values = fields[8].split(":"), fields[sample_index].split(":")
    if len(set(names)) != len(names) or "GT" not in names or len(values) != len(names):
        raise ValueError("native record must contain a complete, unique GT FORMAT binding")
    gt = values[names.index("GT")]
    if gt == ".":
        return "./."
    if not re.fullmatch(r"[01.][/|][01.]", gt):
        raise ValueError("native GT must be diploid and valid for its single ALT; convert allele dosage explicitly")
    return gt


class Reference:
    def __init__(self, path: Path, index: Path):
        self.stream = path.open("rb")
        try:
            with index.open() as handle:
                self.index = {f[0]: tuple(map(int, f[1:5])) for f in (line.rstrip().split("\t") for line in handle)}
            if not self.index or any(len(row) != 4 or row[2] <= 0 or row[3] < row[2] for row in self.index.values()):
                raise ValueError("invalid reference FAI")
        except BaseException:
            self.stream.close()
            raise

    def fetch(self, chrom: str, start: int, end: int) -> str:
        if chrom not in self.index:
            raise ValueError("VCF contig is absent from the frozen reference")
        length, offset, bases, width = self.index[chrom]
        if not 0 <= start < end <= length:
            raise ValueError("reference interval is out of bounds")
        first = offset + start // bases * width + start % bases
        last = offset + (end - 1) // bases * width + (end - 1) % bases + 1
        self.stream.seek(first)
        seq = self.stream.read(last - first).replace(b"\n", b"").replace(b"\r", b"").decode("ascii").upper()
        if len(seq) != end - start:
            raise ValueError("reference FASTA/index mismatch")
        return seq

    def verify_ref(self, fields: list[str]) -> None:
        start = int(fields[1]) - 1
        if self.fetch(fields[0], start, start + len(fields[3])) != fields[3].upper():
            raise ValueError("VCF REF disagrees with the frozen reference")

    def equivalent(self, first: list[str], second: list[str]) -> bool:
        if first[0] != second[0]:
            return False
        left = min(int(first[1]), int(second[1])) - 1
        right = max(int(f[1]) - 1 + len(f[3]) for f in (first, second))
        if right - left > 1_000_000:
            return False
        seq = self.fetch(first[0], left, right)
        haplotypes = []
        for f in (first, second):
            offset = int(f[1]) - 1 - left
            if seq[offset:offset + len(f[3])] != f[3].upper():
                return False
            haplotypes.append(seq[:offset] + f[4].upper() + seq[offset + len(f[3]):])
        return haplotypes[0] == haplotypes[1]


def _headers(lines: list[str]) -> list[str]:
    """Merge native-first definitions without silently changing field types."""
    definitions: dict[str, str] = {}
    other: dict[str, None] = {}
    for line in lines:
        if line.startswith(("##fileformat=", "##source=")):
            continue
        match = re.match(r"##(INFO|FORMAT|FILTER|contig)=<ID=([^,>]+)", line)
        if match:
            label = ":".join(match.groups())
            if label in definitions:
                for name in ("Number", "Type", "length"):
                    old = re.search(rf"(?:<|,){name}=([^,>]+)", definitions[label])
                    new = re.search(rf"(?:<|,){name}=([^,>]+)", line)
                    if (old.group(1) if old else None) != (new.group(1) if new else None):
                        raise ValueError("conflicting VCF header definition: " + label)
            else:
                definitions[label] = line
        else:
            other[line] = None
    return list(definitions.values()) + list(other)


def project_all_sites(native: Path, candidates: Path, destination: Path, sample: str,
                      *, reference: Path | None = None, reference_index: Path | None = None,
                      native_policy: str = "official_default_or_recommended_no_extra_genotype_filtering") -> tuple[int, int]:
    if destination.suffix != ".vcf":
        raise ValueError("adapter output must be an uncompressed .vcf")
    if bool(reference) != bool(reference_index):
        raise ValueError("reference and its declared FAI must be supplied together")
    ref = Reference(reference, reference_index) if reference and reference_index else None
    panel: dict[str, list[str]] = {}
    by_key: dict[tuple[str, int, str, str], str] = {}
    by_allele: dict[str, str] = {}
    panel_headers: list[str] = []
    native_headers: list[str] = []
    selected: dict[str, tuple[list[str], str]] = {}
    trace_path = destination.with_suffix(".native-trace.jsonl.gz")
    audit_path = destination.with_suffix(".projection.json")
    destination.parent.mkdir(parents=True, exist_ok=True)
    native_count = outside = 0
    try:
        saw_panel_header = False
        with open_text(candidates, "rt") as stream:
            for line in stream:
                if line.startswith("##"):
                    panel_headers.append(line.rstrip()); continue
                if line.startswith("#CHROM"):
                    saw_panel_header = True; continue
                if line.startswith("#") or not line.strip():
                    continue
                if not saw_panel_header:
                    raise ValueError("candidate VCF lacks #CHROM header")
                f = line.rstrip("\r\n").split("\t")
                validate_allele(f)
                cid = f[2]
                if cid in {"", "."} or cid in panel or key(f) in by_key:
                    raise ValueError("candidate IDs and normalized allele keys must be unique")
                if ref:
                    ref.verify_ref(f)
                panel[cid] = f
                by_key[key(f)] = cid
                allele = info(f[7]).get("PANGENOME_ALLELE_ID")
                if allele:
                    if allele in by_allele:
                        raise ValueError("ambiguous pangenome allele ID")
                    by_allele[allele] = cid
        if not panel:
            raise ValueError("candidate panel contains no records")
        sample_index = None
        column_count = None
        with open_text(native, "rt") as stream, gzip.open(trace_path, "wt", encoding="utf-8") as trace:
            for line_number, line in enumerate(stream, 1):
                if line.startswith("##"):
                    native_headers.append(line.rstrip()); continue
                if line.startswith("#CHROM"):
                    columns = line.rstrip("\r\n").split("\t")
                    if sample_index is not None or columns[9:].count(sample) != 1 or len(set(columns[9:])) != len(columns[9:]):
                        raise ValueError("native VCF must identify the requested sample exactly once")
                    sample_index = columns.index(sample, 9)
                    column_count = len(columns)
                    continue
                if line.startswith("#") or not line.strip():
                    continue
                if sample_index is None:
                    raise ValueError("native VCF lacks the requested sample header")
                f = line.rstrip("\r\n").split("\t")
                if len(f) != column_count:
                    raise ValueError("native VCF record differs from its column header")
                validate_allele(f)
                gt = genotype(f, sample_index)
                native_count += 1
                claimed = set()
                if f[2] in panel:
                    claimed.add(f[2])
                allele = info(f[7]).get("PANGENOME_ALLELE_ID")
                if allele in by_allele:
                    claimed.add(by_allele[allele])
                sequence_hit = by_key.get(key(f))
                # Validate only candidate proposals, not millions of outside
                # records. Their full native evidence is still streamed to trace.
                if ref and (claimed or sequence_hit is not None):
                    ref.verify_ref(f)
                for cid in claimed:
                    if not (ref.equivalent(f, panel[cid]) if ref else key(f) == key(panel[cid])):
                        raise ValueError("claimed candidate ID conflicts with its native allele: " + cid)
                hits = claimed | ({sequence_hit} if sequence_hit is not None else set())
                if len(hits) > 1:
                    raise ValueError("ambiguous native-to-candidate mapping")
                cid = next(iter(hits), None)
                method = "verified_identity" if claimed else "normalized_sequence"
                if cid is None:
                    outside += 1
                    method = "outside_canonical_universe"
                else:
                    if cid in selected:
                        raise ValueError("native output duplicates candidate " + cid)
                    out = f[:8] + [f[8], f[sample_index]]
                    out[2] = cid
                    names, values = out[8].split(":"), out[9].split(":")
                    values[names.index("GT")] = gt
                    out[9] = ":".join(values)
                    if any(k.startswith("PGBENCH_") for k in info(out[7])):
                        raise ValueError("native output uses reserved PGBENCH audit annotations")
                    suffix = f"PGBENCH_NATIVE_RECORD_ID={quote(f[2], safe='')};PGBENCH_PROJECTION_METHOD={method};PGBENCH_NATIVE_STATUS=observed"
                    out[7] = (out[7] + ";" if out[7] not in {"", "."} else "") + suffix
                    selected[cid] = out, gt
                # Stream every native record, including all quality/FORMAT fields.
                trace.write(json.dumps({"native_record_line": line_number, "candidate_id": cid,
                    "match_strategy": method, "native_record": f, "sample_index": sample_index,
                    "native_gt": gt}, separators=(",", ":")) + "\n")
        if sample_index is None:
            raise ValueError("native VCF lacks the requested sample header")
        headers = _headers(native_headers + panel_headers)
        definitions = {line.split(",", 1)[0] for line in headers}
        with open_text(destination, "wt") as stream:
            stream.write("##fileformat=VCFv4.2\n")
            stream.writelines(line + "\n" for line in headers)
            stream.write("##source=PGBench-reviewed-adapter-projection\n")
            if "##FORMAT=<ID=GT" not in definitions:
                stream.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
            for name in ("PGBENCH_NATIVE_RECORD_ID", "PGBENCH_PROJECTION_METHOD", "PGBENCH_NATIVE_STATUS"):
                stream.write(f'##INFO=<ID={name},Number=1,Type=String,Description="Truth-blind native output trace">\n')
            stream.write(f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}\n")
            for cid, candidate in panel.items():
                if cid in selected:
                    out = selected[cid][0]
                else:
                    out = candidate[:8] + ["GT", "./."]
                    out[5:7] = [".", "."]
                    out[7] = (out[7] + ";" if out[7] not in {"", "."} else "") + "PGBENCH_NATIVE_STATUS=not_observed"
                stream.write("\t".join(out) + "\n")
        no_call = len(panel) - len(selected) + sum("." in gt for _, gt in selected.values())
        audit = {"contract": "pgbench-native-projection-v2", "truth_used": False,
            "reference_equivalence_checked": ref is not None, "caller_alleles_preserved": True,
            "reference_equivalence_scope": "mapped_candidate_alleles",
            "quality_and_filter_preserved": True, "candidate_count": len(panel),
            "matched_candidates": len(selected), "no_call_candidates": no_call,
            "not_observed_candidates": len(panel) - len(selected), "native_records": native_count,
            "outside_canonical_universe": outside, "native_vcf": str(native.resolve()),
            "native_sha256": digest(native), "candidate_sha256": digest(candidates),
            "reference_sha256": digest(reference) if reference else None,
            "reference_fai_sha256": digest(reference_index) if reference_index else None,
            "output_sha256": digest(destination), "trace_path": str(trace_path.resolve()),
            "trace_sha256": digest(trace_path), "native_execution_policy": native_policy,
            "unsupported_representation_policy": "fail_entire_run_not_no_call",
            "biological_performance_validated": False}
        audit_path.write_text(json.dumps(audit, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        return len(selected), no_call
    finally:
        if ref:
            ref.stream.close()
