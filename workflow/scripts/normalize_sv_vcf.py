#!/usr/bin/env python3
"""Conservative single-sample SV VCF normalization for PGBench."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import sys
from pathlib import Path
from typing import TextIO, cast

from build_pangenome_manifest import format_info, infer_svtype, parse_info


class VcfNormalizationError(ValueError):
    """Raised when a tool VCF cannot be normalized without guessing."""


INFO_HEADERS = (
    "##INFO=<ID=CANON_ID,Number=1,Type=String,"
    'Description="Stable PGBench canonical record ID">',
    '##INFO=<ID=ORIG_ID,Number=1,Type=String,Description="Original tool record ID">',
    "##INFO=<ID=ORIG_ALT,Number=1,Type=String,"
    'Description="Original ALT representation">',
    '##INFO=<ID=ORIG_SVTYPE,Number=1,Type=String,Description="Original SVTYPE">',
)


def _open_text(path: Path, mode: str) -> TextIO:
    if path.suffix == ".gz":
        return cast(TextIO, gzip.open(path, mode, encoding="utf-8"))
    return cast(TextIO, path.open(mode, encoding="utf-8"))


def reference_contig_order(reference: Path) -> dict[str, int]:
    order: dict[str, int] = {}
    with reference.open("r", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(">"):
                contig = line[1:].strip().split()[0]
                if contig in order:
                    raise VcfNormalizationError(
                        f"duplicate reference contig {contig!r}"
                    )
                order[contig] = len(order)
    if not order:
        raise VcfNormalizationError("reference FASTA contains no contigs")
    return order


def _parse_integer(info: dict[str, str | bool], key: str, *, default: int) -> int:
    value = info.get(key)
    if value is None or value is True:
        return default
    try:
        return int(str(value).split(",", 1)[0])
    except ValueError as exc:
        raise VcfNormalizationError(f"{key} must be integer, found {value!r}") from exc


def _normalized_sample(fields: list[str]) -> tuple[str, str]:
    if len(fields) < 10:
        return "GT", "./."
    format_keys = fields[8].split(":") if fields[8] not in {"", "."} else []
    sample_values = fields[9].split(":") if fields[9] not in {"", "."} else []
    if "GT" not in format_keys:
        format_keys.insert(0, "GT")
        sample_values.insert(0, "./.")
    if len(sample_values) < len(format_keys):
        sample_values.extend(["."] * (len(format_keys) - len(sample_values)))
    return ":".join(format_keys), ":".join(sample_values[: len(format_keys)])


def _canonical_id(
    *,
    chrom: str,
    pos: int,
    ref: str,
    alt: str,
    svtype: str,
    end: int,
    svlen: int,
    original_id: str,
) -> str:
    payload = json.dumps(
        {
            "alt": alt,
            "chrom": chrom,
            "end": end,
            "original_id": original_id,
            "pos": pos,
            "ref": ref,
            "svlen": svlen,
            "svtype": svtype,
        },
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8")
    return f"CANON_{hashlib.sha256(payload).hexdigest()[:20]}"


def normalize_vcf(
    input_vcf: Path,
    reference: Path,
    output_vcf: Path,
    *,
    sample_name: str = "HG002",
) -> int:
    order = reference_contig_order(reference)
    meta_headers: list[str] = []
    records: list[tuple[int, int, str, list[str]]] = []
    saw_column_header = False

    with _open_text(input_vcf, "rt") as handle:
        for raw_line in handle:
            if raw_line.startswith("##"):
                meta_headers.append(raw_line.rstrip("\n"))
                continue
            if raw_line.startswith("#CHROM"):
                saw_column_header = True
                continue
            if not raw_line.strip():
                continue
            if not saw_column_header:
                raise VcfNormalizationError("VCF is missing #CHROM header")
            fields = raw_line.rstrip("\n").split("\t")
            if len(fields) < 8:
                raise VcfNormalizationError("VCF record has fewer than 8 columns")
            chrom, pos_raw, original_id, ref, alt = fields[:5]
            if chrom not in order:
                raise VcfNormalizationError(
                    f"VCF contig {chrom!r} is absent from reference"
                )
            try:
                pos = int(pos_raw)
            except ValueError as exc:
                raise VcfNormalizationError(f"invalid POS {pos_raw!r}") from exc
            info = parse_info(fields[7])
            original_svtype = str(info.get("SVTYPE", "."))
            svtype = infer_svtype(alt, info)
            end = _parse_integer(info, "END", default=pos)
            svlen = _parse_integer(info, "SVLEN", default=len(alt) - len(ref))
            info["SVTYPE"] = svtype
            info["END"] = str(end)
            info["SVLEN"] = str(svlen)
            info["ORIG_ID"] = original_id
            info["ORIG_ALT"] = alt
            info["ORIG_SVTYPE"] = original_svtype
            canon_id = _canonical_id(
                chrom=chrom,
                pos=pos,
                ref=ref,
                alt=alt,
                svtype=svtype,
                end=end,
                svlen=svlen,
                original_id=original_id,
            )
            info["CANON_ID"] = canon_id
            format_field, sample_field = _normalized_sample(fields)
            normalized = fields[:8] + [format_field, sample_field]
            normalized[2] = canon_id
            normalized[7] = format_info(info)
            records.append((order[chrom], pos, canon_id, normalized))

    if not saw_column_header:
        raise VcfNormalizationError("VCF is missing #CHROM header")

    records.sort(key=lambda item: (item[0], item[1], item[2]))
    output_vcf.parent.mkdir(parents=True, exist_ok=True)
    temporary = output_vcf.with_name(f".{output_vcf.name}.tmp")
    with temporary.open("w", encoding="utf-8") as output:
        existing_info_ids = {
            line.split("ID=", 1)[1].split(",", 1)[0]
            for line in meta_headers
            if line.startswith("##INFO=<ID=")
        }
        for line in meta_headers:
            output.write(line + "\n")
        for header in INFO_HEADERS:
            info_id = header.split("ID=", 1)[1].split(",", 1)[0]
            if info_id not in existing_info_ids:
                output.write(header + "\n")
        output.write(
            f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample_name}\n"
        )
        for _, _, _, fields in records:
            output.write("\t".join(fields) + "\n")
    temporary.replace(output_vcf)
    return len(records)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Normalize a single-sample SV VCF.")
    parser.add_argument("--input-vcf", required=True, type=Path)
    parser.add_argument("--reference", required=True, type=Path)
    parser.add_argument("--output-vcf", required=True, type=Path)
    parser.add_argument("--sample-name", default="HG002")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        normalize_vcf(
            args.input_vcf,
            args.reference,
            args.output_vcf,
            sample_name=args.sample_name,
        )
    except (OSError, VcfNormalizationError) as exc:
        print(f"normalize_sv_vcf: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
