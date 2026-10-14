"""Prepare a lossless, site-merged BayesTyper input from blinded candidates."""
from __future__ import annotations

from collections import OrderedDict
import gzip
from pathlib import Path


def prepare_candidates(panel: Path, destination: Path) -> dict[str, int]:
    opener = gzip.open if panel.suffix == ".gz" else open
    groups = OrderedDict()
    identifiers = set()
    headers = []
    closed_contigs = set()
    previous_chrom = None
    previous_pos = 0
    with opener(panel, "rt") as source:
        for line in source:
            if line.startswith("##"):
                headers.append(line)
                continue
            if line.startswith("#") or not line.strip():
                continue
            row = line.rstrip("\r\n").split("\t")
            # Core's frozen blinded VCF has a GT-only, entirely missing sample.
            # Drop only this empty transport field for native sites-only input;
            # a called or partial genotype must never reach native inference.
            if len(row) == 10:
                if row[8:] != ["GT", "./."]:
                    raise ValueError("BayesTyper candidate sample must be a blinded GT-only ./. placeholder")
                row = row[:8]
            if len(row) != 8 or row[2] in {"", "."} or row[2] in identifiers or ";" in row[2] or "," in row[4]:
                raise ValueError("BayesTyper requires a blinded, unique-ID biallelic candidate panel")
            chrom, pos = row[0], int(row[1])
            if chrom != previous_chrom:
                if chrom in closed_contigs:
                    raise ValueError("Candidate panel is not sorted by contig")
                if previous_chrom is not None:
                    closed_contigs.add(previous_chrom)
                previous_chrom, previous_pos = chrom, 0
            if pos < 1 or pos < previous_pos:
                raise ValueError("Candidate panel is not sorted by position")
            previous_pos = pos
            identifiers.add(row[2])
            groups.setdefault((chrom, pos), []).append(row)
    if not groups:
        raise ValueError("Empty candidate panel")

    with destination.open("wb") as file:
        with gzip.GzipFile(filename="", mode="wb", fileobj=file, mtime=0) as raw:
            for line in headers:
                raw.write(line.encode())
            raw.write(b"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
            for rows in groups.values():
                if len(rows) == 1:
                    row = rows[0]
                else:
                    ref = max((r[3] for r in rows), key=len)
                    alts = []
                    for r in rows:
                        if not ref.startswith(r[3]):
                            raise ValueError("Inconsistent REF sequences at the same candidate position")
                        alt = r[4] + ref[len(r[3]):]
                        if alt not in alts:
                            alts.append(alt)
                    row = [rows[0][0], rows[0][1], ";".join(r[2] for r in rows),
                           ref, ",".join(alts), ".", ".", "."]
                raw.write(("\t".join(row) + "\n").encode())
    return {"candidate_count": len(identifiers), "native_site_count": len(groups),
            "merged_sites": sum(len(rows) > 1 for rows in groups.values())}
