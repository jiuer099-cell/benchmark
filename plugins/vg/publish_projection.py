"""Publish a validated native projection without altering its source attempt."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil


def sha(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--source-vcf", type=Path, required=True)
    p.add_argument("--destination-vcf", type=Path, required=True)
    a = p.parse_args()
    audit = a.source_vcf.with_suffix(".projection.json")
    trace = a.source_vcf.with_suffix(".native-trace.jsonl.gz")
    data = json.loads(audit.read_text())
    if data["output_sha256"] != sha(a.source_vcf) or data["trace_sha256"] != sha(trace):
        raise ValueError("projection publication source hash mismatch")
    for src, dst in ((a.source_vcf, a.destination_vcf), (audit, a.destination_vcf.with_suffix(".projection.json")), (trace, a.destination_vcf.with_suffix(".native-trace.jsonl.gz"))):
        dst.parent.mkdir(parents=True, exist_ok=True)
        temp = dst.with_name(".publish." + dst.name)
        shutil.copyfile(src, temp)
        if sha(src) != sha(temp): raise ValueError("projection copy hash mismatch")
        temp.replace(dst)


if __name__ == "__main__": main()
