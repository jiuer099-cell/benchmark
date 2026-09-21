#!/usr/bin/env python3
"""Make a representation-preserving evaluator query with stable unit traces.

Unlike the all-sites ledger, records are not rewritten into canonical truth
alleles.  The Core only adds a trace ID after deterministic candidate mapping;
allele equivalence remains exclusively the evaluators' decision.
"""
from __future__ import annotations
import argparse,sys
from pathlib import Path
from pgf1 import (
    PGF1Error,
    canonical_alt_dosage,
    json_dump,
    load_units,
    open_text,
    sha256_file,
)


def info_fields(value: str) -> dict[str, str]:
    parsed: dict[str, str] = {}
    for item in value.split(";"):
        if "=" in item:
            key, field_value = item.split("=", 1)
            parsed[key] = field_value
    return parsed


def main(argv:list[str]|None=None)->int:
    p=argparse.ArgumentParser();p.add_argument("--linked-vcf",required=True,type=Path);p.add_argument("--canonical-units",required=True,type=Path);p.add_argument("--output-vcf",required=True,type=Path);p.add_argument("--audit-json",required=True,type=Path);a=p.parse_args(argv);seen=set()
    try:
        units=load_units(a.canonical_units)
        source_records=submitted_records=reference_records=no_call_records=0
        a.output_vcf.parent.mkdir(parents=True,exist_ok=True)
        with open_text(a.linked_vcf,"rt") as src,open_text(a.output_vcf,"wt") as dst:
            for line in src:
                if line.startswith("##"):dst.write(line);continue
                if line.startswith("#CHROM"):
                    dst.write('##PGBENCH_EVALUATOR_QUERY=representation_preserving_trace_v1\n')
                    dst.write('##INFO=<ID=PGBENCH_UNIT_TRACE,Number=1,Type=String,Description="ALT-specific canonical scoring unit">\n')
                    dst.write('##INFO=<ID=PGBENCH_NATIVE_RECORD_ID,Number=1,Type=String,Description="Caller record ID before trace relabeling">\n')
                    dst.write(line);continue
                if not line.strip():continue
                source_records+=1
                f=line.rstrip("\r\n").split("\t")
                if len(f)<10 or f[2] in {"","."}:raise PGF1Error("CANONICAL_MAPPING_FAILURE")
                gt=f[9].split(":")[f[8].split(":").index("GT")] if "GT" in f[8].split(":") else "./."
                dosage=canonical_alt_dosage(gt)
                if dosage is None:no_call_records+=1;continue
                if dosage==0:reference_records+=1;continue
                info=info_fields(f[7]); candidate=info.get("PANGENOME_LINKED_ID","")
                if not candidate or info.get("PANGENOME_CANDIDATE_COUNT")!="1" or candidate in seen:raise PGF1Error("CANONICAL_MAPPING_FAILURE")
                unit_id=f"{candidate}:A1"
                if unit_id not in units:raise PGF1Error("UNIT_SET_MISMATCH")
                native_record_id=f[2]; seen.add(candidate); f[2]=candidate
                suffix=f"PGBENCH_UNIT_TRACE={unit_id};PGBENCH_NATIVE_RECORD_ID={native_record_id}"
                f[7]=(f[7]+";" if f[7] not in {"","."} else "")+suffix;dst.write("\t".join(f)+"\n")
                submitted_records+=1
        json_dump(a.audit_json,{
            "audit_contract":"pgbench_representation_preserving_evaluator_query_v1",
            "linked_vcf_sha256":sha256_file(a.linked_vcf),
            "canonical_unit_set_sha256":sha256_file(a.canonical_units),
            "evaluator_query_sha256":sha256_file(a.output_vcf),
            "source_record_count":source_records,
            "submitted_record_count":submitted_records,
            "reference_record_count":reference_records,
            "no_call_record_count":no_call_records,
            "unit_trace_complete":True,
            "native_record_trace_complete":True,
            "caller_representation_rewritten":False,
        })
    except (OSError,PGF1Error,ValueError) as e:print(f"materialize_pgf1_evaluator_query: {e}",file=sys.stderr);return 2
    return 0
if __name__=="__main__":raise SystemExit(main())
