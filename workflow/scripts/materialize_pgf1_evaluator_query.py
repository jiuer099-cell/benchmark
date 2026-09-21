#!/usr/bin/env python3
"""Make a representation-preserving evaluator query with stable unit traces.

Unlike the all-sites ledger, records are not rewritten into canonical truth
alleles.  The Core only adds a trace ID after deterministic candidate mapping;
allele equivalence remains exclusively the evaluators' decision.
"""
from __future__ import annotations
import argparse,sys
from pathlib import Path
from pgf1 import PGF1Error,open_text
def main(argv:list[str]|None=None)->int:
    p=argparse.ArgumentParser();p.add_argument("--linked-vcf",required=True,type=Path);p.add_argument("--output-vcf",required=True,type=Path);a=p.parse_args(argv);seen=set()
    try:
        a.output_vcf.parent.mkdir(parents=True,exist_ok=True)
        with open_text(a.linked_vcf,"rt") as src,open_text(a.output_vcf,"wt") as dst:
            for line in src:
                if line.startswith("##"):dst.write(line);continue
                if line.startswith("#CHROM"):dst.write('##PGBENCH_EVALUATOR_QUERY=representation_preserving_trace_v1\n');dst.write(line);continue
                if not line.strip():continue
                f=line.rstrip("\r\n").split("\t")
                if len(f)<10 or f[2] in {"","."} or f[2] in seen:raise PGF1Error("CANONICAL_MAPPING_FAILURE")
                gt=f[9].split(":")[f[8].split(":").index("GT")] if "GT" in f[8].split(":") else "./."
                if gt.replace("|","/") in {"0/0","./."}:continue
                seen.add(f[2]); f[7]=(f[7]+";" if f[7] not in {"","."} else "")+"PGBENCH_UNIT_TRACE="+f[2]+":A1";dst.write("\t".join(f)+"\n")
    except (OSError,PGF1Error,ValueError) as e:print(f"materialize_pgf1_evaluator_query: {e}",file=sys.stderr);return 2
    return 0
if __name__=="__main__":raise SystemExit(main())
