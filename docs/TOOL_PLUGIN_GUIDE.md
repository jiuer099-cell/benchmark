# Tool adapter contract

Each adapter lives under `plugins/<tool>/` and contains `tool.yaml`, `run.py`,
an environment definition, and a rule registry. The runner receives only
allowlisted `PGBENCH_*` variables and writes only its assigned output tree.

For the formal track, `supported_modes.end_to_end_from_reads.required_inputs`
must contain `short_fastq_r1`, `short_fastq_r2`, `candidate_panel`, reference,
and panel provenance. Additional indexes or graph bundles must be registered as
`tool_index` or `graph_assets`; their complete directory hashes enter the
information contract.

The required output is an all-sites VCF over the supplied candidate IDs. A
missing native record becomes `./.` plus `missing_output`; it must never be
invented as 0/0. Conversion, unsupported representation, index failure, and
linking failure must remain distinct. The core then canonicalizes, validates,
links alleles, produces the addressability audit, and submits one identical
query VCF to Truvari, Aardvark-GT, and vcfdist.

`information_contract` must explicitly state that HG002 truth, HG002 assembly,
family genotypes, and target-specific external calls were not used. Use either
an official configuration frozen before the test or tuning performed only on
an independent validation sample/region.

An adapter is bundled only when its runner and environment are maintained and
reviewed in this repository. Community adapters are executed through the same
interface and receive no scoring exceptions.
