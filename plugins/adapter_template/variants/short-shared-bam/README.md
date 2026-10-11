# Variant: SHORT reads via the frozen shared BAM (illumina_pe)

Use when your tool consumes an alignment file directly (Paragraph,
GraphTyper2 style). The benchmark maps the Illumina reads once, freezes the
BAM, and every tool genotypes the identical alignment — coverage, mapping
bias, and read selection are held fixed across tools.

## tool.yaml differences

- `capabilities.technology: [illumina_pe]`, `read_class: short`
- `information_inputs` must include `benchmark_supplied_shared_alignment`
- `supported_modes.end_to_end_from_reads.required_inputs`:
  `[shared_shortread_alignment, shared_shortread_alignment_index, reference,
  reference_index, pangenome_manifest, pangenome_panel, candidate_panel]`
- `billable_stages` typically `[genotype, postprocess]` (mapping is shared
  and not re-billed).

## Environment variables injected by the core

| variable | meaning |
| --- | --- |
| `PGBENCH_SHARED_ALIGNMENT_BAM` | frozen shared BAM |
| `PGBENCH_SHARED_ALIGNMENT_BAI` | frozen BAI index |
| `PGBENCH_REFERENCE_FASTA` | GRCh38 reference |
| `PGBENCH_CANDIDATE_VCF` | frozen canonical candidate panel (all-sites input) |
| `PGBENCH_SAMPLE_ID` / `PGBENCH_THREADS` | sample and thread budget |
| `PGBENCH_OUTPUT_DIR` / `PGBENCH_OUTPUT_VCF` | where the projected VCF goes |

Core owns generic BAM/BAI integrity validation. The template only confirms
that Core injected non-empty read-only paths. If a tool requires a matching
`@RG`/`SM`, add that *tool-private prerequisite* to the copied adapter; do not
replace or bypass the Core validation contract.

## Pitfalls learned in this repo

- Choose a reviewed method-valid memory strategy; chunking must preserve
  complete evidence and be declared as a resource exception before execution.
- Freeze official defaults/recommendations. Do not introduce extra genotype
  filtering/rejection or select among multiple scores; one formal PG-F1 only.
- Do not re-map or re-filter the shared reads; if your tool needs different
  alignment parameters, that is a different benchmark question and must be
  declared, not improvised.
