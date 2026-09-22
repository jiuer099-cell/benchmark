# Variant: SHORT reads via raw FASTQ (illumina_pe)

Use when your tool performs its own mapping (vg, BayesTyper, varigraph
style). Your adapter maps the benchmark-supplied reads inside the billed run.

## tool.yaml differences

- `capabilities.technology: [illumina_pe]`, `read_class: short`
- `information_inputs` includes `target_short_reads`
- `required_inputs`: `[short_fastq_r1, short_fastq_r2, reference,
  pangenome_manifest, pangenome_panel, candidate_panel]`
- `optional_inputs` may include `reference_index`, `tool_index`,
  `graph_assets` (declare `accepted_formats`/`accepted_profiles` if used).
- `billable_stages` typically `[index, map, genotype, postprocess]`.

## Environment variables injected by the core

| variable | meaning |
| --- | --- |
| `PGBENCH_INPUT_FASTQ_R1` / `PGBENCH_INPUT_FASTQ_R2` | paired reads |
| `PGBENCH_REFERENCE_FASTA` (and `PGBENCH_REFERENCE_INDEX` if declared optional) | reference |
| `PGBENCH_CANDIDATE_VCF` / `PGBENCH_PANGENOME_MANIFEST` / `PGBENCH_PANEL_VCF` | panel inputs |
| `PGBENCH_GRAPH_DIR` (only if `graph_assets` declared) | benchmark-supplied graph |
| `PGBENCH_INDEX_DIR` (only if `tool_index` declared) | benchmark-supplied tool index |
| `PGBENCH_SAMPLE_ID` / `PGBENCH_THREADS` / `PGBENCH_OUTPUT_DIR` / `PGBENCH_OUTPUT_VCF` | run plumbing |

Core validates FASTQ content once through its immutable shared cache before
the adapter starts. The template only confirms that Core injected a non-empty
read-only path; it deliberately does not reimplement gzip/SHA256/read-count
validation.

## Pitfalls

- Map with your tool's OFFICIAL recommended parameters; record them under
  `parameter_contract` (explicit/implicit) — tuning modes are audited.
- If you build a graph index per run, keep it in `tool-work/` and account for
  the build time in `billable_stages` (it is part of the honest cost).
