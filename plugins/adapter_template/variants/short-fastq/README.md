# Variant: SHORT reads via raw FASTQ (illumina_pe)

Use when your tool consumes raw reads: it may map them (graph tools) or use
k-mer inference. Bill the actual stages; do not invent a mapping step for a
k-mer-only method. External adapters may live outside Core via their manifest.

## tool.yaml differences

- `capabilities.technology: [illumina_pe]`, `read_class: short`
- `information_inputs` includes `target_short_reads`
- `required_inputs`: `[short_fastq_r1, short_fastq_r2, reference,
  reference_index, pangenome_manifest, pangenome_panel, candidate_panel]`
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
- Freeze official defaults/recommendations before running; only documented
  input-compatibility/resource/necessary-conversion exceptions are permitted.
  Add no extra genotype filtering/rejection. One formal PG-F1 per tool/context.
- If you build a graph index per run, keep it in `tool-work/` and account for
  the build time in `billable_stages` (it is part of the honest cost).
