# Variant: LONG reads — PacBio CLR (lr_clr track)

Track: `lr_clr` (config/track_registry.yaml), release `PGBench-LR-CLR-v1.0`.

## tool.yaml differences

- `paradigm: long_read_pangenome_genotyping`
- `capabilities.read_class: long`, `capabilities.technology: [pacbio_clr]`
- `required_inputs`: `[long_reads_fastq, reference, pangenome_manifest,
  pangenome_panel, candidate_panel]`

## Environment variables

Same as the HiFi variant: the frozen CLR FASTQ arrives as
`PGBENCH_INPUT_LONG_READS_FASTQ`.

## CLR-specific notes

- Align with the CLR error model, NOT the HiFi preset:
  `minimap2 -x map-pb` (or pbmm2 `--preset CLR`). Using a HiFi/Oxford preset
  on CLR reads is a silent correctness bug — the alignment looks plausible
  while indel handling is wrong.
- CLR alignments are larger and noisier; if your adapter indexes per-run,
  declare the index stage as billable and size `tool_memory_mb` from a
  measured pilot, not from documentation.
- Do not quality-trim or subset the reads inside the adapter: read selection
  belongs to the benchmark, and an adapter-selected subset is exactly the
  failure mode the ONT policy (see variants/long-ont) was written to prevent.
