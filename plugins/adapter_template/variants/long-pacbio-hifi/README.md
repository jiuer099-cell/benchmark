# Variant: LONG reads — PacBio HiFi (lr_hifi track)

Track: `lr_hifi` (config/track_registry.yaml), release `PGBench-LR-HiFi-v1.0`.

## tool.yaml differences

- `paradigm: long_read_pangenome_genotyping`
- `capabilities.read_class: long`, `capabilities.technology: [pacbio_hifi]`
- `information_inputs` includes `target_long_reads`
- `required_inputs`: `[long_reads_fastq, reference, pangenome_manifest,
  pangenome_panel, candidate_panel]`
- `billable_stages` typically `[index, map, call, genotype, postprocess]`.

## Environment variables

| variable | meaning |
| --- | --- |
| `PGBENCH_INPUT_LONG_READS_FASTQ` | the frozen HiFi FASTQ (single file) |
| others | same plumbing as the short-read variants |

## HiFi-specific notes

- HiFi callers (e.g. pbsv, DeepVariant) want a per-run alignment; align the
  frozen FASTQ yourself with the official preset (`minimap2 -x map-hifi` or
  pbmm2 `--preset HIFI`) inside the billed run.
- Keep the read selection untouched: the benchmark supplies the reads, and
  adapters must not subset them.
- SV genotypers that consume a VCF of calls (not reads) should benchmark on
  the `genotyping_only` paradigm and declare the callset provenance honestly
  (undeclared external callsets are forbidden by the information contract).
