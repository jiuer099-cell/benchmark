# Variant: LONG reads — ONT R9.4.1 Guppy5 SUP pass (lr_ont track)

Track: `lr_ont` (config/track_registry.yaml), release `PGBench-LR-ONT-v1.0`.
This track has the strictest input policy in the benchmark.

## Frozen dataset (binding)

Every ONT adapter consumes this exact release, declared in
`config/track_registry.yaml`:

- `dataset_id: HG002_ONT_R9.4.1_Guppy5.0.6_SUP_pass_v1`
- `object_name: basecalls.fastq.gz`
- `read_selection_policy: official_qscore_pass_only`
- `prohibited_read_selection: [pass_plus_fail, fail_only,
  adapter_selected_subset]`

## tool.yaml differences

- `paradigm: long_read_pangenome_genotyping`
- `capabilities.read_class: long`, `capabilities.technology: [ont]`
- `required_inputs`: `[long_reads_fastq, reference, pangenome_manifest,
  pangenome_panel, candidate_panel]`

## Environment variables

| variable | meaning |
| --- | --- |
| `PGBENCH_INPUT_LONG_READS_FASTQ` | frozen `HG002.ONT...chr1-22.fastq.gz` |
| others | same plumbing as the other long-read variants |

## ONT-specific notes

- NEVER re-filter, trim, or subset the pass reads inside the adapter. The
  official fail partition is a benchmark-level decision; merging pass+fail or
  picking "high quality" subsets is prohibited (`adapter_selected_subset`).
- If your tool consumes a BAM, align the frozen FASTQ yourself with the ONT
  preset (`minimap2 -x map-ont`) inside the billed run and write a proper
  `@RG` header (`SM: HG002`) — see the frozen BAM reheader incident for why
  missing read groups fail closed.
- The SUP basecalling model produces different indel error profiles than
  FAST/HC; use callers/configs validated for SUP data, and record them under
  `parameter_contract`.
