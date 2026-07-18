# HG002 GRCh37 Pangenome SV Benchmark

This repository implements the confirmed design for a pangenome-first,
Snakemake-based structural-variant detection and genotyping benchmark.

The benchmark runs built-in or user-provided tools on HG002 inputs. It does
not accept a precomputed VCF as a formal external-tool submission. Every
eligible tool/run tuple receives one `PGBenchScore` in `[0, 100]`; the report
does not rank tools.

## Current implementation status

Phase 1 is now executable and verified with a synthetic HG002-like fixture:

- 14 executable workflow rules plus the aggregation target complete
  successfully.
- The benchmark launches an external genotyper package itself; no submitted
  VCF is accepted.
- Every executed rule writes a content-hashed manifest, log, and Snakemake
  JSONL benchmark record.
- The pre-score audit verifies a complete eight-job hash lineage before the
  score is calculated.
- The finalizer replays that audit, verifies score/audit manifests, and
  atomically seals a score package with an 11-job final lineage.
- Top-level `score.tsv`, `point_breakdown.tsv`, `metrics.long.tsv`, and
  `metrics.json` preserve tool/config order and never calculate ranks.
- Incremental reruns archive a previously validated VCF/index and manifest;
  failed attempts restore the last trusted outputs and cannot silently reuse
  stale artifacts.
- The whole external-plugin package is a DAG input, and each plugin may
  declare a nested/custom raw VCF path rather than relying on a fixed layout.
- Score finalization rejects mixed run/config/profile/truth contexts and binds
  every report and top-level summary to the sealed package by content hash.
- The regression suite contains 177 tests, including custom plugin layout,
  ratio semantics, context-mixing, report-binding, and failed-rerun recovery.
- The synthetic result is deliberately marked `provisional`; it is a contract
  test, not an HG002 biological result.

Remaining phases:

- Phase 2: real Truvari/Aardvark/vcfdist wrappers plus the first built-in
  callers/genotypers and small-region GRCh37 validation.
- Phase 3: frozen HG002/hs37d5/truth/pangenome assets, locked HPC environments,
  and whole-genome acceptance.

The default production configuration is fail-closed: until Phase 2 assets and
evaluators are supplied, `rule all` performs configuration validation only and
cannot emit a formal numeric score.

The authoritative design is:

`docs/superpowers/specs/2026-07-16-hg002-grch37-sv-benchmark-design.md`

## Environment

Create the workflow environment:

```bash
conda env create -f environment.yaml
conda activate pgbench-sv
```

Run unit tests:

```bash
python -m pytest -q
```

Validate and run the Phase 1 DAG with:

```bash
snakemake \
  --snakefile Snakefile \
  --configfile tests/fixtures/synthetic/config.yaml \
  --cores 1 \
  --dry-run

snakemake \
  --snakefile Snakefile \
  --configfile tests/fixtures/synthetic/config.yaml \
  --cores 2 \
  --rerun-incomplete
```

The environment pins `snakemake-minimal>=9.23.1,<10`. This baseline supports
dynamically named modules, JSONL benchmark records, and update-output failure
recovery; the exact runtime version is recorded in provenance.

## Formal modes

- `caller_only_shared_alignment`: runs the caller/genotyper on a registered
  shared BAM/GAF/GAM. The benchmark has already paid for alignment, so this
  mode isolates calling/genotyping behavior and is diagnostic by default.
- `end_to_end_from_reads`: starts from canonical reads and includes
  tool-specific indexing, mapping/assembly, calling/genotyping, and required
  postprocessing. It answers how the complete user-facing pipeline performs
  and is the default official score mode.

External plugins receive only the inputs allowed by the selected mode.

## External tools

An external tool is provided as executable code under `plugins/`, with a
validated `tool.yaml` and runner script. The benchmark invokes it during the
current run and validates the resulting VCF. Formal untrusted execution
requires a declared bwrap or Apptainer sandbox; the included unsandboxed mock
is development-only.

See `plugins/example_genotyper/` for the executable interface.

## Score

Each applicable tuple
`(run_id, sample, tool, official_score_mode, primary_truth_profile)` receives
one score and no rank:

- evaluator accuracy: 70 points
  (`Truvari:Aardvark:vcfdist = 28:24.5:17.5`, equivalent to `40:35:25`);
- pangenome robustness: 15 points;
- resource efficiency: 10 points;
- traceability: 5 points.

The exact metric-to-point dictionary is frozen in
`config/score_weights.yaml`. Required F1 metrics, genotype concordance,
allele linkage, context strata, resources, and provenance gates are all
reported in the point breakdown.

Phase 1 writes the following primary artifacts:

```text
results/summary/score.tsv
results/summary/point_breakdown.tsv
results/summary/metrics.long.tsv
results/summary/metrics.json
results/summary/<tool>/score-package.json
results/provenance/<tool>/rule-lineage.{json,tsv}
results/provenance/<tool>/provenance-audit.json
results/report/index.html
```

## Obsidian synchronization

The bound canonical design can be synchronized and checked without accessing
the HG002 BAM:

```bash
python workflow/scripts/sync_obsidian_design.py --code-state partial
python workflow/scripts/check_obsidian_sync.py \
  --code-state partial \
  --output results/provenance/obsidian-sync-check.json
```

Equivalent Snakemake targets are `sync_obsidian_design` and
`check_obsidian_sync`.
