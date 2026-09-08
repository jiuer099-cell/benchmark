# Design-document implementation checklist

This checklist maps every substantive requirement in
`PGBench_short_read_pangenome_SV_benchmark_design.md` to executable code or an
explicitly bounded release status. A requirement is not considered completed
merely because a similarly named configuration value exists.

## Benchmark question and candidate universe (sections 1–12)

- **Completed:** one `short_read_fixed_panel_genotyping` track; no long-read or
  caller-only execution modes remain in the active tree.
- **Completed:** tool-neutral `HG002_LOO_HPRC_GRCh38_SV_v1` namespace; canonical
  PGSV IDs and allele ledger are generated independently of tool results.
- **Completed:** panel construction rejects HG002 and aliases plus HG003/HG004
  and their aliases, records all excluded relatives/reasons, and never reads
  HG002 truth when constructing the population panel.
- **Completed:** main-track filters are frozen before scoring: chr1–22,
  biallelic, sequence-resolved simple DEL/INS, 50–10,000 bp.
- **Completed:** multiallelic, overlapping, nested, complex, or otherwise
  excluded source candidates are written to the scope ledger/challenging VCF;
  they are not silently deleted.
- **Completed:** truth preparation and linking create panel-addressability
  evidence. Unmatched/ambiguous truth is `UNSCORABLE`, not an implicit 0/0.
- **Release prerequisite:** the production population VCF, source cohort list,
  extraction method, and SHA-256 values must exist on the server before its
  panel hash can be declared frozen.

## Input and coverage fairness (sections 13–15, 47–48, 61)

- **Completed:** every adapter requires the same `short_fastq_r1` and
  `short_fastq_r2` input names; resolved inputs record hashes, read counts,
  bases, technology, library, and evidence identity.
- **Implementation completed; production run pending:** deterministic paired-read hashing produces byte-identical
  10x/20x/30x subsets for every tool at seeds 1701/1702/1703. The manifest also
  records full depth, counts, bases, paths, and hashes.
- **Implementation completed; production run pending:** `run_coverage_matrix.py` creates exactly nine downsampled run
  configs plus one full-depth config from that shared manifest.
- **Implementation completed; production measurements pending:**
  `summarize_resource_stages.py` refuses combined-only timing, separately
  aggregates one-time build and per-sample wall/CPU/peak-RAM measurements,
  measures index disk bytes, and reports amortized N=1/10/100/1000 costs.
  The release gate remains closed until both scopes have production records;
  a single wrapper runtime cannot be mislabeled as a fair split.

## Adapter and output contract (sections 16–18, 32–45, 66, 68–69)

- **Completed:** plugin core has no tool-name scoring branches. Every plugin
  declares its inputs, runner, environment, modes, capabilities, and output.
- **Implementation completed; real-tool validation pending:** bundled adapters are PanGenie, vg, and Paragraph. External
  adapters are GraphTyper2, Varigraph, and BayesTyper.
- **Implementation completed; per-tool production audits pending:** adapters project native output onto the entire canonical panel.
  The normalization ledger distinguishes `addressable_called`,
  `explicit_no_call`, `missing_output`, `linking_failure`,
  `unsupported_representation`, `ambiguous_mapping`,
  `adapter_conversion_failure`, and `index_build_failure`.
- **Implementation completed; per-tool production audits pending:** addressability manifests retain every canonical candidate,
  tool-native ID, status, and failure reason; the primary denominator is never
  reduced to a tool-specific representable subset.
- **Not applicable yet:** The Great Genotyper has no frozen released adapter
  contract and remains a documented future integration, not fake executable
  code.
- **Planned after the verified adapter track:** a convenience `pgbench score
  --query` wrapper. The underlying canonicalize/evaluate/score modules already
  accept an external query VCF, but the one-command UX is not a release gate
  for the server adapter experiment.

## Evaluation and ME-F1 (sections 19–31, 46, 50–51, 54, 67)

- **Completed:** formal evaluators are exactly Truvari, Aardvark GT, and
  vcfdist. hap.py/RTG are not part of the primary score.
- **Completed:** core deterministically derives `evaluation-query.vcf.gz` from
  `all-sites.vcf.gz` and `candidate-status.tsv`; only 0/1 and 1/1 are emitted.
  Every emitted row is checked against one frozen canonical-panel allele, no
  adapter QUAL/GQ filtering is permitted at this step, and all source/output
  hashes enter provenance. All three evaluators consume that same per-tool
  evaluation-query, prepared truth, BED, reference, and frozen profile.
- **Contract-fixture validation completed; execution with the three installed evaluator binaries pending:** the eight semantic cases cover exact heterozygous calls, both
  heterozygous/homozygous-alt mismatches, ALT↔0/0, representation equivalence,
  explicit no-call, and nearby non-equivalent variants.
- **Completed:** native outputs are converted to benchmark-defined genotype
  TP/FP/FN. Precision, recall, and F1 are recomputed from those counts, so the
  three F1 values have the same biological semantics.
- **Completed:** score-contract v1.0 permits only
  `(F1_Truvari + F1_Aardvark-GT + F1_vcfdist) / 3`. All three inputs are
  required, missing inputs cannot be renormalized, plugins cannot provide the
  score, and strata are explanatory only. The score-contract block and full
  profile have separately verified SHA-256 identities.
- **Completed:** candidate agreement and whole-truth recovery are diagnostics
  only; neither contributes to ME-F1 or determines ranking.
- **Completed:** DEL/INS, five length bins, block bootstrap CI, and
  leave-one-evaluator-out sensitivity are represented in the formal metrics
  contract.
- **Implementation completed; production assets/run pending:** context and
  population-AF strata publish candidate/truth counts, all three evaluator
  TP/FP/FN/P/R/F1 values, ME-F1, call/no-call/failure counts, and
  addressability. The formal scorer fails closed unless exactly the six frozen
  context BEDs and a usable panel AF are present.
- **Implementation completed; production run pending:** paired genomic-block
  bootstrap stores ME-F1/evaluator replicates and 95% CIs; the run aggregator
  computes paired between-tool ME-F1 differences from common resamples.
- **Implementation completed; production run pending:** the repeated-coverage
  summarizer enforces the same three seeds at 10x/20x/30x plus one full-depth
  result and reports mean, sample SD, Student-t 95% CI, minimum, and maximum
  ME-F1. The single full-depth run correctly reports its across-seed CI as
  undefined rather than inventing uncertainty.

## Information, provenance, and release gates (sections 44, 47, 54–56, 63–65)

- **Implementation completed; production manifests pending:** each run freezes `allowed_inputs.json`,
  `parameter_manifest.json`, `external_resource_hashes.json`, and
  `training_or_tuning_status.json`.
- **Completed:** target truth, target assembly, target-family genotypes,
  target-specific calls, and undeclared population resources are forbidden
  during genotyping.
- **Completed:** manifests hash code, config, reference, panel, evaluator
  profile, environments, inputs, outputs, logs, benchmarks, and upstream
  manifests. Incomplete core provenance invalidates the score.
- **Completed:** semantic validation, addressability, and information/tuning are
  mandatory score-validity gates. Synthetic runs cannot become `valid`.
- **Completed:** score-only output and fully reproduced/sealed output are
  distinguishable by status and provenance completeness.

## Generalization and difficult regions (sections 57–60, 62, 70–72)

- **Completed:** the first release is autosomal only and declares HG002 results
  `sample_and_configuration_specific`.
- **Not completed because no additional data were supplied:** independent
  truth samples. No cross-population claim is allowed while
  `independent_samples_status: not_available`.
- **Implementation completed; production asset verification/run pending:** the
  official GIAB GRCh38 v3.1 catalogue freezes source URLs and published MD5s
  for non-repeat, tandem-repeat, segmental-duplication, low-complexity,
  low-mappability, and other-difficult BEDs. The preparation step verifies BED
  syntax and records MD5, SHA-256, size, and interval count; the formal scorer
  records asset hashes and computes per-evaluator/ME-F1 strata. No production
  completion claim is permitted until those exact assets are fetched and
  exercised on the server.

## Required execution order

1. Pass local unit tests and synthetic end-to-end test.
2. Push the exact tested commit to GitHub.
3. Inventory and consolidate old server benchmark material into one managed
   top-level workspace without deleting evidence.
4. Deploy the tested commit into that workspace.
5. Validate real input hashes and build/freeze the production family-aware LOO
   panel.
6. Run evaluator semantic validation in the server environments.
7. Run full-depth and repeated coverage matrices.
8. Publish only sealed scores whose quality gates pass; otherwise publish the
   failure/status evidence instead of a numeric leaderboard value.
