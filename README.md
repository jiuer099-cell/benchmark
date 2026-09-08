# PGBench short-read pangenome SV genotyping benchmark

PGBench compares fixed-panel structural-variant genotypers under one frozen
contract: the same paired Illumina FASTQs, the same tool-neutral family-aware
leave-one-out panel, the same GRCh38 reference, the same panel-addressable
canonical-panel truth denominator, the same autosomal benchmark regions, and the same three genotype-aware
evaluators.

The only primary score is:

```text
ME-F1 = (F1_Truvari + F1_Aardvark-GT + F1_vcfdist) / 3
```

This formula is frozen in score-contract version 1.0. All three evaluator F1
values are mandatory, missing evaluators are never renormalized, plugins cannot
provide ME-F1, and strata never change the primary score. The canonicalized
score-contract block and the complete scoring profile have separate SHA-256
identities in every sealed result.

Each evaluator's TP, FP, FN, precision, recall, and F1 is published alongside
the mean. Evaluator range and standard deviation are disagreement diagnostics;
they never change the ranking score. Candidate-vote agreement and whole-truth
recovery are diagnostics only.

## Frozen main-track scope

- Illumina paired-end short reads only.
- Fixed-panel genotyping, not de novo discovery.
- GRCh38 chromosomes chr1–chr22.
- Biallelic, sequence-resolved, simple DEL/INS of 50–10,000 bp.
- HG002 and known first-degree relatives excluded from panel sources.
- HG002 truth is introduced only after the candidate panel is frozen.
- Truth without a reliable panel allele is `UNSCORABLE`, never inferred as 0/0.
- Every adapter emits every canonical candidate with a genotype or a distinct
  failure/no-call status; unsupported candidates cannot disappear silently.

## Deterministic evaluation layers

```text
native VCF -> normalized VCF -> panel allele linking -> all-sites.vcf.gz
           -> deterministic core projection -> evaluation-query.vcf.gz
           -> Truvari / Aardvark-GT / vcfdist -> unified TP/FP/FN/GT layer
           -> ME-F1
```

`all-sites.vcf.gz` and `candidate-status.tsv` retain every candidate and keep
`addressable_called`, `explicit_no_call`, `unsupported_representation`,
`linking_failure`, `missing_output`, and `ambiguous_mapping` distinct.
`evaluation-query.vcf.gz` contains only canonical 0/1 and 1/1 predictions.
The core verifies that each query row maps to one frozen candidate allele,
applies no extra QUAL/GQ filter, and seals the all-sites, status, panel, and
query hashes. A tool cannot select a reserved core output path.

## Tool ecosystem

Bundled, reviewed reference adapters:

- PanGenie — mapping-free k-mer/population-haplotype genotyping.
- vg giraffe + vg call — whole-pangenome graph mapping/genotyping.
- Paragraph — local SV graph realignment.

Community/external adapters:

- GraphTyper2
- Varigraph
- BayesTyper

The Great Genotyper remains a future adapter until a reproducible released
implementation and frozen interface are available. It is not represented by a
placeholder that could be mistaken for tested software.

## Run locally

```powershell
python -m pytest -q
python -m snakemake -n --cores 1 --configfile tests/fixtures/synthetic/config.yaml
python -m snakemake --cores 1 --configfile tests/fixtures/synthetic/config.yaml
```

Production runs use a tool-specific or multi-tool config whose registrations
all point to the same sample FASTQs and canonical panel source. The workflow
first emits the nine shared 10x/20x/30x subsets (three frozen seeds per depth)
and the full-depth entry. Generate the complete run matrix with:

```text
python workflow/scripts/run_coverage_matrix.py \
  --base-config config/production.yaml \
  --coverage-manifest results/<base-run>/coverage/coverage-manifest.json \
  --output-dir run-configs --execute --cores 32
```

## Validity

A score is `valid` only when evaluator semantics, complete candidate status,
evaluation-query invariants, independent panel provenance, allowed information,
frozen tuning, required context/AF strata, and core provenance pass. Synthetic
scores are always provisional. The main result table starts with ME-F1, then
the three evaluator F1 values, followed by Addressability and call/failure
diagnostics. HG002-only conclusions are explicitly
sample-and-configuration-specific. See `docs/IMPLEMENTATION_CHECKLIST.md` for
the original design crosswalk, `docs/CURRENT_IMPROVEMENTS_CROSSWALK.md` for the
current-improvements crosswalk, and `docs/TOOL_PLUGIN_GUIDE.md` for adapters.
