# PGBench short-read pangenome SV genotyping benchmark

PGBench compares fixed-panel structural-variant genotypers under one frozen
contract: the same paired Illumina FASTQs, the same tool-neutral family-aware
leave-one-out panel, the same GRCh38 reference, the same panel-addressable
truth, the same autosomal benchmark regions, and the same three genotype-aware
evaluators.

The only primary score is:

```text
ME-F1 = (Truvari F1 + Aardvark-GT F1 + vcfdist F1) / 3
BenchmarkScore = 100 * ME-F1
```

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

A score is `valid` only when evaluator semantics, candidate addressability,
allowed information, frozen tuning, and core provenance pass. Synthetic scores
are always provisional. HG002-only conclusions are explicitly
sample-and-configuration-specific. See `docs/IMPLEMENTATION_CHECKLIST.md` for
the design-document crosswalk and `docs/TOOL_PLUGIN_GUIDE.md` for adapters.
