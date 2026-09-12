# PGBench dual-channel pangenome SV genotyping benchmark

PGBench compares fixed-panel structural-variant genotypers under a frozen,
channel-specific evidence contract: SR adapters receive the same HG002
Illumina paired-end reads and HiFi adapters receive the same HG002 PacBio HiFi
reads.  Both channels use the same tool-neutral family-aware leave-one-out
panel, GRCh38 reference, frozen canonical candidate universe, tool-independent
scoring denominator, autosomal benchmark regions, and the same three frozen
evaluators under one unified genotype-aware judgement layer.

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

## Frozen channel scope

- **SR channel:** HG002 Illumina paired-end reads only.
- **HiFi channel:** HG002 PacBio HiFi reads only; CLR and ONT are not
  accepted.
- SR and HiFi are separate comparison tracks and receive separate rankings;
  no cross-channel score, rank, or aggregate is generated.
- Fixed-panel genotyping, not de novo discovery.
- GRCh38 chromosomes chr1–chr22.
- Biallelic, sequence-resolved, simple DEL/INS of 50–10,000 bp.
- HG002 and known first-degree relatives excluded from panel sources.
- HG002 truth is introduced only after the candidate panel is frozen.
- Truth without a reliable panel allele is `UNSCORABLE`, never inferred as 0/0.
- `UNSCORABLE` is decided before any adapter runs, only when a truth event
  cannot be unambiguously linked to the frozen canonical panel.  It is never
  an adapter-specific state.  `unsupported_representation`, `linking_failure`,
  and `missing_output` remain all-sites ledger diagnostics.
- Addressability is diagnostic only and never changes the primary-score
  denominator.
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

## External adapter ecosystem

Core has no tool-specific command, resource, or native-asset logic.  Every
tool, including the five current tools, is registered through the same external
adapter manifest and runner contract.  Adding a tool means adding an adapter
directory and one registry entry; Core source must not change.

- PanGenie — mapping-free k-mer/population-haplotype genotyping.
- vg giraffe + vg call — whole-pangenome graph mapping/genotyping.
- Paragraph — local SV graph realignment.
- GraphTyper2
- Varigraph

BayesTyper and future tools use the same interface.  An adapter that consumes
population or pangenome information must declare the configured Frozen
Haplotype Source Bundle; Core rejects native assets derived from any other
cohort, release, reference, or family-exclusion state.
The frozen biological source is identical; tool-native representations may
differ.  Adapters that do not consume population haplotypes are not given a
population prior merely for symmetry.

| Tool | Adapter | Channel | Evidence | Population context | Current status |
| --- | --- | --- | --- | --- |
| PanGenie | external | SR | Illumina PE | frozen native context | P0 passed; production pending |
| vg Giraffe + call | external | SR | Illumina PE | Frozen Bundle-derived graph | production pending |
| Paragraph | external | SR | shared BAM | none | production pending |
| GraphTyper2 | external | SR | shared BAM | none | production pending |
| Varigraph | external | SR | Illumina PE | Frozen Bundle-derived native assets | production pending |

Adapter status is explicit: `eligible`, `unsupported`, and `adapter_failure`
are distinct.  Unsupported means that an adapter cannot run in the selected
channel; it is never converted to ME-F1 = 0.

## Production isolation and release status

Formal adapters run only in a `bwrap` or Apptainer sandbox with network
disabled, immutable declared inputs, and an adapter-private writable work
directory.  Truth VCFs, GIAB files, and evaluator results are not adapter
inputs and are therefore not mounted into that sandbox.  Parameters are frozen
before the HG002 test run; truth-guided tuning and post-score parameter changes
are forbidden.

PanGenie P0 native-context validation has passed on the server.  The next
release gate is a complete five-adapter production rerun followed by Core
scoring.  Releases are immutable and channel-specific (`PGBench-SR-Illumina`
and `PGBench-LR-HiFi`); changing panel, truth, BED, reference, source bundle,
adapter API, scoring contract, or evaluator version creates a new release.

Every adapter must pass the Adapter Conformance Suite before production.  The
suite covers DEL/INS 0/0, 0/1, 1/1 and no-call cases, unsupported and ambiguous
representations, and canonical biallelic projection.

The Great Genotyper remains a future adapter until a reproducible released
implementation and frozen interface are available. It is not represented by a
placeholder that could be mistaken for tested software.

## Run locally

```powershell
python -m pytest -q
python -m snakemake -n --cores 1 --configfile tests/fixtures/synthetic/config.yaml
python -m snakemake --cores 1 --configfile tests/fixtures/synthetic/config.yaml
```

Production runs use either `config/config.unified-tools.example.yaml` for SR
or `config/config.hifi.example.yaml` for HiFi. Their registrations point to
channel-specific evidence and the same frozen haplotype source. The workflow
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
