# Paper-readiness plan

This document separates implemented benchmark mechanics from the experiments
needed to support a full methods paper. HG002/GRCh38 remains the primary
benchmark and is not replaced by another dataset.

## Implemented methodological core

- one Snakemake DAG with run-isolated results, logs, benchmarks, and provenance;
- one extensible plugin contract for alignment-based, graph-based, and
  panel-genotyping tools, without task-specific rankings;
- technology compatibility plus full streaming FASTQ validation before
  execution (gzip EOF/CRC, four-line structure, sequence/quality lengths,
  paired record counts, and normalized R1/R2 read-name correspondence);
- plugin-specific graph-resource profiles and dependency resolution;
- formal network-isolated execution through `bwrap` or an immutable Apptainer
  container;
- a frozen biallelic DEL/INS 50-to-10,000-bp universe from the GIAB HG002
  GRCh38 T2T-Q100 v0.9 whole-genome draft SV benchmark, with PASS/unfiltered
  truth and query calls, deterministic duplicate-truth removal, and full
  containment in the paired benchmark BED;
- deterministic global one-to-one tolerant matching, so one truth event cannot
  earn duplicate credit;
- separate detection and hidden candidate-site genotype/no-call semantics,
  including a BED/type/size-scoped truth denominator, no-call-as-incorrect
  primary accuracy, called-only concordance, declared handling of omitted
  candidate records, GT confusion/balanced metrics, and safe graph-allele
  linkage diagnostics;
- equal Truvari, Aardvark, and vcfdist detection votes, with no evaluator
  weights and no ordinal ranking;
- a primary blinded-panel genotype macro-F1 score plus a secondary fixed-truth
  global recovery score, paired fixed-genomic-block confidence interval
  that jointly resamples truth/query/FP contributions, type/length stratified
  output, repeated resource measurements, and an HTML report;
- opt-in bounded real-resource preflight and independent formal KanPIG,
  PanGenie, and vg adapter execution gates.
- a frozen per-run evidence profile recording the actual technology, library
  and source evidence IDs, BAM/FASTQ hashes and sizes, coverage/downsampling
  metadata, and FASTQ-derived read/base counts. Suite comparison refuses to
  compare two tools in the same technology/mode cell when those evidence
  identities differ.

These features establish an auditable benchmark implementation. They do not,
by themselves, establish biological validity or population generality.

## Non-negotiable checks before a paper result is frozen

1. **Reference identity**

   Compare FASTA contig names, lengths, and sequence MD5 values with BAM
   `@SQ/SN/LN/M5`, truth VCF contigs, panel contigs, and the graph reference
   path. A shared `GRCh38` label is not enough.

2. **Leave-one-out leakage audit**

   Prove that HG002, NA24385, and all configured aliases are absent from the
   graph and PanGenie panel. For family-independent claims, also exclude the
   relevant HG003/HG004 parental haplotypes or quantify the effect in a
   dedicated ablation.

3. **Immutable software environments**

   Resolve Conda environments to exact package/build lock files. Record
   container image digests rather than mutable tags. Re-run the version and
   parameter fingerprint checks on the production server.

4. **Real end-to-end executions**

   Run the bounded real-region gate first, then whole-genome KanPIG, PanGenie,
   and vg executions. Archive `.complete.json`, rule manifests, logs, resource
   records, normalized VCFs, evaluator ledgers, scores, and HTML reports.
   Execute all three bounded adapter gates on the production Linux server
   before claiming that the bundled adapters were experimentally validated.

5. **External-plugin validation**

   Integrate at least two independently maintained SV/genotyping tools through
   the public plugin contract without editing core workflow code. The tools
   should exercise different input contracts. This is stronger evidence of
   extensibility than another in-repository example adapter.

6. **Stratification assets**

   Materialize and checksum the GRCh38 difficulty, repeat, low-mappability,
   sequence-context, medically relevant, population-frequency, in/out-panel,
   and graph-complexity annotations. Predeclare how boundary-crossing events
   and overlapping strata are handled.

7. **Statistical analysis**

   Retain the implemented deterministic 10 Mb genomic-block interval as the
   per-tool uncertainty estimate. It jointly resamples truth, matched query,
   and unmatched-FP contributions, with a tool-independent block plan. For the
   paper comparison, additionally report the paired distribution and interval
   of the score difference between tools run on the same sample and evidence,
   and repeat the analysis over preregistered alternative block sizes.

8. **Resource-timing separation**

   The built-in repeated measurement is the auditable end-to-end tool rule and
   includes input/output hashing. It creates a fresh attempt-local HOME/XDG
   cache but does not flush the operating-system page cache. For paper resource
   claims, additionally instrument the inner tool command, predeclare warm/cold
   page-cache policy, and report both algorithm-only and audited end-to-end
   cost.

9. **Metric-sensitivity analysis**

   Repeat scoring over a preregistered grid of breakpoint, size-similarity,
   sequence-similarity, and evaluator-credit thresholds. The main table must
   use the frozen primary profile; sensitivity results belong in supplementary
   material.

10. **Truth-uncertainty and discordance review**

    Predeclare a blinded review protocol for a fixed sample of all-three,
    two-of-three, one-of-three, and zero-of-three evaluator disagreements.
    Where feasible, inspect read support and local graph/haplotype
    reconstruction without changing the frozen primary score. Report how much
    observed disagreement is caused by representation, truth incompleteness,
    or genuine tool error.

11. **Reproducibility release**

    Archive the exact source commit, immutable environment/container locks,
    configuration snapshots, resource manifests and checksums, synthetic and
    bounded real-region fixtures, machine-readable result tables, and the
    report-generating commands in a versioned DOI-backed release. Large
    controlled or public inputs may be referenced by verified source URI and
    checksum rather than redistributed.

## Minimum experiment matrix

The primary whole-genome matrix should include:

| Factor | Minimum levels |
| --- | --- |
| Tool | KanPIG, PanGenie, vg, plus two external plugins |
| Evidence | original HG002 PacBio BAM/derived FASTQ and paired HG002 Illumina reads |
| Coverage | full depth plus at least 5x, 10x, 20x, and 30x downsampling where technically valid |
| Replication | at least three fixed downsampling seeds per non-full-depth level |
| Graph/panel | production leave-one-out resource plus panel-size or ancestry-composition ablations |
| Evaluator | frozen primary profile plus the preregistered sensitivity grid |
| Context | overall, SVTYPE, size, difficult/repeat/mappability, AF, in/out-panel, and graph complexity |
| Resources | three measured attempts per full-depth production configuration |

Predeclare primary and secondary research questions before the whole-genome
run. A defensible minimum is: primary ComparableScore comparison under the
frozen contract; secondary candidate-GT/no-call behavior, evaluator
disagreement, difficult-context robustness, graph/panel ablation, coverage
response, and audited versus inner-command resource cost. Confidence intervals
and effect sizes should be shown even when no ordinal ranking is produced.

The original HG002 BAM remains the source of the long-read experiment. Derived
FASTQ and downsampled BAM/FASTQ files must retain a provenance link to it.

For a broad claim that the benchmark produces stable conclusions across
individuals or ancestries, add independently held-out, truth-characterized
samples. If no additional sample is permitted, title and conclusions must
present the study as an HG002 methods case study; the benchmark software may be
general, but population generality has not been experimentally demonstrated.

## Primary score and interpretation

For each submitted in-scope query event, the three frozen evaluators cast equal
binary detection votes. Let `n3`, `n2`, `n1`, and `n0` be the numbers accepted
by three, two, one, or zero evaluators:

```text
softTP = (3*n3 + 2*n2 + n1) / 3
Q = n3 + n2 + n1 + n0
T = number of events in the frozen primary truth universe
ComparableScore = 2 * softTP / (Q + T) * 100
```

Extra calls increase `Q`; missed truth events remain in `T`; duplicated calls
cannot reuse a truth event. The ledger must satisfy `softTP <= T`; violating
that invariant makes the run invalid rather than silently clipping its score.
Genotype correctness and no-call behavior are
reported separately and never change detection credit. The score compares the
observed performance of complete pipelines on one frozen benchmark contract.
It is not a technology-free measure of intrinsic algorithm quality.

## Suggested paper structure

1. benchmark contract and threat model;
2. plugin architecture, sandboxing, and provenance;
3. frozen truth universe and global event matching;
4. equal-evaluator consensus and statistical analysis;
5. HG002 whole-genome results and difficult-context stratification;
6. coverage, graph/panel, and evaluator-sensitivity ablations;
7. external-plugin extensibility experiments;
8. resource/Pareto analysis, limitations, and reproducibility package.

Relevant primary starting points include the
[Genome in a Bottle program](https://www.nist.gov/programs-projects/genome-bottle),
[PanGenie study](https://www.nature.com/articles/s41588-022-01043-w),
[Truvari study](https://pubmed.ncbi.nlm.nih.gov/36575487/), and
[vcfdist study](https://pmc.ncbi.nlm.nih.gov/articles/PMC11446017/).
