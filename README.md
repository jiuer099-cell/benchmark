# HG002 GRCh38 Pangenome SV Benchmark

Snakemake workflow for benchmarking structural-variant detection and genotyping
on HG002 against GRCh38. The workflow runs each tool itself, normalizes its VCF,
and evaluates the same result universe independently with Truvari, Aardvark, and
vcfdist.

## Reference and truth data

- reference: GRCh38 no-alt plus hs38d1 decoy (`GRCh38_no_alt_plus_hs38d1_analysis_set.fasta`), matching the GIAB PacBio CLR GRCh38 BAM;
- sample: HG002 / NA24385;
- primary truth profile: `giab_hg002_grch38_t2tq100_v0_9_sv` (GIAB
  assembly-based whole-genome draft SV benchmark);
- pangenome backbone: GRCh38;
- production assets are fail-closed until paths and SHA-256 values are frozen.

`config/config.example.yaml` is the BAM-based KanPIG caller-only layout.
`config/config.vg.example.yaml` is the PacBio CLR end-to-end graph layout using
GraphAligner, `vg pack`, and `vg call`.
`config/config.pangenie.example.yaml` is the HG002 Illumina short-read
PanGenie layout. The small files in
`tests/fixtures/synthetic/` are contract fixtures, not biological benchmark
data.

## Benchmark architecture

The workflow has a fixed benchmark shell around a replaceable tool middle:

1. common GRCh38/HG002 preparation, pangenome freeze, leakage checks, stable
   allele IDs, and blinded challenge construction;
2. tool-specific processing beginning at the authorized reads/alignment input;
3. common VCF validation, normalization, allele linking, three-evaluator
   consensus scoring, provenance audit, and report generation.

The middle is deliberately not one universal caller rule. A plugin declares a
family (`graph_pangenome`, `discovery_integration`, `assembly_based`,
`mapping_based`, `genotyping_only`, or
`short_read_pangenome_genotyping`) and its own task stages. The included
examples cover GraphAligner/vg, KanPIG, and PanGenie. User-provided genotypers
use the same checked plugin contract; see
[`docs/TOOL_PLUGIN_GUIDE.md`](docs/TOOL_PLUGIN_GUIDE.md).

The bundled production scope is intentionally pangenome-specific. Linear
reference long-read callers are not planned built-ins; they may only be added
as user plugins when a study explicitly needs a linear-reference baseline.

## Resource preflight

Installing callers such as Kanpig is not enough to run the benchmark. Check the
server-side data paths before starting Snakemake:

```bash
python workflow/scripts/check_resources.py \
  --config config/config.example.yaml \
  --repo-root .
```

The command is read-only. It reports every configured path, whether it is
required for the selected execution mode, its status, and its size. It exits `0` only
when all required files are present, regular files, readable, and non-empty;
an incomplete resource set exits `1`. Use `--json` for machine-readable output.
The legacy `--track` option name in this preflight command selects only the
required input set. It does not split tools into task categories or scoring
groups: every formal run uses the same unified ComparableScore contract.

For the two supplied production configurations:

```bash
python workflow/scripts/check_resources.py \
  --config config/config.example.yaml \
  --repo-root . \
  --track caller-only

python workflow/scripts/check_resources.py \
  --config config/config.vg.example.yaml \
  --repo-root . \
  --track end-to-end

python workflow/scripts/check_resources.py \
  --config config/config.pangenie.example.yaml \
  --repo-root . \
  --track end-to-end
```

Both input modes require the same GRCh38 FASTA, FASTA index (`.fai`), sequence
dictionary (`.dict`), GIAB HG002 truth VCF plus Tabix index, benchmark BED,
HG002-free population SV VCF, and only the graph assets required by the
selected plugin profile. The mode-specific inputs are:

| Input mode | Required primary sample input | Not required as a primary input |
| --- | --- | --- |
| `caller_only_shared_alignment` | shared GRCh38 BAM and BAI | FASTQ |
| `end_to_end_from_reads` | registered HG002 FASTQ input: one long-read file or a complete short-read R1/R2 pair, as required by the plugin | shared BAM and BAI |

The technology still has to match the plugin. The vg example expects PacBio
CLR FASTQ. PanGenie requires an explicit HG002 **Illumina paired-end R1/R2**
FASTQ pair and cannot use the PacBio CLR BAM/FASTQ.

Before a plugin starts, every FASTQ is streamed to EOF. The executor validates
gzip integrity, FASTQ record structure, sequence/quality lengths, and (for
Illumina) equal R1/R2 record counts and normalized read-name pairing. Derived
read and base counts are frozen in `resolved_inputs.json`; conflicting
user-entered counts invalidate formal metric materialization.

When both modes are used in separate runs, the corresponding input sets must be
available. A BAM index
may use either `sample.bam.bai` or `sample.bai`; the preflight recognizes both.
The preflight does not download, index, alter, or checksum resources and does
not replace formal configuration/provenance validation.

### Data still needed for a production run

Place these files at the paths in `config/config.example.yaml`, or update that
configuration to the paths already used on the server:

- GRCh38 no-alt plus hs38d1 FASTA, its `.fai`, and its `.dict`;
- for caller-only, the matching HG002 PacBio CLR GRCh38 BAM and BAI;
- for vg end-to-end, the canonical HG002 PacBio CLR FASTQ;
- for PanGenie end-to-end, the canonical HG002 Illumina paired R1/R2 FASTQ
  files;
- GIAB HG002 GRCh38 T2T-Q100 v0.9 whole-genome draft SV truth VCF, its
  `.tbi`, and paired benchmark BED;
- a sequence-resolved GRCh38 population SV VCF plus `.tbi` that excludes
  HG002/NA24385;
- `graph-assets.lock.yaml`, GBZ, minimizer (`.min`), distance (`.dist`), and
  graph sample-list files for the modern vg profile; XG is required only by
  `vg_legacy_xg`;
- the required GRCh38 stratification BED assets named in
  `config/stratifications.yaml` before stratified production reporting.

Kanpig and the other callers are software dependencies, not substitutes for
these data. The evaluator environments must also provide `bcftools`, `tabix`,
Truvari, Aardvark, and vcfdist as configured under `evaluation.commands`.

The population VCF header and graph `samples.txt` are checked fail-closed for
both `HG002` and `NA24385`. A public HPRC graph that still contains the HG002
assembly is therefore not an acceptable benchmark graph: prepare a leave-one-out
graph or another GRCh38 graph whose construction cohort excludes HG002.

Formal execution also refuses a truth catalog whose three data hashes disagree
with the files. `config/truthsets.yaml` contains the hashes already verified
for the named GIAB T2T-Q100 v0.9 files. Recalculate them on every server after transfer
and confirm that all three match before running:

```bash
cd /home/luzhiting/hg002-grch38-pangenome-sv-benchmark

TRUTH_DIR=resources/truth/giab_hg002_grch38_t2tq100_v0_9_sv
TRUTH_S3=s3://giab/data/AshkenazimTrio/analysis/NIST_HG002_DraftBenchmark_defrabbV0.011-20230725
mkdir -p "$TRUTH_DIR"

for FILE in \
  GRCh38_HG002-T2TQ100-V0.9_stvar.vcf.gz \
  GRCh38_HG002-T2TQ100-V0.9_stvar.vcf.gz.tbi \
  GRCh38_HG002-T2TQ100-V0.9_stvar.benchmark.bed \
  checksum.md5 \
  README.md
do
  aws s3 cp --no-sign-request "$TRUTH_S3/$FILE" "$TRUTH_DIR/$FILE"
done
```

The VCF and BED are a required pair; do not combine either one with a BED or
VCF from another release.

```bash
cd /home/luzhiting/hg002-grch38-pangenome-sv-benchmark

sha256sum \
  resources/truth/giab_hg002_grch38_t2tq100_v0_9_sv/GRCh38_HG002-T2TQ100-V0.9_stvar.vcf.gz \
  resources/truth/giab_hg002_grch38_t2tq100_v0_9_sv/GRCh38_HG002-T2TQ100-V0.9_stvar.vcf.gz.tbi \
  resources/truth/giab_hg002_grch38_t2tq100_v0_9_sv/GRCh38_HG002-T2TQ100-V0.9_stvar.benchmark.bed
```

If any value differs, stop and verify the release and transfer instead of
editing the catalog to accept an unexplained file. The workflow checks these
declarations, then materializes one indexed PASS-or-unfiltered,
fully-contained-in-BED, biallelic DEL/INS, 50-to-10,000-bp truth VCF.
Named failing FILTER records and multiallelic events are excluded under a
frozen policy and counted in the truth audit. Duplicate representations of the
same eligible truth event are deterministically collapsed and separately
counted, so every evaluator and the final truth denominator consume the same
unique generated artifact. This GIAB release is a draft and must be described
as such in publications.

The blinded candidate VCF is restricted to that same frozen DEL/INS, size,
FILTER, biallelic, and fully-contained BED universe before any genotyper sees
it. Its hidden ledger contains exactly one row per emitted candidate, while the
challenge audit records the original panel count and mutually exclusive
exclusion counts. This prevents unsupported small or multiallelic records from
changing tool runtime or being silently omitted from scoring.

PanGenie does not need GBZ/XG/min/dist files, but it needs a stricter
PanGenie-ready VCF graph: multi-sample, fully phased, non-overlapping,
sequence-resolved, and with HG002/NA24385 removed. That panel is a separate
biological resource from the ordinary candidate SV VCF.

### Reuse the BAM reads for graph mode

The vg execution mode needs FASTQ, but a second large read download is not mandatory.
After the BAM passes `samtools quickcheck`, extract one primary sequence per
read and exclude secondary/supplementary alignments:

```bash
cd /home/luzhiting/hg002-grch38-pangenome-sv-benchmark
mkdir -p resources/reads

mamba run -n pgbench-bio samtools view \
  -@ 16 -u -F 0x900 \
  resources/alignments/HG002.GRCh38.bam \
  | mamba run -n pgbench-bio samtools fastq -@ 16 -n - \
  | pigz -p 16 \
  > resources/reads/HG002.PacBio_CLR.fastq.gz.tmp

mv resources/reads/HG002.PacBio_CLR.fastq.gz.tmp \
  resources/reads/HG002.PacBio_CLR.fastq.gz
```

Do not run this against a partial BAM or a BAM that still has an `.aria2`
sidecar.

### Additional HG002 reads for PanGenie

PanGenie is a short-read method. The PacBio CLR BAM cannot be converted into
valid Illumina input. Download the GIAB HG002 Illumina 2x250 FASTQ collection
into a separate directory:

```bash
cd /home/luzhiting/hg002-grch38-pangenome-sv-benchmark
mkdir -p resources/reads/HG002_Illumina_2x250/raw

tmux new-session -d -s hg002_illumina_download \
  "cd /home/luzhiting/hg002-grch38-pangenome-sv-benchmark && \
   aws s3 cp --no-sign-request --recursive \
     s3://giab/data/AshkenazimTrio/HG002_NA24385_son/NIST_Illumina_2x250bps/reads/ \
     resources/reads/HG002_Illumina_2x250/raw/ \
     --exclude '*' --include '*.fastq.gz'"

tmux capture-pane -pt hg002_illumina_download | tail -n 30
```

After that tmux job has finished successfully, preserve the mates separately.
The following example combines split lanes within each mate while keeping R1
and R2 distinct (adjust the filename patterns if GIAB used a different lane
naming convention):

```bash
cd /home/luzhiting/hg002-grch38-pangenome-sv-benchmark

find resources/reads/HG002_Illumina_2x250/raw \
  -type f -iname '*R1*.fastq.gz' -print0 \
  | sort -z \
  | xargs -0 cat \
  > resources/reads/HG002.Illumina_2x250.R1.fastq.gz.tmp

find resources/reads/HG002_Illumina_2x250/raw \
  -type f -iname '*R2*.fastq.gz' -print0 \
  | sort -z \
  | xargs -0 cat \
  > resources/reads/HG002.Illumina_2x250.R2.fastq.gz.tmp

gzip -t resources/reads/HG002.Illumina_2x250.R1.fastq.gz.tmp && \
gzip -t resources/reads/HG002.Illumina_2x250.R2.fastq.gz.tmp && \
mv resources/reads/HG002.Illumina_2x250.R1.fastq.gz.tmp \
   resources/reads/HG002.Illumina_2x250.R1.fastq.gz && \
mv resources/reads/HG002.Illumina_2x250.R2.fastq.gz.tmp \
   resources/reads/HG002.Illumina_2x250.R2.fastq.gz
```

This is additional data only for the PanGenie execution mode. The PanGenie-ready
leave-one-out panel VCF is still required at the path in
`config/config.pangenie.example.yaml`; do not substitute an HPRC panel whose
allele construction included HG002.

## Fair scoring on a fixed truth universe

There are no evaluator weights and no resource, pangenome, or provenance points.
The frozen `config/evaluator_profile.yaml` defines the submitted-query filter,
SVTYPE, frozen 50-to-10,000-bp size envelope, biallelic shape, fully-contained benchmark-region,
breakpoint, resolved-sequence, genotype, and no-call semantics. The formal
universe is biallelic DEL/INS from 50 through 10,000 bp: truth and query calls
must be `PASS` or explicitly unfiltered (`.`). Named failing FILTER records remain auditable
but cannot enter the score denominator as submitted detections. Formal matching
is deterministic and one-to-one: a truth event can be assigned to at most one
query event, so duplicated calls cannot earn duplicate truth credit.

Detection, genotype correctness, and no-call are separate semantics. The
unified ComparableScore uses only the detection vote. Candidate-site GT
diagnostics are computed later from the benchmark-private hidden ledger:
returned `0/0`, non-reference GT, and `./.` are compared site by site; a
missing site is interpreted only through the plugin's frozen `hom_ref` or
`no_call` declaration. Only candidate alleles that are fully contained in the
benchmark BED and satisfy the same biallelic DEL/INS 50-to-10,000-bp universe enter the
candidate-GT denominator; all other panel records remain auditable as
`truth_scorable=0` and are never assumed to be `0/0`. Under the frozen primary
profile, a no-call is incorrect in the primary candidate-GT accuracy, while a
separate called-only concordance is also reported. To prevent abundant `0/0`
sites from inflating that diagnostic, the report also includes the full GT
confusion matrix, non-reference precision/sensitivity, specificity, balanced
accuracy, macro-F1, and no-call counts by truth GT class. Candidate-derived outputs
map directly through their blinded candidate ID; discovery outputs may map
through a unique `PANGENOME_LINKED_ID` only when the linker status is
`in_panel_exact` or `in_panel_equivalent`. Conflicts and rejected/ambiguous
links are counted explicitly. These hidden labels are never mounted into the
plugin sandbox. GT concordance and no-call counts cannot silently change
detection credit. For every eligible normalized detection event, each
evaluator casts one equal binary vote. The report records:

- `all_three_correct`: accepted by all three evaluators;
- `exactly_two_correct`: accepted by exactly two;
- `exactly_one_correct`: accepted by exactly one;
- `none_correct`: rejected by all three.

The query-only diagnostic score remains:

```text
ConsensusScore = (3*n3 + 2*n2 + n1) / (3*N) * 100
```

`ConsensusScore` alone is not used for cross-tool comparison: a tool could emit
only a few high-confidence calls and avoid false-negative penalties. Production
runs therefore count the fixed primary-truth VCF records fully contained in the
benchmark BED (`T`) and calculate:

```text
softTP = (3*n3 + 2*n2 + n1) / 3
Q = n3 + n2 + n1 + n0
ComparableScore = 2 * softTP / (Q + T) * 100
```

The one-to-one ledger contract requires `softTP <= T`; a violation invalidates
the run instead of being silently clipped. This value is reported as the
secondary **Global End-to-End SV Recovery Score**. False-positive query calls
increase `Q`; missed truth records increase `T` without increasing `softTP`.
The score is comparable only when sample, exact reference/truth/BED content
hashes, the stable pangenome-manifest contract hash, the hidden challenge-ledger
hash, truth profile, frozen score profile, frozen evaluator-profile SHA-256,
and evaluator version fingerprints are identical. A graph-lock hash is also
frozen for every graph-consuming plugin; it must remain identical when the same
tool and mode are compared across runs, while different tools may legitimately
declare different graph dependencies. Across PacBio and Illumina
it compares the practical accuracy of the complete pipeline, including the
sequencing evidence; the report does not create task-specific lanes or claim a
platform-independent pure-algorithm rank.

The primary **Pangenome Genotyping Score** is `100 * genotype macro-F1` on the
same frozen, blinded population-panel candidate universe for every applicable
plugin. The report separately shows non-reference F1, no-call count, and Panel
Coverage (`truth-positive panel events / all eligible GIAB truth events`).
Coverage is not multiplied into genotype quality: this prevents a small panel
from masquerading as whole-genome discovery while avoiding the former
single-digit ceiling for a correct panel genotyper.

vcfdist receives a benchmark-owned detection-only copy in which unphased `0/1`
is deterministically represented as `0|1`; the original GT is retained for all
genotype scoring. Native execution is content-addressed separately from event
mapping, so parser/report changes reuse the frozen native artifacts. Event
mapping coverage is reported and a consensus with more than 1% unresolved
events is marked invalid rather than converting unresolved rows into errors.

The raw category counts, fixed truth count, query count, soft true-positive
count, comparable precision/recall, and both scores are retained. A
fixed 10 Mb genomic-block bootstrap (1000 deterministic replicates) jointly
resamples truth records, matched queries, and unmatched false positives to
provide the 95% ComparableScore confidence interval. The sampling plan depends
only on the frozen reference/truth/BED/profile assets, so tools evaluated under
the same contract use paired blocks. SVTYPE/length-bin tables expose important
strata. Every formal score also records the actual input technology,
library/source IDs, BAM or FASTQ hashes and sizes, coverage/downsampling
metadata, and validated read/base counts. The unified HTML uses this actual
evidence rather than the plugin's supported-technology list and rejects
same-technology/same-mode comparisons made from different evidence files.
Evaluator version output, parameters, thresholds, and hashes
are sealed into each run. Provenance is
only a validity gate: incomplete non-core provenance makes a result provisional,
while core lineage failure makes it invalid and suppresses numeric scores.
Formal tool execution is measured three times with a fresh, isolated
attempt-local HOME/XDG cache (`isolated_empty_tool_cache`); this does not claim
to flush the operating-system page cache. Median and IQR wall-time/RSS
statistics are reported but never contribute points. The current timed rule is
the auditable end-to-end tool execution, including benchmark hashing overhead;
pure algorithm timing is a separate paper experiment. Every formal plugin runs
without network access inside `bwrap` or Apptainer. The benchmark does not
assign ordinal ranks.

The frozen contracts are `config/consensus_scoring.yaml` and
`config/evaluator_profile.yaml`.

The distinction between implemented benchmark mechanics and the production
experiments required for a full methods paper is maintained in
[`docs/PAPER_READINESS.md`](docs/PAPER_READINESS.md).

## Run

```bash
conda env create -f environment.yaml
conda activate pgbench-sv
python -m pytest -q

snakemake --snakefile Snakefile \
  --configfile tests/fixtures/synthetic/config.yaml \
  --cores 2 --rerun-incomplete
```

Server production commands:

```bash
cd /home/luzhiting/hg002-grch38-pangenome-sv-benchmark

# BAM-based pangenome-panel genotyping with KanPIG.
snakemake --snakefile Snakefile \
  --configfile config/config.example.yaml \
  --cores 16 --use-conda --rerun-incomplete --keep-going

# True read-to-graph mapping/calling for PacBio CLR.
snakemake --snakefile Snakefile \
  --configfile config/config.vg.example.yaml \
  --cores 16 --use-conda --rerun-incomplete --keep-going

# Short-read pangenome genotyping with PanGenie.
snakemake --snakefile Snakefile \
  --configfile config/config.pangenie.example.yaml \
  --cores 16 --use-conda --rerun-incomplete --keep-going
```

The graph adapter environment installs `graphaligner` and `vg>=1.63`; the
KanPIG adapter pins KanPIG 2.0.2, and the PanGenie adapter pins PanGenie 4.2.1.
KanPIG receives the workflow's declared thread allocation through its native
`--threads` option.
`pigz` is needed only for the optional BAM to FASTQ command above. Snakemake
does not download biological resources.

Primary outputs:

```text
results/<run_id>/summary/score.tsv
results/<run_id>/summary/point_breakdown.tsv
results/<run_id>/summary/metrics.long.tsv
results/<run_id>/summary/metrics.json
results/<run_id>/summary/<tool>/score.json
results/<run_id>/summary/<tool>/score-package.json
results/<run_id>/report/index.html
logs/<run_id>/...
benchmarks/<run_id>/...
```

`point_breakdown.tsv` is retained as a compatibility filename; its rows now
contain consensus category counts, never weighted points.

### Combine long- and short-read runs in one HTML

Each Snakemake configuration produces a finalized `score.json` under its
run-id-isolated result directory. After the three example runs finish, build
one unified presentation directly from those current paths:

```bash
cd /home/luzhiting/hg002-grch38-pangenome-sv-benchmark

mkdir -p results/combined_reports

mamba run -n pgbench-bio python workflow/scripts/render_suite_report.py \
  --entry results/hg002_clr_graph_panel_kanpig/summary/kanpig/score.json plugins/kanpig/tool.yaml \
  --entry results/hg002_clr_graphaligner_vg/summary/vg/score.json plugins/vg/tool.yaml \
  --entry results/hg002_illumina_pangenie/summary/pangenie/score.json plugins/pangenie/tool.yaml \
  --output results/combined_reports/pangenome-sv-suite.html
```

The suite renderer refuses to combine scores with different samples, truth
profiles, score-profile hashes, frozen evaluator-profile hashes, shared
benchmark-resource hashes, pangenome-manifest contracts, or hidden challenge
ledgers. It also rejects a changed graph lock for repeated runs of the same
tool/mode, without forcing unrelated plugins to share one graph format. It
preserves the configured display order and shows tool paradigm and sequencing
technology next to both scores. It is one comparison table, not a collection
of task-specific rankings.
