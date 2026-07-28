# HG002 GRCh38 Pangenome SV Benchmark

Snakemake workflow for benchmarking structural-variant detection and genotyping
on HG002 against GRCh38. The workflow runs each tool itself, normalizes its VCF,
and evaluates the same result universe independently with Truvari, Aardvark, and
vcfdist.

## Reference and truth data

- reference: GRCh38 no-alt plus hs38d1 decoy (`GRCh38_no_alt_plus_hs38d1_analysis_set.fasta`), matching the GIAB PacBio CLR GRCh38 BAM;
- sample: HG002 / NA24385;
- primary truth profile: `giab_hg002_grch38_v5_0q`;
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

## Resource preflight

Installing callers such as Kanpig is not enough to run the benchmark. Check the
server-side data paths before starting Snakemake:

```bash
python workflow/scripts/check_resources.py \
  --config config/config.example.yaml \
  --repo-root .
```

The command is read-only. It reports every configured path, whether it is
required for the selected track, its status, and its size. It exits `0` only
when all required files are present, regular files, readable, and non-empty;
an incomplete resource set exits `1`. Use `--json` for machine-readable output.
By default it checks every track in `execution.tracks`; a track can be checked
independently with `--track caller-only` or `--track end-to-end`.

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

Both tracks require the same GRCh38 FASTA, FASTA index (`.fai`), sequence
dictionary (`.dict`), GIAB HG002 truth VCF plus Tabix index, benchmark BED,
HG002-free population SV VCF, and all graph assets enabled by
`pangenome.build_graph_assets`. The track-specific inputs are:

| Track | Required primary sample input | Not required as a primary input |
| --- | --- | --- |
| `caller_only_shared_alignment` | shared GRCh38 BAM and BAI | FASTQ |
| `end_to_end_from_reads` | canonical HG002 FASTQ | shared BAM and BAI |

The technology still has to match the plugin. The vg example expects PacBio
CLR FASTQ. PanGenie expects a single HG002 **Illumina short-read** FASTA/FASTQ
stream and cannot use the PacBio CLR BAM/FASTQ.

When both tracks are configured, both input sets must be available. A BAM index
may use either `sample.bam.bai` or `sample.bai`; the preflight recognizes both.
The preflight does not download, index, alter, or checksum resources and does
not replace formal configuration/provenance validation.

### Data still needed for a production run

Place these files at the paths in `config/config.example.yaml`, or update that
configuration to the paths already used on the server:

- GRCh38 no-alt plus hs38d1 FASTA, its `.fai`, and its `.dict`;
- for caller-only, the matching HG002 PacBio CLR GRCh38 BAM and BAI;
- for end-to-end, the canonical HG002 PacBio CLR FASTQ;
- GIAB HG002 GRCh38 v5.0q SV truth VCF, its `.tbi`, and benchmark BED;
- a sequence-resolved GRCh38 population SV VCF plus `.tbi` that excludes
  HG002/NA24385;
- `graph-assets.lock.yaml`, GBZ, XG, minimizer (`.min`), distance (`.dist`),
  and graph sample-list files for the configured pangenome;
- the required GRCh38 stratification BED assets named in
  `config/stratifications.yaml` before stratified production reporting.

Kanpig and the other callers are software dependencies, not substitutes for
these data. The evaluator environments must also provide `bcftools`, `tabix`,
Truvari, Aardvark, and vcfdist as configured under `evaluation.commands`.

The population VCF header and graph `samples.txt` are checked fail-closed for
both `HG002` and `NA24385`. A public HPRC graph that still contains the HG002
assembly is therefore not an acceptable benchmark graph: prepare a leave-one-out
graph or another GRCh38 graph whose construction cohort excludes HG002.

PanGenie does not need GBZ/XG/min/dist files, but it needs a stricter
PanGenie-ready VCF graph: multi-sample, fully phased, non-overlapping,
sequence-resolved, and with HG002/NA24385 removed. That panel is a separate
biological resource from the ordinary candidate SV VCF.

### Reuse the BAM reads for graph mode

The graph track needs FASTQ, but a second large read download is not mandatory.
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

After that tmux job has finished successfully, create the one FASTQ stream
required by PanGenie. Concatenated gzip members remain a valid gzip stream:

```bash
cd /home/luzhiting/hg002-grch38-pangenome-sv-benchmark

find resources/reads/HG002_Illumina_2x250/raw \
  -maxdepth 1 -type f -name '*.fastq.gz' -print0 \
  | sort -z \
  | xargs -0 cat \
  > resources/reads/HG002.Illumina_2x250.all.fastq.gz.tmp

gzip -t resources/reads/HG002.Illumina_2x250.all.fastq.gz.tmp && \
mv resources/reads/HG002.Illumina_2x250.all.fastq.gz.tmp \
   resources/reads/HG002.Illumina_2x250.all.fastq.gz
```

This is additional data only for the PanGenie track. The PanGenie-ready
leave-one-out panel VCF is still required at the path in
`config/config.pangenie.example.yaml`; do not substitute an HPRC panel whose
allele construction included HG002.

## Fair scoring on a fixed truth universe

There are no evaluator weights and no resource, pangenome, or provenance points.
For every normalized query result, each evaluator casts one equal binary vote.
The report records:

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
runs therefore count the fixed primary-truth VCF records overlapping the
benchmark BED (`T`) and calculate:

```text
softTP = (3*n3 + 2*n2 + n1) / 3
Q = n3 + n2 + n1 + n0
ComparableScore = 2 * min(softTP, T) / (Q + T) * 100
```

`ComparableScore` is the primary end-to-end score. False-positive query calls
increase `Q`; missed truth records increase `T` without increasing `softTP`.
The score is comparable only when sample, reference, truth profile, benchmark
regions, and frozen score-profile SHA-256 are identical. Across PacBio and
Illumina it compares the practical accuracy of the complete pipeline, including
the sequencing evidence. Within one technology track it is the stricter
algorithm comparison.

PanGenie additionally requires interpretation of its panel-limited task:
`ConsensusScore` and panel-stratified genotype metrics describe in-panel
genotyping, while `ComparableScore` retains the complete GIAB truth denominator
and therefore exposes out-of-panel coverage limitations.

The raw category counts, fixed truth count, query count, soft true-positive
count, comparable precision/recall, and both scores are retained. Provenance is
only a validity gate: incomplete non-core provenance makes a result provisional,
while core lineage failure makes it invalid and suppresses numeric scores.
Resource measurements are reported but never contribute points. The benchmark
does not assign ordinal ranks.

The frozen contract is `config/consensus_scoring.yaml`.

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
`pigz` is needed only for the optional BAM to FASTQ command above. Snakemake
does not download biological resources.

Primary outputs:

```text
results/summary/score.tsv
results/summary/point_breakdown.tsv
results/summary/metrics.long.tsv
results/summary/metrics.json
results/summary/<tool>/score-package.json
results/report/index.html
```

`point_breakdown.tsv` is retained as a compatibility filename; its rows now
contain consensus category counts, never weighted points.

### Combine long- and short-read runs in one HTML

Each Snakemake configuration produces a finalized `score.json`. Preserve the
result directories from the PacBio and Illumina runs, then build one
cross-track presentation:

```bash
cd /home/luzhiting/hg002-grch38-pangenome-sv-benchmark

mamba run -n pgbench-bio python workflow/scripts/render_suite_report.py \
  --entry archived_results/kanpig/score.json plugins/kanpig/tool.yaml \
  --entry archived_results/vg/score.json plugins/vg/tool.yaml \
  --entry archived_results/pangenie/score.json plugins/pangenie/tool.yaml \
  --output archived_results/pangenome-sv-suite.html
```

The suite renderer refuses to combine scores with different samples, truth
profiles, or score-profile hashes. It preserves the configured display order
and shows tool paradigm and sequencing technology next to both scores.
