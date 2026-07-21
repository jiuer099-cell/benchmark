# HG002 GRCh38 Pangenome SV Benchmark

Snakemake workflow for benchmarking structural-variant detection and genotyping
on HG002 against GRCh38. The workflow runs each tool itself, normalizes its VCF,
and evaluates the same result universe independently with Truvari, Aardvark, and
vcfdist.

## Reference and truth data

- reference: GRCh38 (`GRCh38_no_alt_analysis_set.fasta`);
- sample: HG002 / NA24385;
- primary truth profile: `giab_hg002_grch38_v5_0q`;
- pangenome backbone: GRCh38;
- production assets are fail-closed until paths and SHA-256 values are frozen.

`config/config.example.yaml` contains the production layout. The small files in
`tests/fixtures/synthetic/` are contract fixtures, not biological benchmark data.

## Unweighted consensus evaluation

There are no evaluator weights and no resource, pangenome, or provenance points.
For every normalized query result, each evaluator casts one equal binary vote.
The report records:

- `all_three_correct`: accepted by all three evaluators;
- `exactly_two_correct`: accepted by exactly two;
- `exactly_one_correct`: accepted by exactly one;
- `none_correct`: rejected by all three.

With `N` evaluated results, the comprehensive score is:

```text
ConsensusScore = (3*n3 + 2*n2 + n1) / (3*N) * 100
```

The raw category counts, total, unanimous-correct rate, majority-correct rate,
and score are all retained. Provenance is only a validity gate: incomplete
non-core provenance makes a result provisional, while core lineage failure makes
it invalid and suppresses the numeric score. The benchmark does not rank tools.

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
