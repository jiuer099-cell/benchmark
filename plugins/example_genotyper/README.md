# Example external genotyper

This plugin is executable code, not a submitted VCF. PGBench launches
`run.py` during the current Snakemake job and supplies only the `PGBENCH_*`
variables allowed by the selected execution mode.

The mock implementation reads the blinded candidate VCF and emits exactly one
deterministic HG002 genotype for every candidate. It never reads truth data and
is intended only for Phase 1 contract and DAG tests. Its
`sandbox_backend: none` setting is a development-fixture exception; formal
execution rejects it. A real external plugin must use `bwrap` or Apptainer.
