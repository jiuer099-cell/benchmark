# Small-region real-tool integration gate

`tests/integration/test_real_tools_region.py` provides two independent,
opt-in checks for the production server.

The preflight gate verifies:

- the KanPIG 2.0.2, PanGenie 4.2.1, PanGenie-index, and vg 1.63+
  command-line programs;
- the `bwrap` or `apptainer` backend declared by every bundled production
  plugin;
- a complete HG002/GRCh38 small-region resource directory shared by the three
  plugin execution modes.

Three independent execution gates run the KanPIG, PanGenie, or vg adapter through
`workflow/scripts/pgbench_exec.py::execute_tool`, with
`execution_purpose=formal`, `sandbox_backend=bwrap`,
the plugin's declared mode, and a hard timeout. Each validates the attempt
record, sandbox status, resolved inputs, logs, and non-empty VCF. The vg gate
also runs vg Giraffe; none of these bounded adapter gates runs an evaluator,
scoring, or the full Snakemake workflow.

Ordinary `pytest` runs skip every test in this module. Setting only the
preflight switch never invokes KanPIG. For `bwrap`, preflight does launch
`/bin/true` in a network-isolated namespace so a server that forbids the
namespace operations required by formal execution fails before any biological
job starts.

## Prepare the bounded resource directory

Create a separate directory outside the repository and copy
`tests/fixtures/real_region/region.example.yaml` to `region.yaml` inside it.
Materialize every path named by the manifest from one interval no longer than
2,000,000 bp. The complete set of unique files must be no larger than 2 GiB.

The directory is deliberately comprehensive:

```text
region.yaml
reference/       regional GRCh38 FASTA, FAI, and dictionary
truth/           regional GIAB VCF/TBI and benchmark BED
alignments/      regional HG002 PacBio BAM/BAI
reads/           regional long-read and paired Illumina R1/R2 FASTQ streams
pangenome/       regional manifest, population, challenge, and PanGenie VCFs
graph/           matching lock, GBZ, XG, MIN, DIST, and sample list
```

All paths in `region.yaml` must be relative and must resolve inside this
directory. Empty files, missing files, files outside the directory, `.aria2`
sidecars, intervals over 2 Mb, and resource sets over 2 GiB fail the test.
Never point this manifest at the full 93 GiB BAM or a whole-genome graph.

Preflight only checks the directory contract. Prepare regional resources with
the appropriate bioinformatics tools and validate their biological coordinate,
sample-exclusion, compression, and index consistency before opting into the
real caller smoke run. `short_fastq_r1` and `short_fastq_r2` must be the two
mates of the same read set; an interleaved or concatenated file is not accepted.

## Select the installed binaries

The test searches `PATH` by default. Because the plugins may be installed in
different Conda/Mamba environments, absolute paths are recommended:

```bash
export PGBENCH_KANPIG_BIN=/path/to/conda/envs/kanpig/bin/kanpig
export PGBENCH_PANGENIE_BIN=/path/to/PanGenie
export PGBENCH_PANGENIE_INDEX_BIN=/path/to/PanGenie-index
export PGBENCH_VG_BIN=/path/to/vg
export PGBENCH_BWRAP_BIN=/usr/bin/bwrap
```

For an Apptainer-backed plugin, use
`PGBENCH_APPTAINER_BIN=/path/to/apptainer`. The plugin manifest must also
declare its immutable container.

For the execution gate, `PGBENCH_KANPIG_BIN` must resolve to
`<conda-prefix>/bin/kanpig`. The test adds that `bin` directory to `PATH` and
sets `CONDA_PREFIX` so `execute_tool` maps the complete runtime prefix read-only
into the bwrap namespace. The regional BAM index and FASTA index are declared
inputs and are also mapped read-only; access does not depend on undeclared
adjacent files.

## Run preflight only

```bash
cd /home/luzhiting/hg002-grch38-pangenome-sv-benchmark

export PGBENCH_REAL_TOOL_TESTS=1
export PGBENCH_REAL_REGION_DIR=/absolute/path/to/hg002-real-region

python -m pytest -q tests/integration/test_real_tools_region.py
```

This switch never executes a production caller.

## Run the formal KanPIG small-region test

Use a separate shell invocation so the execution intent is unambiguous:

```bash
cd /home/luzhiting/hg002-grch38-pangenome-sv-benchmark

unset PGBENCH_REAL_TOOL_TESTS
export PGBENCH_REAL_TOOL_EXECUTION=1
export PGBENCH_REAL_REGION_DIR=/absolute/path/to/hg002-real-region
export PGBENCH_KANPIG_BIN=/path/to/conda/envs/kanpig/bin/kanpig
export PGBENCH_BWRAP_BIN=/usr/bin/bwrap

python -m pytest -q \
  tests/integration/test_real_tools_region.py::test_kanpig_formal_small_region_execution
```

Only the exact value `PGBENCH_REAL_TOOL_EXECUTION=1` authorizes the KanPIG
process. The executor terminates it after 600 seconds. This switch is
independent of `PGBENCH_REAL_TOOL_TESTS`.

## Run the formal PanGenie small-region test

PanGenie and PanGenie-index must be installed in the same environment:

```bash
cd /home/luzhiting/hg002-grch38-pangenome-sv-benchmark

export PGBENCH_REAL_REGION_DIR=/absolute/path/to/hg002-real-region
export PGBENCH_PANGENIE_BIN=/path/to/pangenie_env/bin/PanGenie
export PGBENCH_PANGENIE_INDEX_BIN=/path/to/pangenie_env/bin/PanGenie-index
export PGBENCH_BWRAP_BIN=/usr/bin/bwrap
export PGBENCH_REAL_PANGENIE_EXECUTION=1

python -m pytest -q \
  tests/integration/test_real_tools_region.py::test_pangenie_formal_small_region_execution
```

## Run the formal vg small-region test

For the bounded gate, vg Giraffe and its graph indexes must be available:

```bash
cd /home/luzhiting/hg002-grch38-pangenome-sv-benchmark

export PGBENCH_REAL_REGION_DIR=/absolute/path/to/hg002-real-region
export PGBENCH_VG_BIN=/path/to/vg_env/bin/vg
export PGBENCH_BWRAP_BIN=/usr/bin/bwrap
export PGBENCH_REAL_VG_EXECUTION=1

python -m pytest -q \
  tests/integration/test_real_tools_region.py::test_vg_formal_small_region_execution
```

With neither exact opt-in value, the module remains skipped:

```bash
unset PGBENCH_REAL_TOOL_TESTS
unset PGBENCH_REAL_TOOL_EXECUTION
unset PGBENCH_REAL_PANGENIE_EXECUTION
unset PGBENCH_REAL_VG_EXECUTION
python -m pytest -q tests/integration/test_real_tools_region.py
```

Passing preflight means the executables, declared sandbox, and bounded resource
layout are available. Passing each execution gate additionally means that
specific formal regional adapter run completed and produced a contract-valid
VCF. It is not evidence that the whole-genome benchmark or its biological
conclusions are correct.
