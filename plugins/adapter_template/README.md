# PGBench adapter kit (template)

This directory is a **template, not a runnable tool**. Copy it, rename it,
fill in your tool, and pass the gates described below. The benchmark core
never learns about your tool by name; it only sees this adapter's manifest.

## Directory layout

```
adapter_template/
├── README.md              <- this guide
├── tool.template.yaml     <- copy to <your_tool>/tool.yaml and fill
│                             (schema: workflow/schemas/tool.schema.yaml)
├── run.py                 <- generic scaffold; implement run_tool() only
├── rule-registry.yaml     <- your billable chain declaration
├── envs/environment.yaml  <- pinned conda environment
└── variants/              <- per read-class / technology guidance
    ├── short-shared-bam/  <- short reads, benchmark's frozen BAM (Paragraph-style)
    ├── short-fastq/       <- short reads, you map yourself (vg/BayesTyper-style)
    ├── long-pacbio-hifi/  <- lr_hifi track
    ├── long-pacbio-clr/   <- lr_clr track
    └── long-ont/          <- lr_ont track (frozen dataset binding, strict policy)
```

The manifest is intentionally named `tool.template.yaml` so CI adapter
discovery (`plugins/*/tool.yaml`) never mistakes the scaffold for a formal
tool; rename it to `tool.yaml` in your copy.

## Track × input matrix

| track (config/track_registry.yaml) | read_class | technology | reads arrive as | reference adapters in this repo |
| --- | --- | --- | --- | --- |
| `sr_illumina` | short | `illumina_pe` | frozen shared BAM **or** raw R1/R2 FASTQ | paragraph, graphtyper2 (BAM); vg, bayestyper, varigraph (FASTQ) |
| `lr_hifi` | long | `pacbio_hifi` | frozen long-read FASTQ | example_hifi_adapter (fixture only) |
| `lr_clr` | long | `pacbio_clr` | frozen long-read FASTQ | — |
| `lr_ont` | long | `ont` | frozen `HG002_ONT_R9.4.1_Guppy5.0.6_SUP_pass_v1` FASTQ | — |

One adapter benchmarks one track. Cross-track or cross-technology ranking is
disabled by policy (`track_registry.yaml: policy`).

## The five steps to a formal run

1. **Copy & rename** this directory to `plugins/<your_tool>/`; rename
   `tool.template.yaml` to `tool.yaml` and replace the id `adapter_template`
   everywhere.
2. **Fill `tool.yaml`** using the variant guide that matches your input mode.
   The manifest is schema-validated; keep `information_contract` untouched
   (they are the honesty contract, not boilerplate).
3. **Implement `run_tool()`** in `run.py`. You get input *injection* checks
   and the all-sites projection for free; only the "invoke my tool" part is
   yours. FASTQ/BAM integrity, SHA256, read counts and dataset identity are
   Core-managed immutable validation, performed once before every adapter.
   Do not reimplement or bypass them in an adapter.
   If your tool builds genome-wide graphs in memory, chunk your input like
   `plugins/paragraph/run.py` does — a single whole-genome invocation of that
   kind OOM-killed at 389 GB on the reference server.
4. **Pass the gates, in order**: config `validate` → synthetic/conformance
   run → dry run (only your chain + the shared pre-score chain may execute)
   → production with a fresh `run_id`.
5. **Honest declarations**: `tool_memory_mb` must reflect a measured pilot
   peak, `billable_stages` must include every stage you run, and
   `parameter_contract` must name the exact tool release and parameters.

## Core → adapter environment contract

The core injects exactly the inputs you declare in
`supported_modes...required_inputs` (mapping in
`workflow/scripts/pgbench_exec.py: INPUT_ENVIRONMENT`):

| declared input | environment variable |
| --- | --- |
| `short_fastq_r1` / `short_fastq_r2` | `PGBENCH_INPUT_FASTQ_R1` / `R2` |
| `long_reads_fastq` | `PGBENCH_INPUT_LONG_READS_FASTQ` |
| `shared_shortread_alignment` / `_index` | `PGBENCH_SHARED_ALIGNMENT_BAM` / `BAI` |
| `reference` / `reference_index` | `PGBENCH_REFERENCE_FASTA` / `PGBENCH_REFERENCE_INDEX` |
| `candidate_panel` | `PGBENCH_CANDIDATE_VCF` |
| `pangenome_manifest` / `pangenome_panel` | `PGBENCH_PANGENOME_MANIFEST` / `PGBENCH_PANEL_VCF` |
| `graph_assets` | `PGBENCH_GRAPH_DIR` |
| `tool_index` | `PGBENCH_INDEX_DIR` |
| `adapter_asset.<name>` | declared per-adapter in the config |

Plus always: `PGBENCH_SAMPLE_ID`, `PGBENCH_THREADS`, `PGBENCH_OUTPUT_DIR`,
`PGBENCH_OUTPUT_VCF`.

The template only consumes already validated, read-only public inputs. The
complete declared set is available to `run_tool()` as `context.inputs`, keyed
by contract token (including `adapter_asset.<name>`); do not invent an
environment variable or read an undeclared path. A
tool-specific prerequisite (for example an RG requirement) may be checked in
the copied tool adapter, but it must not replace Core validation or silently
reinterpret input identity.

## Output contract (fixed)

- `outputs.vcf` must cover **every** candidate in `PGBENCH_CANDIDATE_VCF`
  (`candidate_output_contract: all_sites`).
- A candidate your tool did not genotype is emitted as `./.`
  (`absence_semantics: no_call`) — **never** materialize `0/0` for silence.
- `project_all_sites()` in `run.py` implements this; use it.
- Keep the native VCF in `tool-work/`. The template writes the deterministic
  `native-to-canonical-projection.tsv` trace there. A tool-native record
  outside the fixed canonical universe is retained with
  `outside_canonical_universe` (never silently counted); duplicate mapping to
  one canonical candidate fails. PG-F1 evaluator queries preserve the native
  representation; the all-sites VCF is only the canonical GT/accounting ledger.

## What will get your run invalidated

- Reading truth VCFs, target assemblies, relatives' genotypes, or any
  undeclared external callset (`information_contract` is fail-closed).
- Subsetting/re-filtering benchmark reads (worst on the ONT track:
  `adapter_selected_subset` is explicitly prohibited).
- Violating the sandbox: writing outside `PGBENCH_OUTPUT_DIR`, network
  access, or executing undeclared binaries.
- Hardcoding tool identity into the core (CI enforces: core sources must not
  mention any adapter id — see `tests/unit/test_adapter_boundary.py`).
- Declaring a memory ceiling you never measured; the scheduler on the local
  executor does not enforce it, the kernel does (exit 137).
