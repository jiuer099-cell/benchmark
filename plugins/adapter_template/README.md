# PGBench adapter kit (template)

This directory is a **template, not a runnable tool**. Copy it, rename it,
fill in your tool, and pass the gates described below. The benchmark core
never learns about your tool by name; it only sees this adapter's manifest.

## Frozen formal primary policy

所有工具采用预先冻结的官方默认或推荐运行配置。除输入兼容、资源限制和必要的格式转换外，
不添加额外的分型过滤或拒判参数。所有工具只产生一个正式 PG-F1 分数。

Freeze the official source release, explicit flags, implicit output defaults and
any justified exception **before** production. Do not add quality/coverage
filters or no-call switches beyond the preselected official configuration.
Native tool-internal default filters remain part of that configuration; do not
disable or alter them to improve a score. No post-score tuning or selection of
the higher run is allowed. Superseded runs are historical archives, not extra
formal scores. Keep Core panel, truth, evaluator and PG-F1 rules unchanged.

## Directory layout

```
adapter_template/
├── README.md              <- this guide
├── tool.template.yaml     <- copy to <your_tool>/tool.yaml and fill
│                             (schema: workflow/schemas/tool.schema.yaml)
├── run.py                 <- implement and review run_tool(); record invocations
├── adapter_projection.py  <- strict native-preserving projection (do not bypass)
├── gate-plan.template.json <- fill with real per-tool test commands/evidence
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
| `sr_illumina` | short | `illumina_pe` | frozen shared BAM **or** raw R1/R2 FASTQ | paragraph, graphtyper2 (BAM); pangenie, vg, bayestyper (FASTQ); varigraph is an external adapter outside Core |
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
   Fill the mandatory `parameter_contract`, replace every source placeholder,
   and preserve the frozen formal-primary policy. Pin the actual tool runtime
   and lock file, not just the adapter's own version number.
3. **Implement and review `run_tool()`** in `run.py`. You get input *injection* checks
   and the all-sites projection for free; only the "invoke my tool" part is
   yours. FASTQ/BAM integrity, SHA256, read counts and dataset identity are
   Core-managed immutable validation, performed once before every adapter.
   Do not reimplement or bypass them in an adapter.
   Use `run_native(context, argv)` for every invocation, build argv from
   `context.parameters`, and keep original VCFs and conversion inputs in
   `tool-work/`. For native Python extensions, launch `native_python(tool)`;
   do not import them into Core's interpreter. A resource exception such as
   chunking needs a method-valid, evidence-complete plan; per-contig splitting
   is not automatically appropriate for every tool.
4. **Pass actual gates**: config `validate` → real runtime/mini-panel/conformance
   tests → dry run → production with a fresh `run_id`. Copy and fill the gate
   plan, then execute:

   ```sh
   python workflow/scripts/validate_adapter_conformance.py \
     --tool-manifest plugins/<tool>/tool.yaml \
     --checks plugins/<tool>/conformance/gate-plan.json \
     --output plugins/<tool>/conformance/report.json --execute
   ```

   Six real commands must exit successfully and produce fresh nonempty evidence
   under `conformance/`; a no-op command reusing old files is rejected.
   The plan is frozen before execution and checked before/after every gate.
   Plans, receipts, logs and proofs must be distinct, non-symlink paths under
   `conformance/`; proofs cannot overwrite a plan or another gate's evidence.
   Each command, timeout, log, evidence hash, source inventory, parameters and policy
   is bound to the receipt. Hand-written `PASS` labels are rejected. The
   generic executor verifies this receipt before native execution. A changed
   source/env/manifest/evidence invalidates it; re-run gates, never edit a
   PASS receipt. For production, review and set `trusted_reviewed` and a real
   sandbox **before generating the final receipt**. Conformance execution is
   a development operation; run untrusted tool probes in an isolated test
   environment. Receipt validation does not replace human semantic review.
   Receipts also bind the validator implementation; a validator upgrade
   requires rerunning conformance gates, not rerunning existing sealed scores.
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
`PGBENCH_OUTPUT_VCF`. This reviewed kit also requires the Core-verified
`PGBENCH_CONFORMANCE_REPORT` and declared `reference_index`.

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
- Matched records preserve **native CHROM/POS/REF/ALT, QUAL, FILTER, INFO and
  FORMAT**. Only the candidate ID, target sample selection and trace annotations
  are added. No-call placeholders alone use the candidate representation;
  they never inherit candidate QUAL/PASS as caller evidence.
- Select the target sample by its header name, not the first sample column.
  Claimed IDs must match the sequence or reference-verified haplotype;
  ambiguous/duplicate mappings fail the entire run, never choose by quality.
- Without a stable candidate ID, minimal sequence edits are left-aligned
  against the frozen reference for exact lookup, including repeat-shifted
  insertions/deletions. This only determines linkage: the emitted native
  POS/REF/ALT and genotype remain unchanged. Multiple equivalent candidates
  are ambiguous and fail even if a record claims one of their IDs; do not
  match by distance, overlap, quality, truth or evaluator score.
- The generic projection accepts **biallelic, sequence-resolved diploid VCF**.
  Multi-ALT, symbolic/breakend and invalid-GT inputs fail explicitly. Supply a
  declared official dosage-aware conversion before projection and retain its
  original input, output and command. This structural safety check is not a
  biological quality filter or a successful no-call substitution.
- Keep the native VCF in `tool-work/`. The template streams a deterministic
  `raw/calls.native-trace.jsonl.gz` containing full original records. Core
  publishes it with `raw/calls.projection.json`. A tool-native record
  outside the fixed canonical universe is retained with
  `outside_canonical_universe` (never silently counted); duplicate mapping to
  one canonical candidate fails. PG-F1 evaluator queries preserve the native
  representation; the all-sites VCF is only the canonical GT/accounting ledger.
- Declare **plain** `raw/calls.vcf`; ordinary gzip is not BGZF and is not a
  formal adapter output shortcut. Core owns compressed/indexed downstream data.
- Keep both audit flags enabled. `calls.native-audit.json` binds native artifacts,
  commands and eligibility/parameter receipt; the projection audit binds the
  native trace, frozen reference/panel and output. The final Core seal includes
  these artifacts in `score-package.json`, `provenance-audit.json`, `seal.json`.
  Never delete retained `tool-work/` while a release depends on its hashes.
  Disposable temporary indexes may be removed by the adapter before returning;
  all remaining files are retained and hash-bound, so keep only required proof.
- Trace memory is bounded by the candidate set, not the full native record
  count. `no_call` counts both absent records and observed missing GTs.

## Minimum test coverage for a new tool

The six gates must include real version/dependency probes, a fixed mini-panel
with exact expected outputs, native-to-candidate mapping assertions, the
unchanged three-evaluator PG-F1 path, and a bounded real-input pilot. Include:
wrong/missing/duplicate samples, ID/allele conflicts, ALT dosage conversion,
equivalent indel representations, native QUAL/FILTER/FORMAT retention,
missing-call accounting, ambiguous mappings, failed subprocesses, trace/hash
tampering and eligibility rejection. Do not make all six commands a dummy
successful process; Core checks execution and hashes, not the scientific
completeness of arbitrary test programs. Review the actual assertions.

Run native probes with `--execution-purpose development_only` so a new adapter
can generate conformance evidence before it has an eligibility receipt. Core
injects its declared tool ID, frozen parameters and policy; development outputs
are marked unsealed probes, not formal results. Formal mode requires the verified
receipt. This avoids a circular "receipt required to run its own gates" workflow.

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
