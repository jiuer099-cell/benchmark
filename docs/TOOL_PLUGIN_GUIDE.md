# Tool plugin contract

PGBench has a fixed benchmark shell and a replaceable tool middle. This avoids
forcing graph mappers, segment-discovery/integration pipelines, assemblers, and
genotypers through the same algorithmic steps.

## Three benchmark phases

1. **Common preparation**
   validates GRCh38/HG002 inputs, freezes the HG002-excluded pangenome panel,
   creates stable allele IDs and a blinded challenge universe, and records
   content provenance.
2. **Tool-owned middle**
   starts from the input authorized by the selected track. The adapter may map,
   discover segments, assemble, integrate calls, or genotype a frozen panel.
   Its `tool.yaml` declares the paradigm, tasks, accepted read technology,
   required inputs, output semantics, and execution isolation. Its
   `rule-registry.yaml` describes the tool-specific internal stages.
3. **Common postprocessing**
   validates the tool VCF, normalizes it, links it to stable pangenome alleles,
   evaluates it with Truvari/Aardvark/vcfdist against the same fixed GIAB truth
   universe, and produces the unweighted consensus and cross-track comparable
   reports.

The formal benchmark measures the entire tool-owned middle as one isolated,
reproducible execution unit. Internal stages remain tool-specific and are
listed in the plugin rule registry for auditability.

## Built-in adapters

| Adapter | Family | Starts from | Tool-owned stages |
| --- | --- | --- | --- |
| `kanpig` | long-read candidate genotyping | shared GRCh38 BAM | genotype, postprocess |
| `vg` | graph mapping/calling | PacBio/ONT FASTQ | GraphAligner map, vg pack, vg call |
| `pangenie` | short-read pangenome genotyping | Illumina FASTQ | prepare, index, k-mer genotype |

PanGenie and KanPIG are re-genotypers: they cannot discover an allele absent
from their candidate/pangenome panel. The vg adapter has `variant_sites`
semantics and can report novel graph-supported sites.

## Fairness contract

- All official runs must use HG002, the same GRCh38 reference build, primary
  truth profile, benchmark BED, and frozen score-profile hash.
- Tool-native read technology is allowed: PacBio tools use the registered
  PacBio evidence and PanGenie uses the registered HG002 Illumina evidence.
- A graph or population panel must exclude HG002 and NA24385.
- The tool-owned middle cannot read truth assets.
- `ConsensusScore` describes agreement on query results.
- `ComparableScore` uses the fixed eligible truth count as well as the query
  count, so false positives and false negatives both lower the score.
- Cross-technology scores compare complete pipeline utility. Claims about
  algorithm-only superiority must be restricted to one technology/task track.
- Runtime, memory, disk, and provenance are reported or used as validity gates;
  none is an accuracy-score weight.

## Add a user-provided genotyper

Copy `plugins/example_genotyper/` to a new plugin directory and edit:

- `tool.yaml`: unique ID/version, `paradigm`, tasks, read technology, mode input
  contract, output contract, runner/environment, and isolation policy;
- `rule-registry.yaml`: the actual tool-specific stages and their artifacts;
- `envs/environment.yaml`: exact dependencies;
- `run.py` (or another executable runner): consume only the `PGBENCH_*`
  environment variables and create exactly `PGBENCH_OUTPUT_VCF`.

The most useful runner variables are:

```text
PGBENCH_INPUT_FASTQ
PGBENCH_SHARED_ALIGNMENT
PGBENCH_REFERENCE_FASTA
PGBENCH_PANEL_VCF
PGBENCH_CANDIDATE_VCF
PGBENCH_GRAPH_DIR
PGBENCH_SAMPLE_ID
PGBENCH_THREADS
PGBENCH_MEMORY_MB
PGBENCH_OUTPUT_VCF
```

Register the plugin in a copied configuration:

```yaml
external_plugins:
  - id: my_genotyper
    manifest: plugins/my_genotyper/tool.yaml
```

For a BAM-based genotyper, implement `caller_only_shared_alignment`; for a
FASTQ-based method, implement `end_to_end_from_reads`. A genotyping-only tool
must emit every candidate record (`all_sites`, including `0/0` and `./.`).
A discovery/calling tool emits only detected variants (`variant_sites`).

Validate before a production run:

```bash
cd /home/luzhiting/hg002-grch38-pangenome-sv-benchmark

python workflow/scripts/validate_config.py \
  --config config/my-tool.yaml \
  --repo-root .

snakemake --snakefile Snakefile \
  --configfile config/my-tool.yaml \
  --cores 1 --dry-run
```

Unreviewed user code must use `bwrap` or `apptainer` in formal mode. Set
`trust_level: trusted_reviewed` with `sandbox_backend: none` only after a local
code review. The executor fingerprints plugin code and inputs, rejects input
mutation/path escape, validates the fresh VCF, applies a timeout, and records
attempt/log/provenance artifacts.
