# BayesTyper fixed-panel adapter

This adapter runs BayesTyper 1.5 from the benchmark-supplied paired Illumina
FASTQs and the frozen canonical candidate panel. The plugin registration must
set `tool_index` to a frozen directory containing `canon.fa` and `decoy.fa`
from the GRCh38 BayesTyper reference bundle. Those files are hashed into the
resolved-input and allowed-information manifests.

The adapter uses the upstream project's documented k=55 KMC, `makeBloom`,
`cluster`, and `genotype` workflow. It then projects native calls onto every
canonical candidate; candidates without a resolved native call are emitted as
`./.` instead of being dropped.

## Runtime and declared assets

The server example is `config/config.bayestyper.example.yaml`. It declares
`adapter_asset.binaries` pointing to the self-contained official static
`bayesTyper`, `bayesTyperTools`, and KMC 3.2.4 binaries. The whole asset
directory is fingerprinted and mounted read-only by Core. This avoids the
Bioconda BayesTyper build that read a valid KMC database as zero k-mers and
crashed. Do not replace the asset with a symlink to an undeclared download.
The Conda YAML supplies Python and bcftools; it deliberately omits BayesTyper.

KMC uses `-k55 -ci1 -cs255 -sm`, injected threads, and the memory budget
minus 1024 MB. Genotyping units run sequentially. Cluster and genotype receive
the injected random seed. Native Gibbs iterations and filtering stay at their
official defaults; smoke-test iteration reductions are not production options.
The samples TSV uses `M`, accepted by the upstream parser. The supported
comparison universe is chr1–chr22, where M/F have identical diploid ploidy.

BayesTyper 1.5 reads sample counts into unsigned bytes and rejects values above
255. KMC therefore saturates counts at 255 while retaining singleton and
high-frequency k-mers. `-cx255` would remove high-frequency k-mers and must not
replace `-cs255`. Databases produced with the former `-cs65535` recipe are not
compatible inputs and must be regenerated or converted and verified separately.

## Output and regression protections

`raw/calls.vcf` is plain VCF accepted by the formal executor's current output
contract. Native outputs remain compressed in the attempt directory. Input
VCF is converted from BGZF to ordinary gzip for BayesTyper, without changing
records. Native QUAL/FILTER and GQ survive projection; no extra quality filter
or truth-dependent selection is applied. Candidate IDs and graph paths must
be corroborated by reference-haplotype equivalence. Multiallelic GT becomes
ALT-specific diploid dosage. Missing/partial calls remain no-calls; conflicting
mappings get `ambiguous_mapping`, never a quality-based winner.

The native VCF, projection trace and audit sidecars bind every conversion to
its source. Core owns candidate status, evaluator query and scoring outputs.
No complete reference/FASTQ hashing pass is repeated inside the adapter.
Fresh attempt directories prevent stale k-mer or cluster reuse.

BayesTyper v1.5 uses all non-nested variant types for sequential noise
estimation and accepts few or no SNVs. Fewer than 100,000 eligible variants
triggers a warning, not a requirement to add 100,000 isolated SNVs. The adapter
does not add a variation prior from another cohort or target-specific calls.
The current configuration remains experimental until real-input conformance
and scoring gates pass. Synthetic checks do not establish biological accuracy.

## Noise stability and the explicit joint-noise configuration

Every completed unit's noise table and an independent diagnostic JSON are
hashed in the native audit. Complete chains are screened with classical split
R-hat after discarding the first 100 iterations and retaining 250 draws from
each of 20 chains. The local threshold 1.05 flags review; it is not an upstream
acceptance rule. `NO_SCREEN_FLAG` is not proof of convergence or accuracy.
Diagnostic flags do not rewrite GT, filtering, scoring or the seal.

The log reports `unmapped_candidates`, `gt_called` and `gt_no_call` separately.
An addressable variant can still have a native missing genotype.

`tool.yaml` keeps the sequential official defaults. Use the separate
`tool.noise-genotyping.yaml` and `config/config.bayestyper.noise-genotyping.example.yaml`
for the documented optional `--noise-genotyping` policy, under a new frozen
run ID. The latter estimates noise and genotypes jointly, including nested
variants. It can be slower and use more memory. Its noise table intentionally
has no final 0/0 summary row; the diagnostic accepts that upstream format only
for the explicitly recorded joint policy. No other Gibbs or filter parameter
changes. Source: https://github.com/bioinformatics-centre/BayesTyper/releases/tag/v1.5

Historical runs can be screened without deleting failed attempts or replacing
executed-source records (Python 3.12):

```bash
python3.12 -B plugins/bayestyper/stability.py \
  --acceptance /absolute/path/to/acceptance \
  --attempt SUCCESSFUL_ATTEMPT_ID \
  --output /absolute/path/to/acceptance/supplemental-noise-stability/NEW_REPORT
```

The selected successful attempt must match published VCF and native hashes,
the configuration and preserved executed-source hashes. Historical runs that
did not hash noise tables at execution are explicitly identified as such;
an independent report cannot manufacture a missing historical hash binding.
The command writes only a new supplemental report. Old scores and seals remain
bound to their original configuration; comparison runs are never chosen by
the larger score.

Review fixes on 2026-10-13:

- `reference_index` is required and must be the declared `.fai` beside the
  reference FASTA. This is checked before counting reads.
- Candidates are checked for ordering before KMC. Same-position records are
  padded to a common REF and merged for BayesTyper, retaining their IDs. The
  original biallelic panel remains the projection and scoring universe.
- A called `*` or other non-sequence placeholder cannot become a complete
  diploid call such as `0/0`; it is represented as `./.` with the original
  genotype retained in the native trace. Ordinary sequence ALT dosage is
  unchanged. Empty/headerless native VCFs are rejected.
- The reference paths in the example include the server's `.v1` directory.
  Read/panel/provenance entries are still example values and must be bound to
  the selected frozen run before a biological benchmark.

The 2026-10-13 Core/bwrap integration check used seven synthetic candidates,
including two at the same position. All IDs survived and all seven genotypes
matched the simulated values with the native default inference parameters.

Validation on 2026-10-12: actual default-parameter tool chain in Core's
network-disabled bwrap sandbox, six synthetic DEL/INS spanning 0/0, 0/1 and
1/1; all six observed genotypes matched the simulated values. Projection and
resource-limit regressions are in `tests/unit/test_bayestyper_adapter.py`.
