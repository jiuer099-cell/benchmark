# GATK-SV 1.1.1: frozen-BAM compatible discovery adapter

Production entry: `run_rg_compatible.py --config <frozen-config.json> --execute`.
`run.py` contains shared official-WDL/backend helpers; it is not the production
entry for the known RG-incomplete HG002 BAM. This is an independent SV-discovery
experiment, not an eligible fixed-panel PG-F1 plugin or leaderboard entry.

## Frozen input and permitted compatibility transformation

The existing shared HG002 BAM and BAI remain immutable and are the sole original
alignment inputs. This adapter creates an attempt-private intermediate BAM with
official samtools 1.10 `addreplacerg -m orphan_only -R <existing-header-RG-ID>`.
This only fills missing RG tags. It requires one source-declared read group with
the requested sample; it must not collapse multiple libraries/read groups.

Before any biological workflow starts, every transformed SAM record is traversed.
The auditor requires the declared RG on every read and removes only RG before
hashing the complete remaining records in order. This hash must equal the
pre-existing frozen full-record SHA256. Coordinates, CIGAR, flags, sequences,
qualities, MAPQ, pairing, non-RG auxiliary fields and record order cannot change.
The frozen source BAM/BAI identities are checked before/after; intermediate BAM,
index, commands and hashes are recorded in `compatibility/`.
The source byte SHA256 is inherited explicitly from its pinned freeze manifest,
not misrepresented as newly rehashed; full SAM-record equivalence is checked anew.

## Reference compatibility and official workflow stages

1. Run the unmodified pinned official `Manta.wdl` with the exact original alignment
   reference FASTA/FAI, whose contig names, order and lengths must match the BAM.
2. Run unmodified official `GatherSampleEvidence.wdl` with its original Broad
   resources and the precomputed native Manta VCF. Shared sequences are checked
   against the frozen alignment dictionary MD5. Default biological filters and
   thresholds are unchanged. MELT and optional metrics remain disabled as in the
   original single-sample configuration. The official `run_localize_reads=false`
   switch avoids an unnecessary copy on the shared filesystem; it changes I/O,
   not read or genotype selection.
3. Require nonempty PE/SR/SD, positive coverage and nonzero accepted reads before
   downstream inference. Failed attempt 001's empty evidence is never reused.
4. Supply these native outputs through official `case_*` inputs to the unmodified
   `GATKSVPipelineSingleSample.wdl`. Preserve native SV types, genotypes and filters;
   publish the native final VCF only on workflow success.

No WellformedReadFilter bypass, RG-check exemption, truth access, panel projection,
quality-filter relaxation or score-based parameter selection is permitted.

## Provenance and isolation

Build the record auditor on the target Linux host before running the RG policy
or its unit tests: `gcc -O2 -o plugins/gatk-sv/rg_record_audit
plugins/gatk-sv/rg_record_audit.c -lcrypto` (one shell command).

`resume.py --config CONFIG [--execute]` supports the authenticated attempt-002
gzip-validation checkpoint. Its manifest must bind completed Manta/evidence
metadata and outputs, the alignment identities and the original configuration.
It rejects changed checkpoints and settings, mounts the prior attempt read-only,
and writes a new attempt. It does not resume arbitrary failed Cromwell tasks.
The repaired evidence gate handles both gzip and plain counts.

Cromwell localizes inputs with `copy` only. Official gCNV tasks run `gunzip`
in place: symbolic links produce ELOOP, while hard links may let task writes
affect upstream data. Each call and retry therefore receives regular files.
This increases disk usage; budget space for task input copies. Biological WDL,
filters and resource versions remain pinned. See
`tests/integration/gatk_sv_localization_probe.py` for a real two-call
Cromwell/Apptainer regression with source-integrity checks.

Use a new attempt directory. Keep original input identities, conversion report,
program and official WDL hashes, pinned SIF manifest, stage inputs, metadata and
logs. Each stage is network-disabled, receives immutable resources read-only,
and has a private writable execution root. Earlier stages are mounted read-only
for subsequent stages. Task concurrency is one; CPU and declared task memory are
bounded by the frozen configuration. This is not OS-level RSS enforcement.

Validation consists of unit tests, official small-data RG conversion with complete
record comparison, positive default-filter GATK coverage/SR checks, original-reference
Manta configuration and official WDL/input validation. It is not a claim that the
full HG002 biological workflow has already succeeded.
