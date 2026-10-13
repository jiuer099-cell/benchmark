# vg and Paragraph native validation repair — 2026-10-13

Paragraph could silently turn a missing chunk result into all-sites no-calls,
and its projector read the first native sample without checking its identity.
The adapter now rejects missing/empty chunk VCFs, missing or duplicate headers,
wrong or multiple samples, records preceding the header, and records without
GT. Both fresh/count-replay projection and authenticated frozen-VCF replay use
the same single-sample checks. A valid header-only VCF retains the existing
unobserved-site no-call semantics.

Paragraph's reference index is now a required input. Its adjacency to the
reference, nonempty content and FAI parsing are checked before native commands.
The depth/read acquisition policy, candidate projection and 10,000-character
evaluator length protection remain unchanged.

vg now rejects empty/headerless native VCFs instead of publishing a successful
all-no-call projection. GT values referencing non-sequence ALT alleles such as
`*` produce no-calls. Uncalled placeholders leave sequence dosages unchanged,
and the full original native genotype remains in the trace. Pack/call recipes
and mapping checkpoint logic are unchanged.

## Verification

- 47 focused regression tests pass with Python 3.12 on the server in an
  isolated checkout, including existing Paragraph and vg projection tests.
- Local new/old comparisons confirm unchanged valid Paragraph output and
  unchanged valid vg projection/decoded native trace. Gzip filename/mtime
  transport hashes are excluded from that comparison.
- Explicit tests cover wrong/multiple samples, missing headers/chunk outputs,
  reference preflight, called and uncalled placeholders, normal multiallelic
  sequence dosage, frozen replay, and preservation of the length guard.
- The deployment script checks baseline hashes before any changes, backs up
  previous files, and supports applying `vg` or `paragraph` separately.

## Production deployment boundary

At inspection, the existing Paragraph vcfdist rule was still running and listed
`plugins/paragraph/tool.yaml` as a guarded input. `pgbench_rule_exec.py` verifies
that input's content, size and mtime again after the evaluator exits. The
Paragraph production checkout must therefore retain its current files until
the rule successfully completes. A change during this run would invalidate the
guard even though calling has already finished.

Publishing this commit in a separate checkout does not modify that running
rule's inputs. The vg production files can be updated after confirming no vg
caller is active. Existing raw VCFs, evaluator results and sealed scores are
not rewritten by this code repair. Future runs must record the new adapter
hashes. GitHub publication is separate from the deferred Paragraph production
deployment; this commit does not launch an automatic updater or new analysis.
