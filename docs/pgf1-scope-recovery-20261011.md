# PG-F1 native query-scope recovery (2026-10-11)

The VG native projection preserves QUAL/FILTER. One submitted positive query
(`PGSV_80390f246bbe8fe8b906`, FILTER `lowad`) is correctly excluded by the
frozen `pass_or_unfiltered` evaluator policy. All three raw evaluator ledgers
contain its explicit resolved, scope-ineligible, non-event, zero-vote row.
The v1 normalizer incorrectly rejected that row as incomplete evidence.

The v2 normalizer independently recomputes scope from the actual evaluator
query and frozen profile, reconciles IDs and ALT dosage with all-sites output,
and permits only corroborated resolved zero-vote exclusions. Missing rows,
unresolved positives, inconsistent scope, or acceptance of non-events remain
formal failures. The excluded call keeps its genotype in PG-F1 accounting;
the correction does not turn it into a no-call or remove it from the panel.
Raw evaluator ledgers and their source provenance are not rewritten.

Recovery uses explicit seal targets and only these allowed/forced rules:
`normalize_pgf1_evidence`, `materialize_pg_f1`, `materialize_pgf1_invariants`,
`validate_pgf1_leaderboard_admission`, `finalize_pgf1_release`.
VG and Paragraph retain their frozen queries and completed evaluator outputs.
Panel generation, adapter execution, mapping, packing, genotyping and all
three evaluator executions are outside the recovery rule whitelist.

Both server checkouts are advanced to the same tested recovery commit before
scoring. Historical run contexts and evaluator manifests retain their original
versions; the finalization manifest records the live recovery code version.
VG retains its legacy Q0 pack and native calls, as declared in its projection
audit; this recovery does not retroactively claim the new Q5 recipe was used.
PanGenie and GT2 are not modified or rerun by this recovery.
