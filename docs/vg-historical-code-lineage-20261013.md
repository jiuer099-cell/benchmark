# VG historical result and adapter versions — 2026-10-13

The completed `vg_sr_pgf1_production_v1_20260922` score must be interpreted using
its recorded stage versions, not the latest files in `code/plugins/vg`.
`vg-pgf1-code/plugins/vg/run.py` and `projection.py` are symlinks to the current
plugin directory; those historical-looking paths are therefore mutable.

## Verified stage mapping

| Stage | Code commit | Verification |
|---|---|---|
| Original native calling | `f7f3709054421c8f1b11b44bed0425d1d87b8fbb` | Runner SHA `d598b6a1aa8f5c1c7f153a3250cf73d8a8c5966847f14b67dd2916ffcd6a1122` and tool manifest SHA `8006cdf7ddbd2048f56242712949a004c9ecaf4fb77c94f5e4e833ec77642738` match the successful attempt |
| Native reprojection | `6a12f14213a18c20a130b1201f3cc97f221c9c48` | Projection SHA `463bdbbf57e72ed6190c8bc84ffe87dd854090e83eaf147b4f992f8959274662` and tool manifest match `projection-v2-20261010/producer.json` |
| Reprojection publication | `6a12f14213a18c20a130b1201f3cc97f221c9c48` | Publisher SHA `9bca41f2d851a9f7c43015146d226de4392897f7ce3fdc6bb6fe9258f36c0446` matches the publication rule manifest |
| Repaired adapter for future executions | `310141c7818b7f5767b0e61bcdd2dc9adf99301a` | Current validation repair; it did not generate the historical score |

The old score's execution policy remains `legacy_pack_Q0_native_replay`. Its
raw scoring VCF SHA is
`d424bd29c38db585f1cf25d199d3038248c382d54538d062f510444655aca4d0`.
The raw hash matches both the projection audit and publication rule manifest.

## Historical code archive

An additive archive is stored at:

```text
/home/luzhiting/benchmark/recovery-assets/vg-run-version-ledger-20261013/
```

It contains a separate `plugins/vg` Git archive for native calling, reprojection,
publication and the future repaired adapter, plus copied historical evidence and
`version-ledger.json` with per-file and archive hashes. Archive files are set
read-only. Each historical script is byte-verified against the recorded SHA,
rather than assumed to match merely because its Git date or directory fits.

The attempt, producer/publication manifests, raw VCF, projection audit, score
and seal hashes were checked before and after creating this archive and remain
unchanged. No historical rule path or sealed metadata was rewritten. The live
plugin layout still has a single development entry point; the archives are
historical evidence, not additional active adapters.

## Usage

- To interpret or reproduce an old stage, use its recorded commit/archive and
  input/output hashes. Do not substitute current symlink targets.
- A future call or projection with the repaired adapter needs a new run identity
  and its own code/output manifests. Do not overwrite the historical VCF or
  assign its score to the repaired adapter.
- No re-run is needed merely to add this version ledger. It does not establish
  a new Q5 score or change the old Q0 policy.
