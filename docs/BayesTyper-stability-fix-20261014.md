# BayesTyper stability repairs (2026-10-14)

The completed default run `bayestyper_sr_pgf1_default_v1_20261013` retains
PG-F1 42.01973302379571 and its original seal. All 29 protected artifacts
were independently verified unchanged when the repair was deployed.
Its noise stability is `REVIEW_REQUIRED_CHAIN_DISAGREEMENT`: classical
split R-hat 4.696478865989162, 20 complete chains, 250 retained draws each.
The chain means range from 0.177451476 to 0.208186076. The diagnostic is a
screening flag, not a new upstream acceptance rule or a revised score.

## Changes

- Log unmapped candidate counts separately from actual missing GT. The
  default run has 5666 missing GT out of 18164 candidates; its historical
  `no_call=0` log label actually meant no unmapped candidates.
- Bind every unit's complete noise table and a standalone diagnostic JSON
  by SHA256 in new native audits. Record diagnostic state without altering GT,
  QUAL, FILTER, GQ or scoring.
- Provide `stability.py --acceptance ... --attempt ... --output ...` for
  explicit successful-attempt selection. Verify native/published VCF,
  configuration and archived executed-source bindings. Multiple historical
  attempts are retained. Historical noise tables with no execution-time hash
  are explicitly reported as such; no historical binding is fabricated.
- Correct the old SNV caveat: v1.5 accepts few or no SNVs, with low eligible
  variant counts producing a warning. The default run's native log records
  18079 eligible variants. See the [official v1.5 release](https://github.com/bioinformatics-centre/BayesTyper/releases/tag/v1.5).
- Keep sequential defaults in `tool.yaml`. Add the explicitly documented
  optional mode in `tool.noise-genotyping.yaml` and its own example config.
  Joint mode's upstream table has no final `0/0` mean row; the parser accepts
  this format only when the joint policy is declared.

## Verification

30 tests passed on the server, including real KMC count saturation and
singleton retention, native genotype/quality preservation, joint/default
commands, whole adapter audit generation, truncated chains, chain
disagreement, archived sources, failed attempts and tampered files.

The server's old completed-run supplemental diagnostic is stored under
`acceptance/bayestyper_pgf1_default_v1_20261013/supplemental-noise-stability/completed-review-20261014`.
The old status now references this diagnostic. The original native audit,
score, seal and frozen implementation manifest are preserved.

## Separate joint-mode experiment

The prepared run is `bayestyper_sr_pgf1_joint_noise_v1_20261014`, with a frozen
independent Core checkout and a fresh native cache. The full 17-job dry-run
passes. Preparation preserved all 29 old protected artifacts.
The documented optional mode was selected to address observed chain
disagreement. It is an experimental diagnostic run and does not automatically
replace the default primary result or select the larger score.

Joint-mode stability must be evaluated from its completed output. Code fixes
and a successful dry-run do not establish biological stability or performance.
Existing KMC/Bloom/cluster files were not imported because their reuse lacks a
separate reviewed cache contract. The new run regenerates its native intermediates.
