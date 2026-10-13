# VG Q5 primary run plan — 2026-10-13

The new run is selected before native execution, without inspecting its score:
`vg_sr_pgf1_q5_v2_20261013`.

- Code: repaired VG adapter, with Q5 pack and all-reference-snarls calling
  (`vg call -a`). The prior Q0 score remains historical.
- Inputs: the same frozen HG002 FASTQs, GRCh38 reference, family-excluded graph
  and 18,164-candidate panel contract as the previous benchmark.
- Mapping: 16 threads, VG 1.73.0, exactly the legacy Giraffe recipe and fixed
  fragment parameters. Eight GAM checkpoint metadata files identify this
  matching mapping identity. The new run receives hardlinks to those GAMs;
  the adapter verifies sizes and SHA-256 before reuse.
- Recomputed stages: Q5 pack, call, native projection, evaluator inputs, all
  three evaluators and final score/seal. No Q0 pack, native VCF, evaluator
  evidence or old score is imported.
- Execution: an isolated code checkout; Snakemake 9.23.1, 16 cores, shared
  existing Conda environments and 262,144 MB declared native memory. Existing
  VG memory gates and abort thresholds remain in force. vcfdist isolated
  recovery uses two workers as a scheduling setting.

The 23-job end-to-end DAG passes dry-run checks. Required native and Core Conda
environments already exist, including VG 1.73.0. This plan does not claim the
new score is complete or that all checkpoint contents have already passed
their full SHA scan.

Server locations:

```text
Acceptance: /home/luzhiting/benchmark/acceptance/vg_pgf1_q5_v2_20261013
Results:    /home/luzhiting/benchmark/runs/vg_sr_pgf1_q5_v2_20261013
Code:       /home/luzhiting/benchmark/adapter-development/vg-paragraph-validation-20261013
Controller: /home/luzhiting/benchmark/recovery-assets/run-vg-q5-20261013.py
```

Paragraph is assessed read-only for this request. Its completed native calls
are not affected by the validation-only repair, and its current vcfdist
evaluation should continue. The new run does not stop, replace or restart that
evaluation, nor deploy Paragraph's deferred production changes.
