# PG-F1 实现与验收说明

早期 ME-F1 产物仅为 provenance 留存，不进入最终 PG-F1 leaderboard。
PanGenie 已有 genotyping native output 可以复用，但 evaluator、normalizer、
PG-F1 与 admission 链必须按本合同重新回放。

版本身份分为两层：

- benchmark release：`PGBench-SR-Illumina-PGF1-v1.0`；
- score contract：`pgf1_v1`（代码中的 profile ID 为
  `pgbench_pgf1_v1`）。

以后若只修改评分合同，应升级为 `pgf1_v2`；track/release 身份单独管理，
不会把二者混成同一个版本号。

## 正式链路

```text
frozen 18,164 panel + hidden truth ledger
  → canonical_units.tsv.gz (candidate × ALT)
  → representation-preserving evaluator query
  → Truvari / Aardvark / vcfdist raw per-unit trace
  → normalized evidence (NONE / LM / AM → 0 / 0 / 1)
  → fixed 2-of-3 vote + Core exact unphased GT
  → pgf1_evidence.tsv.gz → PG-F1 → leaderboard admission
```

`all-sites.vcf.gz` 只承担 canonical GT、no-call 与 TP/FP/FN accounting；
它不是 evaluator 的表示重写输入。任何 1:N、N:1、N:M 或 vcfdist 部分
component mapping 以 `AMBIGUOUS_COMPLEX_MAPPING` fail-closed，不会静默删
unit、改分母或扩大 `UNSCORABLE`。

## 方案落实点

- `config/pg_f1_scoring.yaml` 与 `workflow/scripts/pgbench_pgf1_contract.py`
  冻结 PG-F1、phase/order-insensitive dosage、固定 2-of-3、失败状态。
- `materialize_canonical_units.py` 冻结稳定 `candidate_id:A<n>` unit ID、
  ALT-specific truth GT 和 unit-set SHA256；不重新选择 18,164 universe。
- `normalize_evaluator_evidence.py` 要求每个 unit 都有可追踪 raw ledger
  hash、normalizer version/hash、one-to-one mapping；禁止 aggregate→unit
  inference。缺任意 evaluator 为 `CONSENSUS_EVIDENCE_INCOMPLETE`，绝不作
  2/2 vote。
- `materialize_pgf1_evidence.py` 只由 Core 一次性检查 exact GT，并按冻结
  TP/FP/FN 表计分。caller/evaluator/adapter 故障会抛出 formal failure，
  不会变为 PG-F1=0。
- `validate_leaderboard_admission.py` 绑定 track/release、panel/unit、truth/
  BED/reference、PG-F1 contract、evaluator bundle、query construction 和
  normalizer contract；未 PASS 的 score 不可入榜。
- `validate_adapter_conformance.py` 定义 schema、真实 runtime、mini-panel、
  canonical mapping、PG-F1 evidence、real input 六项全部 PASS 才可 production。
- `tests/conformance/fixed-mini-panel/` 覆盖 C01–C15，且固定 voter、GT、
  missing-evidence 和 complex-mapping 失败语义。

## Adapter template 审阅修复

通用模板已改为仅消费 Core 已验证的只读输入。它不再逐 adapter 重做 FASTQ
gzip / SHA256 / read-count 或共享 BAM 的通用完整性检查；这保证所有 adapter
复用同一份 Core immutable input-validation cache。模板还修复了实际注入环境
变量名，并拒绝 duplicate candidate/allele mapping，避免静默选一个匹配。
