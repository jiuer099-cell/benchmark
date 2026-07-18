---
project_id: hg002-grch37-pangenome-sv-benchmark
project_slug: hg002-grch37-pangenome-sv-benchmark
repo_root: /Users/wanghao/Desktop/lzt/hg002-grch37-pangenome-sv-benchmark
vault_root: /Users/wanghao/Desktop/risefl_mvp/memory/Research/hg002-grch37-pangenome-sv-benchmark
hub_note: Research/hg002-grch37-pangenome-sv-benchmark/00-Hub.md
language: zh
last_sync_at: 2026-07-19T00:00:00Z
status: active
auto_sync: true
code_state: partial
canonical_note: Research/hg002-grch37-pangenome-sv-benchmark/Writing/HG002-GRCh37泛基因组SV-Benchmark设计规范.md
---

# Project Memory: hg002-grch37-pangenome-sv-benchmark

## Current Focus

- 以 hs37d5 backbone 和 population SV panel 构建泛基因组 benchmark。
- 固定 Truvari:Aardvark:vcfdist 为 40:35:25；每个适用工具 tuple 输出且只输出一个 PGBenchScore，不排名。
- 让内置与外部 genotyper 都由 benchmark 实际运行并产生 VCF。
- 为每个 Snakemake rule 建立规范 registry、manifest、hash lineage 和资源记录。
- 默认正式模式为 end-to-end；caller-only 作为共享 alignment 上的诊断模式。

## Active Tasks

- 冻结 hs37d5、GIAB truth、population SV panel 和 stratification checksums。
- 实现真实 Truvari、Aardvark、vcfdist wrapper/parser。
- 接入首批内置 caller/genotyper，并在 Linux/HPC 上验证正式隔离环境。
- 在目标服务器运行 HG002 CLR 全基因组验收。

## Recent Sync Summary

- 2026-07-16：建立 repo 与 Obsidian 项目绑定；当前代码状态为 design-only。
- 2026-07-16：明确 caller-only/end-to-end、100-point point dictionary、pangenome allele linkage 和全 rule registry。
- 2026-07-17：Phase 1 synthetic MVP 已贯通外部工具执行、泛基因组 panel、标准化、allele linkage、评分、报告和 rule lineage；代码状态更新为 partial。
- 2026-07-17：91 个本地测试通过；7/7 pre-score jobs 可追溯，synthetic score 为 provisional。
- 2026-07-19：完成 Phase 1 最终封存链路、标准 metrics 与无排名长表；177 个测试通过，8/8 pre-score jobs 及 11-job final lineage 可追溯。
- 2026-07-19：补齐失败重跑恢复、插件完整代码指纹、跨 manifest 上下文一致性、ratio 语义和 finalizer-to-report hash 门禁；准备以干净 commit 重跑并同步 Obsidian。
