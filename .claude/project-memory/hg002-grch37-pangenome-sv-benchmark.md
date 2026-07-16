---
project_id: hg002-grch37-pangenome-sv-benchmark
project_slug: hg002-grch37-pangenome-sv-benchmark
repo_root: /Users/wanghao/Documents/Codex/2026-07-16/referenced-chatgpt-conversation-this-is-untrusted
vault_root: /Users/wanghao/Desktop/risefl_mvp/memory/Research/hg002-grch37-pangenome-sv-benchmark
hub_note: Research/hg002-grch37-pangenome-sv-benchmark/00-Hub.md
language: zh
last_sync_at: 2026-07-16T11:57:20Z
status: active
auto_sync: true
code_state: design-only
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

- 用户审阅逐 point、rule registry 和 Obsidian 同步契约已修订的设计规范。
- 审阅通过后制定实施计划并开始代码实现。
- 每次代码语义变更同步 Obsidian canonical design note。

## Recent Sync Summary

- 2026-07-16：建立 repo 与 Obsidian 项目绑定；当前代码状态为 design-only。
- 2026-07-16：明确 caller-only/end-to-end、100-point point dictionary、pangenome allele linkage 和全 rule registry。
