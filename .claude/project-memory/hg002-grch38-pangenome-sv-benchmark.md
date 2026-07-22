---
project_id: hg002-grch38-pangenome-sv-benchmark
project_slug: hg002-grch38-pangenome-sv-benchmark
repo_root: .
vault_root: memory/Research/hg002-grch38-pangenome-sv-benchmark
hub_note: Research/hg002-grch38-pangenome-sv-benchmark/00-Hub.md
language: zh
last_sync_at: 2026-07-21T00:00:00Z
status: active
auto_sync: true
code_state: partial
canonical_note: Research/hg002-grch38-pangenome-sv-benchmark/Writing/HG002-GRCh38泛基因组SV-Benchmark设计规范.md
---

# Project Memory: hg002-grch38-pangenome-sv-benchmark

## Current Focus

- 使用 GRCh38 参考基因组和 GIAB HG002 GRCh38 v5.0q truth set。
- 使用 Truvari、Aardvark 和 vcfdist 的一致性结果进行综合测评，不使用权重打分。
- 保持路径相对于仓库或运行配置可解析，支持复制到 Windows Server 后运行。
- 由 Snakemake 记录规则、清单、哈希谱系和资源使用情况。

## Active Tasks

- 在 Windows Server 冻结真实 GRCh38、HG002 和泛基因组资源及其校验和。
- 接入并验证真实 Truvari、Aardvark、vcfdist wrapper/parser。
- 使用真实数据完成全流程验收并生成最终报告。
