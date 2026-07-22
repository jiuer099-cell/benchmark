# HG002 GRCh37 泛基因组 SV Benchmark 实施计划

日期：2026-07-17

状态：Phase 1 已完成；Phase 2/3 待真实资产与 HPC

设计基线：`docs/superpowers/specs/2026-07-16-hg002-grch37-sv-benchmark-design.md`

## 2026-07-19 Phase 1 最终验收

已完成：

- 配置、tool、rule registry、run manifest、pangenome manifest 和 metrics schemas；
- Snakemake 9 动态外部工具模块；
- benchmark 当前 job 实际运行外部 genotyper，而非读取用户提交 VCF；
- synthetic pangenome panel、稳定 allele ID、盲化 challenge panel；
- VCF raw gate、canonicalization 和 pangenome allele linkage；
- 标准 metrics schema/字典语义、100-point `pgbench_v1` 计分、逐 point TSV、no-ranking HTML；
- 14 个执行 rule 的 manifest/log/benchmark；
- 8 个 pre-score rule 的固定作业清单审计和完整 hash lineage；
- 最终 score package 原子封存、pre-score audit 重放和 11-job final provenance audit；
- `score.tsv`、`point_breakdown.tsv`、`metrics.long.tsv`、`metrics.json` 顶层无排名汇总；
- 安全增量重跑：旧 VCF/index 与旧 manifest 进入可恢复 archive，失败 run 不复用旧结果；
- Snakemake update-output 失败恢复已实测；archive 身份使用内容 hash 与 size，不依赖 Snakemake 会更新的 mtime；
- 外部插件完整代码包进入 DAG 指纹，manifest 可声明任意嵌套 raw VCF 路径；
- 执行前后冻结并复核输入、代码、配置、profile、reference 和上游 manifest 指纹；
- run/config/profile/truth/context 跨 manifest 一致性和 ratio 分子/分母语义均 fail-closed；
- HTML 报告与顶层汇总显式绑定 finalizer manifest、score package、score 和 metrics hash；
- Obsidian canonical design sync/check 脚本和 Snakemake targets。

当前本地回归测试为 177 个；synthetic DAG 为 14 个执行 rule 加 1 个聚合 target。
synthetic fixture 的 `PGBenchScore=78.06`、状态 `provisional`，环境锁缺失导致
traceability 为 4/5。该 fixture 只验证契约；真实 Truvari、Aardvark、vcfdist parser/wrapper
和 HG002 全基因组数据仍在 Phase 2/3。

## 实施原则

1. 每个活动 Snakemake rule 在首次落地时就必须登记 registry、manifest、log 和 benchmark，不能事后补 provenance。
2. `config/score_weights.yaml` 和 `config/metric_dictionary.yaml` 是评分事实源；Python 代码只负责加载、校验和计算。
3. Synthetic evaluator 只能产生 `evaluation_mode=synthetic_smoke`、`score_status=provisional`，不能伪装成正式 `pgbench_v1 valid`。
4. 外部工具必须由当前 benchmark job 实际执行，不接受预计算 VCF。
5. Obsidian 同步属于 post-score workflow，不影响工具分，但每个语义 commit 后必须更新实现状态与哈希。

## Phase 0：冻结确认状态

- 将设计状态改为“已确认，进入实现阶段”。
- 建立独立实现分支。
- 保存本实施计划。
- 验证：`git diff --check`。

## Phase 1：本地可运行 MVP

MVP 不访问真实 HG002 BAM，也不下载全量 graph。它使用 synthetic reference、population panel、truth
和 mock external genotyper，贯通配置校验、实际工具运行、VCF、三 evaluator 契约、100-point
smoke score、rule lineage 和 HTML report。

### P1.1 环境与骨架

文件：

- `environment.yaml`
- `pyproject.toml`
- `README.md`
- `Snakefile`
- `profiles/local/config.yaml`
- `profiles/slurm/config.yaml`
- `workflow/envs/core.yaml`

验证：

```bash
python -m pytest --version
snakemake --version
```

### P1.2 配置、Schema 与 Rule Registry

文件：

- `config/config.example.yaml`
- `config/config.schema.yaml`
- `config/truthsets.yaml`
- `config/stratifications.yaml`
- `config/score_weights.yaml`
- `config/metric_dictionary.yaml`
- `workflow/schemas/*.yaml`
- `workflow/rule-registry.yaml`
- `tests/unit/test_schemas.py`
- `tests/unit/test_rule_registry.py`

要求：

- JSON Schema draft 2020-12；
- official mode、外部工具模式契约和 no-ranking gate 可验证；
- 全部 normative rules 登记；未实现 rule 标记 `planned`。

### P1.3 Provenance 基础设施

文件：

- `workflow/scripts/pgbench_provenance.py`
- `workflow/scripts/pgbench_exec.py`
- `workflow/scripts/write_rule_manifest.py`
- `workflow/scripts/build_rule_lineage.py`
- `workflow/scripts/audit_provenance.py`
- `workflow/rules/provenance.smk`
- `tests/unit/test_provenance.py`

要求：

- 文件、目录、symlink 的确定性 SHA-256；
- stable manifest ID；
- 原子 JSON 写入；
- 成功/失败 attempt record；
- lineage 和核心缺失审计。

### P1.4 Synthetic 泛基因组资产

文件：

- `tests/fixtures/synthetic/`
- `workflow/scripts/build_pangenome_manifest.py`
- `workflow/scripts/build_challenge_panel.py`
- `workflow/rules/{preflight,reference,truth,preprocessing,pangenome,panel}.smk`

覆盖：

- in-panel 和 out-of-panel truth；
- 正/负 genotype candidates；
- stable `PANGENOME_ALLELE_ID`；
- truth leakage audit；
- blinded candidate IDs。

### P1.5 外部工具与 Mock Genotyper

文件：

- `workflow/modules/generic_external/Snakefile`
- `workflow/scripts/pgbench_exec.py`
- `workflow/scripts/validate_tool_output.py`
- `plugins/example_genotyper/`
- `tests/unit/test_tool_contract.py`

要求：

- benchmark 实际启动工具；
- caller-only/end-to-end 挂载隔离；
- `resolved_inputs.json` 固定资产路径与 hash；
- 当前 run 的 mtime、run ID、命令和输出 hash 可审计。

### P1.6 标准化与泛基因组等位链接

文件：

- `workflow/rules/normalization.smk`
- `workflow/scripts/normalize_sv_vcf.py`
- `workflow/scripts/link_pangenome_alleles.py`
- `tests/unit/test_normalization.py`
- `tests/unit/test_allele_linkage.py`

必须覆盖：

```text
in_panel_exact
in_panel_equivalent
out_of_panel
ambiguous
unresolved
invalid_allele_id
wrong_link
```

### P1.7 Evaluator 契约与确定性计分

文件：

- `workflow/scripts/parse_{truvari,aardvark,vcfdist}.py`
- `workflow/scripts/weight_evaluators.py`
- `workflow/scripts/pgbench_scoring.py`
- `workflow/scripts/score_tools.py`
- `workflow/rules/{evaluation,evaluator_weighting,scoring}.smk`
- `tests/unit/test_scoring.py`

当前检查点只完成 evaluator 输入契约、固定 fixture 和确定性融合计分；
三个真实 evaluator wrapper/parser 顺延到 Phase 2。

门禁：

- 28 + 24.5 + 17.5 + 15 + 10 + 5 = 100；
- evaluator 权重严格为 40:35:25；
- 无相位工具的 vcfdist 固定 0，不重分配；
- infrastructure failure 不等于工具零分；
- 输出不存在 rank/leaderboard 字段；
- smoke evaluator 的状态强制为 provisional。

### P1.8 Finalization、报告与 Obsidian

文件：

- `workflow/scripts/render_report.py`
- `workflow/scripts/check_obsidian_sync.py`
- `workflow/scripts/sync_obsidian_design.py`
- `workflow/rules/report.smk`
- `workflow/report/templates/index.html.j2`

输出：

```text
results/summary/metrics.long.tsv
results/summary/score.tsv
results/summary/point_breakdown.tsv
results/summary/<tool>/score-package.json
results/provenance/<tool>/rule-lineage.json
results/provenance/<tool>/provenance-audit.json
results/report/index.html
```

### P1.9 MVP 验收

```bash
python -m pytest -q
snakemake --lint
snakemake --dry-run --cores 1 \
  --configfile tests/fixtures/synthetic/config.yaml
snakemake --cores 2 \
  --configfile tests/fixtures/synthetic/config.yaml \
  results/report/index.html
```

完成条件：

- mock external genotyper 确实被当前 DAG 调用；
- 只产生一个 smoke/provisional PGBenchScore；
- 所有强制 F1 和 point 明细存在；
- 所有活动 rule 可追溯；
- lineage 到达 synthetic inputs、panel、truth、配置和代码；
- Obsidian `code_state`、commit 和实现哈希一致。

## Phase 2：真实 Evaluator 与首批工具

- Truvari、Aardvark、vcfdist 的固定版本 wrapper；
- Sniffles2、Sniffles2 force、kanpig、SVJedi-graph；
- 一个 mapping baseline 和一个正式 graph/pangenome 主路径；
- GRCh37 小区域真实数据测试；
- `pgbench_v1 valid` gate。

## Phase 3：HG002 全基因组与 HPC

- 冻结 hs37d5、GIAB v5.0q/v0.6、GIAB stratifications 和 1000G SV checksums；
- production pangenome manifest 和 graph recipes；
- SLURM/Conda/Apptainer 环境；
- 在目标路径运行 HG002 CLR 全基因组；
- 发布逐工具唯一 PGBenchScore、无排名报告与完整 provenance。

## 当前环境约束

- 本机没有真实 HG002 BAM；
- 初始环境缺少 Snakemake、samtools、bcftools、bgzip、tabix 和 Apptainer；
- 全基因组与真实 evaluator 验收必须在 Linux/HPC 完成；
- 本地 Phase 1 使用隔离环境与 synthetic fixtures，不对外声称正式生物学结果。
