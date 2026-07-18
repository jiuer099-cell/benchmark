# HG002 GRCh37 泛基因组结构变异检测与基因分型 Benchmark 设计规范

日期：2026-07-16

状态：已确认，进入实现阶段

工作流：Snakemake 9

默认样本：HG002

默认输入 BAM：
`/home/zhj/Experiment/SVDF-main/data/HG002.Sequel.10kb.pbmm2.hs37d5.whatshap.haplotag.RTG.10x.trio.bam`

## 0. 当前实现快照

截至 2026-07-19，仓库已完成 Phase 1 synthetic MVP，并通过本地可执行验证：

- `workflow/rule-registry.yaml` 登记 53 个规范 rule/template；
- synthetic DAG 实际执行 14 个带命令的 rule 和 1 个聚合 target；
- benchmark 在当前 Snakemake job 内实际启动
  `plugins/example_genotyper/run.py`，不接收预计算 VCF；
- 8 个 pre-score job 以及评分/封存链路的 manifest/log/benchmark 完整率为 100%，hash lineage 无缺失、无环；
- 已实现 pangenome panel、稳定 `PGSV_*` allele ID、盲化 challenge panel、VCF 标准化、
  allele linkage、标准 metrics 合约、固定 100-point 评分、原子 score seal、无排名长表/报告和 Obsidian canonical-note sync/check；
- 外部插件完整代码包进入 DAG 指纹且可声明自定义嵌套 raw VCF；失败重跑会恢复可信输出并按内容 hash/size 归档；
- provenance 会拒绝混合 run/config/profile/truth/context，ratio 会复核分子/分母语义，报告与汇总由 hash 绑定最终封存 manifest；
- 当前 synthetic evaluator 使用固定 fixture，只用于接口与计分回归测试。示例得分为
  `78.06 / 100`、状态 `provisional`；177 个本地测试通过。其中环境锁尚未冻结，因此 traceability 为 `4 / 5`。
  该数字不是 HG002 生物学结论，也不能代替正式 Truvari/Aardvark/vcfdist 运行；
- 真实 HG002 BAM preflight、冻结 truth/panel 下载、三个真实 evaluator、内置工具模块、
  Conda lock/container digest 与 Linux/HPC 全基因组验收仍属于后续阶段。

当前实现事实源包括：

```text
Snakefile
workflow/rules/*.smk
workflow/modules/generic_external/Snakefile
workflow/scripts/*.py
workflow/rule-registry.yaml
config/score_weights.yaml
plugins/example_genotyper/
tests/fixtures/synthetic/
```

## 1. 目标

构建一个以泛基因组候选等位、图或 panel 资产为核心，可扩展、可复现、可用于论文评测的结构变异 benchmark。框架需要：

1. 对 HG002 PacBio Sequel/CLR、hs37d5 BAM 执行统一前处理。
2. 在 hs37d5 backbone 上构建带稳定等位 ID 的 population pangenome panel，并为图工具和线性工具提供语义一致的候选等位。
3. 支持从 canonical reads 开始分流的完整端到端评测，以及共享 BAM/graph alignment 上的 caller-only 诊断评测。
4. 兼容以下分析范式：
   - mapping-based；
   - assembly-based；
   - graph/pangenome-based；
   - fragment/signature discovery and integration；
   - genotyping-only。
5. 内置若干结构变异检测与基因分型工具。
6. 允许使用者提供自己的可执行基因分型工具，由 benchmark 实际运行该工具分析 HG002，并生成待评分 VCF。
7. 统一校验和标准化工具输出，但保留复杂 SV 的原始表示。
8. 使用多个 SV 比较工具，按照固定 evaluator 权重和明确 point 明细生成唯一最终得分。
9. 必须报告 precision、recall 和 F1；F1 是最终得分的核心组成。
10. 覆盖检测、基因型、相位、类型、长度、区域、泛基因组上下文、断点、等位序列和资源消耗。
11. 所有 Snakemake rule 必须输出可追溯 manifest，并形成从输入到最终得分的 rule lineage。
12. 输出可审计的长表指标、单工具 score card、逐变异判定表和 HTML 报告，不生成排名。
13. 对每个适用工具在正式模式下定义且只定义一个 `PGBenchScore`，范围为 0–100。

## 2. 非目标

第一版不承诺：

- 在本地完整运行用户服务器上的全基因组 BAM；
- 自动下载并默认启用所有大型 pangenome graph 和 population panel；
- 把短读长专用工具强行用于当前 PacBio CLR 数据；
- 用单个数字掩盖评测器分歧或 truth set 局限；
- 接受使用者预先生成的 VCF 作为自研工具正式评测结果。
- 按最终得分对工具排序或输出 leaderboard/rank。

预计算 VCF 可以用于开发阶段的格式调试，但不得进入正式得分。正式结果必须由当前 benchmark run 调用工具后生成，并具有完整 provenance。

## 3. 已确认的核心原则

### 3.1 工具无关的评分层

评分层只依赖标准结果契约，不依赖内置工具白名单。内置工具和外部工具必须经过相同的：

```text
HG002 输入准备
  -> 工具执行
  -> 输出校验
  -> canonical VCF
  -> evaluator-specific VCF
  -> 多评测器
  -> 分层指标
  -> PGBenchScore
```

### 3.2 从 reads 开始的端到端公平性

用户提供的是已比对 BAM。框架从 BAM 一次性提取 canonical FASTQ，作为端到端赛道所有工具的共同起点。

同时保留 caller-only 模式，在共享 BAM 或共享 graph alignment 上诊断只负责 calling/genotyping 的部分：

- `caller_only_shared_alignment`
- `end_to_end_from_reads`

每次正式 run 只能指定一个 `official_score_mode`。默认值是
`end_to_end_from_reads`，它为每个适用工具 tuple 产生唯一正式 `PGBenchScore`；
caller-only 默认只输出诊断指标和资源消耗，不产生第二个正式总分。若研究问题明确只比较 caller，
可以把 `official_score_mode` 切换为 `caller_only_shared_alignment`，但每个工具仍只保留一个正式总分。

“唯一得分”的精确粒度是：

```text
(run_id, sample, tool, official_score_mode, primary_truth_profile)
```

每个适用 tuple 恰好一条 `PGBenchScore`。一次 run 可以测多个工具，因此 `score.tsv` 可以有多行工具分；
benchmark 不把这些工具分再次聚合成 run-level 总分，也不据此生成排名。

### 3.3 适用性不是零分

工具因输入技术、候选 panel 或图资产根本无法进入所选 official mode 时，状态为 `not_applicable`，不生成正式得分。工具一旦声明支持并进入正式 PGBenchScore run，固定 point 不再重新分配；缺少相位、out-of-panel discovery 或合法 allele linkage 等能力时，对应 point 为 0。这样可以区分“不能运行该 benchmark”和“能运行但能力不完整”。
SVTYPE 也按此原则处理：只有工具对 official profile 的全部计分 SVTYPE 都不支持时，整个 tuple 才是
`not_applicable`；若只支持其中一部分，仍进入正式评分，未支持的 truth-positive SVTYPE strata F1 为 0。

### 3.4 原始表示必须保留

复杂 SV 可能存在多个等价 VCF 表示。框架生成 canonical VCF，但不对所有工具全局执行 BND 分解、强制多等位拆分或不可逆的表示转换。每个评测器使用独立视图。

### 3.5 泛基因组优先

主评测对象是基于 population pangenome panel 或 graph 的 HG002 SV 检测与基因分型流程。默认 pangenome 不直接使用 HG002 truth 构图，而是使用独立 population SV 资源在 hs37d5 backbone 上构建：

```text
hs37d5 backbone
  + population SV alleles
  + stable PANGENOME_ALLELE_ID
  -> candidate panel VCF
  -> graph/index assets
```

图工具读取 graph/index；线性 genotyper 读取从同一 pangenome manifest 导出的 candidate VCF。mapping-based 和 assembly-based 工具作为泛基因组流程的可比基线存在，但所有结果最终映射到相同的 pangenome allele ledger 与 GRCh37 truth coordinate system。

## 4. 总体架构

主 Snakefile 负责：

1. 加载并校验配置；
2. 构建公共输入；
3. 动态加载已启用工具模块；
4. 汇总模块输出；
5. 执行标准化、评测器加权、评分和报告。

每个工具是独立 Snakemake 模块。简单外部工具可通过通用 runner 接入；多阶段外部工具可提供自己的模块。

```text
                          +----------------------+
                          | config + JSON schema |
                          +----------+-----------+
                                     |
                           +---------v---------+
                           | preflight / truth |
                           +---------+---------+
                                     |
                         +-----------v------------+
                         | pangenome manifest     |
                         | panel + graph + hashes |
                         +-----------+------------+
                                     |
              +----------------------+----------------------+
              |                                             |
      +-------v--------+                            +-------v--------+
      | shared BAM/GAF |                            | canonical reads |
      | caller-only    |                            | end-to-end      |
      +-------+--------+                            +-------+--------+
              |                                             |
              +----------------------+----------------------+
                                     |
                         +-----------v------------+
                         | dynamic tool modules   |
                         +-----------+------------+
                                     |
                         +-----------v------------+
                         | output validation      |
                         +-----------+------------+
                                     |
                    +----------------v----------------+
                    | canonical + evaluator views    |
                    +----------------+----------------+
                                     |
             +-----------------------+-----------------------+
             |                       |                       |
       +-----v-----+           +-----v-----+           +-----v-----+
       | Truvari  |           | Aardvark |           | vcfdist   |
       +-----+-----+           +-----+-----+           +-----+-----+
             +-----------------------+-----------------------+
                                     |
                          +----------v-----------+
                          | weight / stratify    |
                          +----------+-----------+
                                     |
                          +----------v-----------+
                          | one score / HTML     |
                          +----------------------+
```

Snakemake 9 动态模块用于按配置加载工具，并通过命名空间避免规则冲突。参考：
<https://snakemake.readthedocs.io/en/stable/snakefiles/modularization.html>

## 5. 建议目录结构

```text
.
├── Snakefile
├── README.md
├── config/
│   ├── config.example.yaml
│   ├── config.schema.yaml
│   ├── samples.example.tsv
│   ├── truthsets.yaml
│   ├── stratifications.yaml
│   ├── score_weights.yaml
│   └── metric_dictionary.yaml
├── profiles/
│   ├── local/
│   └── slurm/
├── workflow/
│   ├── rule-registry.yaml
│   ├── rules/
│   │   ├── common.smk
│   │   ├── preflight.smk
│   │   ├── reference.smk
│   │   ├── truth.smk
│   │   ├── pangenome.smk
│   │   ├── preprocessing.smk
│   │   ├── panel.smk
│   │   ├── normalization.smk
│   │   ├── evaluation.smk
│   │   ├── evaluator_weighting.smk
│   │   ├── scoring.smk
│   │   ├── provenance.smk
│   │   └── report.smk
│   ├── modules/
│   │   ├── sniffles2/
│   │   ├── cutesv/
│   │   ├── pbsv/
│   │   ├── svim/
│   │   ├── flye_dipcall/
│   │   ├── svim_asm/
│   │   ├── svjedi_graph/
│   │   ├── kanpig/
│   │   ├── vg/
│   │   ├── pangenie/
│   │   └── generic_external/
│   ├── schemas/
│   │   ├── tool.schema.yaml
│   │   ├── rule-registry.schema.yaml
│   │   ├── pangenome-manifest.schema.yaml
│   │   ├── run-manifest.schema.yaml
│   │   └── metrics.schema.yaml
│   ├── scripts/
│   │   ├── inspect_bam.py
│   │   ├── build_pangenome_manifest.py
│   │   ├── build_challenge_panel.py
│   │   ├── validate_tool_output.py
│   │   ├── normalize_sv_vcf.py
│   │   ├── link_pangenome_alleles.py
│   │   ├── pgbench_exec.py
│   │   ├── write_rule_manifest.py
│   │   ├── build_rule_lineage.py
│   │   ├── audit_provenance.py
│   │   ├── parse_truvari.py
│   │   ├── parse_aardvark.py
│   │   ├── parse_vcfdist.py
│   │   ├── weight_evaluators.py
│   │   ├── score_tools.py
│   │   ├── check_obsidian_sync.py
│   │   ├── sync_obsidian_design.py
│   │   └── render_report.py
│   └── envs/
├── plugins/
│   └── example_genotyper/
├── resources/
│   ├── references/
│   ├── truth/
│   ├── pangenome/
│   ├── panels/
│   └── stratifications/
├── results/
├── logs/
├── benchmarks/
├── tests/
│   ├── unit/
│   ├── integration/
│   ├── fixtures/
│   └── mock_tools/
└── docs/
```

`results/`、`logs/`、`benchmarks/` 和下载资源不进入版本控制。

## 6. 配置模型

主配置至少包含：

```yaml
project:
  id: hg002_grch37_sv_benchmark
  run_id: hg002_clr_baseline

sample:
  id: HG002
  sex: male
  technology: pacbio_clr
  bam: /home/zhj/Experiment/SVDF-main/data/HG002.Sequel.10kb.pbmm2.hs37d5.whatshap.haplotag.RTG.10x.trio.bam

reference:
  id: hs37d5
  fasta: resources/references/hs37d5.fa

pangenome:
  id: hs37d5_1kg_phase3_sv_v1
  backbone: hs37d5
  population_panel: 1000g_phase3_grch37_sv
  truth_leakage_policy: exclude_hg002_truth
  build_graph_assets: true

truth:
  primary: giab_hg002_grch37_v5_0q
  secondary:
    - giab_hg002_grch37_v0_6_legacy

execution:
  tracks:
    - caller_only_shared_alignment
    - end_to_end_from_reads
  official_score_mode: end_to_end_from_reads
  use_conda: true
  use_apptainer: false
  keep_going: true
  benchmark_repeats: 1

caller_only:
  alignment_kind: bam
  shared_bam: /home/zhj/Experiment/SVDF-main/data/HG002.Sequel.10kb.pbmm2.hs37d5.whatshap.haplotag.RTG.10x.trio.bam
  shared_graph_alignments: {}

score:
  profile: pgbench_v1
  final_score_name: PGBenchScore
  produce_ranking: false

obsidian:
  enabled: true
  vault_root: null
  vault_root_env: PGBENCH_OBSIDIAN_VAULT
  project_relpath: Research/hg002-grch37-pangenome-sv-benchmark
  canonical_note: Writing/HG002-GRCh37泛基因组SV-Benchmark设计规范.md
  binding_file: .claude/project-memory/hg002-grch37-pangenome-sv-benchmark.md

tools:
  enabled:
    - sniffles2
    - cutesv
    - pbsv
    - svim
    - sniffles2_force
    - kanpig
    - svjedi_graph

external_plugins: []
```

配置必须由 JSON Schema 校验。无效枚举、缺少 reference、工具 ID 冲突和不兼容输入在 DAG 构建前报错。
`execution.tracks` 表示本次允许执行和展示的模式；`official_score_mode` 必须是其中之一，并且只允许一个。
因此两个模式可以同时产生结果，但正式 `PGBenchScore` 只来自 `official_score_mode` 指定的模式。

Obsidian vault 按 `--obsidian-vault` CLI、`PGBENCH_OBSIDIAN_VAULT` 环境变量、本地 binding、
`obsidian.vault_root` 的顺序解析。公共示例配置不硬编码个人 vault 绝对路径；本机 binding 可以记录实际路径。
vault 不可用只会使 post-score sync/check 失败，不得改变工具分。

## 7. 统一前处理

### 7.1 BAM preflight

执行：

- 文件可读性与 `samtools quickcheck`；
- BAM/BAI 是否存在且一致；
- `@SQ` 的 `SN`、`LN` 和可用的 `M5` 与 hs37d5 比较；
- truth VCF contig 与 BAM/reference 兼容性检查；
- `@RG SM` 是否为 HG002；
- primary、secondary、supplementary、unmapped 记录比例；
- hard-clipping 比例；
- read length、MAPQ、coverage 概览；
- `HP`、`PS` 标签覆盖率；
- contig naming 风格；
- 输入文件大小和校验和。

以下情况硬失败：

- reference 与 BAM 长度或 M5 明显不匹配；
- BAM 损坏；
- truth 与 reference 版本不一致；
- HG002 样本名不匹配且未显式允许 override。

### 7.2 Canonical reads

从 BAM 一次性生成：

```text
prepared/HG002/all.fastq.gz
prepared/HG002/hp1.fastq.gz
prepared/HG002/hp2.fastq.gz
prepared/HG002/unphased.fastq.gz
```

所有端到端工具使用 `all.fastq.gz` 或明确声明的 HP 分组文件。提取规则记录 read 数、base 数、hard-clipped 风险和 BAM 到 FASTQ 的可逆性限制。

### 7.3 共享比对

caller-only 默认使用用户提供的 pbmm2 hs37d5 BAM。

框架还提供可选的共享重比对：

```text
canonical FASTQ -> minimap2 -x map-pb -> sorted/indexed BAM
```

共享重比对用于分离“caller 差异”和“原始 BAM 特性”。工具自带 mapper 的完整流程只进入端到端赛道。

### 7.4 Caller-only 与 End-to-end 的定义

#### Caller-only

Caller-only 用来回答：

> 在上游 reads、reference/pangenome panel 和比对结果完全相同的情况下，哪个 caller 或 genotyper 本身表现如何？

正式 mode 名为 `caller_only_shared_alignment`。固定输入是 benchmark 已准备好的共享 BAM、GAF/GAM
或其他在 tool manifest 中声明的标准 alignment，以及同一 candidate panel。`alignment_kind` 必须显式记录为
`bam|cram|gaf|gam|other`，不能让工具自行选择未登记的上游比对。计时从 caller/genotyper 进程启动开始，
到原始 VCF 写出结束。它不计入：

- BAM/graph alignment 生成；
- 公共 pangenome panel 下载与构建；
- 公共 truth 和 stratification 准备。

它适合定位误差来自 mapping 还是 calling，也适合开发者快速调参和调试自研 genotyper。因为它没有覆盖完整用户流程，默认只生成诊断指标，不生成正式总分。

共享只在输入契约兼容的 `alignment_group_id` 内成立，例如线性 BAM 工具共用一个组，读取同一 vg graph
的 GAF/GAM 工具共用另一个组。`alignment_group_id` 写入 score tuple 的 provenance；不同组之间不声称
上游完全等价，也不做排名。

#### End-to-end

End-to-end 用来回答：

> 使用者只有统一 HG002 reads 和同一 pangenome manifest 时，这个完整工具流程最终能产生多准确、耗费多少资源的 VCF？

固定起点是 canonical FASTQ、hs37d5 backbone 和不包含 HG002 truth 标签的 pangenome manifest。计时与资源统计包括：

- 从公共 pangenome 语义资产转换出的工具专用 graph/index；
- reads mapping 或 graph alignment；
- signature/fragment discovery；
- assembly，若该流程需要；
- SV calling 或 genotyping；
- 工具自身后处理；
- 原始 VCF 生成。

benchmark 公共 truth 下载、公共报告生成和 evaluator 执行不计入被测工具资源分。默认正式 `PGBenchScore` 来自 end-to-end，因为它最接近真实使用成本和最终效果。

资源边界固定如下：

- 不计入工具分：hs37d5、population source、canonical panel VCF、allele FASTA、stable allele ledger、
  pangenome manifest、backend-neutral graph recipe、truth 和 stratifications 的一次性公共准备；
- 若 benchmark 发布了所有工具都可读取的 canonical graph，它作为正式输入资产，不重复计建图成本；
- 计入 end-to-end：工具为了读取上述公共语义资产而创建的专用 graph、index、database、cache，
  以及 mapping/assembly/calling/genotyping/postprocess；
- 计入 caller-only：caller/genotyper 与其必要 postprocess；共享 alignment 和共享 candidate panel 不计；
- evaluator、计分、report 和 Obsidian 同步不计入任何工具资源 point。

正式资源测量使用干净的工具工作目录。工具专用 index 不得通过预热 cache 规避计费；公共输入资产可以复用，
但其路径、hash 和 cache policy 必须写入 manifest。

如果同时运行两个模式，报告将 caller-only 放在“诊断”页，把 end-to-end 放在“正式得分”页，不对两个模式互相排序。

### 7.5 泛基因组资产

默认 pangenome asset 使用 hs37d5 作为 backbone，并从独立 population SV 资源构建候选等位。第一版默认 population source 是 1000 Genomes Phase 3 GRCh37 integrated SV map；HG002 v5.0q truth 只用于评测和 truth leakage 排除，不能作为被测图或 panel 的训练标签。

统一 manifest 至少包含：

```yaml
pangenome_id: hs37d5_1kg_phase3_sv_v1
backbone_reference: hs37d5
population_sources:
  - 1000g_phase3_grch37_sv
allele_id_namespace: PGSV
truth_samples_excluded:
  - HG002
panel_vcf_sha256: ...
graph_build_recipe_sha256: ...
```

预期公共产物：

```text
resources/pangenome/<id>/panel.vcf.gz
resources/pangenome/<id>/alleles.fasta.gz
resources/pangenome/<id>/manifest.yaml
resources/pangenome/<id>/graphs/svjedi/
resources/pangenome/<id>/graphs/vg/
resources/pangenome/<id>/indexes/
```

每个 pangenome allele 都有稳定 `PANGENOME_ALLELE_ID`，用于连接：

- population panel；
- graph path/bubble；
- 工具候选记录；
- HG002 输出 VCF；
- truth match；
- evaluator 判定；
- 最终 point 明细。

## 8. 工具范式与内置模块

| 范式 | 标准入口 | 第一版状态 | 代表模块 |
|---|---|---|---|
| Graph/pangenome | FASTQ + pangenome graph/panel | 主路径 | vg、SVJedi-graph、外部图 genotyper |
| Mapping-based | BAM 或 canonical FASTQ + pangenome candidate panel | 基线与兼容路径 | Sniffles2、cuteSV、pbsv、SVIM |
| Assembly-based | canonical/HP FASTQ | 可启用 | Flye + dipcall、SVIM-asm |
| Discovery + integration | BAM/FASTQ -> signatures -> integrated VCF | 通用接口与示例 | generic two-stage adapter |
| Genotyping-only | BAM/FASTQ + blinded candidate VCF | 核心可运行 | Sniffles2 force、kanpig、SVJedi-graph |
| Short-read pangenome genotyping | short FASTQ + panel | 默认不适用 | PanGenie |

PanGenie 官方要求短读长输入，不能把当前 CLR reads 当作等价输入：
<https://github.com/eblerjana/pangenie>

第一版保留 PanGenie 模块和完整资产检查，但只有在用户另行配置短读长及 panel 时才进入 DAG。

## 9. 工具适配器契约

### 9.1 Manifest

每个工具必须声明：

```yaml
interface_version: 1
id: kanpig
source: builtin
paradigm: genotyping_only

tasks:
  - genotyping

inputs:
  reads: bam
  reference: required
  candidates: required
  pangenome_manifest: required
  pangenome_panel: required
  graph_assets:
    required: false
    accepted_formats:
      - vg
      - xg_gbwt
      - svjedi
  shared_alignment:
    required: false
    accepted_formats:
      - bam
      - gaf
      - gam

capabilities:
  technology:
    - pacbio_clr
    - pacbio_hifi
    - ont
  svtypes:
    - DEL
    - INS
    - DUP
    - INV
  phasing: false
  novel_discovery: false
  allele_linkage_output: none
  allele_namespace:
    - PGSV

execution:
  module: workflow/modules/kanpig/Snakefile
  environment: workflow/envs/kanpig.yaml

supported_modes:
  caller_only_shared_alignment:
    required_inputs:
      - shared_alignment
      - reference
      - pangenome_manifest
      - pangenome_panel
    forbidden_inputs:
      - canonical_fastq
      - original_input_bam
    billable_stages:
      - genotype
      - postprocess
  end_to_end_from_reads:
    required_inputs:
      - canonical_fastq
      - reference
      - pangenome_manifest
      - pangenome_panel
    forbidden_inputs:
      - original_input_bam
      - shared_alignment
    billable_stages:
      - prepare_assets
      - index
      - map
      - genotype
      - postprocess

outputs:
  vcf: raw/calls.vcf.gz
  candidate_output_contract: all_sites
  absence_semantics: no_call
```

`genotyping_only` 工具必须使用 `candidate_output_contract: all_sites`：每个 challenge candidate
均输出一个 GT，缺行或 `./.` 是 no-call，不能暗中解释为 `0/0`。discovery/calling 工具可以声明
`variant_sites`；此时未出现的 candidate 只在 `absence_semantics: hom_ref` 被 schema 允许且提前冻结时
解释为隐式 `0/0`，否则仍是 no-call。该语义不能在看到结果后修改。

### 9.2 统一输出

每个工具运行必须生成：

```text
results/HG002/<track>/<tool>/raw/calls.vcf.gz
results/HG002/<track>/<tool>/raw/calls.vcf.gz.tbi
results/HG002/<track>/<tool>/meta/run_manifest.json
results/HG002/<track>/<tool>/meta/status.json
logs/<tool>/HG002.log
benchmarks/<tool>/HG002.jsonl
```

`run_manifest.json` 包含：

- benchmark run ID；
- 工具 ID、版本与来源；
- 完整命令；
- 参数；
- container digest 或 Conda environment export；
- 输入路径、大小和校验和；
- reference ID；
- pangenome manifest、panel、graph/index asset ID 与校验和；
- allele namespace 与 candidate panel ID；
- official mode、alignment kind 与计费 stage；
- threads 和资源请求；
- 开始、结束时间；
- 退出码；
- 输出校验和。

### 9.3 所有 Rule 的可追溯契约

不只工具 rule，reference、pangenome、preprocessing、normalization、evaluation、evaluator weighting、scoring 和 report 的每个 Snakemake rule 都必须生成一个 rule manifest：

```text
results/provenance/rules/<rule_name>/<job_key>.json
logs/rules/<rule_name>/<job_key>.log
benchmarks/rules/<rule_name>/<job_key>.jsonl
```

每个 manifest 必须包含：

```text
manifest_schema_version
manifest_id
attempt_id
run_id
rule_name
job_key
wildcards
module_or_tool_id
snakefile_path
rule_source_sha256
script_or_wrapper_path
script_or_wrapper_sha256
command
params
threads
requested_resources
input_paths
input_sha256
input_size
input_mtime
output_paths
output_sha256
config_snapshot_sha256
pangenome_manifest_sha256
reference_sha256
truth_profile
conda_lock_sha256
container_uri
container_digest
git_head
git_dirty
git_diff_sha256
snakemake_version
execution_profile
hardware_fingerprint_sha256
random_seed
upstream_manifest_ids
started_at
finished_at
exit_code
status
```

规则要求：

1. rule 的输入、输出和日志必须显式命名，不允许关键文件只存在于未声明的临时目录。
2. 负责计时的规则优先使用 `shell`、`script` 或 `wrapper`，避免无法可靠 benchmark 的复杂 `run` block。
3. 任何脚本或 wrapper 的内容哈希发生变化，manifest 必须变化。
4. 所有 output checksum 在 rule 成功后计算。
5. evaluator weighting 和 scoring manifest 必须记录所使用 evaluator 版本、权重和每个 point 的来源指标。
6. `input_sha256` 和 `output_sha256` 是 path 到 hash 的映射，不是无顺序列表。目录使用按相对路径排序后计算的
   Merkle-style tree hash；symlink 同时记录 link target 和 target content hash。
7. manifest 自身、临时文件和日志不进入自身 `output_sha256`，避免自哈希递归。
8. 不适用字段必须写显式 `null` 和 `not_applicable_reason`，不得用空字符串混淆“未知”。
9. `manifest_id` 由 schema version、run ID、rule name、job key、attempt ID 和 rule source hash 确定，
   同一 attempt 稳定且全局唯一。
10. 最终生成：

```text
results/provenance/rule-lineage.json
results/provenance/rule-lineage.tsv
results/provenance/provenance-audit.json
```

`rule-lineage` 从最终 `PGBenchScore` 反向链接到 evaluator metrics、canonical VCF、tool raw VCF、pangenome assets、reads/BAM 和原始配置。任何缺失链路都会扣除 traceability points；核心链路缺失则最终分状态为 `invalid`。

每个 job 通过统一 `pgbench-exec` runner 启动。runner 先在 `log:` 路径原子写入
`logs/provenance_attempts/<rule_name>/<job_key>.<attempt_id>.json`，并用退出 trap 更新状态，因此即使
Snakemake 删除失败 output，失败 attempt 记录仍保留。成功时再以临时文件加 `os.replace` 写正式 manifest。
Python `script:`、shell、wrapper 和外部插件都必须经过同一 runner，不能绕过。

### 9.4 Normative Rule Registry

仓库必须维护 `workflow/rule-registry.yaml`。下表是第一版规范化 rule ID；实现可以拆分内部 helper，
但任何进入 DAG 的可执行 rule 都必须登记，且不能只有代码而没有 registry 条目。

| Rule ID 或模板 | 主要输入 | 主要输出 | 直接上游/用途 |
|---|---|---|---|
| `validate_config` | config、schemas | validated config | DAG 前配置校验 |
| `snapshot_run_context` | config、Git、profile、hardware | run context JSON | 冻结 score profile、seed、代码和硬件 |
| `inspect_input_bam` | HG002 BAM/BAI、reference dict | BAM QC JSON | 输入完整性与样本检查 |
| `prepare_reference` | hs37d5 source | FASTA/FAI/DICT | 统一 backbone |
| `prepare_primary_truth` | GIAB v5.0q source | truth VCF/BED | 正式 truth |
| `prepare_legacy_truth` | GIAB v0.6 source | legacy VCF/BED | 仅历史报告 |
| `prepare_stratifications` | GIAB GRCh37 BEDs | normalized BEDs | 固定 context universe |
| `extract_canonical_reads` | input BAM | all FASTQ | end-to-end 共同起点 |
| `split_haplotype_reads` | input BAM | HP1/HP2/unphased FASTQ | 可选 haplotype 输入 |
| `prepare_population_sv_source` | population SV source | source VCF | 独立 population alleles |
| `normalize_population_panel` | source VCF、hs37d5 | normalized panel VCF | 统一表示与 contig |
| `assign_pangenome_allele_ids` | normalized panel | ID panel、allele ledger | 分配稳定 ID |
| `build_pangenome_manifest` | ID panel、alleles、recipes | manifest YAML | 冻结 pangenome 事实源 |
| `audit_truth_leakage` | manifest、source metadata、HG002 truth | audit JSON | 禁止 truth 注入 |
| `build_canonical_graph_recipe` | manifest、panel、alleles | backend-neutral recipe | 图构建共同语义 |
| `build_shared_alignment__{kind}` | reads、registered graph/reference | BAM/GAF/GAM | caller-only 共享 alignment |
| `match_hidden_panel_truth` | panel、primary truth、BED | hidden label ledger | 构建正/负 genotype labels |
| `build_blinded_challenge_panel` | panel、hidden ledger、seed | blinded VCF | 工具可见候选集 |
| `audit_challenge_panel` | blinded VCF、hidden ledger | leakage/QC JSON | 检查负位点与盲化 |
| `tool__{tool}__preflight` | tool manifest、available assets | compatibility JSON | 工具适用性 |
| `tool__{tool}__execute` | validated manifest、mode-authorized assets | raw VCF、attempt、resolved inputs | 通用外部 runner 的具体执行 rule |
| `tool__{tool}__prepare_assets` | pangenome manifest | tool asset recipe | 工具特有资产准备 |
| `tool__{tool}__index` | reference/graph/panel | tool index | end-to-end 计费 stage |
| `tool__{tool}__map` | canonical reads、tool index | tool alignment | end-to-end mapping |
| `tool__{tool}__discover` | reads/alignment | signatures/fragments | discovery stage |
| `tool__{tool}__assemble` | canonical/HP reads | assembly | assembly stage |
| `tool__{tool}__integrate` | fragments/signatures | integrated candidates | integration stage |
| `tool__{tool}__call` | declared tool inputs | raw calls | discovery/calling |
| `tool__{tool}__genotype` | alignment、blinded panel | raw genotypes | genotyping-only |
| `tool__{tool}__postprocess` | native output | raw VCF | 工具自身必要转换 |
| `tool__{tool}__validate_raw` | raw VCF | validation JSON | 输出契约门禁 |
| `canonicalize_vcf` | raw VCF、reference | canonical VCF | 公共可逆标准化 |
| `link_pangenome_alleles` | canonical VCF、allele ledger | linked VCF/TSV | 稳定 allele linkage |
| `build_rule_lineage` | pre-score manifests | lineage JSON/TSV | 评分前哈希 lineage |
| `make_truvari_view` | linked VCF | Truvari view | evaluator 专用视图 |
| `make_aardvark_view` | linked VCF | Aardvark view | evaluator 专用视图 |
| `make_vcfdist_view` | linked VCF | vcfdist view | evaluator 专用视图 |
| `evaluate_truvari` | Truvari view、truth | raw evaluator artifacts | bench + refine |
| `evaluate_aardvark` | Aardvark view、truth | raw evaluator artifacts | haplotype evaluation |
| `evaluate_vcfdist` | vcfdist view、truth | raw evaluator artifacts | representation/phase |
| `parse_truvari_metrics` | Truvari artifacts | standard metrics/ledger | 固定 parser schema |
| `parse_aardvark_metrics` | Aardvark artifacts | standard metrics/ledger | 固定 parser schema |
| `parse_vcfdist_metrics` | vcfdist artifacts | standard metrics/ledger | 固定 parser schema |
| `stratify_evaluator_metrics` | evaluator ledgers、BED/allele bins | stratum metrics | type/length/context/pangenome |
| `fuse_evaluator_metrics` | three standard metrics | weighted metrics | 固定 40:35:25 |
| `collect_tool_resources` | stage benchmark records | resource metrics | official-mode 边界求和 |
| `audit_score_inputs` | score ancestor manifests | pre-score audit | 计算 traceability component |
| `compute_pgbench_score` | weighted metrics、resource、audit | score/point breakdown | 唯一数值总分 |
| `finalize_score_provenance` | score、score manifest、ancestor lineage | final lineage/audit | 封装正式 score package |
| `render_report` | finalized score package | HTML/TSV/JSON | 展示，不反向影响得分 |
| `render_report_index` | tool HTML cards | no-ranking index | 多工具导航，不生成排名 |
| `check_obsidian_sync` | repo/KB hashes | sync check JSON | 只读、独立于 BAM |
| `sync_obsidian_design` | repo design、binding、旧 note hash | sync record（含新 note hash） | 原子更新 canonical note，不进入得分 |

工具只实例化其 manifest 声明的 stage，但生成后的
`results/provenance/resolved-rule-registry.{json,tsv}` 必须列出每个 concrete rule/job、active/inactive
原因、输入输出契约、上游 rule、环境、资源 class、source hash 和 manifest path。自研工具的自定义
Snakefile 若增加 rule，也必须在插件 `rule_registry` 中登记，并在加载时合并校验。

### 9.5 Score DAG 边界

`pre-score audited rules` 定义为 `audit_score_inputs` 的所有传递祖先，不包含
`audit_score_inputs` 自身。该集合不含 `compute_pgbench_score`、finalization、report 或 sync，
因此 pre-score audit 不会审计自己，也不会形成 scoring→report→scoring 循环。

`compute_pgbench_score` 完成后，`finalize_score_provenance` 再校验 `audit_score_inputs` 和 scoring rule
自身的 manifest，并封装最终 score package。`audit_score_inputs` 与 `compute_pgbench_score` 的 manifest
属于 formal validity gate，但不反向进入已计算的 2-point completeness 分母。
`render_report`、`check_obsidian_sync` 和 `sync_obsidian_design` 属于 post-score
rules：仍然必须有完整 manifest/log/benchmark，但不进入当前工具的数值分，避免 HTML 或外部 vault
故障改变工具准确性得分。若 `finalize_score_provenance` 失败，数值只能标记 `invalid`，不能发布为正式分。

## 10. 外部自研基因分型工具

### 10.1 正式评测入口

使用者提供可执行工具包，而不是 VCF：

```text
plugins/my_genotyper/
├── tool.yaml
├── run.sh
├── envs/environment.yaml
├── scripts/normalize.py   # 原生输出不是 VCF 时可选
└── tests/smoke/
```

以下二选一：

1. `run.sh`：适合单步或少量步骤工具；
2. `workflow/Snakefile`：适合多阶段工具。

`normalize.py` 仅在工具原生输出不是合规 VCF 时使用。它仍然必须在同一个 benchmark run 内执行。

### 10.2 通用 runner 环境

benchmark 将以下变量传给外部工具：

```text
PGBENCH_RUN_ID
PGBENCH_SAMPLE_ID
PGBENCH_INPUT_FASTQ
PGBENCH_HP1_FASTQ
PGBENCH_HP2_FASTQ
PGBENCH_REFERENCE_FASTA
PGBENCH_CANDIDATE_VCF
PGBENCH_PANGENOME_MANIFEST
PGBENCH_PANEL_VCF
PGBENCH_ALLELES_FASTA
PGBENCH_GRAPH_DIR
PGBENCH_INDEX_DIR
PGBENCH_ALLELE_NAMESPACE
PGBENCH_SHARED_ALIGNMENT
PGBENCH_ALIGNMENT_KIND
PGBENCH_RESOLVED_INPUTS
PGBENCH_OUTPUT_DIR
PGBENCH_OUTPUT_VCF
PGBENCH_THREADS
PGBENCH_MEMORY_MB
```

用户的 `run.sh` 负责把统一变量翻译为自研工具 CLI。不得读取评分 truth 路径。
未声明需要的 graph/index/alignment 变量为未设置，而不是指向其他工具资产。所有输入默认只读挂载；
只有 `PGBENCH_OUTPUT_DIR` 和显式 tool workdir 可写。接口前缀固定为 `PGBENCH_`。

runner 在每个工具 job 前生成 `resolved_inputs.json`，记录本模式允许的 asset ID、格式、只读挂载路径、
content hash 和 billable stage。caller-only 不暴露 canonical FASTQ 或原始 BAM 语义入口，只暴露
已登记的 shared alignment；end-to-end 不暴露原始 BAM 或 shared alignment。环境变量也只为允许输入设置，
不把整个公共资源目录挂载给外部工具，禁止插件通过目录扫描自行选择其他 graph/index。
若 caller-only 工具确实只为补取 read sequence 需要 FASTQ，必须声明
`read_sequences_auxiliary: true`，FASTQ 以无 reference/index 的只读挂载提供，runner 审计不得产生新 alignment；
否则状态为 `mode_contract_violation`，不生成正式分。

### 10.3 禁止预计算结果替代运行

正式评测要求：

- VCF 位于本次 Snakemake job 分配的输出目录；
- mtime、run ID 和输出校验和写入 manifest；
- 生成 VCF 的工具命令由当前 job 执行；
- job 退出前执行 VCF 校验；
- 缺少本次运行 provenance 时状态为 `provenance_failed`；
- `source=builtin|external` 只用于审计，不影响评分权重。

## 11. HG002 基因分型 Challenge Panel

### 11.1 目的

genotyping-only 工具必须同时面对：

- HG002 中存在的 SV：`0/1`、`1/1` 或适当的单倍体基因型；
- HG002 中不存在的候选 SV：`0/0`；
- 缺失或无法可靠判定的位点。

只用 truth-positive 位点会高估基因分型性能，因此正式 `panel-genotyping` 赛道必须包含负位点。

### 11.2 正位点

正位点不能直接把 HG002 truth records 注入 candidate panel。正式流程是：

1. 从独立 population pangenome manifest 读取候选等位；
2. 保留完全处于 HG002 benchmark BED 的候选；
3. 使用高灵敏、多表示 matcher 与 GIAB v5.0q truth 匹配；
4. 匹配成功的 pangenome candidate 获得隐藏的 HG002 `0/1`、`1/1` 或单倍体 truth label；
5. 工具只看到候选等位，不看到 HG002 truth genotype。

HG002 truth 中存在但 population pangenome 不包含的 SV 不会被偷偷加入 genotyping panel；它们计入 `pangenome allele coverage` 和 end-to-end novel discovery 指标。

### 11.3 负位点

默认候选来源为统一 pangenome panel，其第一版 population source 是 1000 Genomes Phase 3 GRCh37 integrated SV map：
<https://www.internationalgenome.org/phase-3-structural-variant-dataset/>

负位点必须：

1. 完全处于 v5.0q benchmark BED；
2. 在所有适用的高灵敏匹配视图中均不与 HG002 truth 等价；
3. 不位于 truth 复杂 cluster 的安全排除窗口；
4. 不属于 evaluator-discordant、无法解析的 BND 或不完整等位；
5. 具有可供被测 genotyper 使用的合法候选表示；
6. 按 SVTYPE、长度和 genome context 与正位点分层抽样，避免简单位置占主导。

为降低把等价表示误标为 `0/0` 的风险，若任一高灵敏 matcher 找到合理等价关系，该候选从负集排除，而不是通过多数票保留。

### 11.4 盲化

工具接收的 challenge VCF：

- 删除 HG002 truth GT；
- 删除 truth/negative 标签；
- 删除来源数据集标识；
- 使用固定随机种子生成中性 candidate ID；
- 保留工具完成基因分型所需的 REF、ALT、SVTYPE、END、SVLEN 和等位序列；
- truth labels 存放在独立的 scoring artifact 中，不传给工具 job。

盲化用于减少直接读取 truth 标签的风险，但不是对恶意执行环境的密码学隔离。正式对比推荐使用只挂载 reads、reference、candidate panel 和输出目录的容器执行方式。

## 12. Truth Sets 与参考资源

### 12.1 主 truth

主 truth：

```text
HG002_GRCh37_v5.0q_stvar.vcf.gz
HG002_GRCh37_v5.0q_stvar.vcf.gz.tbi
HG002_GRCh37_v5.0q_stvar.benchmark.bed
```

来源：
<https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/v5.0q/>

v5.0q 为 assembly-based draft benchmark。其 README 建议使用能够处理复杂表示的 Truvari refine、hap-eval、vcfdist 或 Aardvark：
<https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/v5.0q/NIST_HG002_v5.0q_variant-benchmarksets_README.md>

报告必须明确显示 `truth_status=draft`，并保留评测器分歧与人工复核候选。

### 12.2 历史 truth

GIAB v0.6 只作为 legacy profile：

```text
HG002_SVs_Tier1_v0.6.vcf.gz
HG002_SVs_Tier1_v0.6.bed
```

它主要覆盖 isolated、sequence-resolved DEL/INS。必须仅使用 PASS 且位于 Tier1 BED 的变异。v0.6 指标不得与 v5.0q 合并计分，只用于复现历史论文。

参考论文：
<https://www.nature.com/articles/s41587-020-0538-8>

### 12.3 Reference

默认 reference 为 GIAB 分发的 hs37d5：
<https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/references/GRCh37/>

reference FASTA、FAI、DICT 和所有工具索引均由统一 reference ID 管理。禁止在同一次 run 中混用 UCSC `chr1` 风格 hg19 和 hs37d5 `1` 风格而不进行显式映射。

### 12.4 Genome stratifications

默认使用 GIAB v3.6 GRCh37 stratifications。至少包含：

- all difficult / easy；
- segmental duplications；
- tandem repeats，按 repeat 长度分层；
- homopolymers；
- low complexity；
- low mappability；
- GC extremes；
- coding/exonic regions；
- medically relevant regions；
- autosomes、X、Y 和 PAR/non-PAR。

## 13. VCF 校验与标准化

### 13.1 Canonical VCF

公共标准化执行：

- 单样本提取；
- 样本名改为 HG002；
- contig 映射与 header 修复；
- 按 reference 顺序排序；
- bgzip 和 tabix；
- VCF 语法校验；
- 对可可靠推断的记录补充 `SVTYPE`、`SVLEN`、`END`；
- 为每条记录增加稳定 `CANON_ID`；
- 对 panel-derived 记录保留 `PANGENOME_ALLELE_ID`；
- 保留 `ORIG_ID`、`ORIG_SVTYPE`、`ORIG_ALT`；
- 记录 symbolic allele、missing sequence、multi-allelic 和 unpaired BND 状态。

不会默认：

- 删除所有小于 50 bp 的邻近变异；
- 全局拆分复杂多等位记录；
- 全局把所有 SV 分解为 BND；
- 猜测缺失的插入序列；
- 把 DUP 无条件重写为 INS。

### 13.2 Evaluator-specific views

每个 evaluator 声明自己的输入转换：

- Truvari raw bench view；
- Truvari refine view；
- Aardvark haplotype view；
- vcfdist phased、eligible-size view；
- 可选 hap-eval/Sansa view。

所有转换生成 manifest，确保结果可追溯到 canonical 和 original record。

### 13.3 泛基因组等位链接

`link_pangenome_alleles` rule 在 canonical VCF 与 evaluator view 之间执行。它不能只信任工具写入的
`PANGENOME_ALLELE_ID`，而是按以下固定顺序验证或推断：

1. 校验工具给出的 stable allele ID 是否存在于当前 manifest，且 REF/ALT 或 graph path 一致；
2. 尝试 canonical coordinate、SVTYPE、SVLEN 和完整等位序列的 exact match；
3. 尝试局部单倍型重构后的 representation-equivalent match；
4. 若没有 panel candidate 与调用等价，标记为 `out_of_panel`；
5. 若多个候选无法唯一消歧，保留 `ambiguous`，不选择对工具最有利的 ID；
6. 若工具声称的 ID 与 benchmark 唯一推断出的 panel allele 冲突，标记 `wrong_link` 并同时保存推断 ID；
7. 若序列或配对 BND 信息不足，标记 `unresolved`。

每条调用必须得到以下一个状态：

```text
in_panel_exact
in_panel_equivalent
out_of_panel
ambiguous
unresolved
invalid_allele_id
wrong_link
```

输出：

```text
results/HG002/<track>/<tool>/canonical/allele_links.tsv
results/HG002/<track>/<tool>/canonical/linked.calls.vcf.gz
results/HG002/<track>/<tool>/canonical/linked.calls.vcf.gz.tbi
```

`allele_links.tsv` 至少包含 `CANON_ID`、工具原始 allele ID、最终
`PANGENOME_ALLELE_ID`、link status、候选数、匹配方法、相似度、panel AF bin、graph-complexity bin
和 hidden truth consistency。该表是 pangenome point、逐变异复核和 provenance lineage 的共同事实源。

## 14. 多评测器设计

### 14.1 默认评测器

1. **Truvari bench + refine**
   - event precision、recall、F1；
   - genotype concordance 和 GT matrix；
   - breakpoint、size、sequence similarity；
   - complex representation refinement。

2. **Aardvark**
   - haplotype sequence context；
   - variant-type agnostic comparison；
   - exact genotype assessment；
   - inexact representation partial credit。

3. **vcfdist**
   - phased、sequence-resolved eligible variants；
   - complex variant representation；
   - switch 和 flip error；
   - 按固定规则计算 evaluator score。

### 14.2 可选评测器

- hap-eval：pre-release，作为敏感性分析；
- Sansa：breakpoint-centric sensitivity check；
- Wittyer：已归档，仅保留 legacy compatibility，不进入默认评分。

### 14.3 指标长表

每个 evaluator 输出统一长表：

```text
sample
track
tool
tool_source
truth_profile
evaluator
metric
value
numerator
denominator
svtype
length_bin
region
context
eligible
status
```

## 15. 多评测器加权

三个默认 evaluator 都产生一个归一化到 `[0, 100]` 的 evaluator score。固定权重为：

```yaml
truvari: 0.40
aardvark: 0.35
vcfdist: 0.25
```

评测器准确性总分：

```text
EvaluatorAccuracy =
  0.40 * TruvariScore
  + 0.35 * AardvarkScore
  + 0.25 * VcfdistScore
```

`EvaluatorAccuracy` 范围为 0–100，最终贡献 70 个 point。

每个 evaluator 的原始 precision、recall、F1、GT、phase 和 representation metrics 始终保留。权重写入版本化 score profile，不能由单次运行临时修改。

状态处理：

- evaluator 因安装、运行环境或 benchmark 自身错误失败：本次正式得分为 `invalid`，不能把它当作工具零分；
- evaluator 因被测工具缺少必须输出而不可运行，例如完全没有相位信息导致 vcfdist 不适用：对应 evaluator
  score 为 0，并记录 `completed_with_tool_ineligible_zero`；其所有计分字段仍写数值 0、固定 denominator、
  `tool_ineligible_reason` 和 capability evidence，不能只删除该 evaluator 行；
- 可选 evaluator 只用于敏感性分析，不改变固定 40:35:25 权重；
- evaluator disagreement 单独报告，不通过删除不利 evaluator 来提高最终得分。

所有 event F1 使用：

```text
precision = TP / (TP + FP)
recall = TP / (TP + FN)
F1 = 2 * precision * recall / (precision + recall)
```

当 precision 与 recall 分母均为 0 时，该 stratum 标记为 `undefined`；预定义 truth-positive stratum 中工具无有效调用时，F1 记为 0。

逐变异 ledger 包含：

```text
CANON_ID
original_id
tool
truvari_status
aardvark_status
vcfdist_status
consensus_status
review_reason
```

至少两个适用 evaluator 判定为匹配时标记 `consensus_TP`。其他情况保留为 discordant，不强行变成 FP 或 TP。

## 16. 评测维度

### 16.1 检测

- TP、FP、FN；
- precision、recall、F1；
- PR curve；
- call count 和 callable fraction；
- evaluator consensus 与 disagreement。

### 16.2 基因型

- exact GT concordance；
- `0/0`、`0/1`、`1/1` 三分类 confusion matrix；
- macro-F1 和 balanced accuracy；
- no-call rate；
- heterozygous/homozygous concordance；
- allele dosage error；
- 可选 GQ calibration。

### 16.3 相位

- phased genotype concordance；
- phase call rate；
- switch error；
- flip error；
- phase-set continuity；
- truth 允许时的 haplotype assignment consistency。

### 16.4 类型

至少：

```text
DEL INS DUP INV BND CNV CPX OTHER
```

truth 或工具不支持的类型明确标记 eligibility。

### 16.5 长度

默认 bins：

```text
50-99
100-499
500-999
1,000-4,999
5,000-9,999
10,000-49,999
>=50,000
```

可增加 `<50` 作为大 indel/小 SV 边界分析，但不混入标准 `>=50 bp` 总分。

### 16.6 区域与上下文

- autosomes、X、Y、PAR/non-PAR；
- easy/difficult；
- segmental duplication；
- tandem repeat；
- homopolymer；
- low complexity；
- low mappability；
- GC extreme；
- coding/exonic；
- medically relevant。

### 16.7 表示与断点

- start/end distance；
- length similarity；
- inserted allele sequence similarity；
- reciprocal overlap；
- haplotype edit-distance partial credit；
- symbolic versus sequence-resolved rate。

### 16.8 资源

Snakemake `benchmark` 与 `--benchmark-extended` 记录：

- wall time；
- CPU time 和 mean load；
- peak RSS、PSS、USS、VMS；
- I/O；
- input size；
- threads 和 requested resources。

额外记录：

- 临时文件峰值；
- 最终磁盘占用；
- container/image 大小；
- scheduler queue 和实际运行时间，若 profile 可提供。

## 17. 唯一最终得分：PGBenchScore

每个适用的
`(run_id, sample, tool, official_score_mode, primary_truth_profile)` tuple
只输出一个 `PGBenchScore`，范围为 0–100。一次 run 可以有多个工具 tuple，但不存在第二个 run-level
聚合分。报告不生成 `rank`、`leaderboard_position` 或按分数排序的结论。

默认 score profile 为 `pgbench_v1`。所有 point 固定，不因工具缺少某项能力而重新分配。工具没有相位、无法恢复 out-of-panel truth 或没有合法 pangenome allele linkage 时，相应 point 为 0；benchmark 基础设施故障则整个 score 状态为 `invalid`，不能伪装成工具低分。

### 17.1 总体 point 构成

| 大项 | Point |
|---|---:|
| 多评测器准确性 | 70 |
| 泛基因组与困难区域稳健性 | 15 |
| 资源效率 | 10 |
| Rule 可追溯与可复现性 | 5 |
| **合计** | **100** |

```text
PGBenchScore =
  EvaluatorPoints
  + PangenomeRobustnessPoints
  + ResourcePoints
  + TraceabilityPoints
```

内部计算保留完整浮点精度，最终展示四舍五入到小数点后两位。

### 17.2 Truvari：28 point

Truvari 占 70-point evaluator layer 的 40%，即最多 28 point。

| Truvari 指标 | Point | 计算 |
|---|---:|---|
| Overall event F1 | 10 | `10 * F1` |
| SVTYPE macro-F1 | 3 | `3 * macro_F1_by_svtype` |
| Length-bin macro-F1 | 3 | `3 * macro_F1_by_length` |
| Exact GT concordance | 5 | `5 * concordance` |
| Breakpoint/size/sequence fidelity | 4 | `4 * mean_similarity` |
| Difficult-context F1 | 3 | `3 * difficult_macro_F1` |

其中：

```text
mean_similarity =
  0.35 * sequence_similarity
  + 0.25 * size_similarity
  + 0.20 * start_accuracy
  + 0.20 * end_accuracy
```

```text
size_similarity =
  min(abs(query_SVLEN), abs(truth_SVLEN))
  / max(abs(query_SVLEN), abs(truth_SVLEN))

breakpoint_scale =
  max(100 bp, 0.10 * abs(truth_SVLEN))

start_accuracy =
  max(0, 1 - abs(query_POS - truth_POS) / breakpoint_scale)

end_accuracy =
  max(0, 1 - abs(query_END - truth_END) / breakpoint_scale)
```

`sequence_similarity` 使用 evaluator 输出的归一化等位序列相似度；缺失插入序列时该子项为 0，不把权重转移给其他子项。

### 17.3 Aardvark：24.5 point

Aardvark 占 evaluator layer 的 35%，即最多 24.5 point。

| Aardvark 指标 | Point | 计算 |
|---|---:|---|
| Exact haplotype event F1 | 8 | `8 * exact_F1` |
| Partial-credit F1 | 6 | `6 * partial_F1` |
| Exact GT concordance | 4 | `4 * concordance` |
| Complex-cluster recovery F1 | 4 | `4 * complex_F1` |
| Allele-sequence fidelity | 2.5 | `2.5 * sequence_similarity` |

partial-credit F1 使用 evaluator 给出的匹配权重：

```text
weighted_precision = sum(match_weight) / query_count
weighted_recall = sum(match_weight) / truth_count
partial_F1 =
  2 * weighted_precision * weighted_recall
  / (weighted_precision + weighted_recall)
```

### 17.4 vcfdist：17.5 point

vcfdist 占 evaluator layer 的 25%，即最多 17.5 point。

| vcfdist 指标 | Point | 计算 |
|---|---:|---|
| Representation-aware F1 | 7 | `7 * F1` |
| GT concordance | 3.5 | `3.5 * concordance` |
| Phase accuracy | 3 | `3 * phase_accuracy` |
| Complex representation consistency | 2 | `2 * consistency` |
| Eligible-call coverage | 2 | `2 * coverage` |

```text
phase_accuracy =
  0.50 * phased_GT_concordance
  + 0.25 * (1 - switch_error_rate)
  + 0.25 * (1 - flip_error_rate)
```

```text
consistency =
  representation_equivalent_clusters
  / evaluated_complex_clusters

coverage =
  evaluated_eligible_truth_clusters
  / predefined_vcfdist_eligible_truth_clusters
```

若工具不输出 vcfdist 所需的相位信息，vcfdist evaluator status 为
`completed_with_tool_ineligible_zero`，所有 vcfdist 计分 metric 明确写 0 并保留固定 denominator 和原因；
其 17.5 point 为 0，权重不转移给其他 evaluator。score card 仍展示这些 F1/GT/phase 字段及状态。

### 17.5 Evaluator 权重核对

```text
Truvari:  28.0 / 70 = 40%
Aardvark: 24.5 / 70 = 35%
vcfdist:  17.5 / 70 = 25%
```

等价公式：

```text
EvaluatorPoints =
  0.70 * (
    0.40 * TruvariScore
    + 0.35 * AardvarkScore
    + 0.25 * VcfdistScore
  )
```

三个 evaluator score 均先归一化到 0–100。

```text
TruvariScore  = 100 * TruvariPoints  / 28.0
AardvarkScore = 100 * AardvarkPoints / 24.5
VcfdistScore  = 100 * VcfdistPoints  / 17.5
```

因此按 evaluator score 加权和直接把各 evaluator point 相加严格等价，不存在隐藏的第四套融合公式。

### 17.6 泛基因组与困难区域稳健性：15 point

| 指标 | Point | 计算 |
|---|---:|---|
| In-panel exact genotype macro-F1 | 4 | `4 * fused_in_panel_GT_macro_F1` |
| Pangenome allele-link precision | 2 | `2 * allele_link_precision` |
| Pangenome allele-link coverage | 2 | `2 * allele_link_coverage` |
| Population-AF robustness macro-F1 | 2 | `2 * macro_F1_by_population_AF` |
| Graph-complexity robustness macro-F1 | 1 | `1 * macro_F1_by_graph_complexity` |
| Difficult-context macro-F1 | 2 | `2 * macro_F1_by_context` |
| Out-of-panel truth recovery F1 | 2 | `2 * novel_truth_F1` |

这些指标只读取冻结后的 pangenome allele ledger 和 hidden scoring labels。除两个 allele-link 指标外，
F1 均按相同 evaluator 权重融合：

```text
fused_F1(stratum) =
  0.40 * Truvari_F1(stratum)
  + 0.35 * Aardvark_F1(stratum)
  + 0.25 * vcfdist_F1(stratum)
```

in-panel genotype 使用 candidate-level 事实表：

```text
results/HG002/<track>/<tool>/evaluation/pangenome_gt_ledger.tsv
```

每个 `PANGENOME_ALLELE_ID × evaluator` 至少记录 hidden truth GT、原始/标准化 tool GT、
explicit/implicit/no-call、evaluator assignment、truth genotype class、predicted class 以及
one-vs-rest TP/FP/FN。三个 evaluator 分别输出：

```text
truvari.pangenome.in_panel.gt.macro_f1
aardvark.pangenome.in_panel.gt.macro_f1
vcfdist.pangenome.in_panel.gt.macro_f1
```

最终：

```text
fused_in_panel_GT_macro_F1 =
  0.40 * truvari.pangenome.in_panel.gt.macro_f1
  + 0.35 * aardvark.pangenome.in_panel.gt.macro_f1
  + 0.25 * vcfdist.pangenome.in_panel.gt.macro_f1
```

对每个 genotype class，正确预测计 TP；错误 GT 同时造成 truth class FN 和 predicted class FP；
no-call 造成 truth class FN；只有 manifest 预先允许 `variant_sites + hom_ref` 时，缺失 candidate 才是
显式审计过的隐式 `0/0`。

权重不因某工具无相位或某 evaluator 对该工具不可运行而重新分配；工具导致的 evaluator
`completed_with_tool_ineligible_zero` 对应 F1 为 0。profile 本身不支持的 stratum 在运行前从全部工具的共同 universe
中排除，并在 `score_weights.yaml` 中版本化，不能运行后选择。

```text
allele_link_precision =
  truth_consistent_in_panel_exact_or_equivalent_calls
  / calls_with_link_status_in {
      in_panel_exact,
      in_panel_equivalent,
      ambiguous,
      invalid_allele_id,
      wrong_link
    }

allele_link_coverage =
  truth_positive_in_panel_alleles_with_a_correct_linked_call
  / all_eligible_truth_positive_in_panel_alleles
```

当 precision 分母为 0 但固定 universe 中存在 eligible in-panel truth 时，precision 为 0。
`invalid_allele_id`、`ambiguous` 和 `wrong_link` 进入 precision 分母但不进入分子；
`out_of_panel` 与 `unresolved` 不进入 link precision 分母，后者会降低 coverage。

population AF bins 固定为：

```text
ultrarare:       0 < AF < 0.001
rare:        0.001 <= AF < 0.01
low_frequency: 0.01 <= AF < 0.05
common:       0.05 <= AF <= 1
```

graph-complexity bins 固定为：

```text
simple_biallelic
multiallelic
nested_or_overlapping
repeat_ambiguous
```

`graph_complexity_v1` 在工具运行前为每个 panel allele 冻结唯一 class。分类优先级从高到低：

1. `repeat_ambiguous`：SV reference span 与 tandem-repeat/segmental-duplication union 的重叠比例
   `>=0.50`；INS 使用 `POS±50 bp` 作为 span；
2. `nested_or_overlapping`：同一 graph connected component 中存在另一个 candidate，二者 reference span
   至少重叠 1 bp，且不只是相同 anchor 的多 ALT；
3. `multiallelic`：归一化后相同 reference anchor/locus 有至少两个不同 ALT；
4. `simple_biallelic`：其余 allele。

分类算法版本、annotation hashes、connected-component recipe 和最终 class 写入 allele ledger。
若上级 class 命中，不再落入下级 class，保证每个 allele 只计一次。

macro-F1 对 profile 预定义且 truth count 大于 0 的 strata 等权平均。工具在某个 truth-positive
stratum 无调用时，该 stratum F1 为 0，不能从 macro 中删除。

`out-of-panel truth recovery F1` 只在 population pangenome 不含、但 HG002 truth 中存在的 cluster
上计算。纯 genotyping-only 工具不能发现此类 SV 时，该项为 0。

对每个 evaluator：

```text
novel_truth_universe =
  primary-truth clusters inside benchmark BED
  with no exact/equivalent population-panel allele

novel_query_universe =
  canonical query calls with link_status = out_of_panel

novel_TP = evaluator-matched pairs between the two universes
novel_FP = unmatched calls in novel_query_universe
novel_FN = unmatched clusters in novel_truth_universe
novel_truth_F1 = fixed 40:35:25 fusion of evaluator-specific novel F1
```

`ambiguous`、`unresolved`、`invalid_allele_id` 和 `wrong_link` 不冒充 novel TP；它们通过 link
precision/coverage 被扣分并单独报告。`novel_truth_universe` 在 profile freeze 时必须非空，否则
该 score profile 无效，不能把 2 point 动态移给其他项。

### 17.7 资源效率：10 point

资源 point 使用固定 `target` 和 `limit`，不依赖其他工具的表现，因此不会形成隐式排名。

| 资源 | Point |
|---|---:|
| Official-mode wall time | 4 |
| Peak RSS memory | 3 |
| Peak temporary + final disk | 2 |
| CPU-hour efficiency | 1 |

单项归一化：

```text
resource_fraction(x) = 1, x <= target
resource_fraction(x) = 0, x >= limit
resource_fraction(x) =
  (log(limit) - log(x))
  / (log(limit) - log(target)),
  target < x < limit
```

单项 point 等于 `max_point * resource_fraction(x)`。`target`、`limit`、硬件指纹和并发策略写入 `pgbench_v1` profile 与最终 manifest。

official-mode 资源聚合固定为：

```text
wall_time   = 最早计费 stage started_at 到最晚计费 stage finished_at 的 elapsed time
peak_RSS    = 同一时间点全部计费子进程 RSS 的最大和
disk_peak   = tool-scoped work/output directory 相对开始时基线的最大增量
CPU_hours   = 全部计费进程 user+system CPU seconds 之和 / 3600
```

若 `benchmark_repeats > 1`，四项均取成功重复的中位数；所有工具必须使用相同 repeat 数和 cache policy。
任一正式重复运行失败时工具状态为 `runtime_failed`，不能只挑成功重复计分。

`pgbench_v1` 默认预算：

| Official mode | 资源 | Target | Limit |
|---|---|---:|---:|
| end-to-end | wall time | 12 h | 72 h |
| end-to-end | peak RSS | 64 GiB | 512 GiB |
| end-to-end | disk | 250 GiB | 2 TiB |
| end-to-end | CPU-hours | 96 | 2000 |
| caller-only | wall time | 4 h | 24 h |
| caller-only | peak RSS | 64 GiB | 256 GiB |
| caller-only | disk | 100 GiB | 500 GiB |
| caller-only | CPU-hours | 32 | 512 |

如果正式实验平台确实无法满足这些预算，只能通过发布新的版本化 score profile 修改，不能在看到工具结果后调整。

### 17.8 Rule 可追溯与可复现性：5 point

| 要求 | Point |
|---|---:|
| 所有 pre-score audited rule/job 均有合法 manifest、log 和 benchmark record | 2 |
| 最终 score 可通过 input/output hashes 追溯到原始 BAM/FASTQ、pangenome、truth 和 raw VCF | 1 |
| Conda lock 或 container digest 完整 | 1 |
| config snapshot、随机种子、代码 commit 和 evaluator weight profile 完整 | 1 |

前 2 point 的精确公式：

```text
manifest_completeness =
  audited_jobs_with_valid_manifest_log_benchmark
  / all_executed_pre_score_audited_jobs

manifest_points = 2 * manifest_completeness
```

其余三个 1-point 条件为 binary gate：全部字段和 hash 校验通过得 1，否则得 0，不给部分分。
若工具执行、canonicalization、allele linkage、任一默认 evaluator、evaluator fusion、resource collection、
pre-score audit 或 scoring 的核心 manifest 缺失，最终 score 状态直接为 `invalid`。
非核心 score ancestor 缺失时数值按上式扣分且状态至少为 `provisional`；只有
`manifest_completeness=1` 且三个 binary gate 全通过，才允许 `valid`。

### 17.9 F1 的强制要求

最终 score card 必须明确列出：

- Truvari overall F1；
- Truvari SVTYPE macro-F1；
- Truvari length-bin macro-F1；
- Aardvark exact F1；
- Aardvark partial-credit F1；
- vcfdist representation-aware F1；
- pangenome in-panel genotype macro-F1；
- population-AF、graph-complexity 和 difficult-context macro-F1；
- out-of-panel truth recovery F1。

没有上述 F1 明细的报告不能标记为正式 `valid` score。

### 17.10 置信区间与 score 状态

precision、recall、F1、GT concordance 和主要 point component 使用 variant-cluster 分层 bootstrap 生成 95% CI。

最终状态：

```text
valid       三个默认 evaluator 均成功，或按 profile 合法记录 completed_with_tool_ineligible_zero；point 明细和核心 provenance 完整
provisional 指标可计算，但存在明确的非核心缺失；不得用于正式结论
invalid     evaluator 基础设施失败、核心 rule lineage 缺失或 truth/reference 不一致
```

置信区间不参与 point 求和。bootstrap 使用 score profile 中固定 seed，并以 truth variant-cluster 为
resampling unit，避免同一复杂 cluster 的多个表示被当成独立样本。

### 17.11 Metric Dictionary 与缺失值规则

仓库必须提供版本化 `config/metric_dictionary.yaml`。三个 evaluator parser 只能输出该 schema
登记的标准字段；计分脚本不能直接读取未登记的临时列。

| Point component | 标准字段 | 固定 universe / 分母 |
|---|---|---|
| Truvari overall event F1 | `truvari.event.overall.f1` | primary BED 内 profile-eligible truth clusters 与 query calls |
| Truvari SVTYPE macro-F1 | `truvari.event.svtype.macro_f1` | 预定义 truth-positive SVTYPE strata |
| Truvari length macro-F1 | `truvari.event.length.macro_f1` | 预定义 truth-positive length bins |
| Truvari exact GT concordance | `truvari.gt.exact_concordance` | Truvari matcher 链接到的固定 challenge genotype universe；no-call 为不一致 |
| Truvari fidelity | `truvari.match.{sequence,size,start,end}_similarity` | matched truth clusters；缺 query sequence 的 sequence 值为 0 |
| Truvari difficult F1 | `truvari.context.difficult.macro_f1` | score profile 列出的 difficult BED strata |
| Aardvark exact F1 | `aardvark.haplotype.exact_f1` | primary truth local-haplotype clusters |
| Aardvark partial F1 | `aardvark.haplotype.partial_f1` | 同一 clusters；match weight 限定在 `[0,1]` |
| Aardvark exact GT concordance | `aardvark.gt.exact_concordance` | Aardvark linkage 后的固定 challenge genotype universe |
| Aardvark complex F1 | `aardvark.complex.exact_f1` | truth 标记 CPX，或同一 local window 中至少两个相互作用 truth allele 的 clusters |
| Aardvark sequence fidelity | `aardvark.haplotype.sequence_similarity` | 具有 truth sequence 的 clusters；missing query sequence 为 0 |
| vcfdist representation F1 | `vcfdist.representation.f1` | 运行前冻结的 vcfdist eligible truth clusters |
| vcfdist GT concordance | `vcfdist.gt.exact_concordance` | vcfdist eligible challenge genotypes |
| vcfdist phase accuracy | `vcfdist.phase.accuracy` | truth phase-eligible heterozygous clusters/phase sets |
| vcfdist consistency | `vcfdist.complex.consistency` | predefined eligible complex clusters |
| vcfdist coverage | `vcfdist.eligible.coverage` | predefined vcfdist eligible truth clusters |
| Pangenome in-panel GT macro-F1 | `{truvari,aardvark,vcfdist}.pangenome.in_panel.gt.macro_f1` 与 `pangenome.in_panel.gt.fused_macro_f1` | candidate-level ledger 中 truth count > 0 的 `0/0`,`0/1`,`1/1` 和适用 haploid classes |
| Allele-link precision/coverage | `pangenome.link.{precision,coverage}` | 17.6 的固定 allele ledger 分母 |
| AF/graph/context/novel F1 | `pangenome.{af,graph,context,novel}.*_f1` | manifest/score profile 在工具运行前冻结的 strata |

统一规则：

1. 进入 point 公式的所有比例、相似度、concordance 和 F1 必须在 `[0,1]`。仅允许对
   `1e-12` 以内的浮点越界 clamp；更大越界视为 parser/infrastructure error，score 为 `invalid`。
2. event F1 的 TP/FP/FN 来自各 evaluator 自己的 matching ledger；不能把三套 evaluator 的 TP
   先投票后再冒充某个 evaluator F1。
3. `macro_F1 = mean(F1_s)`，其中 `s` 只包括 score profile 预先声明且 truth count 大于 0 的 strata。
   truth-positive stratum 中工具无调用时 `F1_s=0`。
4. exact genotype macro-F1 使用 one-vs-rest 方式分别计算每个预定义 genotype class 的 F1 后等权平均；
   no-call 计为该 truth class 的 FN。基因型 concordance 为 exact matches / 全部 eligible genotypes。
5. matched-pair similarity 对固定 truth cluster 求平均；一个 truth cluster 有多个 query match 时使用
   evaluator 决定的 primary assignment，其余调用仍按 evaluator 规则计 FP，不能选最高相似度。
6. 工具缺字段、缺相位、缺 sequence 或输出 no-call 时，对应 metric 按上文规则记 0；固定 point 不转移。
7. evaluator 安装失败、parser schema 不匹配、truth/reference 错误或固定 universe 无法构建属于
   infrastructure error，整个 tuple 为 `invalid`，不能把 metric 写成 0。
8. 若某 stratum 在冻结 universe 中 truth count 为 0，记 `undefined` 并从所有工具的同一 macro
   中排除；不能按工具输出动态决定是否排除。
9. 所有 numerator、denominator、eligible count、undefined reason 和 parser source field
   必须写入 `metrics.long.tsv` 与 point breakdown。

## 18. 报告

主要文件：

```text
results/summary/metrics.long.tsv
results/summary/metrics.json
results/summary/score.tsv
results/summary/point_breakdown.tsv
results/summary/evaluator_disagreement.tsv
results/summary/discordant_review.tsv
results/summary/failed_or_ineligible_tools.tsv
results/summary/run_provenance.json
results/report/index.html
results/report/tool_cards/<tool>.html
```

HTML 至少包含：

- 执行概览和输入 QC；
- truth/profile 说明；
- 每个适用工具 tuple 的唯一 `PGBenchScore`、状态和逐 point 明细；
- SVTYPE/length/context heatmap；
- GT confusion matrix；
- switch/flip errors；
- breakpoint、size、sequence similarity；
- evaluator agreement/upset plot；
- runtime-memory-disk scatter（不标记名次或“最佳”）；
- 失败与不适用工具；
- 可下载的逐变异 review 表。

多工具运行时，每个工具获得独立 score card。HTML 按配置顺序或工具 ID 展示，不按得分排序，不显示名次，也不生成“最佳工具”结论。

报告首页明确显示：

- 主 truth 为 v5.0q draft；
- 当前 hardware fingerprint；
- 评分 profile；
- evaluator 版本；
- 外部工具来源；
- 是否为 caller-only 或 end-to-end。

## 19. 失败处理

### 19.1 硬失败

- BAM/reference/truth 不兼容；
- 配置 schema 无效；
- 工具 ID 冲突；
- 必需 reference 或 truth 校验和失败；
- challenge panel truth 泄漏检查失败。

### 19.2 工具级失败

单个工具失败时：

- 写入 `status.json`；
- 保存 stdout/stderr 和退出码；
- 不生成伪造的零分结果；
- 使用 `snakemake --keep-going` 继续其他工具；
- 报告状态为 `runtime_failed`。

### 19.3 输出失败

- VCF 语法无效：`invalid_output`；
- 无样本或无 GT：不能进入 genotyping track；
- sample/reference 不符：`incompatible_output`；
- 空 VCF：保留事实并由适用赛道评分；若工具宣称运行成功但格式不完整，则为 invalid；
- 缺失 provenance：`provenance_failed`。

### 19.4 不适用

- 缺 short reads；
- 缺 graph/panel；
- 工具不支持 CLR；
- 工具对 official profile 的全部计分 SVTYPE 均不支持；

以上状态为 `not_applicable`，不等于失败或零分。

工具能够完成 official mode 但无相位输出时，不把整个工具标为 `not_applicable`；其 vcfdist/phase 相关固定 point 为 0，并继续计算其余 point。
工具只支持部分 SVTYPE 时同理：支持类型照常评测，未支持类型的固定 strata F1 为 0，不重新分配 point。

## 20. 测试

### 20.1 单元测试

- 配置和 tool manifest schema；
- BAM/reference 兼容性逻辑；
- candidate panel 正负标签生成；
- truth 信息盲化；
- VCF canonical normalization；
- evaluator parser；
- evaluator 40:35:25 加权与 point 公式；
- PGBenchScore 及 point 明细；
- tool-ineligible、基础设施失败和 invalid score 边界条件；
- 每个 rule manifest 与 rule-lineage 审计。

### 20.2 Mock external tool

测试必须包含一个自研 genotyper mock：

1. benchmark 创建 challenge VCF；
2. Snakemake 调用 mock 工具；
3. mock 工具读取 BAM/候选 panel；
4. mock 工具在当前 run 内生成 VCF；
5. VCF 通过校验；
6. 进入多评测器和 score；
7. provenance 能证明不是预计算 VCF。

### 20.3 集成测试

使用小型 synthetic reference、BAM、truth 和 population candidates：

- `snakemake -n`；
- DAG 构建；
- local profile 执行；
- 至少一个 mapping caller；
- 至少一个内置 genotyper；
- 一个 external genotyper；
- evaluator、加权、PGBenchScore 和 HTML report。

### 20.4 全基因组验收

由于当前本地环境无法访问用户服务器 BAM，全基因组验收在目标 Linux/HPC 上进行：

1. preflight；
2. 构建并校验 hs37d5 population pangenome manifest、panel 和 graph assets；
3. `caller_only_shared_alignment` 诊断；
4. 一个正式 end-to-end pangenome/graph tool；
5. Sniffles2 force、kanpig 和 SVJedi-graph；
6. 一个 external mock 或用户工具；
7. v5.0q 三评测器评测与 40:35:25 加权；
8. PGBenchScore、rule-lineage、资源和 HTML 报告。

## 21. 第一版交付范围

### 21.1 完整实现

- 项目骨架和 schema；
- BAM/reference/truth preflight；
- canonical FASTQ 和 HP FASTQ；
- hs37d5 population pangenome manifest、panel、stable allele IDs 和基础 graph assets；
- v5.0q、v0.6 和 v3.6 资源配置；
- Sniffles2、cuteSV、pbsv、SVIM；
- Sniffles2 force、kanpig、SVJedi-graph；
- generic external tool runner；
- blinded genotyping challenge panel；
- canonical 和 evaluator-specific normalization；
- Truvari、Aardvark、vcfdist；
- evaluator weighting；
- 单一 PGBenchScore 与 point breakdown；
- 全 DAG rule manifests、rule-lineage 和 provenance audit；
- Obsidian design sync script、binding 校验和 KB lint；
- TSV/JSON/HTML report；
- local 与 SLURM profile 示例；
- mock 和小型集成测试。

### 21.2 可启用实现

- Flye + dipcall；
- SVIM-asm；
- vg；
- generic fragment/signature integration；
- PanGenie，要求短读长和 panel；
- hap-eval、Sansa、Wittyer compatibility。

这些模块具有完整 manifest、输入检查和输出契约，但默认不因缺少大型资产而进入当前 HG002 CLR DAG。

## 22. 验收标准

设计完成后的实现必须满足：

1. 示例配置可通过 JSON Schema。
2. `snakemake -n` 能构建核心 DAG。
3. 缺少用户服务器 BAM 时给出明确 preflight 信息，而不是 Python traceback。
4. 工具模块输出统一 artifacts。
5. 外部 genotyper 由 benchmark 实际执行，不接受预先提交的 VCF。
6. benchmark 为外部 genotyper 提供 HG002、reference 和 blinded candidate panel。
7. 外部与内置工具使用相同 evaluator、固定权重和 PGBenchScore。
8. `not_applicable`、`runtime_failed` 和 `invalid_output` 被区分。
9. v5.0q 和 v0.6 结果独立报告。
10. 三个默认 evaluator 均按 40:35:25 固定权重进入正式得分。
11. 逐变异 disagreement 可追溯到原始 VCF 记录。
12. 资源指标与 hardware fingerprint 一起报告。
13. 小型 integration test 在无全基因组数据时可执行。
14. README 提供新增自研 genotyper 的最小示例。
15. 每个适用工具 tuple 只有一个 PGBenchScore，不包含 run-level 聚合分、rank 或 leaderboard。
16. score card 明确列出所有强制 F1、evaluator point 和最终 point 求和。
17. 每个正式 rule 都能通过 manifest 追溯到输入、代码、环境、上游 rule 和输出哈希。
18. Obsidian canonical design note 的 repo commit、设计哈希和实现状态与仓库一致。

## 23. 关键风险与缓解

### v5.0q 仍为 draft

缓解：保留 evaluator 原始结果、固定权重、分歧报告、人工复核表和 truth version。

### 复杂 SV 表示差异

缓解：canonical 与 evaluator-specific views 分离；至少使用两个 haplotype/representation-aware evaluator。

### `0/0` 负位点误标

缓解：只在 benchmark BED 内选择；任一高灵敏 matcher 命中即排除；排除复杂 cluster 和 evaluator-discordant 位点；报告 panel 构建审计表。

### 从 BAM 提取 reads 的信息损失

缓解：报告 hard clipping 和 tag loss；caller-only 使用原 BAM；HP reads 单独提取。

### 不同范式资源不可比

缓解：每次 run 只指定一个 official score mode；使用固定 target/limit budgets；同时展示原始资源。

### 外部工具任意代码

缓解：推荐固定 digest 容器、受限挂载、独立输出目录和完整 provenance。框架不把外部代码视为可信代码。

## 24. Obsidian 设计—代码同步契约

本项目绑定到：

```text
/Users/wanghao/Desktop/risefl_mvp/memory/
  Research/hg002-grch37-pangenome-sv-benchmark/
```

最新设计的 Obsidian canonical note：

```text
Writing/HG002-GRCh37泛基因组SV-Benchmark设计规范.md
```

该 note 必须包含：

```yaml
repo_root: /Users/wanghao/Documents/Codex/2026-07-16/referenced-chatgpt-conversation-this-is-untrusted
repo_design_path: docs/superpowers/specs/2026-07-16-hg002-grch37-sv-benchmark-design.md
repo_commit: <current commit>
git_head: <current commit>
git_dirty: false
git_diff_sha256: <sha256 of current git diff; clean tree uses empty-diff sha256>
git_status_sha256: <sha256 of porcelain status bytes>
untracked_files_sha256: <sha256 of sorted untracked path/content inventory>
repo_state_sha256: <sha256 combining HEAD, status, diff and untracked inventory>
tracked_tree_sha256: <sha256 of sorted tracked file inventory and blob ids>
design_sha256: <current design file sha256>
code_commit: <commit represented by current code_state>
code_state: design-only|partial|implemented|verified
score_profile: pgbench_v1
score_profile_sha256: <hash or null when design-only>
pangenome_id: hs37d5_1kg_phase3_sv_v1
rule_registry_sha256: <hash or null when design-only>
tool_registry_sha256: <hash or null when design-only>
sync_status: current|stale
```

同步规则：

1. 当前仓库中的 rule、tool adapter、score profile、pangenome manifest schema 或 evaluator parser 发生语义变化时，必须同时检查设计规范。
2. 设计规范变化后，必须同步覆盖 Obsidian canonical note，而不是新增多个版本造成事实源分叉。
3. 每次同步更新 Git 状态、`design_sha256`、`code_commit`、`code_state`、registry/profile hashes 和 `updated`。
4. 更新 `00-Hub.md`、`01-Plan.md`、`02-Index.md` 与 `_system/registry.md`。
5. 执行 registry、index、wikilink 和项目 lint。
6. 若 hash 或 commit 不一致，note 必须标记 `sync_status: stale`，不能声称是最新成果。
7. 实现完成后，Obsidian note 的工具矩阵、rule 清单、score point 和代码路径必须与实际代码一致。
8. canonical note 的同步正文必须位于
   `<!-- PGBENCH:BEGIN GENERATED DESIGN -->` 与
   `<!-- PGBENCH:END GENERATED DESIGN -->` 之间；同步只替换该区域，保留区域外人工注释。
9. 同步先写同目录临时文件，完成 frontmatter、generated markers、hash 和 wikilink 校验后再原子替换。
   失败时保留旧 canonical note，并将 check 结果写为 stale。
10. `check_obsidian_sync` 是只读检查，可在没有 HG002 BAM 时独立执行。它比较设计 hash、Git 状态、
    score profile、pangenome manifest schema、resolved rule registry、tool registry 与实现 binding。
11. `repo_commit`/`git_head` 表示同步时仓库 HEAD；`code_commit` 表示 `code_state` 所描述的实现基线。
    在 `design-only` 阶段二者可以相同，但 registry/profile hash 必须显式为 null，不能假装代码已实现。
12. Obsidian sync/check 属于 post-score workflow；vault 不可用不得影响 `PGBenchScore`。
13. 只有 `git status --porcelain=v1 -z --untracked-files=all` 输出为空时，canonical note 才能标记
    `sync_status: current` 和 `git_dirty: false`。存在 staged、unstaged 或未跟踪实现文件时只能同步为
    `stale`；不得漏掉 untracked 文件。被 `.gitignore` 明确排除的 runtime outputs 不属于实现状态。

Git hash 算法固定为：

```text
git_status_bytes =
  exact bytes from:
  git status --porcelain=v1 -z --untracked-files=all

git_diff_bytes =
  exact bytes from:
  git diff --binary --no-ext-diff HEAD --

tracked_tree_bytes =
  exact NUL-delimited bytes from:
  git ls-files -s -z

untracked_inventory_bytes =
  for each untracked path sorted by UTF-8 byte order:
  NUL + relative_path + NUL + sha256(file_content_or_tree_manifest)

git_status_sha256     = sha256(git_status_bytes)
git_diff_sha256       = sha256(git_diff_bytes)
tracked_tree_sha256   = sha256(tracked_tree_bytes)
untracked_files_sha256 = sha256(untracked_inventory_bytes)
repo_state_sha256 =
  sha256(
    git_head + NUL
    + git_status_sha256 + NUL
    + git_diff_sha256 + NUL
    + tracked_tree_sha256 + NUL
    + untracked_files_sha256
  )
```

命令必须在 repo root、`LC_ALL=C` 下运行；hash 为小写十六进制 SHA-256。clean tree 的 status、diff
和 untracked inventory 都使用空字节串的 SHA-256。

第一版代码应提供：

```text
workflow/scripts/sync_obsidian_design.py
workflow/scripts/check_obsidian_sync.py
snakemake sync_obsidian_design
snakemake check_obsidian_sync
```

同步目标读取 repo 设计文件和 Git 状态，更新 canonical note 元数据，并运行 KB lint。只读 check
输出 `results/provenance/obsidian-sync-check.json`。vault path 由 CLI、环境变量、本地 binding 或配置解析，
不硬编码进可移植公共配置。

## 25. 参考

- NIST GIAB v5.0q：
  <https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/v5.0q/>
- NIST GIAB 项目：
  <https://www.nist.gov/programs-projects/genome-bottle>
- GIAB v0.6 paper：
  <https://www.nature.com/articles/s41587-020-0538-8>
- GIAB stratifications：
  <https://pmc.ncbi.nlm.nih.gov/articles/PMC11489684/>
- 1000 Genomes Phase 3 SV：
  <https://www.internationalgenome.org/phase-3-structural-variant-dataset/>
- Snakemake modularization：
  <https://snakemake.readthedocs.io/en/stable/snakefiles/modularization.html>
- Snakemake benchmarking：
  <https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html>
- Truvari：
  <https://github.com/ACEnglish/Truvari>
- Aardvark：
  <https://github.com/PacificBiosciences/aardvark>
- vcfdist：
  <https://github.com/TimD1/vcfdist>
- Sniffles2：
  <https://github.com/fritzsedlazeck/Sniffles>
- kanpig：
  <https://github.com/ACEnglish/kanpig>
- SVJedi-graph：
  <https://github.com/SandraLouise/SVJedi-graph>
- PanGenie：
  <https://github.com/eblerjana/pangenie>
