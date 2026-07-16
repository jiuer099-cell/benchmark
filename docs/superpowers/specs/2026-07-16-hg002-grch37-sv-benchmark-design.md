# HG002 GRCh37 结构变异检测与基因分型 Benchmark 设计规范

日期：2026-07-16

状态：已确认设计

工作流：Snakemake 9

默认样本：HG002

默认输入 BAM：
`/home/zhj/Experiment/SVDF-main/data/HG002.Sequel.10kb.pbmm2.hs37d5.whatshap.haplotag.RTG.10x.trio.bam`

## 1. 目标

构建一个可扩展、可复现、可用于论文评测的结构变异 benchmark。框架需要：

1. 对 HG002 PacBio Sequel/CLR、hs37d5 BAM 执行统一前处理。
2. 支持从 canonical reads 开始分流的完整端到端评测，以及共享 BAM 上的 caller-only 评测。
3. 兼容以下分析范式：
   - mapping-based；
   - assembly-based；
   - graph/pangenome-based；
   - fragment/signature discovery and integration；
   - genotyping-only。
4. 内置若干结构变异检测与基因分型工具。
5. 允许使用者提供自己的可执行基因分型工具，由 benchmark 实际运行该工具分析 HG002，并生成待评分 VCF。
6. 统一校验和标准化工具输出，但保留复杂 SV 的原始表示。
7. 使用多个 SV 比较工具，分别保留原始指标并生成带分歧惩罚的融合指标。
8. 覆盖检测、基因型、相位、类型、长度、区域、基因组上下文、断点、等位序列和资源消耗。
9. 输出可审计的长表指标、排行榜、逐变异判定表和 HTML 报告。
10. 定义独立的 `SVBench-100` 评分体系。

## 2. 非目标

第一版不承诺：

- 在本地完整运行用户服务器上的全基因组 BAM；
- 自动下载并默认启用所有大型 pangenome graph 和 population panel；
- 把短读长专用工具强行用于当前 PacBio CLR 数据；
- 用单个数字掩盖评测器分歧或 truth set 局限；
- 接受使用者预先生成的 VCF 作为自研工具正式评测结果。

预计算 VCF 可以用于开发阶段的格式调试，但不得进入正式排行榜。正式结果必须由当前 benchmark run 调用工具后生成，并具有完整 provenance。

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
  -> SVBench-100
```

### 3.2 从 reads 开始的端到端公平性

用户提供的是已比对 BAM。框架从 BAM 一次性提取 canonical FASTQ，作为端到端赛道所有工具的共同起点。

同时保留 caller-only 赛道，在共享 BAM 上比较只负责 calling/genotyping 的工具。两个赛道不能混合排名：

- `caller_only_shared_bam`
- `end_to_end_from_reads`

### 3.3 适用性不是零分

工具因输入技术、SV 类型、候选 panel、相位或图资产不满足要求时，状态为 `not_applicable`。只有声明支持且完成运行的维度才进入对应赛道。

### 3.4 原始表示必须保留

复杂 SV 可能存在多个等价 VCF 表示。框架生成 canonical VCF，但不对所有工具全局执行 BND 分解、强制多等位拆分或不可逆的表示转换。每个评测器使用独立视图。

## 4. 总体架构

主 Snakefile 负责：

1. 加载并校验配置；
2. 构建公共输入；
3. 动态加载已启用工具模块；
4. 汇总模块输出；
5. 执行标准化、评测、融合、评分和报告。

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
              +----------------------+----------------------+
              |                                             |
      +-------v--------+                            +-------v--------+
      | shared BAM     |                            | canonical reads |
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
                          | fusion / stratify    |
                          +----------+-----------+
                                     |
                          +----------v-----------+
                          | scores / HTML report |
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
│   └── score_weights.yaml
├── profiles/
│   ├── local/
│   └── slurm/
├── workflow/
│   ├── rules/
│   │   ├── common.smk
│   │   ├── preflight.smk
│   │   ├── reference.smk
│   │   ├── truth.smk
│   │   ├── preprocessing.smk
│   │   ├── panel.smk
│   │   ├── normalization.smk
│   │   ├── evaluation.smk
│   │   ├── fusion.smk
│   │   ├── scoring.smk
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
│   │   ├── run-manifest.schema.yaml
│   │   └── metrics.schema.yaml
│   ├── scripts/
│   │   ├── inspect_bam.py
│   │   ├── build_challenge_panel.py
│   │   ├── validate_tool_output.py
│   │   ├── normalize_sv_vcf.py
│   │   ├── parse_truvari.py
│   │   ├── parse_aardvark.py
│   │   ├── parse_vcfdist.py
│   │   ├── fuse_metrics.py
│   │   ├── score_tools.py
│   │   └── render_report.py
│   └── envs/
├── plugins/
│   └── example_genotyper/
├── resources/
│   ├── references/
│   ├── truth/
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

truth:
  primary: giab_hg002_grch37_v5_0q
  secondary:
    - giab_hg002_grch37_v0_6_legacy

execution:
  tracks:
    - caller_only_shared_bam
    - end_to_end_from_reads
  use_conda: true
  use_apptainer: false
  keep_going: true
  benchmark_repeats: 1

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

## 8. 工具范式与内置模块

| 范式 | 标准入口 | 第一版状态 | 代表模块 |
|---|---|---|---|
| Mapping-based | BAM 或 canonical FASTQ | 核心可运行 | Sniffles2、cuteSV、pbsv、SVIM |
| Assembly-based | canonical/HP FASTQ | 可启用 | Flye + dipcall、SVIM-asm |
| Graph/pangenome | FASTQ + graph/panel | 可启用 | vg、SVJedi-graph |
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

execution:
  module: workflow/modules/kanpig/Snakefile
  environment: workflow/envs/kanpig.yaml

outputs:
  vcf: raw/calls.vcf.gz
```

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
- candidate panel ID；
- threads 和资源请求；
- 开始、结束时间；
- 退出码；
- 输出校验和。

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
SVBENCH_RUN_ID
SVBENCH_SAMPLE_ID
SVBENCH_INPUT_BAM
SVBENCH_INPUT_FASTQ
SVBENCH_HP1_FASTQ
SVBENCH_HP2_FASTQ
SVBENCH_REFERENCE_FASTA
SVBENCH_CANDIDATE_VCF
SVBENCH_OUTPUT_DIR
SVBENCH_OUTPUT_VCF
SVBENCH_THREADS
SVBENCH_MEMORY_MB
```

用户的 `run.sh` 负责把统一变量翻译为自研工具 CLI。不得读取评分 truth 路径。

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

正位点来自：

- GIAB HG002 GRCh37 v5.0q structural-variant benchmark VCF；
- 完全处于对应 benchmark BED；
- 经过 engine-specific truth preparation；
- 保留用于复杂表示匹配所需的邻近小变异。

### 11.3 负位点

默认候选来源为 1000 Genomes Phase 3 GRCh37 integrated SV map：
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

它主要覆盖 isolated、sequence-resolved DEL/INS。必须仅使用 PASS 且位于 Tier1 BED 的变异。v0.6 指标不得与 v5.0q 直接融合，只用于复现历史论文。

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
   - 仅在其支持范围内进入融合。

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

## 15. 多评测器融合

只融合语义相同的指标。默认 evaluator 权重：

```yaml
truvari: 0.4
aardvark: 0.4
vcfdist_or_other_eligible: 0.2
```

对值域为 `[0, 1]` 的指标：

```text
central = evaluator 加权中位数
disagreement = evaluator 加权 MAD
fused = clip(central - 0.25 * disagreement, 0, 1)
```

要求：

- 原始 evaluator 指标始终保留；
- evaluator 不适用时不填 0；
- 只有一个 evaluator 可用时标记 `single_engine`，不伪装为多引擎融合；
- 核心 discovery 和 genotype 排名原则上要求至少两个独立 evaluator；
- evaluator disagreement 作为报告的重要结果，而不是只当噪声删除。

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

## 17. SVBench-100

不同范式按赛道排名，不建立跨赛道的伪公平总榜。

### 17.1 权重

| 分项 | Discovery | Genotyping phase-aware | End-to-end phase-aware |
|---|---:|---:|---:|
| SV 检测或候选位点覆盖 | 45 | 15 | 30 |
| 基因型准确性 | 0 | 45 | 20 |
| 相位一致性 | 0 | 10 | 10 |
| 类型/长度/区域稳健性 | 25 | 15 | 15 |
| 断点与等位序列质量 | 15 | 5 | 10 |
| 资源效率 | 10 | 5 | 10 |
| 可靠性与可复现性 | 5 | 5 | 5 |

不输出相位的 genotyper 进入单独的 `genotyping-core` 排行榜，其权重重新归一化并在 score profile 中固定。它不能与 `genotyping phase-aware` 榜混排。

### 17.2 检测分

```text
Detection =
  100 * (
    0.60 * fused_overall_F1
    + 0.40 * macro_F1_across_type_and_length
  )
```

### 17.3 基因型分

正式 panel-genotyping：

```text
Genotype =
  100 * (
    0.50 * exact_GT_concordance
    + 0.35 * three_class_macro_F1
    + 0.15 * genotype_call_rate
  )
```

`positive-sites` 结果只作为补充，不进入完整 genotyping 排行榜。

### 17.4 相位分

```text
Phase =
  100 * (
    0.50 * phased_GT_concordance
    + 0.25 * (1 - switch_error_rate)
    + 0.25 * (1 - flip_error_rate)
  )
```

只在足够 phased truth 与 phased calls 时计算。

### 17.5 稳健性分

对预定义 context strata 的 fused F1 或 genotype macro-F1：

```text
Robustness =
  100 * (
    0.70 * macro_mean
    + 0.30 * lower_decile
  )
```

低分位项防止总体平均掩盖困难区域严重失败。

### 17.6 断点与等位序列分

```text
Resolution =
  100 * (
    0.35 * sequence_similarity
    + 0.25 * size_similarity
    + 0.20 * breakpoint_accuracy
    + 0.20 * haplotype_partial_credit
  )
```

不适用于某一 SV 类型的子项在预定义 type-specific profile 中重新归一化，不能由工具运行后临时选择。

### 17.7 资源分

每个赛道定义 `target` 和 `limit` 资源预算。单项资源得分：

```text
resource_score(x) = 100, x <= target
resource_score(x) = 0,   x >= limit
```

在 target 与 limit 之间按 `log(x)` 线性下降。

```text
Efficiency =
  0.40 * runtime_score
  + 0.30 * peak_memory_score
  + 0.20 * disk_score
  + 0.10 * cpu_hour_score
```

资源分只能在同一硬件指纹、同一数据版本和同一赛道内比较。报告同时保留原始值，避免只看归一化分数。

### 17.8 可靠性分

```text
Reliability =
  0.40 * output_validity
  + 0.30 * successful_completion
  + 0.20 * provenance_completeness
  + 0.10 * repeat_stability
```

一次正式 run 时，`repeat_stability` 只报告、不进入总分；进行重复 benchmark 时才启用该权重，并使用固定 profile。单次运行 profile 会把前三项按 `40:30:20` 重新归一化为 100%，避免可靠性分只累计到 90%。

### 17.9 置信区间

precision、recall、F1、GT concordance 和主要 score component 使用分层 bootstrap 生成 95% CI。bootstrap 单元为 variant cluster，避免把邻近复杂表示当作完全独立事件。

## 18. 报告

主要文件：

```text
results/summary/metrics.long.tsv
results/summary/metrics.json
results/summary/scores.tsv
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
- 总体 PR 和 score；
- SVTYPE/length/context heatmap；
- GT confusion matrix；
- switch/flip errors；
- breakpoint、size、sequence similarity；
- evaluator agreement/upset plot；
- runtime-memory-disk Pareto；
- 失败与不适用工具；
- 可下载的逐变异 review 表。

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
- 工具不支持目标 SVTYPE；
- 无相位输出。

以上状态为 `not_applicable`，不等于失败或零分。

## 20. 测试

### 20.1 单元测试

- 配置和 tool manifest schema；
- BAM/reference 兼容性逻辑；
- candidate panel 正负标签生成；
- truth 信息盲化；
- VCF canonical normalization；
- evaluator parser；
- 融合公式；
- SVBench-100；
- N/A、失败和单引擎边界条件。

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
- evaluator、fusion、score 和 HTML report。

### 20.4 全基因组验收

由于当前本地环境无法访问用户服务器 BAM，全基因组验收在目标 Linux/HPC 上进行：

1. preflight；
2. `caller_only_shared_bam`；
3. 一个端到端 mapping tool；
4. Sniffles2 force、kanpig 和 SVJedi-graph；
5. 一个 external mock 或用户工具；
6. v5.0q 多引擎评测；
7. 资源和 HTML 报告。

## 21. 第一版交付范围

### 21.1 完整实现

- 项目骨架和 schema；
- BAM/reference/truth preflight；
- canonical FASTQ 和 HP FASTQ；
- v5.0q、v0.6 和 v3.6 资源配置；
- Sniffles2、cuteSV、pbsv、SVIM；
- Sniffles2 force、kanpig、SVJedi-graph；
- generic external tool runner；
- blinded genotyping challenge panel；
- canonical 和 evaluator-specific normalization；
- Truvari、Aardvark、vcfdist；
- metric fusion；
- track-specific SVBench-100；
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
7. 外部与内置工具使用相同 evaluator、fusion 和 score。
8. `not_applicable`、`runtime_failed` 和 `invalid_output` 被区分。
9. v5.0q 和 v0.6 结果独立报告。
10. 至少两个 evaluator 支撑核心融合指标。
11. 逐变异 disagreement 可追溯到原始 VCF 记录。
12. 资源指标与 hardware fingerprint 一起报告。
13. 小型 integration test 在无全基因组数据时可执行。
14. README 提供新增自研 genotyper 的最小示例。

## 23. 关键风险与缓解

### v5.0q 仍为 draft

缓解：保留 evaluator 原始结果、分歧惩罚、人工复核表和 truth version。

### 复杂 SV 表示差异

缓解：canonical 与 evaluator-specific views 分离；至少使用两个 haplotype/representation-aware evaluator。

### `0/0` 负位点误标

缓解：只在 benchmark BED 内选择；任一高灵敏 matcher 命中即排除；排除复杂 cluster 和 evaluator-discordant 位点；报告 panel 构建审计表。

### 从 BAM 提取 reads 的信息损失

缓解：报告 hard clipping 和 tag loss；caller-only 使用原 BAM；HP reads 单独提取。

### 不同范式资源不可比

缓解：分赛道排名；使用 track-specific budgets；同时展示原始资源。

### 外部工具任意代码

缓解：推荐固定 digest 容器、受限挂载、独立输出目录和完整 provenance。框架不把外部代码视为可信代码。

## 24. 参考

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
