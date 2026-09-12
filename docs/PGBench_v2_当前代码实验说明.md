# PGBench v2 当前代码实验说明

> 文档状态：实现说明，不是已发布的 benchmark 结果。  
> 更新：2026-09-12。  
> 当前阶段：外接 Adapter 与 SR/HiFi 双通道代码已实现；PanGenie P0 已在服务器通过；五工具 production 重跑和不可变 release 尚未发布。

## 1. 实验目的

PGBench v2 用统一、可追溯的流程比较结构变异（SV）固定 panel
genotyping 工具。比较对象是候选 panel 中的基因型判断，而不是 de novo
SV discovery。所有评分由 PGBench Core 完成；工具 adapter 只能产生其原生
调用结果，不能自行产生 ME-F1、排名或修改评分规则。

实验的核心目标是保证：

1. 同一通道内的工具使用相同证据数据、参考基因组、候选 panel 和评分规则；
2. 工具无法表达某个候选位点时，该位点不会从分母消失；
3. 使用 population/pangenome 信息的工具均从同一冻结生物学来源派生原生资产；
4. 短读长与长读长分别评测、分别发布，不进行混合排名；
5. 新工具通过新增外接 adapter 接入，而不是向 Core 增加工具专用分支。

## 2. Benchmark 范围

### 2.1 固定 SV 范围

当前配置将正式比较范围冻结为：

| 项目 | 约束 |
| --- | --- |
| 参考基因组 | GRCh38 |
| 染色体 | chr1–chr22 |
| SV 类型 | DEL、INS |
| 长度 | 50–10,000 bp |
| 位点形式 | biallelic、sequence-resolved |
| 比较任务 | fixed-panel genotyping |
| 候选 universe | 冻结的 **18,164 canonical scoring universe**；所有工具须确定性投影到该评分题库 |
| 主评分 | ME-F1 |

候选 panel 在工具运行前冻结。候选 ID、panel allele、reference、truth 和
evaluation BED 都进入 provenance；工具不能替换 panel 或选择自己的评分分母。

### 2.2 两条独立通道

| 通道 | benchmark track | 固定测序证据 | 配置 |
| --- | --- | --- | --- |
| SR | `short_read_fixed_panel_genotyping` | HG002 Illumina PE FASTQ | `config/config.unified-tools.example.yaml` |
| HiFi | `long_read_fixed_panel_genotyping` | HG002 PacBio HiFi FASTQ | `config/config.hifi.example.yaml` |

`benchmark_contract.cross_track_ranking` 固定为 `false`。每次运行只选择一条
track；adapter 的 `capabilities.read_class` 与该 track 不匹配时标为
`unsupported_for_selected_track`，不进入 DAG、评分集合或零分统计。

CLR 与 ONT 不属于当前 HiFi 通道。现有五个生产 adapter 目前均声明为 SR；
HiFi 配置中的 `example_hifi_adapter` 是接口/流程测试 fixture，不是可以用于
正式排行榜的科学工具。正式 HiFi release 需要接入真实、可复现且原生支持
PacBio HiFi 的 adapter。

## 3. 数据和冻结生物学来源

### 3.1 样本与测序数据

- 样本：HG002。
- SR：同一份 HG002 Illumina paired-end R1/R2 FASTQ 提供给 SR adapter。
- HiFi：同一份 HG002 PacBio HiFi FASTQ 提供给 HiFi adapter。
- 正式输入会记录 FASTQ 的内容哈希、read count、read bases、coverage、
  downsampling seed 和 library/source evidence ID。
- 当前仓库的示例配置只提供资源路径，不包含完整 production FASTQ、GRCh38
  和 population VCF；因此示例配置通过 schema/semantic validation 不等于
  production 已经执行。

### 3.2 Frozen Haplotype Source Bundle

所有需要 population 或 pangenome 信息的 adapter 必须声明：

```yaml
native_assets:
  uses_population_or_pangenome_information: true
  frozen_haplotype_source_bundle: HG002_LOO_HPRC_GRCh38_SV_v1
```

对应配置中的 bundle 固定：

```yaml
id: HG002_LOO_HPRC_GRCh38_SV_v1
source_cohort_id: HG002_LOO_HPRC_GRCh38_SV_v1
release: HPRC_GRCh38_SV_v1
reference: grch38
target_family_excluded: true
derivation_policy: native_assets_must_be_derived_from_this_bundle
content_lock: frozen-haplotype-source.lock.yaml
```

HG002、其父母及别名样本（HG002/NA24385、HG003/NA24149、HG004/NA24143）从
population source 中排除。这里冻结的是**生物学来源**：cohort、release、
reference、haplotype roster 和 family exclusion 必须相同。`content_lock` 还封存
GFA/GBZ、population VCF、sample roster、haplotype roster、family-exclusion manifest
和 reference 的 SHA256，以及构建这些资产的软件名、版本与 SHA256；不同工具的原生表示
可以不同，例如 GBZ、phased VCF、PGIN 或 tool-native index 不要求字节相同。

不使用 population haplotype prior 的 adapter（如 shared-BAM 型工具）不应被
强行给予 population prior。

## 4. Adapter 架构与当前工具

### 4.1 统一 Adapter API

每个 adapter 在 `plugins/<adapter>/tool.yaml` 中声明：

- `interface_version: 2`；
- `source: external`；
- read class 与可支持 technology；
- 所需输入、可选输入、billable stages；
- 是否消费 Frozen Haplotype Source Bundle；
- runner、环境、rule registry、sandbox backend；
- 固定 all-sites 输出契约。

Core 的泛化执行器为：

```text
workflow/modules/generic_external/Snakefile
workflow/scripts/pgbench_exec.py
```

工具私有资产由注册项的 `adapter_assets` 传入，键会映射为
`adapter_asset.<name>`。Core 只验证、哈希、只读挂载和传递这些资产名称，不应
理解具体工具语义。PanGenie 的 phased panel、biallelic panel 和 converter
目前即通过这种机制传递。

关键边界由 `tests/unit/test_adapter_boundary.py` 自动检查：核心 Snakefile、
generic external module、config validator 和 executor 中出现当前五工具的具体
身份即失败。该门禁防止“接入第六工具时修改 Core”的回归。

### 4.2 当前工具表

| 工具 | adapter ID | 通道 | 主要输入 | population/pangenome 原生信息 | 状态 |
| --- | --- | --- | --- | --- | --- |
| PanGenie | `pangenie` | SR | Illumina PE | Frozen Bundle 派生的 private phased/biallelic context | P0 已在服务器通过；production 待重跑 |
| vg Giraffe + call | `vg` | SR | Illumina PE + graph assets | Frozen Bundle 派生图资产 | production 待跑 |
| Paragraph | `paragraph` | SR | shared BAM/BAI | 不使用 population prior | production 待跑 |
| GraphTyper2 | `graphtyper2` | SR | shared BAM/BAI | 不使用 population prior | production 待跑 |
| Varigraph | `varigraph` | SR | Illumina PE | Frozen Bundle 派生 native assets | production 待跑 |

BayesTyper 也保留为 external adapter 示例，但不属于用户要求的五工具正式
重跑集合。

### 4.3 Adapter eligibility

adapter 状态与分数严格分离：

| 状态 | 含义 | 是否计为 ME-F1=0 |
| --- | --- | --- |
| `eligible` | 与本次 track、technology、mode 兼容并可调度 | 否，正常评分 |
| `unsupported_for_selected_track` | 工具本身不支持本通道 | 否，不进入评分集合 |
| `adapter_failure` | 本应可运行但执行或输出契约失败 | 否，不伪装为工具性能分数；保留失败审计 |

## 5. 工具运行期的 truth 隔离

正式运行 (`execution_purpose=formal`) 拒绝 `sandbox_backend: none`，只能使用
`bwrap` 或 Apptainer。执行器为正式 attempt 写入如下 isolation contract：

```yaml
network: disabled
inputs_read_only: true
truth_visible: false
evaluator_results_visible: false
```

在 bwrap 下，runner 可看到的路径只包括：自身 adapter 目录、已声明的只读
输入、必要系统运行时和自身可写 attempt work directory；使用 `--unshare-net`。
Apptainer 使用 `--containall --cleanenv --no-home --network none`，只绑定 adapter
目录与声明输入。truth VCF、GIAB 文件、hidden truth ledger 和 evaluator 输出均
不是 adapter 输入，不会被挂载。

adapter 的信息契约还要求：

```yaml
target_truth_used_for_calling: false
target_assembly_used: false
target_family_genotypes_used: false
target_specific_external_callset: false
```

## 6. 调用、全位点账本和固定分母

adapter 产生原始 VCF 后，Core 完成以下步骤：

```text
adapter raw VCF
  -> normalization
  -> canonical allele linking
  -> all-sites VCF + candidate-status TSV
  -> evaluation-query VCF
  -> evaluator ledgers
  -> Core ME-F1
```

`all-sites.vcf.gz` 和 `candidate-status.tsv` 保留所有 canonical candidate。常见
状态包括：

- `addressable_called`
- `explicit_no_call`
- `unsupported_representation`
- `linking_failure`
- `missing_output`
- `ambiguous_mapping`

其中 addressability 只是诊断指标。工具不能表示某个 canonical candidate 时，
该候选仍留在 all-sites ledger；不会因为某个工具“只 address 了较小子集”而缩小
主评分分母。

`UNSCORABLE` 只在 tool-independent 的 panel/truth 链接阶段确定：当 truth event
在任何工具运行前就无法无歧义映射到冻结 canonical panel 时，才可标为
`UNSCORABLE`。某个特定工具的表达能力不足必须保留为 adapter/candidate status，
不能转换为 `UNSCORABLE`。

## 7. Truth、评分和 ME-F1

### 7.1 Truth 使用时机

truth 只在 Core 的 challenge construction、hidden ledger 和 evaluation 阶段使用，
不提供给 adapter calling 阶段。候选 panel 对 adapter 是 blinded；truth label 不在
candidate VCF 中暴露。

### 7.2 三评估器统一判断层

Core 固定使用 Truvari、Aardvark-GT、vcfdist 的结果，并将其转换到统一的
TP/FP/FN/GT judgement layer。主评分为：

```text
ME-F1 = (F1_Truvari + F1_Aardvark-GT + F1_vcfdist) / 3
```

三项均必须存在，不对缺失 evaluator 重新归一化。rank 不在 adapter 或中间 metrics
中生成；SR 和 HiFi 的 release 内可各自展示排序，但没有跨 track 总排名。

### 7.3 Tuning 隔离

`config/tuning_policy.yaml` 明确禁止：

```yaml
tune_on_hg002_truth: true
tune_on_final_me_f1: true
choose_parameters_after_viewing_final_test_score: true
```

允许的 parameter source 是官方默认、官方文档推荐配置，或在独立 validation
sample 上完成并在 HG002 前冻结的参数。每次运行冻结 parameter hash、tool
manifest hash、resolved-input hash 与 score profile hash。

## 8. 运行资源和性能记录

每个 adapter attempt 记录：tool version、runner hash、plugin root hash、完整命令、
threads、memory、timeout、cache policy、sandbox backend、开始/结束时间、输入哈希
和输出哈希。Core 的 provenance 还固定 config snapshot、random seed、环境定义、
上游 rule manifest、reference、pangenome manifest 和 evaluator profile。

性能比较在 production 时应区分：

| 类别 | 需记录指标 |
| --- | --- |
| 一次性 build/index | wall time、CPU time、peak RAM、磁盘占用 |
| per-sample calling | wall time、CPU time、peak RAM、threads、磁盘占用 |

当前代码已经记录 rule/tool benchmark 与资源配置；生产 release 前仍需在目标服务器
上固定 CPU/hardware fingerprint、容器或 conda lock、线程与内存上限，并审计
build/index 与 per-sample 阶段的分离。

## 9. Adapter Conformance Suite

synthetic fixture 位于：

```text
tests/fixtures/synthetic/
```

adapter 接入前至少必须通过：

1. DEL 与 INS 的 0/0、0/1、1/1；
2. `./.` no-call；
3. unsupported representation；
4. ambiguous mapping；
5. multi-allelic 到 canonical biallelic projection；
6. all-sites completeness；
7. 无 truth 泄露、无 adapter score/rank；
8. SR/HiFi capability 与输入证据匹配。

当前本地回归命令：

```powershell
.\.venv\Scripts\python.exe -m pytest -q
python workflow/scripts/validate_config.py --config config/config.unified-tools.example.yaml --repo-root .
python workflow/scripts/validate_config.py --config config/config.hifi.example.yaml --repo-root .
```

最近一次代码回归结果为 `205 passed, 4 skipped`。skip 原因是 Windows symlink
权限和未安装的可选 `edlib`，并非测试失败。

## 10. 当前状态与下一步

### 已完成

- 外接 Adapter API v2；
- 五个现有工具迁移为 external adapter；
- SR/HiFi 独立 track、独立评分边界；
- Frozen Haplotype Source Bundle 约束；
- Core all-sites、candidate ledger、统一 ME-F1；
- 正式 sandbox truth isolation；
- adapter/Core 边界的自动测试；
- PanGenie P0 已在服务器通过。

### 尚未发布为正式结果

- 五工具在完整正式 SR 资源上的 production 重跑；
- 真实 HiFi genotyper adapter 的接入与 HiFi production run；
- SR/HiFi 各自的不可变 release manifest 和正式 leaderboard；
- production 硬件、container/conda lock、runtime fairness 的最终封存。

因此当前代码可以用于完成正式运行准备、adapter 接入和 smoke/contract 验证，但
不得把现有示例配置或 synthetic 输出表述为五工具正式 benchmark 结果。

## 11. 建议的正式执行顺序

1. 在服务器确认冻结 SR/HiFi FASTQ、GRCh38、population source、truth、BED 和
   原生资产全部存在且哈希正确；
2. 以已通过 P0 的 PanGenie native context 为基线，执行五个 SR adapter；
3. 对每个 adapter 执行 validate、synthetic conformance、production calling、
   all-sites/ledger 审计和统一 Core 评分；
4. 固定 runtime/hardware/container provenance；
5. 封存 `PGBench-SR-Illumina-v1.0`，不覆盖旧结果；
6. 接入真实 HiFi adapter 后重复同一流程，封存独立的
   `PGBench-LR-HiFi-v1.0`；
7. 仅在各自 release 内展示排名，禁止生成 SR+HiFi 总排名。
