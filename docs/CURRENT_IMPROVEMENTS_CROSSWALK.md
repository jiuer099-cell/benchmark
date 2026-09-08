# 当前版本改进文档逐项落实表

本表对应 `PGBench_当前版本仍可提升的地方.md`。状态区分“代码已落实”和
“生产数据/正式运行待验证”，防止用尚未取得的服务器资产伪造完成状态。

| 文档项 | 代码状态 | 正式发布状态 |
|---|---|---|
| 2 Addressability | 已实现完整 candidate-status、六类核心状态、完整 panel 分母、ME-F1_ADDRESSABLE 诊断值和主表列 | 等正式五工具运行验证 |
| 3 all-sites / evaluation-query 分离 | 已实现确定性两层 VCF；query 仅含 0/1、1/1；panel allele、唯一 ID、无额外过滤和全部 SHA-256 均硬校验 | 等正式产物验证 |
| 4 简化 ME-F1 | 已冻结 v1.0 唯一公式，三 evaluator 缺一即失败，不重归一，不接受插件分数，分层不参与主分 | 已完成代码门禁 |
| 5 non-reference 解释 | README/报告明确 ME-F1 是 genotype-aware non-reference SV F1；同时报告 P/R、Call Rate、No-call Rate、Addressability；all-site GT macro-F1 仅诊断 | 等正式结果 |
| 6 panel 来源独立性 | 已增加独立 provenance catalogue、source/population/sample/extraction/merge/normalization/filter 字段、source/canonical hash 和 fail-closed gate | 上游 assembly/extraction 元数据尚未提供，因此诚实保持 `unverified_information` |
| 7 Allowed-information | 已增加全局 allow/forbid policy、插件声明、resolved input 映射、未声明/禁用资源拒绝及 hash | 等正式 resolved-input 审计 |
| 8 Tuning | 已增加三种合法模式、parameter source/hash、冻结时点和禁止查看 HG002 truth/final ME-F1 的门禁 | 等正式参数清单 |
| 9 evaluator disagreement | 已报告 evaluator SD/range、三种 leave-one-evaluator-out sensitivity | 等五工具结果 |
| 10 runtime 拆分 | 已新增严格 stage summarizer；必须同时提供 one-time build 与 per-sample 的 wall/CPU/peak RAM，另计 index disk 和 N=1/10/100/1000 摊销；combined-only 会拒绝 | 各工具阶段级 benchmark 尚待服务器正式测量，release gate 保持关闭 |
| 11 coverage uncertainty | 10x/20x/30x 强制相同至少三个 seed；输出 mean、sample SD、Student-t 95% CI、min/max；full 单次运行 CI 明确为 undefined | 等 coverage matrix 正式运行 |
| 12 paired comparison | 已使用相同 genomic-block bootstrap draws 输出任意工具对 delta ME-F1 和 95% CI | 等五工具共同结果 |
| 13 多样本 | 按文档属于 v2；v1 明确限定为 HG002 sample/configuration-specific | 未扩张当前主轨 |
| 14 complex/sex/population tracks | 按文档属于后续独立 tracks；当前代码继续只允许 chr1-22、simple biallelic DEL/INS、50–10,000 bp、短读长 | 未混入 v1 主榜 |
| 15 结果结构 | 跨工具主表顺序为 Tool、ME-F1、三 evaluator F1，再列 Addressability/Call Rate/各失败计数；单工具 TSV 同样把 ME-F1 放在所有分数首列 | 等正式渲染 |
| 16 发布门槛 | release gate 已覆盖 evaluator semantics、addressability、evaluation-query、panel provenance、family LOO、information/tuning、三 evaluator、strata/coverage/CI、runtime split、provenance | 任何生产项未通过均不得发布 verified leaderboard |

主分契约只有：

```text
ME-F1 = (F1_Truvari + F1_Aardvark-GT + F1_vcfdist) / 3
```

`ME-F1_DEL`、`ME-F1_INS`、`ME-F1_LEN_*`、`ME-F1_AF_*`、
`ME-F1_<REGION>`、`ME-F1_10X/20X/30X/FULL` 均携带
`role: explanatory` 和 `affects_primary_score: false`。可表示子集只使用
`ME-F1_ADDRESSABLE`，其角色为 `diagnostic`。
