# F39亚组与缺失统计复现说明

F39脚本生成补充表S37所需的亚组统计及S38所需的逐变量缺失统计。它读取既有派生数据并应用冻结预测公式，不重新拟合预测模型；校准回归仅用于性能评价。输出为CSV、可选Excel及Markdown报告，不包括Word稿件表格排版。

## 来源与适配

代码来自2026年10月17日本地F39补充分析。本次仅替换临时目录链接、调整开头说明及文本编码/换行，统计表达式保持不变。中文命名的主脚本与 `F39_TRIPOD.R` 原始字节相同，因此仅收录一份。

适配前本地文件SHA-256：

| 文件 | SHA-256 |
|---|---|
| `F39_TRIPOD.R` | `bb40c3ef854077f045029ff12aa483199528d3deebc1fa69ed5fc54cfeed8d94` |
| `F39_check_calibration.R` | `851c8d095cbb32c38c4c6b157fffb0540c9f15f2bf2fc35c4cd4a3ab2cc18241` |

仓库适配后的校验值见 `source-manifest.csv`。本次不收录本地生成结果、日志、模型对象或参与者数据。中文分支沿用中文项目子目录，不能直接混用英文分支准备后的目录布局。

## 输入准备

按[复现说明](reproduction.md)运行 `tools/prepare_workspace.py`。两个脚本将复制到 `<project_root>/16.论文修订第三轮/_核查/脚本/`，项目根路径与R库路径由配置迁移，无需创建旧的 `D:/xu_*` 目录链接。

配置后的项目根目录下须具备：

| 输入 | 来源及要求 |
|---|---|
| `output/02_五波数据与质控/waves_cognition_2026-09-13.rds` | 脚本03；标识已核验，含性别、城乡、认知、预测及调查变量 |
| `output/10_筛查指数/screening_index_memory_v3_2026-09-26.rds` | 脚本19；冻结 `model_u` 及开发数据 |
| `external_data/HRS/analysis/HRS_analysis_2026-09-26.rds` | 脚本23，HRS的S0风险集 |
| `external_data/ELSA/analysis/ELSA_analysis_2026-09-26.rds` | 脚本23，ELSA的S0风险集 |
| `external_data/SHARE/analysis/SHARE_analysis_2026-09-26.rds` | 脚本23，SHARE的S0风险集 |

主脚本核对固定系数与保存模型、开发参与者标识、整数评分和历史人数：开发4441人、658事件，宽样本完整病例8296人。这些断言用于复现稿件版本，并非面向任意新数据的通用分析入口。若脚本23通过 `A11_OUT` 将结果写到其他位置，须先按本表准备对应输入。

## 执行顺序

以下为配置后的示例路径，请替换为实际项目根目录：

```powershell
Rscript "E:/memory-reproduction/project/16.论文修订第三轮/_核查/脚本/F39_TRIPOD.R"
Rscript "E:/memory-reproduction/project/16.论文修订第三轮/_核查/脚本/F39_check_calibration.R"
```

输出目录为 `<project_root>/16.论文修订第三轮/F39_TRIPOD补充分析_2026-10-17/`。重复运行会覆盖该目录内同名输出，应使用独立复现工作区。

| 对应内容 | 输出文件 |
|---|---|
| S37开发亚组 | `F39_23a_CHARLS_dev_subgroups.csv` |
| S37外部S0亚组 | `F39_23a_external_S0_subgroups.csv` |
| 亚组合并结果 | `F39_23a_all_subgroups.csv` |
| S38 CHARLS逐阶段缺失 | `F39_20b_CHARLS_missing_by_var.csv` |
| S38外部风险集缺失 | `F39_20b_external_missing_by_var.csv` |
| 校准核查 | `F39_check_calibration_O_vs_E.csv` |
| 汇总与运行记录 | `F39_TRIPOD_report_zh.md`、`F39_TRIPOD_tables.xlsx`（已安装openxlsx时生成）、`F39_run_log.txt` |

校准核查脚本应在主脚本之后运行，会向开发与合并亚组CSV补充 `cal_note` 字段，解释开发表观教育层offset截距接近零的原因，不修改冻结模型。

## 样本与指标口径

CHARLS覆盖性别、教育、年龄段、65岁分组及城乡。既有外部分析文件没有性别和城乡字段，外部亚组仅分析可用的教育和年龄分组；字段不可用与缺失率为零须区分。

S37外部连续模型口径使用 `*_S0_modelComplete` 行；CSV另含 `*_S0_scoreFull` 全评分风险集及整数评分AUC。校准指标来自连续模型预测概率，不是整数评分的校准。全评分风险集行的概率指标会排除缺失预测，而 `n` 仍是全评分样本人数，不能把该人数当作连续模型指标分母。

CSV包含比稿件S37更多的行和指标。选取、排版时须保持队列、亚组、模型/评分及完整病例口径一致。源码补齐与软件检查不等于重新完成端到端统计复现；仍需获授权的数据、正确派生输入及历史版本一致性检查。
