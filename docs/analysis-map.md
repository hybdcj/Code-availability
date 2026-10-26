# 分析入口与版本优先级

文件名前缀是项目历史编号，以下中文说明用于解释其用途。编号不是严格的全流程执行顺序；不得遍历目录并批量运行全部 R 文件。

| 模块 | 主要脚本（相对于仓库根目录） | 依赖/说明 |
|---|---|---|
| 2020 数据清理 | [01_build_2020.R](../analysis/core/01_build_2020.R) | CHARLS 原始模块；生成 2020 基线派生数据 |
| 五波合并 | [03_waves_build.R](../analysis/core/03_waves_build.R) | 2011—2018 原始数据及 2020 派生数据 |
| 个体标识修复 | [33_fix_waves_cognition_idn.R](../analysis/core/33_fix_waves_cognition_idn.R) | 运行在合并之后；临时副本修复，须核验再回填 |
| 开发来源队列 | [04_cohort_change.R](../analysis/core/04_cohort_change.R) | 使用经过标识核验的五波数据 |
| 测量与纵向分析 | `analysis/core/05*`—`16*` | 包括风险因素、轨迹、测量等价、发生率、加权与失访分析；文件名相同数字前缀可对应不同模块 |
| 记忆结局评分开发 | [19_screening_index_memory_v3.R](../analysis/core/19_screening_index_memory_v3.R) | 当前记忆结局主入口；18 是其他认知结局的历史脚本 |
| Bootstrap 内部验证 | [20_internal_validation_bootstrap.R](../analysis/core/20_internal_validation_bootstrap.R)、[21_internal_validation_figures.R](../analysis/core/21_internal_validation_figures.R) | 已构建的评分/模型对象 |
| 外部数据统一 | [23_external_cohort_harmonize.R](../analysis/core/23_external_cohort_harmonize.R) | 需自行取得 HRS/RAND、ELSA、SHARE 数据 |
| 外部验证 | [24_external_validation.R](../analysis/core/24_external_validation.R) | 23 的分析数据与 19 的冻结模型 |
| 外部图表/异质性 | [25_a12_figures_report.R](../analysis/core/25_a12_figures_report.R)、[26_a13_cross_cohort_heterogeneity.R](../analysis/core/26_a13_cross_cohort_heterogeneity.R) | 外部验证输出 |
| 模型冻结与检查 | [27_finalize_model_A14.R](../analysis/core/27_finalize_model_A14.R)、[28_A14_final_freeze_checklist.R](../analysis/core/28_A14_final_freeze_checklist.R) | 开发与验证结果 |
| 敏感性分析 | [29_A15_sensitivity_analyses.R](../analysis/core/29_A15_sensitivity_analyses.R) | 本次收录版本已含行序/插补汇总函数修正；F29/F30 用于更新和核查旧产物 |
| 原始表图导出 | [31_A16_tables.R](../analysis/core/31_A16_tables.R)、[32_A16_figures.R](../analysis/core/32_A16_figures.R) | 属历史导出，须应用后续修正；不能直接当作最终正文图表 |
| 第三轮时段/全人群准备 | `analysis/revisions/F02*`—`F22*` | 保留后续脚本引用的准备表达式；含被替代的结果分支 |
| 筛查指标修正 | [F24_recompute_screening.R](../analysis/revisions/F24_recompute_screening.R) | 引用 F03/F14/F15/F16/F19/F22；脚本参数为第三轮目录；涉及 BDSI/Tu 时还需遵循后续限制 |
| BDSI 修正 | [F26_bdsi_audit.R](../analysis/revisions/F26_bdsi_audit.R) | 需 F22 派生数据和相关原始变量；优先于旧 BDSI 比较 |
| 开发样本不重叠验证 | [F27_nonoverlap_validation.R](../analysis/revisions/F27_nonoverlap_validation.R) | F26 结果、开发对象和队列数据 |
| 结局口径登记 | [F28_outcome_registry.R](../analysis/revisions/F28_outcome_registry.R) | 冻结开发及外部数据 |
| 内部性能行序核查 | [F29_修复内部性能行序.R](../analysis/revisions/F29_修复内部性能行序.R) | 针对既存 A15 输出；不是从零复现的必要前置步骤 |
| 插补区间核查/更新 | [F30_插补性能区间核查.R](../analysis/revisions/F30_插补性能区间核查.R)、[F30_更新插补区间.R](../analysis/revisions/F30_更新插补区间.R) | 先核查生成逐次性能，再更新旧产物；2.5%—97.5% 为经验分位区间 |
| 时段与测量核查 | [F31_时段与测量口径核查.R](../analysis/revisions/F31_时段与测量口径核查.R) | F04 准备表达式及开发模型 |
| 全人群分母修正 | [F32_全人群入组与分母修正.R](../analysis/revisions/F32_全人群入组与分母修正.R) | 引用 F15 准备表达式；整数评分样本不再额外要求 ADL 完整 |
| 修订图表 | [F34_更新正文图件.R](../analysis/revisions/F34_更新正文图件.R)、[F34_补充表构建.R](../analysis/revisions/F34_补充表构建.R)、[F36_纵向图件.R](../analysis/revisions/F36_纵向图件.R) | TIFF 图件及原格式表格；依赖已生成结果和部分排版缓存，见限制说明 |
| TRIPOD亚组与缺失统计（S37/S38） | [F39_TRIPOD.R](../analysis/revisions/F39_TRIPOD.R)，随后运行 [F39_check_calibration.R](../analysis/revisions/F39_check_calibration.R) | 已核验五波数据、冻结开发模型及23生成的外部S0风险集；详见[F39复现说明](f39-reproduction.md) |

路径准备工具会恢复原项目相对布局：`analysis/core/*` → `<project_root>/code/*`；`analysis/revisions/*` → `<project_root>/16.论文修订第三轮/_核查/脚本/*`。这是为了保留原脚本之间的引用关系。
