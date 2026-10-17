# 复现准备与运行

## 1. 环境

整理时检查环境为 Windows、R 4.6.1。R 包版本记录于 `environment/r-packages.csv`，这是整理当日的已安装版本快照，不代表每次历史分析实际使用的版本，也不是锁定全部传递依赖的 renv.lock。`install_dependencies.R` 安装缺少的直接依赖，不强制降级现有包。若需版本一致，请依据快照另建隔离环境。

Python 辅助工具需要 Python 3.10 或以上并使用标准库；测量汇总脚本另需 openpyxl，版本见 requirements.txt。中文路径及 UTF-8 编码需可用；原 R 脚本包含 Windows 中文 locale 设置。TIFF 绘图依赖 ragg，部分历史图件使用 SimSun 字体；其他操作系统需调整 locale/字体。不要未经验证声称跨平台运行一致。

## 2. 配置路径

复制 `config.example.json` 为 `config.local.json`，修改为本机绝对路径。建议 project_root 使用全新的复现工作目录；scratch_root 放中间数据；external_staging_root 放解压的外部原始文件；r_library 为本机实际 R 库目录。

```powershell
python tools/prepare_workspace.py --config config.local.json
```

工具只复制文本脚本和变量映射说明，替换四组已知本机路径前缀，并写入迁移清单。不会复制、下载或上传数据，不执行分析。若目标已有不同内容的源文件则停止，不覆盖原项目。`source-manifest.csv` 校验的是仓库中的原始字节；准备后的脚本因路径替换而有不同哈希。

## 3. 数据与派生输入

参照 `data-access.md` 准备授权数据。`input-references.csv` 列出读取表达式及行号，变量拼接、动态读取和报告缓存仍需查看脚本上下文。源脚本含日期固定文件名，不要随意改名；也不要将不同时期或口径不同的 RDS 替换为相同别名。

主开发路径为 01 → 03 → 标识核验/33 → 04 → 19 → 20/21；其他描述性和纵向模块按输入依赖单独执行。33 使用 scratch 中的 `waves_cognition_orig.rds`，产生 `waves_cognition_fixed.rds`，需检查日志和标识一致性后将修复结果回填至五波数据位置，再执行依赖五波数据的模块。该回填不由本仓库工具自动执行。

外部路径为 23 → 24 → 25/26，再结合已有开发结果执行冻结与敏感性模块。此顺序只是模块导航，不能取代各脚本的输入检查。已有冻结结果的审计与从原始数据重新拟合是两种不同操作，应保存新旧结果的版本标识，避免覆盖既有论文结果。

第三轮脚本还使用 scratch_root/tmp 中的短文件名，常见映射如下：

| 短名 | 来源 |
|---|---|
| `waves_cognition.rds` | 核验过的 `output/02_五波数据与质控/waves_cognition_2026-09-13.rds` |
| `cohort_change.rds` | `output/04_纵向队列/cohort_change_2026-09-13.rds` |
| `a07_memory_v3.rds` | `output/10_筛查指数/screening_index_memory_v3_2026-09-26.rds` |
| `ext_hrs.rds` / `ext_elsa.rds` / `ext_share.rds` | 对应外部队列主分析 RDS（23 的输出） |
| `ext_*_sens.rds` | 对应外部队列敏感性风险集（23 的输出） |
| `F16_base_*.rds` | F16_build 的准备结果 |
| `F22_d.rds` | F22a 的准备结果；其旧比较结果不作为有效原模型验证 |
| `F26_bdsi_audit.rds` | F26 的审计结果 |

原始模块的临时别名还包括 `w1_hs.dta`、`w2_hs.dta`、`w3_hs.dta`、`w1_biomarker.dta`、`w3_biomarker.dta`、`w4_cognition.dta`、`w4_hsf.dta`、`w5_hsf.dta`。这些历史别名的完整自动准备链未收录，须根据变量清单和源波次确认对应模块后在本地准备，不能仅凭近似文件名替代。F29/F30 等更新脚本要求旧输出已存在，首次分析不应直接运行这些历史产物修复步骤。

## 4. 分阶段运行与检查

以下为路径配置后的单模块调用形式，不代表已准备好其输入：

```powershell
Rscript "E:/memory-reproduction/project/code/19_screening_index_memory_v3.R"
Rscript "E:/memory-reproduction/project/16.论文修订第三轮/_核查/脚本/F24_recompute_screening.R" "E:/memory-reproduction/project/16.论文修订第三轮"
```

从项目根目录启动原分析时，先确认输入、缺失规则、输出目录及对应阶段；检查样本量/事件数、模型与评分样本差异、阈值与同分处理。脚本含写出和更新产物的操作，正式结果应在独立目录运行并与冻结版本逐项比较。不要把所有输出或 participant-level 模型对象提交到 GitHub。

## 5. 本次验证边界

只校验源文件完整性、R/Python 语法、路径准备工具和仓库内容，不运行原始统计分析、不修改论文、不声称全新环境端到端复现成功。详见 `validation.md` 及 `known-limitations.md`。
