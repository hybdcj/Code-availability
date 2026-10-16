# 数据来源与本地布局

研究使用 CHARLS、HRS（含 RAND 整理数据）、ELSA 及 SHARE/easySHARE。个体数据须经对应提供方获取并按其条件使用；本仓库不分发原始数据、派生个体数据、含参与者行的模型对象或临时缓存。

| 队列 | 本地分析入口 | 输入布局 |
|---|---|---|
| CHARLS | `analysis/core/01_build_2020.R`、`03_waves_build.R` | `<project_root>/data/raw/` 下各波模块；具体大小写/文件名见脚本 |
| HRS / RAND HRS | `analysis/core/23_external_cohort_harmonize.R` | `<external_staging_root>/hrs/H_HRS_d.dta` 和 `randhrs1992_2022v1.dta` |
| ELSA | 同上 | `<external_staging_root>/elsa/gh_elsa_h.dta` |
| SHARE / easySHARE | 同上 | `<external_staging_root>/share/easySHARE_rel9-0-0.dta` |

文件名记录本研究脚本的预期输入，不表示这些数据可免申请下载，也不表示其他发布版本可以直接替换。实际数据版本、变量字典及权限以所使用提供方资料为准。主分析外部派生数据写入 `<project_root>/external_data/`。

`docs/input-references.csv` 是代码读取表达式的静态索引，便于追踪输入。它不含任何数据行。第三轮临时别名及手工准备环节见 `reproduction.md`。

`.gitignore` 排除了常见数据/模型/图表文件，但仅靠忽略规则不能替代提交前检查；包括 CSV 在内的任意格式都可能含个体信息。仓库当前 CSV 仅用于源文件清单、输入表达式索引和依赖版本记录。
