# 社区中老年记忆低分风险评分：分析代码

Analysis code for a memory low-score risk model in community-dwelling middle-aged and older adults, developed in CHARLS and evaluated in HRS, ELSA and SHARE.

本仓库收录论文相关分析源代码、第三轮修正脚本、依赖版本记录和复现准备工具。当前为投稿前私有整理版；作者计划投稿后公开。本仓库不包含参与者数据、模型对象、论文文件或生成的图表。

**开始前请先读 [分析入口](docs/analysis-map.md) 和 [复现说明](docs/reproduction.md)。** 源代码保留原有文件名与内容，内部编号仅用于追溯，不是论文中的研究分组。仓库不是可直接一键重算全部论文结果的流水线；数据授权、路径配置、派生数据和部分历史中间产物需另行准备。

## 目录

| 路径 | 内容 |
|---|---|
| `analysis/core/` | 数据整理、纵向分析、评分开发、内部与外部验证、敏感性分析 |
| `analysis/revisions/` | 第三轮补充分析、筛查指标和统计口径修正、图表源代码 |
| [docs/analysis-map.md](docs/analysis-map.md) | 主要脚本路径、作用、依赖和版本优先级 |
| [docs/data-access.md](docs/data-access.md) | 数据来源及本地输入布局 |
| [docs/reproduction.md](docs/reproduction.md) | 环境设置、路径迁移和分阶段执行说明 |
| [docs/known-limitations.md](docs/known-limitations.md) | 历史分析限制、已修正内容和未完成的复现环节 |
| [docs/input-references.csv](docs/input-references.csv) | 源代码中读取数据/脚本的位置索引（静态扫描，非完整依赖图） |
| `environment/` | R 直接依赖、版本快照、Python 依赖和安装脚本 |
| `tools/` | 源文件校验、工作目录准备和 R 语法检查工具 |
| [source-manifest.csv](source-manifest.csv) | 源脚本相对原项目路径、文件大小和 SHA-256 |

## 最快开始

```powershell
git clone https://github.com/hybdcj/Data-availability.git
cd Data-availability
python tools/verify_sources.py
Copy-Item config.example.json config.local.json
# 编辑 config.local.json，使用一个专门的复现工作目录和本机实际 R 库路径
python tools/prepare_workspace.py --config config.local.json
Rscript environment/install_dependencies.R
python -m pip install -r environment/requirements.txt
Rscript tools/check_syntax.R .
```

以上命令校验/准备代码，不拟合模型。`prepare_workspace.py` 将脚本复制到配置的工作目录，并替换已知绝对路径；默认拒绝覆盖已有不同内容的文件。之后须根据[复现说明](docs/reproduction.md)准备获授权的数据和派生输入，按模块执行。

## 研究口径

主开发分析使用 CHARLS 2011—2015 年资料；连续模型与整数评分须分别报告。连续模型含年龄、教育、抑郁及 ADL；整数评分的非零条目为年龄、教育和抑郁。外部验证使用冻结模型与评分，并区分不同结局阈值定义。具体纳入、缺失处理及敏感性分析以对应脚本和论文方法为准。

第三轮修正优先于历史导出：筛查指标见 F24，BDSI 实现见 F26，开发样本不重叠验证见 F27，内部性能行序见 F29，插补区间解释见 F30，全人群评分分母见 F32。旧的 Tu 2026 比较保留为历史探索，不能视为原模型的完整外部验证。

## Data and code availability

Individual-level data are not distributed in this repository. Access to CHARLS, HRS/RAND HRS, ELSA and SHARE/easySHARE data must be obtained from the respective data providers under their applicable conditions. This repository contains analysis source code and reproduction documentation. The repository is currently private and is intended to be made public after manuscript submission.

本次整理完成的验证见 [validation.md](docs/validation.md)。尚未在全新环境中重跑整套模型，依赖快照也不等同于原始分析时的完整锁定环境。未擅自选择开放源代码许可证；公开前由作者确定许可和正式论文引用信息。
