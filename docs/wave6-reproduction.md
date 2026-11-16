# 第六期时间迁移与截距更新复现说明

本扩展覆盖2018年至2021/2022/2023年的质控、事件审计、冻结模型时间迁移及
仅截距更新的内部验证。需要已授权的CHARLS数据及既有冻结开发产物；
本流程不从全部历史原始数据重建整篇论文。

## 环境及工作目录

本次运行环境为Windows、Python 3.14及R 4.6.1；R直接依赖见environment目录，
包括pROC。其他Python版本需要满足固定包版本的要求，未在本轮测试。
在仓库根目录按[主复现说明](reproduction.md)设置config.local.json，
将project_root指向仓库以外的专门目录，再运行：

```text
python -m pip install -r environment/requirements-wave6.txt
python tools/verify_sources.py
python tools/prepare_workspace.py --config config.local.json
```

准备工具只复制代码与方案摘要。依赖版本是本次检查快照，不是完整传递依赖锁文件。
在配置的project_root下准备以下获授权输入：

| 输入 | 相对project_root的路径 |
| --- | --- |
| 2018认知与样本模块 | data/raw/wave4_2018/Cognition.dta；data/raw/wave4_2018/Sample_Infor.dta |
| 第六期年度认知与人口模块 | data/raw/wave6_2021_2023/{2021,2022,2023}/Cognition_Depression.dta及Demographics.dta |
| 第六期样本模块 | data/raw/wave6_2021_2023/2023/Sample_Infor.dta |
| 已核验五波认知数据 | output/02_五波数据与质控/waves_cognition_2026-09-13.rds |
| 开发来源队列 | output/04_纵向队列/cohort_change_2026-09-13.rds |
| 冻结模型对象 | output/12_模型定稿/A14_final_model_2026-09-26.rds |

历史RDS文件名在output下必须各有且仅有一份匹配文件。
冻结对象须含$model及其开发数据（4,441人、658例事件）。
事件导出脚本直接从冻结R对象生成全精度系数CSV，无须手工另建系数文件。
原始数据、派生队列和拟合模型对象不随仓库分发。

## 按序执行

先将当前工作目录切换至配置的project_root，随后按以下顺序执行，
任一步失败均应停止并排查，不继续使用旧输出：

```text
Rscript code/export_wave6_qc_inputs_20261113.R
python code/wave6_qc_20261113.py
python code/verify_wave6_qc_20261113.py
Rscript code/export_wave6_event_inputs_20261114.R
python code/wave6_event_audit_20261114.py
python code/validate_wave6_20261114.py
Rscript code/verify_wave6_validation_20261114.R
python code/validate_wave6_intercept_update_20261116.py
Rscript code/_checks/verify_wave6_intercept_update_20261116.R
```

输出目录依次为output/wave6_qc_20261113、output/wave6_events_20261114、
output/wave6_validation_20261114、output/wave6_intercept_validation_20261116。
其中部分文件含个体标识、协变量、预测或抽样成员信息，须保留在本地仓库外。

## 预期结果与解释

不重叠的2021/2022/2023年度组分别4,041/1,188/672人，事件785/269/160例。
四年组1,188人、212个社区，原冻结模型AUC约0.685815，
截距更新量约0.499309、校准斜率约1.033425。

截距校正按社区有放回抽样10,000次（种子20261116），每次固定原预测变量系数、
仅重估截距，再计算原样本测试损失减抽样训练损失，平均后加到表观Brier上。

| 指标 | 预期值 |
| --- | ---: |
| 原冻结模型Brier | 0.166192323801 |
| 更新后表观Brier | 0.160236888712 |
| 乐观度校正增量 | 0.000297193034 |
| 更新后校正Brier | 0.160534081746 |
| 校正更新后减原模型Brier | -0.005658242056 |
| 校正增量Monte Carlo标准误 | 0.000070513239 |

Monte Carlo标准误仅反映有限次模拟误差，不是抽样标准误或置信区间；
本轮未计算校正估计及其差值的抽样区间。
这是时间验证样本中截距更新程序的内部验证，不是更新模型的独立样本验证。
截距平移不改变同一样本预测排序及AUC。
词表、学习次数、生存和应答选择、日期代理等限制仍然存在，
三/五年结果不能视为四年风险估计。

## 来源及中英文对应

9个分析/核验脚本与3份中文方案摘要纳入source-manifest.csv。
[来源记录](wave6-source-provenance.csv)列出原脚本哈希与公开适配。
中文分支与[英文发布版本](https://github.com/hybdcj/Code-availability/commit/342fc25aec1a3fe123ea26ebd249affba696426b)
使用相同统计实现；仅模型系数目录/文件名按中文布局适配，方案摘要按中文提供。
文件名日期标识分析版本；方案摘要的公开发布不构成前瞻性注册。
实际核验范围见[检查记录](validation.md)，完整仓库限制见[限制说明](known-limitations.md)。
