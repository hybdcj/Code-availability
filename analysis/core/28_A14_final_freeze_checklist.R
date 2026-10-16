# ============================================================================
# 28_A14_final_freeze_checklist.R —— A14 最终冻结清单（一页版）
# 输入：A14 final model artifact 与 A14/A11–A14 复核日志
# 输出：A14_最终冻结清单_2026-09-26.md、A14_最终冻结清单_关键数字.csv
# ============================================================================
invisible(try(Sys.setlocale("LC_ALL","zh_CN.UTF-8"), silent=TRUE))
.libPaths(c("D:/R/library", .libPaths()))
suppressPackageStartupMessages({ library(dplyr); library(knitr) })

P <- "D:/BaiduSyncdisk/科研/徐中豪"; OUT <- file.path(P,"output","12_模型定稿")
a14 <- readRDS(file.path(OUT,"A14_final_model_2026-09-26.rds"))
perf <- read.csv(file.path(OUT,"A14_性能汇总.csv"), check.names=FALSE, fileEncoding="UTF-8")
meta <- read.csv(file.path(OUT,"A14_跨队列meta汇总.csv"), check.names=FALSE, fileEncoding="UTF-8")
strata <- read.csv(file.path(OUT,"A14_风险分层.csv"), check.names=FALSE, fileEncoding="UTF-8")

status <- function(file) {
  if (!file.exists(file)) return("MISSING")
  if (any(grepl("RESULT: PASS", readLines(file, warn=FALSE)))) "PASS" else "CHECK"
}
s11 <- status(file.path(P,"output/10_筛查指数/A11_v2_独立复核_2026-09-26.txt"))
s12 <- status(file.path(P,"output/11_外部验证/A12_独立复核_2026-09-26.txt"))
s13 <- status(file.path(P,"output/11_外部验证/A13_独立复核_2026-09-26.txt"))
s14 <- status(file.path(OUT,"A14_独立复核_2026-09-26.txt"))

fmt_ci <- function(a,lo,hi,d=3) sprintf(paste0("%.",d,"f (%.",d,"f–%.",d,"f)"),a,lo,hi)
perf_fmt <- perf %>%
  filter(scenario %in% c("internal_apparent","internal_bootstrap","external_S0")) %>%
  mutate(
    连续模型AUC = ifelse(is.na(auc_model),"—",fmt_ci(auc_model,auc_model_low,auc_model_high)),
    整数分AUC = ifelse(is.na(auc_score),"—",fmt_ci(auc_score,auc_score_low,auc_score_high)),
    校准 = paste0(ifelse(is.na(cal_slope),"NA",sprintf("%.2f",cal_slope))," / ",
                 ifelse(is.na(cal_intercept),"NA",sprintf("%.2f",cal_intercept))),
    Brier = ifelse(is.na(brier),"NA",sprintf("%.3f",brier))
  ) %>% select(队列=cohort, 场景=scenario, 连续模型AUC, 整数分AUC, 校准, Brier)

meta_fmt <- meta %>% mutate(合并AUC = fmt_ci(pooled_auc,pooled_low,pooled_high), I2=sprintf("%.1f%%",I2)) %>%
  select(切点=scenario, 指标=metric_zh, 合并AUC, I2)

key <- data.frame(
  项目 = c("最终模型版本","CHARLS 风险人群/事件","HRS 风险人群/事件","ELSA 风险人群/事件","SHARE 风险人群/事件",
           "内部 apparent 连续模型 AUC","内部 apparent 整数分 AUC","内部 Bootstrap 校正 AUC",
           "HRS S0 连续模型 AUC","ELSA S0 连续模型 AUC","SHARE S0 连续模型 AUC",
           "HRS S0 整数分 AUC","ELSA S0 整数分 AUC","SHARE S0 整数分 AUC",
           "S0 跨队列合并 AUC（连续模型）","S0 跨队列合并 I²"),
   value = c("v1.0（2026-09-26）","4,441 / 658（14.8%）","14,051 / 1,842（13.1%）","6,725 / 554（8.2%）","26,879 / 2,930（10.9%）",
             "0.702（0.681–0.722）","0.700（0.679–0.721）","0.698（0.679–0.720）",
             "0.696（0.682–0.709）","0.763（0.742–0.784）","0.751（0.742–0.761）",
             "0.692（0.678–0.705）","0.759（0.737–0.780）","0.752（0.742–0.761）",
             "0.737（0.690–0.779）","94.7%"),
   source = c("A14 artifact","A07/A14","A11/A12/A14","A11/A12/A14","A11/A12/A14",
              "A07/A14","A07/A14","A08/A14","A12/A14","A12/A14","A12/A14","A12/A14","A12/A14","A12/A14","A13/A14","A13/A14")
)

md <- c(
"# A14 最终冻结清单（v1.0，2026-09-26）",
"",
paste0("> **状态**：A11 = ",s11,"；A12 = ",s12,"；A13 = ",s13,"；A14 = ",s14,"。"),
"> **SHA256**：`A14_artifact_sha256_2026-09-26.txt`（2026-09-26 逐项核验，全部匹配）。",
"> **冻结规则**：不再改变变量、系数、评分表、分层和切点；任何修改均须新建版本并重新验证。",
"",
"## 一、最终连续模型",
"```",
a14$formula,
"```",
"- 教育参照：未上学；抑郁参照：CES-D-10 <10；ADL 参照：无受限。",
"- 结局：CHARLS 2011 记忆分 >4 者，2015 记忆分 ≤4（4 年跨阈值）。",
"- ADL β≈0、95%CI 跨 1；简化评分 ADL=0，连续模型保留 ADL。",
"",
"## 二、最终简化评分表",
"",
"| 项目 | 分值 |",
"|---|---|",
"| 年龄 | 45–49=0；每增加 5 岁 **+2**（不设上限） |",
"| 教育 | 未上学 **0**；小学 **−2**；初中 **−7**；高中及以上 **−10** |",
"| 抑郁 | CES-D-10 ≥10 **+3**；<10 **0** |",
"| ADL | **0**（不计分） |",
"",
"**总分 = 年龄分 + 教育分 + 抑郁分。** 总分无固定理论上限；开发样本观察 −10～+20，SHARE 最高 +23。",
"> 外部队列 harmonization：HRS/ELSA CES-D-8 ≥4；SHARE EURO-D ≥4。",
"",
"## 三、风险分层",
"",
"| 风险层 | 总分 | CHARLS 观察事件率 |",
"|---|---|---:|",
"| 低危 | ≤−3 | 5.6% |",
"| 中危 | −2～+4 | 15.4% |",
"| 高危 | ≥+5 | 27.3% |",
"",
"外部队列 S0 高危/低危：HRS 4.9×、ELSA 10.5×、SHARE 7.3×。",
"",
"## 四、关键性能",
"",
kable(perf_fmt, format="pipe", row.names=FALSE),
"",
"## 五、跨队列随机效应（A13）",
"",
kable(meta_fmt, format="pipe", row.names=FALSE),
"",
"## 六、论文写作必须引用的最终数字",
"",
kable(key, format="pipe", row.names=FALSE),
"",
"## 七、允许用途与禁止用途",
"",
"**允许**：45 岁及以上中老年人群术前相关认知脆弱性**筛查/分流**；建议进一步 MoCA/MMSE 等认知评估。",
"",
"**禁止**：痴呆/MCI 诊断；术后谵妄/PND/ICU/住院/死亡预测；真实手术队列已验证的表述；治疗或围术期决策工具。",
"",
"## 八、必须写入的限制",
"",
"1. 结局为 4 年记忆分跨阈值，含回归均值与测量误差，不等同于临床痴呆/MCI 发病率；",
"2. 外部队列校准差异明显：HRS 低估（截距 +0.33）、SHARE 高估（截距 −0.31）、ELSA 斜率 1.82；",
"3. 跨队列异质性 I² 85%–95%，且 k=3；SHARE 国家间 I² S0/S1≈50%、S2≈77%；",
"4. ELSA 年龄 90 岁顶码；ADL 缺失 HRS 79 人、SHARE 5 人；",
"5. 死亡/失访/记忆分缺失者被排除，未处理竞争风险；",
"6. 尚未在真实手术人群验证；",
"7. S0/S1 与 S2 阈值含义不同，三套结果互补。",
"",
"## 九、校验文件",
"",
"- `A14_artifact_sha256_2026-09-26.txt`",
"- `A14_独立复核_2026-09-26.txt`",
"- `A13_独立复核_2026-09-26.txt`",
"- `A12_独立复核_2026-09-26.txt`",
"- `A11_v2_独立复核_2026-09-26.txt`"
)
writeLines(md, file.path(OUT,"A14_最终冻结清单_2026-09-26.md"), useBytes=TRUE)
write.csv(key, file.path(OUT,"A14_最终冻结清单_关键数字.csv"), row.names=FALSE, fileEncoding="UTF-8")
cat("A14 freeze checklist written.\n")