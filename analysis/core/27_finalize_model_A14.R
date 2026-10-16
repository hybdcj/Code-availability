# ============================================================================
# 27_finalize_model_A14.R —— A14：最终模型与评分表冻结
# 输入：A07 模型 RDS、A08 Bootstrap、A11/A12/A13 外部验证结果
# 输出：output/12_模型定稿/
# ============================================================================
invisible(try(Sys.setlocale("LC_ALL","zh_CN.UTF-8"), silent=TRUE))
.libPaths(c("D:/R/library", .libPaths()))
suppressPackageStartupMessages({ library(dplyr); library(pROC); library(knitr) })

P   <- "D:/BaiduSyncdisk/科研/徐中豪"
OUT <- file.path(P, "output", "12_模型定稿")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

a07 <- readRDS(file.path(P, "output/10_筛查指数/screening_index_memory_v3_2026-09-26.rds"))
a08 <- readRDS(file.path(P, "output/10_筛查指数/A08_bootstrap_2026-09-26.rds"))
a12 <- readRDS(file.path(P, "output/11_外部验证/A12_results_2026-09-26.rds"))
a13 <- readRDS(file.path(P, "output/11_外部验证/A13_results_2026-09-26.rds"))

model_u <- a07$model_u
scenario_lab <- c(baseline_all_mean_sd="S0_baseline_all", complete_fu_mean_sd="S1_complete_fu", fixed_charls_le4="S2_fixed_le4")
model_levels <- model_u$xlevels$edu_group
coef_mat <- summary(model_u)$coefficients
ci_mat <- confint.default(model_u)
term_labels <- c(
  "(Intercept)" = "截距",
  "age" = "年龄（每岁）",
  "edu_groupPrimary school" = "教育：小学（vs 未上学）",
  "edu_groupMiddle school" = "教育：初中（vs 未上学）",
  "edu_groupHigh school or above" = "教育：高中及以上（vs 未上学）",
  "dep" = "抑郁：CES-D-10 ≥10（vs <10）",
  "adl" = "ADL：任一受限（vs 无受限）"
)
coef_table <- data.frame(
  term = rownames(coef_mat),
  term_zh = term_labels[rownames(coef_mat)],
  beta = coef_mat[, "Estimate"],
  se = coef_mat[, "Std. Error"],
  z = coef_mat[, "z value"],
  p = coef_mat[, "Pr(>|z|)"],
  beta_low = ci_mat[, 1],
  beta_high = ci_mat[, 2],
  OR = exp(coef_mat[, "Estimate"]),
  OR_low = exp(ci_mat[, 1]),
  OR_high = exp(ci_mat[, 2]),
  row.names = NULL, check.names = FALSE
)

# 最终连续模型公式
beta <- setNames(coef_table$beta, coef_table$term)
formula_text <- sprintf(
  "logit(P(4年记忆跨阈值)) = %.4f + %.4f×年龄 + (%.4f)×小学 + (%.4f)×初中 + (%.4f)×高中及以上 + %.4f×抑郁 + (%.4f)×ADL",
  beta[["(Intercept)"]], beta[["age"]], beta[["edu_groupPrimary school"]],
  beta[["edu_groupMiddle school"]], beta[["edu_groupHigh school or above"]],
  beta[["dep"]], beta[["adl"]]
)

# 最终评分表
score_table <- data.frame(
  项目 = c("年龄","教育","教育","教育","教育","抑郁","ADL"),
  类别 = c("45–49 岁=0；50–54=+2；每增加 5 岁再 +2，不设上限","未上学（Less than primary）",
           "小学（Primary school）","初中（Middle school）","高中及以上（High school or above）",
           "CES-D-10 ≥10","任一 ADL 受限"),
  分值 = c("+2 / 每 5 岁","0","−2","−7","−10","+3","0"),
  说明 = c("45–49 岁为 0 分，不设上限","参照组","相对未上学","相对未上学","相对未上学","vs CES-D-10 <10","β≈0，不计分；连续模型仍保留"),
  row.names = NULL
)

# 风险分层
strata_table <- a07$strata %>%
  mutate(风险层 = as.character(stratum),
         风险层定义 = c("≤−3","−2 ～ +4","≥+5")[match(stratum, c("低","中","高"))],
         事件率_num = as.numeric(sub("%","",事件率))/100) %>%
  select(风险层, 风险层定义, n, 事件, 事件率, 平均分)

# 内部性能
p_hat <- predict(model_u, type = "response")
r_cont <- roc(a07$data$incident_2015, p_hat, direction = "<", quiet = TRUE)
r_score <- roc(a07$data$incident_2015, a07$data$score, direction = "<", quiet = TRUE)
ci_cont <- as.numeric(ci.auc(r_cont)); ci_score <- as.numeric(ci.auc(r_score))
lp <- qlogis(pmin(pmax(p_hat, 1e-6), 1 - 1e-6))
cal_slope_fit <- glm(a07$data$incident_2015 ~ lp, family = binomial())
cal_slope_in <- coef(cal_slope_fit)[["lp"]]; cal_slope_ci <- confint.default(cal_slope_fit)["lp", ]
cal_int_fit <- glm(a07$data$incident_2015 ~ offset(lp), family = binomial())
cal_int_in <- coef(cal_int_fit)[["(Intercept)"]]; cal_int_ci <- confint.default(cal_int_fit)["(Intercept)", ]
brier_in <- mean((p_hat - a07$data$incident_2015)^2)
internal_apparent <- data.frame(
  cohort = "CHARLS", scenario = "internal_apparent",
  n_total = nrow(a07$data), events_total = sum(a07$data$incident_2015),
  n_model = nrow(a07$data), events_model = sum(a07$data$incident_2015),
  auc_model = as.numeric(auc(r_cont)), auc_model_low = ci_cont[1], auc_model_high = ci_cont[3],
  auc_score = as.numeric(auc(r_score)), auc_score_low = ci_score[1], auc_score_high = ci_score[3],
  cal_slope = cal_slope_in, cal_slope_low = cal_slope_ci[1], cal_slope_high = cal_slope_ci[2],
  cal_intercept = cal_int_in, cal_intercept_low = cal_int_ci[1], cal_intercept_high = cal_int_ci[2],
  brier = brier_in
)
internal_boot <- data.frame(
  cohort = "CHARLS", scenario = "internal_bootstrap",
  n_total = nrow(a07$data), events_total = sum(a07$data$incident_2015),
  n_model = nrow(a07$data), events_model = sum(a07$data$incident_2015),
  auc_model = a08$corrected[["auc"]], auc_model_low = a08$ci_auc[["2.5%"]], auc_model_high = a08$ci_auc[["97.5%"]],
  auc_score = NA_real_, auc_score_low = NA_real_, auc_score_high = NA_real_,
  cal_slope = a08$corrected[["slope"]], cal_slope_low = NA_real_, cal_slope_high = NA_real_,
  cal_intercept = a08$corrected[["intercept"]], cal_intercept_low = NA_real_, cal_intercept_high = NA_real_,
  brier = a08$corrected[["brier"]]
)
external_s0 <- a12$perf %>% filter(cut_definition == "baseline_all_mean_sd") %>%
  transmute(cohort, scenario = "external_S0",
            n_total, events_total, n_model, events_model,
            auc_model, auc_model_low, auc_model_high,
            auc_score = auc_score_full, auc_score_low = auc_score_full_low, auc_score_high = auc_score_full_high,
            cal_slope, cal_slope_low, cal_slope_high,
            cal_intercept, cal_intercept_low, cal_intercept_high, brier)
external_sens <- a12$perf %>% filter(cut_definition != "baseline_all_mean_sd") %>%
  transmute(cohort, scenario = paste0("external_", scenario_lab[cut_definition]),
            n_total, events_total, n_model, events_model,
            auc_model, auc_model_low, auc_model_high,
            auc_score = auc_score_full, auc_score_low = auc_score_full_low, auc_score_high = auc_score_full_high,
            cal_slope, cal_slope_low, cal_slope_high,
            cal_intercept, cal_intercept_low, cal_intercept_high, brier)
performance <- bind_rows(internal_apparent, internal_boot, external_s0, external_sens)

# 外部风险分层 S0
external_strata <- a12$strata %>% filter(cut_definition == "baseline_all_mean_sd") %>%
  transmute(cohort, stratum = as.character(stratum), n, events, rate, pred_mean, oe_model)

# 外部 meta
meta_summary <- a13$cross %>% mutate(scenario = scenario_lab[cut_definition], metric_zh = ifelse(metric=="model","连续模型","整数分"))

# 保存 CSV
write.csv(coef_table, file.path(OUT,"A14_最终模型系数.csv"), row.names=FALSE, fileEncoding="UTF-8")
write.csv(score_table, file.path(OUT,"A14_最终评分表.csv"), row.names=FALSE, fileEncoding="UTF-8")
write.csv(strata_table, file.path(OUT,"A14_风险分层.csv"), row.names=FALSE, fileEncoding="UTF-8")
write.csv(performance, file.path(OUT,"A14_性能汇总.csv"), row.names=FALSE, fileEncoding="UTF-8")
write.csv(external_strata, file.path(OUT,"A14_外部风险分层_S0.csv"), row.names=FALSE, fileEncoding="UTF-8")
write.csv(meta_summary, file.path(OUT,"A14_跨队列meta汇总.csv"), row.names=FALSE, fileEncoding="UTF-8")

# 保存最终模型 artifact
final_model <- list(
  model = model_u,
  coef_table = coef_table,
  formula = formula_text,
  score_table = score_table,
  strata = strata_table,
  performance = performance,
  external_strata = external_strata,
  meta = meta_summary,
  metadata = list(
    version = "v1.0",
    frozen_date = "2026-09-26",
    development = "CHARLS 2011→2015; n=4,441; events=658",
    primary_outcome = "4-year memory cross-threshold (2015 memory score ≤4; baseline >4)",
    external_validation = "HRS w10→w12; ELSA w4→w6; SHARE w4→w6",
    intended_use = "preoperative-relevant cognitive vulnerability screening/triage; not a diagnostic or postoperative outcome prediction tool",
    no_refit = TRUE
  )
)
saveRDS(final_model, file.path(OUT,"A14_final_model_2026-09-26.rds"))
writeLines(capture.output(sessionInfo()), file.path(OUT,"A14_sessionInfo_2026-09-26.txt"), useBytes=TRUE)

# 报告表格式化
fmt <- function(x,d=3) ifelse(is.na(x),"NA",sprintf(paste0("%.",d,"f"),x))
coef_tab <- coef_table %>% mutate(
  项 = term_zh,
  β = fmt(beta,4), SE = fmt(se,4), OR = fmt(OR,3),
  OR_95CI = paste0(fmt(OR_low,3),"–",fmt(OR_high,3)),
  p = ifelse(p<0.001,"<0.001",fmt(p,3))
) %>% mutate(OR = ifelse(term=="(Intercept)","—",OR),
             OR_95CI = ifelse(term=="(Intercept)","—",OR_95CI)) %>%
  select(项, β, SE, OR, OR_95CI, p)
score_tab <- score_table %>% select(项目, 类别, 分值, 说明)
strata_tab <- strata_table
fmt_ci <- function(a,lo,hi,d=3) ifelse(is.na(a), "—", paste0(fmt(a,d), " (", fmt(lo,d), "–", fmt(hi,d), ")"))
perf_tab <- performance %>% filter(scenario %in% c("internal_apparent","internal_bootstrap","external_S0")) %>%
  mutate(
    AUC模型 = fmt_ci(auc_model, auc_model_low, auc_model_high),
    AUC整数分 = fmt_ci(auc_score, auc_score_low, auc_score_high),
    校准斜率 = fmt(cal_slope,2), 校准截距 = fmt(cal_intercept,2), Brier = fmt(brier,3)
  ) %>% select(队列=cohort, 场景=scenario, 全风险n=n_total, 全风险事件=events_total, 模型n=n_model, 事件=events_model, AUC模型, AUC整数分, 校准斜率, 校准截距, Brier)
sens_tab <- performance %>% filter(scenario %in% c("external_S1_complete_fu","external_S2_fixed_le4")) %>%
  mutate(
    AUC模型 = fmt_ci(auc_model, auc_model_low, auc_model_high),
    AUC整数分 = fmt_ci(auc_score, auc_score_low, auc_score_high),
    校准斜率 = fmt(cal_slope,2), 校准截距 = fmt(cal_intercept,2), Brier = fmt(brier,3)
  ) %>% select(队列=cohort, 场景=scenario, 全风险n=n_total, 全风险事件=events_total, 模型n=n_model, 事件=events_model, AUC模型, AUC整数分, 校准斜率, 校准截距, Brier)
meta_tab <- meta_summary %>% mutate(
  合并AUC = paste0(fmt(pooled_auc), " (", fmt(pooled_low), "–", fmt(pooled_high), ")"),
  I2 = sprintf("%.1f%%", I2)
) %>% select(切点=scenario, 指标=metric_zh, 合并AUC, I2, k)

md <- c(
"# A14：最终模型与评分表冻结（v1.0，2026-09-26）",
"",
"> **冻结声明**：自 2026-09-26 起，最终模型、评分表、风险分层和公式冻结。",
"> 外部队列仅用于验证，**未重新拟合、未重新选变量、未改变切点**；后续任何修改须标注为事后分析并重新验证。",
"",
"## 一、最终连续模型",
"",
paste0("**模型公式**：", formula_text),
"",
kable(coef_tab, format="pipe", row.names=FALSE),
"",
"> 参考类别：教育=未上学（Less than primary）；抑郁=CES-D-10 <10；ADL=无受限。",
"> 结局：CHARLS 2011 年记忆分 >4 者，2015 年记忆分 ≤4（4 年跨阈值）。",
"> ADL 系数接近 0、95%CI 跨 1，因此最终简化评分表中 ADL 计 0 分；连续模型保留 ADL 以保持预先设定的完整模型。",
"",
"## 二、最终简化评分表",
"",
kable(score_tab, format="pipe", row.names=FALSE),
"",
"**总分** = 年龄分 + 教育分 + 抑郁分；ADL 不计分。",
"> 年龄分不设上限（每增加 5 岁 +2）。总分不是固定 −10～+20；开发样本观察范围 −10～+20，外部 SHARE 观察最高 +23。
> 外部队列应用时，抑郁切点按预先 harmonization：HRS/ELSA CES-D-8 ≥4；SHARE EURO-D ≥4（详见 A11 v2.2 映射表）。",
"",
"## 三、风险分层",
"",
kable(strata_tab, format="pipe", row.names=FALSE),
"",
"> 分层来自 CHARLS 开发样本整数总分三分位：低 ≤−3；中 −2～+4；高 ≥+5。",
"> 开发样本观察事件率：低 5.6%、中 15.4%、高 27.3%；外部队列 S0 事件率见 A14_外部风险分层_S0.csv。",
"",
"## 四、内部与外部性能汇总",
"",
kable(perf_tab, format="pipe", row.names=FALSE),
"",
"> 内部 Bootstrap 校正后 AUC 0.698（0.679–0.720），校准斜率 0.98，Brier 0.1184。 表中“模型n/事件”为连续模型完整样本；“全风险n/事件”为整数分 AUC 使用的样本。HRS 全风险集 14,051/1,842；ELSA 6,725/554；SHARE 26,879/2,930。",
"> 外部队列 S0 连续模型 AUC：HRS 0.696、ELSA 0.763、SHARE 0.751；校准存在明显队列差异（HRS 低估、SHARE 高估、ELSA 斜率偏大）。",
"",
"## 五、正式敏感性 S1/S2",
"",
kable(sens_tab, format="pipe", row.names=FALSE),
"",
"> S1=complete-fu mean−1SD；S2=固定 CHARLS ≤4。S2 AUC 更高但事件率更低，且依赖跨队列记忆分绝对可比性；主分析仍为 S0。",
"",
"## 六、跨队列随机效应汇总（A13）",
"",
kable(meta_tab, format="pipe", row.names=FALSE),
"",
"> 跨队列 I² 很高（85%–95%），合并值仅作描述；论文应重点报告各队列结果、校准差异和异质性来源。",
"",
"## 七、最终用途与边界",
"",
"**允许的用途**：",
"- 45 岁及以上中老年人群的术前相关认知脆弱性**筛查/分流**；",
"- 识别可能需要进一步 MoCA/MMSE 等认知评估者；",
"- 作为公开数据开发、内部验证和跨队列外部验证的筛查指数。",
"",
"**禁止的用途**：",
"- 不作为痴呆/MCI 诊断工具；",
"- 不预测术后谵妄、PND、ICU、住院或死亡；",
"- 不声称已在真实手术队列验证；",
"- 不作为治疗决策或围术期决策工具。",
"",
"## 八、限制",
"",
"1. 结局是 4 年记忆分跨阈值，含基线位置、测量误差和回归均值成分，不等同于临床痴呆/MCI 发病率；",
"2. 外部队列校准存在队列差异，HRS 低估风险、SHARE 高估风险、ELSA 校准斜率 1.82；",
"3. 跨队列异质性高（I² 85%–95%），且 k=3；",
"4. ELSA 年龄 90 岁顶码；ADL 缺失 HRS 79 人、SHARE 5 人；",
"5. 死亡/失访/记忆分缺失者被排除，未处理竞争风险；",
"6. 尚未在真实手术人群验证；",
"7. S0/S1 使用各队列内 mean−1SD，与 CHARLS 固定 ≤4 不是同一绝对结局；S2 固定 ≤4 假设跨队列记忆分绝对可比。三套结果互补，不能互相替代。",
"",
"## 九、可复现性与文件",
"",
"- 模型 artifact：`A14_final_model_2026-09-26.rds`",
"- 系数：`A14_最终模型系数.csv`",
"- 评分表：`A14_最终评分表.csv`",
"- 风险分层：`A14_风险分层.csv`",
"- 性能汇总：`A14_性能汇总.csv`、`A14_外部风险分层_S0.csv`、`A14_跨队列meta汇总.csv`",
"- Session info：`A14_sessionInfo_2026-09-26.txt`",
"- 冻结校验：`A14_artifact_sha256_2026-09-26.txt`",
"- 脚本：`code/27_finalize_model_A14.R`",
"- 独立复核：`A14_独立复核_2026-09-26.txt`（由 verify_a14 生成）",
"",
"## 十、结论",
"",
"> 最终模型与简化评分表冻结，版本 v1.0（2026-09-26）。",
"> 最终工具是**术前相关认知脆弱性筛查/分流工具**，不是术后结局预测工具。"
)
writeLines(md, file.path(OUT,"A14_最终模型规范_2026-09-26.md"), useBytes=TRUE)

# 一页评分卡
scorecard <- c(
"# 术前相关认知脆弱性筛查评分卡（v1.0，2026-09-26）",
"",
"> 用途：45 岁及以上中老年人群的术前认知评估分流；不作诊断，不预测术后谵妄/PND。",
"",
"## 评分步骤",
"",
"### 1. 年龄分",
"| 年龄 | 分值 |",
"|---|---|",
"| 45–49 岁 | 0 |",
"| 50–54 岁 | +2 |",
"| 55–59 岁 | +4 |",
"| 60–64 岁 | +6 |",
"| 65–69 岁 | +8 |",
"| 70–74 岁 | +10 |",
"| 75–79 岁 | +12 |",
"| 80–84 岁 | +14 |",
"| 85–89 岁 | +16 |",
"| 90–94 岁 | +18 |",
"| 每再增加 5 岁 | 再 +2 |",
"",
"### 2. 教育分",
"| 教育 | 分值 |",
"|---|---|",
"| 未上学 | 0 |",
"| 小学 | −2 |",
"| 初中 | −7 |",
"| 高中及以上 | −10 |",
"",
"### 3. 抑郁分",
"| CES-D-10 | 分值 |",
"|---|---|",
"| <10 | 0 |",
"| ≥10 | +3 |",
"",
"### 4. ADL",
"ADL 不计分（β≈0）；连续模型仍保留 ADL。",
"",
"### 5. 总分与风险层",
"| 风险层 | 总分 | CHARLS 观察事件率 |",
"|---|---|---|",
"| 低危 | ≤−3 | 5.6% |",
"| 中危 | −2 ～ +4 | 15.4% |",
"| 高危 | ≥+5 | 27.3% |",
"",
"**示例**：70 岁、小学、CES-D-10 ≥10 → 10 + (−2) + 3 = **11 分 → 高危**。",
"",
"## 重要提醒",
"- 高风险只表示建议进一步认知评估，不等于诊断；",
"- 外部队列校准存在差异，个体概率解释需谨慎；",
"- 不用于预测术后谵妄、PND 或围术期结局；",
"- 尚未在真实手术队列验证。"
)
writeLines(scorecard, file.path(OUT,"A14_评分卡_2026-09-26.md"), useBytes=TRUE)

cat("A14 finalized.\n")