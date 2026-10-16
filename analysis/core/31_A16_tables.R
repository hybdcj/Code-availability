# ============================================================================
# 31_A16_tables.R —— A16：Table 1–5 中英文投稿版
# 输入：A07/A14/A15、A08/A09、A12/A13 冻结结果
# 输出：output/14_A16_图表/{中文,English,data}/
# 生成：2026-09-27
# ============================================================================
invisible(try(Sys.setlocale("LC_ALL", "zh_CN.UTF-8"), silent = TRUE))
.libPaths(c("D:/R/library", .libPaths()))
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(openxlsx)
  library(officer); library(flextable)
})
P <- "D:/BaiduSyncdisk/科研/徐中豪"
OUT <- file.path(P, "output", "14_A16_图表")
ZH <- file.path(OUT, "中文"); EN <- file.path(OUT, "English"); DATA <- file.path(OUT, "data")
for (d in c(OUT, ZH, EN, DATA)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# 通用格式函数
# ---------------------------------------------------------------------------
fmt_mean_sd <- function(x, d = 1) sprintf(paste0("%.", d, "f (%.", d, "f)"), mean(x, na.rm = TRUE), sd(x, na.rm = TRUE))
fmt_pct <- function(x, d = 1) sprintf(paste0("%.", d, "f%%"), 100 * x)
fmt_n <- function(x) trimws(format(round(x), big.mark = ","))
fmt_n_pct <- function(k, n, d = 1) sprintf("%s (%.1f%%)", fmt_n(k), 100 * k / n)
fmt_ci <- function(est, lo, hi, d = 3) ifelse(is.na(est), "—", sprintf(paste0("%.", d, "f (%.", d, "f–%.", d, "f)"), est, lo, hi))
fmt_num <- function(x, d = 2) ifelse(is.na(x), "—", sprintf(paste0("%.", d, "f"), x))
parse_leading_num <- function(x) suppressWarnings(as.numeric(sub(" .*$", "", trimws(x))))
en_dash <- function(x) gsub("-", "–", x, fixed = TRUE)
clean_est <- function(x) en_dash(sub(", p=.*$", "", x))

safe_write_csv <- function(x, file, ...) {
  for (i in seq_len(10)) {
    ok <- tryCatch({ utils::write.csv(x, file, ...); TRUE }, error = function(e) FALSE)
    if (ok) return(invisible(TRUE))
    Sys.sleep(0.5)
  }
  stop("无法写入文件（重试10次）: ", file)
}
safe_write_xlsx <- function(x, file, ...) {
  for (i in seq_len(10)) {
    ok <- tryCatch({ openxlsx::write.xlsx(x, file = file, ...); TRUE }, error = function(e) FALSE)
    if (ok) return(invisible(TRUE))
    Sys.sleep(0.5)
  }
  stop("无法写入文件（重试10次）: ", file)
}
safe_save_docx <- function(x, path) {
  for (i in seq_len(10)) {
    ok <- tryCatch({ print(x, target = path); TRUE }, error = function(e) { message("DOCX错误: ", conditionMessage(e)); FALSE })
    if (ok) return(invisible(TRUE))
    Sys.sleep(0.5)
  }
  stop("无法写入文件（重试10次）: ", path)
}

# ---------------------------------------------------------------------------
# Table 1：队列基线
# ---------------------------------------------------------------------------
a07 <- readRDS(file.path(P, "output/10_筛查指数/screening_index_memory_v3_2026-09-26.rds"))
ch <- a07$data %>% mutate(event = incident_2015, mem_base = mem2011)
cohort_data <- list(
  CHARLS = list(d = ch, role_zh = "开发队列", role_en = "Development",
                cut_zh = "≤4（固定切点，raw 4.765）", cut_en = "≤4 (fixed, raw 4.765)"),
  HRS = list(d = NULL, role_zh = "外部验证", role_en = "External validation",
             cut_zh = NULL, cut_en = NULL),
  ELSA = list(d = NULL, role_zh = "外部验证", role_en = "External validation",
              cut_zh = NULL, cut_en = NULL),
  SHARE = list(d = NULL, role_zh = "外部验证", role_en = "External validation",
               cut_zh = NULL, cut_en = NULL)
)
for (co in c("HRS", "ELSA", "SHARE")) {
  x <- readRDS(file.path(P, "external_data", co, "analysis",
                         paste0(co, "_sensitivity_risk_sets_2026-09-26.rds")))[[paste0(co, "__baseline_all_mean_sd")]]
  cut_raw <- unique(x$cut)[1]; cut_int <- unique(x$cut_int)[1]
  cohort_data[[co]]$d <- x
  cohort_data[[co]]$cut_zh <- sprintf("≤%d（队列 mean−1SD，raw %.2f）", cut_int, cut_raw)
  cohort_data[[co]]$cut_en <- sprintf("≤%d (cohort mean−1SD, raw %.2f)", cut_int, cut_raw)
}

cohort_rows <- c("role", "n", "age", "age65", "edu_lt", "edu_primary", "edu_middle", "edu_high",
                 "dep", "adl", "mem", "events", "cut", "fu")
cohort_values <- list()
for (co in names(cohort_data)) {
  z <- cohort_data[[co]]; d <- z$d
  edu <- prop.table(table(factor(d$edu_group, levels = c("Less than primary", "Primary school", "Middle school", "High school or above"))))
  cohort_values[[co]] <- c(
    role = z$role_zh,
    n = fmt_n(nrow(d)),
    age = fmt_mean_sd(d$age),
    age65 = fmt_pct(mean(d$age >= 65, na.rm = TRUE)),
    edu_lt = fmt_pct(edu[["Less than primary"]]),
    edu_primary = fmt_pct(edu[["Primary school"]]),
    edu_middle = fmt_pct(edu[["Middle school"]]),
    edu_high = fmt_pct(edu[["High school or above"]]),
    dep = fmt_pct(mean(d$dep, na.rm = TRUE)),
    adl = sprintf("%s (missing %s)", fmt_pct(mean(d$adl, na.rm = TRUE)), fmt_n(sum(is.na(d$adl)))),
    mem = fmt_mean_sd(d$mem_base),
    events = fmt_n_pct(sum(d$event, na.rm = TRUE), nrow(d)),
    cut = z$cut_zh,
    fu = "4 年"
  )
}
cohort_values_en <- list()
for (co in names(cohort_data)) {
  z <- cohort_data[[co]]; d <- z$d
  edu <- prop.table(table(factor(d$edu_group, levels = c("Less than primary", "Primary school", "Middle school", "High school or above"))))
  cohort_values_en[[co]] <- c(
    role = z$role_en,
    n = fmt_n(nrow(d)),
    age = fmt_mean_sd(d$age),
    age65 = fmt_pct(mean(d$age >= 65, na.rm = TRUE)),
    edu_lt = fmt_pct(edu[["Less than primary"]]),
    edu_primary = fmt_pct(edu[["Primary school"]]),
    edu_middle = fmt_pct(edu[["Middle school"]]),
    edu_high = fmt_pct(edu[["High school or above"]]),
    dep = fmt_pct(mean(d$dep, na.rm = TRUE)),
    adl = sprintf("%s (missing %s)", fmt_pct(mean(d$adl, na.rm = TRUE)), fmt_n(sum(is.na(d$adl)))),
    mem = fmt_mean_sd(d$mem_base),
    events = fmt_n_pct(sum(d$event, na.rm = TRUE), nrow(d)),
    cut = z$cut_en,
    fu = "4 years"
  )
}
rows_zh <- c(role = "队列角色", n = "样本量", age = "年龄（岁）", age65 = "年龄≥65 岁",
             edu_lt = "教育：未上学", edu_primary = "教育：小学", edu_middle = "教育：初中",
             edu_high = "教育：高中及以上", dep = "抑郁（队列特异性切点）",
             adl = "ADL 受限", mem = "基线记忆分（0–20）", events = "随访事件，n（%）",
             cut = "认知受损切点", fu = "随访时长")
rows_en <- c(role = "Cohort role", n = "Sample size", age = "Age, years", age65 = "Age ≥65 years",
             edu_lt = "Education: less than primary", edu_primary = "Education: primary school",
             edu_middle = "Education: middle school", edu_high = "Education: high school or above",
             dep = "Depression (cohort-specific threshold)", adl = "ADL limitation",
             mem = "Baseline memory score (0–20)", events = "Follow-up events, n (%)",
             cut = "Impairment cutpoint", fu = "Follow-up duration")

tab1_zh <- data.frame(特征 = rows_zh, CHARLS = cohort_values$CHARLS[names(rows_zh)],
                      HRS = cohort_values$HRS[names(rows_zh)], ELSA = cohort_values$ELSA[names(rows_zh)],
                      SHARE = cohort_values$SHARE[names(rows_zh)], row.names = NULL, check.names = FALSE)
tab1_en <- data.frame(Characteristic = rows_en, CHARLS = cohort_values_en$CHARLS[names(rows_en)],
                      HRS = cohort_values_en$HRS[names(rows_en)], ELSA = cohort_values_en$ELSA[names(rows_en)],
                      SHARE = cohort_values_en$SHARE[names(rows_en)], row.names = NULL, check.names = FALSE)
tab1_note_zh <- "缩写：ADL，日常生活活动能力；CES-D，流行病学研究中心抑郁量表；EURO-D，欧洲抑郁量表。外部队列抑郁切点为 HRS/ELSA CES-D-8 ≥4、SHARE EURO-D ≥4；CHARLS 为 CES-D-10 ≥10。外部队列 ADL 缺失：HRS 79 人、SHARE 5 人。外部队列事件按 S0（各队列基线全样本 mean−1SD）切点定义；CHARLS 为固定切点 ≤4。性别未纳入本表，因外部验证数据未保留性别变量，且性别不是最终指数变量。"
tab1_note_en <- "Abbreviations: ADL, activities of daily living; CES-D, Center for Epidemiologic Studies Depression Scale; EURO-D, European Depression Scale. Depression thresholds: HRS/ELSA CES-D-8 ≥4, SHARE EURO-D ≥4, CHARLS CES-D-10 ≥10. Missing ADL: HRS 79, SHARE 5. External events use the S0 cutpoint (cohort-specific baseline mean−1SD); CHARLS uses fixed ≤4. Sex was not included because it was not retained in the harmonized external validation datasets and is not part of the final index."

# ---------------------------------------------------------------------------
# Table 2：候选变量筛选与效应量
# ---------------------------------------------------------------------------
orm <- read.csv(file.path(P, "output/03_五波横断面/五波分波Logistic_OR矩阵_2026-09-19.csv"),
                check.names = FALSE, fileEncoding = "UTF-8")
risk <- read.csv(file.path(P, "output/04_纵向队列/cohort_risk_2026-09-13.csv"),
                 check.names = FALSE, fileEncoding = "UTF-8")
coef_final <- read.csv(file.path(P, "output/12_模型定稿/A14_最终模型系数.csv"),
                       check.names = FALSE, fileEncoding = "UTF-8")
iadl2 <- read.csv(file.path(P, "output/04_纵向队列/iadl2_longitudinal_models_2026-09-26.csv"),
                  check.names = FALSE, fileEncoding = "UTF-8")
or_range <- function(pattern) {
  idx <- grep(pattern, orm[[1]])
  if (length(idx) == 0) return("—")
  vals <- unlist(orm[idx[1], 3:7], use.names = FALSE)
  nums <- parse_leading_num(vals); nums <- nums[!is.na(nums)]
  if (length(nums) == 0) return("—")
  sprintf("%.2f–%.2f", min(nums), max(nums))
}
risk_est <- function(rowname) {
  idx <- which(risk[[1]] == rowname)
  if (length(idx) == 0) return("—")
  paste0("OR ", clean_est(risk[idx[1], 2]), "; HR ", clean_est(risk[idx[1], 3]))
}
final_or <- function(term) {
  z <- coef_final[coef_final$term == term, ]
  if (nrow(z) == 0) return("—")
  fmt_ci(z$OR, z$OR_low, z$OR_high)
}
final_or5 <- function(term) {
  z <- coef_final[coef_final$term == term, ]
  if (nrow(z) == 0) return("—")
  fmt_ci(exp(5*z$beta), exp(5*z$beta_low), exp(5*z$beta_high))
}
iadl2_str_zh <- function() {
  getv <- function(model) iadl2$unweighted[iadl2$exposure == "iadl2" & iadl2$model == model]
  sprintf("OR %s；HR %s；变化 β %s → %s（基线校正后）",
          en_dash(getv("logistic")), en_dash(getv("cox")), en_dash(getv("change")), en_dash(getv("change_adj")))
}
iadl2_str_en <- function() {
  getv <- function(model) iadl2$unweighted[iadl2$exposure == "iadl2" & iadl2$model == model]
  sprintf("OR %s; HR %s; change β %s → %s after baseline adjustment",
          en_dash(getv("logistic")), en_dash(getv("cox")), en_dash(getv("change")), en_dash(getv("change_adj")))
}
tab2_zh <- data.frame(
  变量 = c("年龄（每 5 岁）", "教育（高中及以上 vs 未上学）", "抑郁（CES-D-10 ≥10）",
           "ADL 受限", "IADL 受限（补充）", "睡眠 <6 h", "睡眠 >8 h", "规律体力活动",
           "卒中", "糖尿病", "高血压", "疼痛", "PLT/WBC/CRP、Neu/NLR"),
  角色 = c("核心", "核心", "核心", "核心", "补充", "增量", "增量", "增量",
           "增量", "增量", "背景/探索", "排除", "排除"),
  `最终模型 OR` = c(final_or5("age"), final_or("edu_groupHigh school or above"), final_or("dep"), final_or("adl"),
                  "—", "—", "—", "—", "—", "—", "—", "—", "—"),
  `横断面 OR 范围` = c("—", "—", or_range("抑郁"), or_range("ADL 受限"), or_range("IADL 受限"),
                 or_range("睡眠 <6"), or_range("睡眠 >8"), or_range("规律体力活动"),
                 or_range("卒中"), or_range("糖尿病"), or_range("高血压"), or_range("躯体疼痛"), "—"),
  纵向估计 = c("—", "—", risk_est("抑郁（CES-D-10 ≥10）"), risk_est("ADL 受限"),
             iadl2_str_zh(), risk_est("睡眠时长（每 1 小时）"), "—", risk_est("规律体力活动"),
             "—", risk_est("糖尿病"), risk_est("高血压"), risk_est("疼痛"), "无稳健关联（Neu/NLR 不可得）"),
  敏感性 = c("最终模型", "最终模型", "加权/IPW 方向一致", "加权/IPW 方向一致",
           "变化模型显示回归均值，不进主评分", "原分析较弱", "原分析较弱",
           "加权后显著性不稳定", "原分析为探索性", "原分析较弱", "不进核心评分",
           "加权/IPW 方向不稳定", "不进入核心工具"),
  最终处理 = c("纳入连续模型与评分表（每 5 岁 +2）", "纳入连续模型与评分表（0/−2/−7/−10）",
           "纳入连续模型与评分表（+3）", "连续模型保留；最终记忆结局 β≈0，评分计 0",
           "仅补充分析", "增量/探索性", "增量/探索性", "增量/探索性", "增量/探索性",
           "增量/探索性", "不作为核心评分", "排除核心评分", "排除核心评分"),
  row.names = NULL, check.names = FALSE
)
tab2_en <- data.frame(
  Variable = c("Age (per 5 years)", "Education (high school or above vs less than primary)",
               "Depression (CES-D-10 ≥10)", "ADL limitation", "IADL limitation (supplementary)",
               "Sleep <6 h", "Sleep >8 h", "Regular physical activity", "Stroke", "Diabetes",
               "Hypertension", "Pain", "PLT/WBC/CRP and Neu/NLR"),
  Role = c("Core", "Core", "Core", "Core", "Supplementary", "Incremental", "Incremental",
           "Incremental", "Incremental", "Incremental", "Background/exploratory", "Excluded", "Excluded"),
  `Final model OR` = c(final_or5("age"), final_or("edu_groupHigh school or above"), final_or("dep"), final_or("adl"),
                      "—", "—", "—", "—", "—", "—", "—", "—", "—"),
  `Cross-sectional OR range` = c("—", "—", or_range("抑郁"), or_range("ADL 受限"), or_range("IADL 受限"),
                 or_range("睡眠 <6"), or_range("睡眠 >8"), or_range("规律体力活动"),
                 or_range("卒中"), or_range("糖尿病"), or_range("高血压"), or_range("躯体疼痛"), "—"),
  `Longitudinal estimate` = c("—", "—", risk_est("抑郁（CES-D-10 ≥10）"), risk_est("ADL 受限"),
             iadl2_str_zh(), risk_est("睡眠时长（每 1 小时）"), "—", risk_est("规律体力活动"),
             "—", risk_est("糖尿病"), risk_est("高血压"), risk_est("疼痛"), "No robust association (Neu/NLR unavailable)"),
  `Sensitivity` = c("Final model", "Final model", "Consistent under weighting/IPW",
                    "Consistent under weighting/IPW", "Change model showed regression to the mean; not in score",
                    "Weak in primary analysis", "Weak in primary analysis", "Significance unstable under weighting",
                    "Exploratory", "Weak in primary analysis", "Not in core score",
                    "Direction unstable under weighting/IPW", "Not in core tool"),
  `Final decision` = c("Included in continuous model and score (+2 per 5 years)",
                       "Included in continuous model and score (0/−2/−7/−10)",
                       "Included in continuous model and score (+3)",
                       "Retained in continuous model; β≈0 for memory endpoint, score 0",
                       "Supplementary only", "Incremental/exploratory", "Incremental/exploratory",
                       "Incremental/exploratory", "Incremental/exploratory", "Incremental/exploratory",
                       "Not in core score", "Excluded from core score", "Excluded from core tool"),
  row.names = NULL, check.names = FALSE
)
tab2_note_zh <- "横断面 OR 范围为 2011–2020 五波 ms10 认知受损 Logistic 模型的 OR 最小值–最大值；纵向估计来自 T06：2020 终点评定 OR 与首次受损 Cox HR（未加权主分析），IADL 另列 iadl2 连续变化校正前后 β。“最终模型 OR”来自 2011→2015 记忆分跨阈值结局的 A14 冻结模型；年龄为每 5 岁 OR，教育为高中及以上 vs 未上学。其余候选变量的横断面/纵向估计来自 T06/T07，仅用于筛选与稳定性展示。"
tab2_note_en <- "Cross-sectional OR ranges are the minimum–maximum ORs across the 2011–2020 five-wave ms10 impairment logistic models. Longitudinal estimates are from T06: OR for the 2020 endpoint and HR for first impairment (unweighted primary analysis); IADL also shows unadjusted and baseline-adjusted change β from the iadl2 analysis. “Final model OR” is from the A14 frozen model for the 2011→2015 memory cross-threshold endpoint; age is per 5 years and education is high school or above vs less than primary. Other candidate estimates are from T06/T07 and are shown for screening and stability only."

# ---------------------------------------------------------------------------
# Table 3：最终评分规则与风险分层
# ---------------------------------------------------------------------------
score_tab <- read.csv(file.path(P, "output/12_模型定稿/A14_最终评分表.csv"),
                      check.names = FALSE, fileEncoding = "UTF-8")
strata_tab <- read.csv(file.path(P, "output/12_模型定稿/A14_风险分层.csv"),
                       check.names = FALSE, fileEncoding = "UTF-8")
tab3a_zh <- score_tab %>% transmute(项目, 类别, 分值, 说明)
tab3a_en <- data.frame(
  Item = c("Age", "Education", "Education", "Education", "Education", "Depression", "ADL"),
  Category = c("45–49 years = 0; +2 per additional 5 years, no upper limit",
               "Less than primary", "Primary school", "Middle school", "High school or above",
               "CES-D-10 ≥10", "Any ADL limitation"),
  Points = c("+2 per 5 years", "0", "−2", "−7", "−10", "+3", "0"),
  Note = c("Reference: 45–49 years; no upper limit", "Reference group", "vs less than primary",
           "vs less than primary", "vs less than primary", "vs CES-D-10 <10",
           "β≈0; score 0; retained in continuous model"),
  row.names = NULL, check.names = FALSE
)
tab3b_zh <- strata_tab %>% transmute(风险层, 风险层定义, n = fmt_n(n), 事件 = fmt_n(事件),
                                     事件率, 平均分)
tab3b_en <- data.frame(
  `Risk stratum` = c("Low", "Intermediate", "High"),
  `Score range` = c("≤−3", "−2 to +4", "≥+5"),
  n = fmt_n(strata_tab$n), Events = fmt_n(strata_tab$事件),
  `Event rate` = strata_tab$事件率, `Mean score` = strata_tab$平均分,
  row.names = NULL, check.names = FALSE
)
tab3_note_zh <- "总分 = 年龄分 + 教育分 + 抑郁分；ADL 计 0 分。风险层来自 CHARLS 开发样本整数总分三分位：低 ≤−3、中 −2～+4、高 ≥+5。外部队列应用时，HRS/ELSA 抑郁为 CES-D-8 ≥4、SHARE 为 EURO-D ≥4。"
tab3_note_en <- "Total score = age score + education score + depression score; ADL contributes 0 points. Risk strata are based on tertiles of the integer score in the CHARLS development sample: low ≤−3, intermediate −2 to +4, high ≥+5. For external application, depression is defined as CES-D-8 ≥4 in HRS/ELSA and EURO-D ≥4 in SHARE."

# ---------------------------------------------------------------------------
# Table 4：内部验证
# ---------------------------------------------------------------------------
perf <- read.csv(file.path(P, "output/12_模型定稿/A14_性能汇总.csv"),
                 check.names = FALSE, fileEncoding = "UTF-8")
a08 <- read.csv(file.path(P, "output/10_筛查指数/A08_bootstrap_optimism_2026-09-26.csv"),
                check.names = FALSE, fileEncoding = "UTF-8")
strata_pred <- read.csv(file.path(P, "output/10_筛查指数/A09_strata_observed_vs_predicted_2026-09-26.csv"),
                        check.names = FALSE, fileEncoding = "UTF-8")
dca <- read.csv(file.path(P, "output/10_筛查指数/A09_dca_2026-09-26.csv"),
                check.names = FALSE, fileEncoding = "UTF-8")
p_app <- perf[perf$scenario == "internal_apparent", ]
p_boot <- perf[perf$scenario == "internal_bootstrap", ]
tab4a_zh <- data.frame(
  指标 = c("AUC（连续模型，表观）", "AUC（整数分，表观）", "AUC（Bootstrap 校正）",
           "校准斜率", "校准截距", "Brier 分数"),
  估计 = c(fmt_ci(p_app$auc_model, p_app$auc_model_low, p_app$auc_model_high),
           fmt_ci(p_app$auc_score, p_app$auc_score_low, p_app$auc_score_high),
           fmt_ci(p_boot$auc_model, p_boot$auc_model_low, p_boot$auc_model_high),
           sprintf("%s（表观 %s）", fmt_num(p_boot$cal_slope), fmt_num(p_app$cal_slope)),
           sprintf("%s（表观 %s）", fmt_num(p_boot$cal_intercept), fmt_num(p_app$cal_intercept)),
           sprintf("%s（表观 %s）", fmt_num(p_boot$brier, 4), fmt_num(p_app$brier, 4))),
   row.names = NULL, check.names = FALSE
)
tab4a_en <- data.frame(
  Metric = c("AUC (continuous model, apparent)", "AUC (integer score, apparent)",
             "AUC (bootstrap-corrected)", "Calibration slope", "Calibration intercept", "Brier score"),
  Estimate = tab4a_zh$估计, row.names = NULL, check.names = FALSE
)
tab4b_zh <- data.frame(
  风险层 = c("低", "中", "高"), 分值范围 = c("≤−3", "−2～+4", "≥+5"),
  n = fmt_n(strata_pred$n), 事件 = fmt_n(strata_pred$事件),
  观测事件率 = fmt_pct(strata_pred$观测率), 模型预测率 = fmt_pct(strata_pred$预测率),
  平均分 = strata_pred$平均分, row.names = NULL, check.names = FALSE
)
tab4b_en <- data.frame(
  `Risk stratum` = c("Low", "Intermediate", "High"), `Score range` = c("≤−3", "−2 to +4", "≥+5"),
  n = fmt_n(strata_pred$n), Events = fmt_n(strata_pred$事件),
  `Observed event rate` = fmt_pct(strata_pred$观测率), `Predicted risk` = fmt_pct(strata_pred$预测率),
  `Mean score` = strata_pred$平均分, row.names = NULL, check.names = FALSE
)
dca_sel <- dca[dca$阈值 %in% c(0.10, 0.15, 0.20, 0.25, 0.30, 0.35), ]
tab4c_zh <- data.frame(
  阈值 = sprintf("%d%%", round(100 * dca_sel$阈值)),
  `模型净获益` = fmt_num(dca_sel$净获益_模型, 4),
  `全部干预净获益` = fmt_num(dca_sel$净获益_全部干预, 4),
  `不干预净获益` = fmt_num(dca_sel$净获益_不干预, 4),
  row.names = NULL, check.names = FALSE
)
tab4c_en <- data.frame(
  Threshold = tab4c_zh$阈值, `Model net benefit` = tab4c_zh$`模型净获益`,
  `Treat-all net benefit` = tab4c_zh$`全部干预净获益`, `Treat-none net benefit` = tab4c_zh$`不干预净获益`,
  row.names = NULL, check.names = FALSE
)
tab4_note_zh <- "内部验证：Bootstrap 1,000 次；AUC 的 95% CI 为 Bootstrap 百分位法；校准斜率/截距在 logit 尺度计算。DCA 使用开发样本表观预测，因此净获益可视为略偏乐观的上限；模型在 5%–35% 阈值内净获益为正且优于全部干预策略。"
tab4_note_en <- "Internal validation: 1,000 bootstrap resamples; 95% CI for AUC by percentile bootstrap; calibration slope/intercept on the logit scale. DCA uses apparent predictions in the development sample and is therefore slightly optimistic; net benefit was positive from 5% to 35% and superior to treat-all."

# ---------------------------------------------------------------------------
# Table 5：外部验证
# ---------------------------------------------------------------------------
ext_s0 <- perf[perf$scenario == "external_S0", ]
ext_s1 <- perf[perf$scenario == "external_S1_complete_fu", ]
ext_s2 <- perf[perf$scenario == "external_S2_fixed_le4", ]
ext_strata <- read.csv(file.path(P, "output/12_模型定稿/A14_外部风险分层_S0.csv"),
                       check.names = FALSE, fileEncoding = "UTF-8")
cohort_order <- c("HRS", "ELSA", "SHARE")
ext_s0 <- ext_s0[match(cohort_order, ext_s0$cohort), ]
ext_s1 <- ext_s1[match(cohort_order, ext_s1$cohort), ]
ext_s2 <- ext_s2[match(cohort_order, ext_s2$cohort), ]
tab5a_zh <- data.frame(
  队列 = ext_s0$cohort, 总样本 = fmt_n(ext_s0$n_total), 总事件 = fmt_n(ext_s0$events_total),
  模型样本 = fmt_n(ext_s0$n_model), 模型事件 = fmt_n(ext_s0$events_model),
  连续模型AUC = fmt_ci(ext_s0$auc_model, ext_s0$auc_model_low, ext_s0$auc_model_high),
  整数分AUC = fmt_ci(ext_s0$auc_score, ext_s0$auc_score_low, ext_s0$auc_score_high),
  校准斜率 = fmt_num(ext_s0$cal_slope), 校准截距 = fmt_num(ext_s0$cal_intercept),
  Brier = fmt_num(ext_s0$brier, 4), row.names = NULL, check.names = FALSE
)
tab5a_en <- data.frame(
  Cohort = ext_s0$cohort, `Total n` = fmt_n(ext_s0$n_total), `Total events` = fmt_n(ext_s0$events_total),
  `Model n` = fmt_n(ext_s0$n_model), `Model events` = fmt_n(ext_s0$events_model),
  `AUC continuous` = fmt_ci(ext_s0$auc_model, ext_s0$auc_model_low, ext_s0$auc_model_high),
  `AUC integer score` = fmt_ci(ext_s0$auc_score, ext_s0$auc_score_low, ext_s0$auc_score_high),
  `Calibration slope` = fmt_num(ext_s0$cal_slope), `Calibration intercept` = fmt_num(ext_s0$cal_intercept),
  Brier = fmt_num(ext_s0$brier, 4), row.names = NULL, check.names = FALSE
)
tab5b_zh <- data.frame(
  队列 = cohort_order,
  S1连续AUC = fmt_ci(ext_s1$auc_model, ext_s1$auc_model_low, ext_s1$auc_model_high),
  S1整数分AUC = fmt_ci(ext_s1$auc_score, ext_s1$auc_score_low, ext_s1$auc_score_high),
  S2连续AUC = fmt_ci(ext_s2$auc_model, ext_s2$auc_model_low, ext_s2$auc_model_high),
  S2整数分AUC = fmt_ci(ext_s2$auc_score, ext_s2$auc_score_low, ext_s2$auc_score_high),
  row.names = NULL, check.names = FALSE
)
tab5b_en <- data.frame(
  Cohort = cohort_order, `S1 AUC continuous` = tab5b_zh$S1连续AUC, `S1 AUC score` = tab5b_zh$S1整数分AUC,
  `S2 AUC continuous` = tab5b_zh$S2连续AUC, `S2 AUC score` = tab5b_zh$S2整数分AUC,
  row.names = NULL, check.names = FALSE
)
ratio_rows <- do.call(rbind, lapply(cohort_order, function(co) {
  z <- ext_strata[ext_strata$cohort == co, ]
  lo <- z$rate[z$stratum == "低"]; hi <- z$rate[z$stratum == "高"]
  data.frame(Cohort = co, `Low-risk event rate` = fmt_pct(lo), `High-risk event rate` = fmt_pct(hi),
             `High/low ratio` = sprintf("%.1f×", hi / lo), row.names = NULL, check.names = FALSE)
}))
tab5c_zh <- data.frame(队列 = ratio_rows$Cohort, 低危事件率 = ratio_rows$`Low-risk event rate`,
                       高危事件率 = ratio_rows$`High-risk event rate`, 高危_低危比 = ratio_rows$`High/low ratio`,
                       row.names = NULL, check.names = FALSE)
tab5c_en <- ratio_rows
tab5_note_zh <- "外部队列：HRS w10→w12、ELSA w4→w6、SHARE w4→w6；均按 A14 固定公式验证，不重新拟合。S0=各队列基线全样本 mean−1SD（主分析）；S1=complete-fu mean−1SD；S2=固定 CHARLS ≤4。连续模型分析使用 ADL 完整病例：HRS 缺失 79 人（其中 9 个事件）、SHARE 缺失 5 人（0 个事件）；整数分 AUC 使用全部风险人群。跨队列 I² 85%–95%，校准异质性明显，合并 AUC 仅作描述。"
tab5_note_en <- "External cohorts: HRS w10→w12, ELSA w4→w6, SHARE w4→w6; the A14 fixed formula was applied without refitting. S0 = cohort-specific baseline-all mean−1SD (primary); S1 = complete-fu mean−1SD; S2 = fixed CHARLS ≤4. Continuous-model analyses used complete ADL cases: HRS missing 79 (9 events), SHARE missing 5 (0 events); integer-score AUC used the full risk set. Cross-cohort I² was 85%–95% with substantial calibration heterogeneity; pooled AUC is descriptive only."

# ---------------------------------------------------------------------------
# 输出：CSV / XLSX / DOCX（中文、English）
# ---------------------------------------------------------------------------
write_panel_csv <- function(data, lang, id) {
  safe_write_csv(data, file.path(DATA, paste0(id, "_", lang, ".csv")),
                 row.names = FALSE, fileEncoding = "UTF-8")
}
make_ft <- function(data, lang) {
  ft <- flextable(data)
  ft <- theme_booktabs(ft)
  ft <- font(ft, fontname = if (lang == "ZH") "SimSun" else "Times New Roman", part = "all")
  ft <- fontsize(ft, size = if (nrow(data) > 12) 7 else 8, part = "all")
  ft <- bold(ft, part = "header")
  ft <- autofit(ft)
  ft
}
write_panel_docx <- function(panels, lang, path) {
  doc <- read_docx()
  for (i in seq_along(panels)) {
    z <- panels[[i]]
    doc <- body_add_par(doc, z$title, style = "Normal")
    doc <- body_add_flextable(doc, make_ft(z$data, lang))
    doc <- body_add_par(doc, z$note, style = "Normal")
    if (i < length(panels)) doc <- body_add_break(doc)
  }
  safe_save_docx(doc, path = path)
}
write_panel_xlsx <- function(panels, path) {
  sheets <- setNames(lapply(panels, `[[`, "data"), vapply(panels, `[[`, "", "sheet"))
  safe_write_xlsx(sheets, file = path, asTable = TRUE, overwrite = TRUE)
}

panels_zh <- list(
  list(sheet = "Table1", title = "表 1. 开发队列与外部验证队列的协调化基线特征", data = tab1_zh, note = tab1_note_zh),
  list(sheet = "Table2", title = "表 2. 候选变量筛选与效应量", data = tab2_zh, note = tab2_note_zh),
  list(sheet = "Table3A", title = "表 3A. 最终筛查指数评分规则", data = tab3a_zh, note = tab3_note_zh),
  list(sheet = "Table3B", title = "表 3B. 风险分层", data = tab3b_zh, note = tab3_note_zh),
  list(sheet = "Table4A", title = "表 4A. 内部验证：区分度与校准", data = tab4a_zh, note = tab4_note_zh),
  list(sheet = "Table4B", title = "表 4B. 内部验证：风险分层", data = tab4b_zh, note = tab4_note_zh),
  list(sheet = "Table4C", title = "表 4C. 内部验证：决策曲线净获益", data = tab4c_zh, note = tab4_note_zh),
  list(sheet = "Table5A", title = "表 5A. 外部验证主分析（S0）", data = tab5a_zh, note = tab5_note_zh),
  list(sheet = "Table5B", title = "表 5B. 外部验证敏感性分析（S1/S2）", data = tab5b_zh, note = tab5_note_zh),
  list(sheet = "Table5C", title = "表 5C. 外部验证风险分层梯度", data = tab5c_zh, note = tab5_note_zh)
)
panels_en <- list(
  list(sheet = "Table1", title = "Table 1. Harmonized baseline characteristics of the development and external validation cohorts", data = tab1_en, note = tab1_note_en),
  list(sheet = "Table2", title = "Table 2. Candidate variable screening and effect estimates", data = tab2_en, note = tab2_note_en),
  list(sheet = "Table3A", title = "Table 3A. Final screening index scoring rule", data = tab3a_en, note = tab3_note_en),
  list(sheet = "Table3B", title = "Table 3B. Risk stratification", data = tab3b_en, note = tab3_note_en),
  list(sheet = "Table4A", title = "Table 4A. Internal validation: discrimination and calibration", data = tab4a_en, note = tab4_note_en),
  list(sheet = "Table4B", title = "Table 4B. Internal validation: risk stratification", data = tab4b_en, note = tab4_note_en),
  list(sheet = "Table4C", title = "Table 4C. Internal validation: decision-curve net benefit", data = tab4c_en, note = tab4_note_en),
  list(sheet = "Table5A", title = "Table 5A. External validation: primary analysis (S0)", data = tab5a_en, note = tab5_note_en),
  list(sheet = "Table5B", title = "Table 5B. External validation: sensitivity analyses (S1/S2)", data = tab5b_en, note = tab5_note_en),
  list(sheet = "Table5C", title = "Table 5C. External validation: risk gradient", data = tab5c_en, note = tab5_note_en)
)
for (i in seq_along(panels_zh)) write_panel_csv(panels_zh[[i]]$data, "ZH", panels_zh[[i]]$sheet)
for (i in seq_along(panels_en)) write_panel_csv(panels_en[[i]]$data, "EN", panels_en[[i]]$sheet)
write_panel_xlsx(panels_zh, file.path(ZH, "A16_Tables_ZH.xlsx"))
write_panel_xlsx(panels_en, file.path(EN, "A16_Tables_EN.xlsx"))
write_panel_docx(panels_zh, "ZH", file.path(ZH, "A16_TableCollection_ZH_v2.docx"))
write_panel_docx(panels_en, "EN", file.path(EN, "A16_TableCollection_EN_v2.docx"))
for (i in seq_along(panels_zh)) {
  write_panel_docx(panels_zh[i], "ZH", file.path(ZH, paste0(panels_zh[[i]]$sheet, "_ZH.docx")))
  write_panel_docx(panels_en[i], "EN", file.path(EN, paste0(panels_en[[i]]$sheet, "_EN.docx")))
}
manifest <- data.frame(
  表格 = vapply(panels_zh, `[[`, "", "sheet"),
  中文标题 = vapply(panels_zh, `[[`, "", "title"),
  English_title = vapply(panels_en, `[[`, "", "title"),
  中文文件 = paste0(ZH, "/", vapply(panels_zh, `[[`, "", "sheet"), "_ZH.docx"),
  English_file = paste0(EN, "/", vapply(panels_en, `[[`, "", "sheet"), "_EN.docx"),
  row.names = NULL, check.names = FALSE
)
safe_write_csv(manifest, file.path(DATA, "A16_Table_Index.csv"), row.names = FALSE, fileEncoding = "UTF-8")
cat("A16 tables written.\n")
