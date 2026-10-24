# ============================================================================
# 31_A16_tables.R —— A16: Table 1–5 Chinese and English submission versions
# input: A07/A14/A15, A08/A09, A12/A13 frozen results
# output: output/14_A16_figures_and_tables/{Chinese,English,data}/
# generate: 2026-09-27
# ============================================================================
invisible(try(Sys.setlocale("LC_ALL", "zh_CN.UTF-8"), silent = TRUE))
.libPaths(c("D:/R/library", .libPaths()))
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(openxlsx)
  library(officer); library(flextable)
})
P <- "D:/BaiduSyncdisk/research/Xu_Zhonghao"
OUT <- file.path(P, "output", "14_A16_figures_and_tables")
ZH <- file.path(OUT, "Chinese"); EN <- file.path(OUT, "English"); DATA <- file.path(OUT, "data")
for (d in c(OUT, ZH, EN, DATA)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# general formatting functions
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
  stop("cannot write file(retry 10 times): ", file)
}
safe_write_xlsx <- function(x, file, ...) {
  for (i in seq_len(10)) {
    ok <- tryCatch({ openxlsx::write.xlsx(x, file = file, ...); TRUE }, error = function(e) FALSE)
    if (ok) return(invisible(TRUE))
    Sys.sleep(0.5)
  }
  stop("cannot write file(retry 10 times): ", file)
}
safe_save_docx <- function(x, path) {
  for (i in seq_len(10)) {
    ok <- tryCatch({ print(x, target = path); TRUE }, error = function(e) { message("DOCX incorrect: ", conditionMessage(e)); FALSE })
    if (ok) return(invisible(TRUE))
    Sys.sleep(0.5)
  }
  stop("cannot write file(retry 10 times): ", path)
}

# ---------------------------------------------------------------------------
# Table 1: cohort baseline
# ---------------------------------------------------------------------------
a07 <- readRDS(file.path(P, "output/10_screening_index/screening_index_memory_v3_2026-09-26.rds"))
ch <- a07$data %>% mutate(event = incident_2015, mem_base = mem2011)
cohort_data <- list(
  CHARLS = list(d = ch, role_zh = "development cohort", role_en = "Development",
                cut_zh = "≤4(fixed cutoff, raw 4.765)", cut_en = "≤4 (fixed, raw 4.765)"),
  HRS = list(d = NULL, role_zh = "external validation", role_en = "External validation",
             cut_zh = NULL, cut_en = NULL),
  ELSA = list(d = NULL, role_zh = "external validation", role_en = "External validation",
              cut_zh = NULL, cut_en = NULL),
  SHARE = list(d = NULL, role_zh = "external validation", role_en = "External validation",
               cut_zh = NULL, cut_en = NULL)
)
for (co in c("HRS", "ELSA", "SHARE")) {
  x <- readRDS(file.path(P, "external_data", co, "analysis",
                         paste0(co, "_sensitivity_risk_sets_2026-09-26.rds")))[[paste0(co, "__baseline_all_mean_sd")]]
  cut_raw <- unique(x$cut)[1]; cut_int <- unique(x$cut_int)[1]
  cohort_data[[co]]$d <- x
  cohort_data[[co]]$cut_zh <- sprintf("≤%d(cohort mean−1SD, raw %.2f)", cut_int, cut_raw)
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
    fu = "4 year"
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
rows_zh <- c(role = "cohort role", n = "sample size", age = "age(years old)", age65 = "age≥65 years old",
             edu_lt = "education: no schooling", edu_primary = "education: primary school", edu_middle = "education: middle school",
             edu_high = "education_high_school_or_above", dep = "depression(cohort specific cutoff)",
             adl = "ADL limited", mem = "baseline memory score(0–20)", events = "follow up events, n(%)",
             cut = "cognitive impairment cutoff", fu = "follow up duration")
rows_en <- c(role = "Cohort role", n = "Sample size", age = "Age, years", age65 = "Age ≥65 years",
             edu_lt = "Education: less than primary", edu_primary = "Education: primary school",
             edu_middle = "Education: middle school", edu_high = "Education: high school or above",
             dep = "Depression (cohort-specific threshold)", adl = "ADL limitation",
             mem = "Baseline memory score (0–20)", events = "Follow-up events, n (%)",
             cut = "Impairment cutpoint", fu = "Follow-up duration")

tab1_zh <- data.frame(characteristics = rows_zh, CHARLS = cohort_values$CHARLS[names(rows_zh)],
                      HRS = cohort_values$HRS[names(rows_zh)], ELSA = cohort_values$ELSA[names(rows_zh)],
                      SHARE = cohort_values$SHARE[names(rows_zh)], row.names = NULL, check.names = FALSE)
tab1_en <- data.frame(Characteristic = rows_en, CHARLS = cohort_values_en$CHARLS[names(rows_en)],
                      HRS = cohort_values_en$HRS[names(rows_en)], ELSA = cohort_values_en$ELSA[names(rows_en)],
                      SHARE = cohort_values_en$SHARE[names(rows_en)], row.names = NULL, check.names = FALSE)
tab1_note_zh <- "abbreviations: ADL, activities of daily living; CES-D, Center for Epidemiologic Studies Depression Scale; EURO-D, European Depression Scale.external depression cutoffs are HRS/ELSA CES-D-8 ≥4, SHARE EURO-D ≥4; CHARLS is CES-D-10 ≥10.external cohort ADL missing: HRS 79 participants, SHARE 5 participants.external events based on S0(full baseline sample in each cohort mean−1SD)cutoff definition; CHARLS is fixed cutoff ≤4.sex omitted from this table, external validation data did not retain sex, sex not a final index variable."
tab1_note_en <- "Abbreviations: ADL, activities of daily living; CES-D, Center for Epidemiologic Studies Depression Scale; EURO-D, European Depression Scale. Depression thresholds: HRS/ELSA CES-D-8 ≥4, SHARE EURO-D ≥4, CHARLS CES-D-10 ≥10. Missing ADL: HRS 79, SHARE 5. External events use the S0 cutpoint (cohort-specific baseline mean−1SD); CHARLS uses fixed ≤4. Sex was not included because it was not retained in the harmonized external validation datasets and is not part of the final index."

# ---------------------------------------------------------------------------
# Table 2: candidate variable screening and effect estimates
# ---------------------------------------------------------------------------
orm <- read.csv(file.path(P, "output/03_five_wave_cross_sectional/separate_wave_analyses_Logistic_OR_matrix_2026_09_19.csv"),
                check.names = FALSE, fileEncoding = "UTF-8")
risk <- read.csv(file.path(P, "output/04_longitudinal_cohort/cohort_risk_2026-09-13.csv"),
                 check.names = FALSE, fileEncoding = "UTF-8")
coef_final <- read.csv(file.path(P, "output/12_final_model/A14_final_model_coefficients.csv"),
                       check.names = FALSE, fileEncoding = "UTF-8")
iadl2 <- read.csv(file.path(P, "output/04_longitudinal_cohort/iadl2_longitudinal_models_2026-09-26.csv"),
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
  sprintf("OR %s; HR %s; change β %s → %s(after baseline adjustment)",
          en_dash(getv("logistic")), en_dash(getv("cox")), en_dash(getv("change")), en_dash(getv("change_adj")))
}
iadl2_str_en <- function() {
  getv <- function(model) iadl2$unweighted[iadl2$exposure == "iadl2" & iadl2$model == model]
  sprintf("OR %s; HR %s; change β %s → %s after baseline adjustment",
          en_dash(getv("logistic")), en_dash(getv("cox")), en_dash(getv("change")), en_dash(getv("change_adj")))
}
tab2_zh <- data.frame(
  variables = c("age(per 5 years old)", "education(high school or above vs no schooling)", "depression(CES-D-10 ≥10)",
           "ADL limited", "IADL limited(supplementary)", "sleep <6 h", "sleep >8 h", "regular physical activity",
           "stroke", "diabetes", "hypertension", "pain", "PLT/WBC/CRP, Neu/NLR"),
  role = c("core", "core", "core", "core", "supplementary", "incremental", "incremental", "incremental",
           "incremental", "incremental", "background/exploration", "exclude", "exclude"),
  `final_model_OR` = c(final_or5("age"), final_or("edu_groupHigh school or above"), final_or("dep"), final_or("adl"),
                  "—", "—", "—", "—", "—", "—", "—", "—", "—"),
  `cross_sectional_OR_range` = c("—", "—", or_range("depression"), or_range("ADL limited"), or_range("IADL limited"),
                 or_range("sleep <6"), or_range("sleep >8"), or_range("regular physical activity"),
                 or_range("stroke"), or_range("diabetes"), or_range("hypertension"), or_range("bodily pain"), "—"),
  longitudinal_estimate = c("—", "—", risk_est("depression(CES-D-10 ≥10)"), risk_est("ADL limited"),
             iadl2_str_zh(), risk_est("sleep duration(per 1 hours)"), "—", risk_est("regular physical activity"),
             "—", risk_est("diabetes"), risk_est("hypertension"), risk_est("pain"), "no robust association(Neu/NLR unavailable)"),
  sensitivity_analysis = c("final model", "final model", "weighted/IPW consistent direction", "weighted/IPW consistent direction",
           "change model indicates regression to mean, excluded from primary score", "weaker in original analysis", "weaker in original analysis",
           "unstable significance after weighting", "original analysis exploratory", "weaker in original analysis", "excluded from core score",
           "weighted/IPW unstable direction", "excluded from core tool"),
  final_handling = c("included in continuous model and score table(per 5 years old +2)", "included in continuous model and score table(0/−2/−7/−10)",
           "included in continuous model and score table(+3)", "retain in continuous model; final memory outcome β≈0, score points 0",
           "supplementary analysis only", "incremental/exploratory", "incremental/exploratory", "incremental/exploratory", "incremental/exploratory",
           "incremental/exploratory", "not part of core score", "excluded from core score", "excluded from core score"),
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
  `Cross-sectional OR range` = c("—", "—", or_range("depression"), or_range("ADL limited"), or_range("IADL limited"),
                 or_range("sleep <6"), or_range("sleep >8"), or_range("regular physical activity"),
                 or_range("stroke"), or_range("diabetes"), or_range("hypertension"), or_range("bodily pain"), "—"),
  `Longitudinal estimate` = c("—", "—", risk_est("depression(CES-D-10 ≥10)"), risk_est("ADL limited"),
             iadl2_str_zh(), risk_est("sleep duration(per 1 hours)"), "—", risk_est("regular physical activity"),
             "—", risk_est("diabetes"), risk_est("hypertension"), risk_est("pain"), "No robust association (Neu/NLR unavailable)"),
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
tab2_note_zh <- "cross sectional OR range is 2011–2020 five waves ms10 cognitive impairment Logistic of model OR minimum–maximum; longitudinal estimates from T06: 2020 endpoint assessment OR and first impairment Cox HR(unweighted primary analysis), IADL listed separately iadl2 continuous change before and after adjustment β.“final model OR”from 2011→2015 memory threshold crossing outcome A14 frozen model; age per 5 years old OR, education high school or above vs no schooling.cross sectional results of other candidates/longitudinal estimates from T06/T07, for screening and stability display only."
tab2_note_en <- "Cross-sectional OR ranges are the minimum–maximum ORs across the 2011–2020 five-wave ms10 impairment logistic models. Longitudinal estimates are from T06: OR for the 2020 endpoint and HR for first impairment (unweighted primary analysis); IADL also shows unadjusted and baseline-adjusted change β from the iadl2 analysis. “Final model OR” is from the A14 frozen model for the 2011→2015 memory cross-threshold endpoint; age is per 5 years and education is high school or above vs less than primary. Other candidate estimates are from T06/T07 and are shown for screening and stability only."

# ---------------------------------------------------------------------------
# Table 3: final scoring rules and risk strata
# ---------------------------------------------------------------------------
score_tab <- read.csv(file.path(P, "output/12_final_model/A14_final_score_table.csv"),
                      check.names = FALSE, fileEncoding = "UTF-8")
strata_tab <- read.csv(file.path(P, "output/12_final_model/A14_risk_stratification.csv"),
                       check.names = FALSE, fileEncoding = "UTF-8")
tab3a_zh <- score_tab %>% transmute(project, class, point_value, notes)
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
tab3b_zh <- strata_tab %>% transmute(risk_stratum, risk_stratum_definition, n = fmt_n(n), events = fmt_n(events),
                                     event_rate, mean_score)
tab3b_en <- data.frame(
  `Risk stratum` = c("Low", "Intermediate", "High"),
  `Score range` = c("≤−3", "−2 to +4", "≥+5"),
  n = fmt_n(strata_tab$n), Events = fmt_n(strata_tab$events),
  `Event rate` = strata_tab$event_rate, `Mean score` = strata_tab$mean_score,
  row.names = NULL, check.names = FALSE
)
tab3_note_zh <- "total score = age points + education points + depression points; ADL score 0 points.risk strata derived from CHARLS tertiles of development sample integer score: low ≤−3, in −2～+4, high ≥+5.when applied to external cohorts, HRS/ELSA depression is CES-D-8 ≥4, SHARE is EURO-D ≥4."
tab3_note_en <- "Total score = age score + education score + depression score; ADL contributes 0 points. Risk strata are based on tertiles of the integer score in the CHARLS development sample: low ≤−3, intermediate −2 to +4, high ≥+5. For external application, depression is defined as CES-D-8 ≥4 in HRS/ELSA and EURO-D ≥4 in SHARE."

# ---------------------------------------------------------------------------
# Table 4: internal validation
# ---------------------------------------------------------------------------
perf <- read.csv(file.path(P, "output/12_final_model/A14_performance_summary.csv"),
                 check.names = FALSE, fileEncoding = "UTF-8")
a08 <- read.csv(file.path(P, "output/10_screening_index/A08_bootstrap_optimism_2026-09-26.csv"),
                check.names = FALSE, fileEncoding = "UTF-8")
strata_pred <- read.csv(file.path(P, "output/10_screening_index/A09_strata_observed_vs_predicted_2026-09-26.csv"),
                        check.names = FALSE, fileEncoding = "UTF-8")
dca <- read.csv(file.path(P, "output/10_screening_index/A09_dca_2026-09-26.csv"),
                check.names = FALSE, fileEncoding = "UTF-8")
p_app <- perf[perf$scenario == "internal_apparent", ]
p_boot <- perf[perf$scenario == "internal_bootstrap", ]
tab4a_zh <- data.frame(
  metrics = c("AUC(continuous model, apparent)", "AUC(integer score, apparent)", "AUC(Bootstrap adjustment)",
           "calibration_slope", "calibration_intercept", "Brier_score"),
  estimate = c(fmt_ci(p_app$auc_model, p_app$auc_model_low, p_app$auc_model_high),
           fmt_ci(p_app$auc_score, p_app$auc_score_low, p_app$auc_score_high),
           fmt_ci(p_boot$auc_model, p_boot$auc_model_low, p_boot$auc_model_high),
           sprintf("%s(apparent %s)", fmt_num(p_boot$cal_slope), fmt_num(p_app$cal_slope)),
           sprintf("%s(apparent %s)", fmt_num(p_boot$cal_intercept), fmt_num(p_app$cal_intercept)),
           sprintf("%s(apparent %s)", fmt_num(p_boot$brier, 4), fmt_num(p_app$brier, 4))),
   row.names = NULL, check.names = FALSE
)
tab4a_en <- data.frame(
  Metric = c("AUC (continuous model, apparent)", "AUC (integer score, apparent)",
             "AUC (bootstrap-corrected)", "Calibration slope", "Calibration intercept", "Brier score"),
  Estimate = tab4a_zh$estimate, row.names = NULL, check.names = FALSE
)
tab4b_zh <- data.frame(
  risk_stratum = c("low", "in", "high"), point_range = c("≤−3", "−2～+4", "≥+5"),
  n = fmt_n(strata_pred$n), events = fmt_n(strata_pred$events),
  observed_event_rate = fmt_pct(strata_pred$observed_rate), model_predicted_rate = fmt_pct(strata_pred$predicted_rate),
  mean_score = strata_pred$mean_score, row.names = NULL, check.names = FALSE
)
tab4b_en <- data.frame(
  `Risk stratum` = c("Low", "Intermediate", "High"), `Score range` = c("≤−3", "−2 to +4", "≥+5"),
  n = fmt_n(strata_pred$n), Events = fmt_n(strata_pred$events),
  `Observed event rate` = fmt_pct(strata_pred$observed_rate), `Predicted risk` = fmt_pct(strata_pred$predicted_rate),
  `Mean score` = strata_pred$mean_score, row.names = NULL, check.names = FALSE
)
dca_sel <- dca[dca$threshold %in% c(0.10, 0.15, 0.20, 0.25, 0.30, 0.35), ]
tab4c_zh <- data.frame(
  threshold = sprintf("%d%%", round(100 * dca_sel$threshold)),
  `model_net_benefit` = fmt_num(dca_sel$net_benefit_model, 4),
  `net_benefit_of_intervening_in_all` = fmt_num(dca_sel$net_benefit_intervene_in_all, 4),
  `net_benefit_of_no_intervention` = fmt_num(dca_sel$net_benefit_no_intervention, 4),
  row.names = NULL, check.names = FALSE
)
tab4c_en <- data.frame(
  Threshold = tab4c_zh$threshold, `Model net benefit` = tab4c_zh$`model_net_benefit`,
  `Treat-all net benefit` = tab4c_zh$`net_benefit_of_intervening_in_all`, `Treat-none net benefit` = tab4c_zh$`net_benefit_of_no_intervention`,
  row.names = NULL, check.names = FALSE
)
tab4_note_zh <- "internal validation: Bootstrap 1,000 times; AUC of 95% CI is Bootstrap percentile method; calibration slope/intercept in logit scale calculation.DCA uses apparent development sample predictions, net benefit is a mildly optimistic upper bound; model within 5%–35% thresholds has positive net benefit exceeding intervene all."
tab4_note_en <- "Internal validation: 1,000 bootstrap resamples; 95% CI for AUC by percentile bootstrap; calibration slope/intercept on the logit scale. DCA uses apparent predictions in the development sample and is therefore slightly optimistic; net benefit was positive from 5% to 35% and superior to treat-all."

# ---------------------------------------------------------------------------
# Table 5: external validation
# ---------------------------------------------------------------------------
ext_s0 <- perf[perf$scenario == "external_S0", ]
ext_s1 <- perf[perf$scenario == "external_S1_complete_fu", ]
ext_s2 <- perf[perf$scenario == "external_S2_fixed_le4", ]
ext_strata <- read.csv(file.path(P, "output/12_final_model/A14_external_risk_stratification_S0.csv"),
                       check.names = FALSE, fileEncoding = "UTF-8")
cohort_order <- c("HRS", "ELSA", "SHARE")
ext_s0 <- ext_s0[match(cohort_order, ext_s0$cohort), ]
ext_s1 <- ext_s1[match(cohort_order, ext_s1$cohort), ]
ext_s2 <- ext_s2[match(cohort_order, ext_s2$cohort), ]
tab5a_zh <- data.frame(
  cohort = ext_s0$cohort, total_sample = fmt_n(ext_s0$n_total), total_events = fmt_n(ext_s0$events_total),
  model_sample = fmt_n(ext_s0$n_model), model_events = fmt_n(ext_s0$events_model),
  continuous_model_AUC = fmt_ci(ext_s0$auc_model, ext_s0$auc_model_low, ext_s0$auc_model_high),
  integer_score_AUC = fmt_ci(ext_s0$auc_score, ext_s0$auc_score_low, ext_s0$auc_score_high),
  calibration_slope = fmt_num(ext_s0$cal_slope), calibration_intercept = fmt_num(ext_s0$cal_intercept),
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
  cohort = cohort_order,
  S1_continuous_AUC = fmt_ci(ext_s1$auc_model, ext_s1$auc_model_low, ext_s1$auc_model_high),
  S1_integer_score_AUC = fmt_ci(ext_s1$auc_score, ext_s1$auc_score_low, ext_s1$auc_score_high),
  S2_continuous_AUC = fmt_ci(ext_s2$auc_model, ext_s2$auc_model_low, ext_s2$auc_model_high),
  S2_integer_score_AUC = fmt_ci(ext_s2$auc_score, ext_s2$auc_score_low, ext_s2$auc_score_high),
  row.names = NULL, check.names = FALSE
)
tab5b_en <- data.frame(
  Cohort = cohort_order, `S1 AUC continuous` = tab5b_zh$S1_continuous_AUC, `S1 AUC score` = tab5b_zh$S1_integer_score_AUC,
  `S2 AUC continuous` = tab5b_zh$S2_continuous_AUC, `S2 AUC score` = tab5b_zh$S2_integer_score_AUC,
  row.names = NULL, check.names = FALSE
)
ratio_rows <- do.call(rbind, lapply(cohort_order, function(co) {
  z <- ext_strata[ext_strata$cohort == co, ]
  lo <- z$rate[z$stratum == "low"]; hi <- z$rate[z$stratum == "high"]
  data.frame(Cohort = co, `Low-risk event rate` = fmt_pct(lo), `High-risk event rate` = fmt_pct(hi),
             `High/low ratio` = sprintf("%.1f×", hi / lo), row.names = NULL, check.names = FALSE)
}))
tab5c_zh <- data.frame(cohort = ratio_rows$Cohort, low_risk_event_rate = ratio_rows$`Low-risk event rate`,
                       high_risk_event_rate = ratio_rows$`High-risk event rate`, high_risk_ratio_to_low_risk = ratio_rows$`High/low ratio`,
                       row.names = NULL, check.names = FALSE)
tab5c_en <- ratio_rows
tab5_note_zh <- "external cohort: HRS w10→w12, ELSA w4→w6, SHARE w4→w6; all according to A14 fixed formula validation, without refitting.S0=full baseline sample in each cohort mean−1SD(primary analysis); S1=complete-fu mean−1SD; S2=fixed CHARLS ≤4.continuous model analysis uses ADL complete cases: HRS missing 79 participants(of which 9 events), SHARE missing 5 participants(0 events); integer score AUC use the full population at risk.cross cohort I² 85%–95%, substantial calibration heterogeneity, combine AUC descriptive only."
tab5_note_en <- "External cohorts: HRS w10→w12, ELSA w4→w6, SHARE w4→w6; the A14 fixed formula was applied without refitting. S0 = cohort-specific baseline-all mean−1SD (primary); S1 = complete-fu mean−1SD; S2 = fixed CHARLS ≤4. Continuous-model analyses used complete ADL cases: HRS missing 79 (9 events), SHARE missing 5 (0 events); integer-score AUC used the full risk set. Cross-cohort I² was 85%–95% with substantial calibration heterogeneity; pooled AUC is descriptive only."

# ---------------------------------------------------------------------------
# output: CSV / XLSX / DOCX(Chinese, English)
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
  list(sheet = "Table1", title = "table 1. harmonized baseline characteristics of development and external cohorts", data = tab1_zh, note = tab1_note_zh),
  list(sheet = "Table2", title = "table 2. candidate variable screening and effect estimates", data = tab2_zh, note = tab2_note_zh),
  list(sheet = "Table3A", title = "table 3A. final screening index scoring rules", data = tab3a_zh, note = tab3_note_zh),
  list(sheet = "Table3B", title = "table 3B. risk stratification", data = tab3b_zh, note = tab3_note_zh),
  list(sheet = "Table4A", title = "table 4A. internal validation: discrimination and calibration", data = tab4a_zh, note = tab4_note_zh),
  list(sheet = "Table4B", title = "table 4B. internal validation: risk stratification", data = tab4b_zh, note = tab4_note_zh),
  list(sheet = "Table4C", title = "table 4C. internal validation: decision curve net benefit", data = tab4c_zh, note = tab4_note_zh),
  list(sheet = "Table5A", title = "table 5A. primary external validation(S0)", data = tab5a_zh, note = tab5_note_zh),
  list(sheet = "Table5B", title = "table 5B. external validation sensitivity(S1/S2)", data = tab5b_zh, note = tab5_note_zh),
  list(sheet = "Table5C", title = "table 5C. external risk stratification gradient", data = tab5c_zh, note = tab5_note_zh)
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
  tables = vapply(panels_zh, `[[`, "", "sheet"),
  Chinese_title = vapply(panels_zh, `[[`, "", "title"),
  English_title = vapply(panels_en, `[[`, "", "title"),
  Chinese_file = paste0(ZH, "/", vapply(panels_zh, `[[`, "", "sheet"), "_ZH.docx"),
  English_file = paste0(EN, "/", vapply(panels_en, `[[`, "", "sheet"), "_EN.docx"),
  row.names = NULL, check.names = FALSE
)
safe_write_csv(manifest, file.path(DATA, "A16_Table_Index.csv"), row.names = FALSE, fileEncoding = "UTF-8")
cat("A16 tables written.\n")
