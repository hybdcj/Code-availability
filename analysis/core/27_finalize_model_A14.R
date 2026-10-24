# ============================================================================
# 27_finalize_model_A14.R —— A14: freeze final model and score table
# input: A07 model RDS, A08 Bootstrap, A11/A12/A13 external validation results
# output: output/12_final_model/
# ============================================================================
invisible(try(Sys.setlocale("LC_ALL","zh_CN.UTF-8"), silent=TRUE))
.libPaths(c("D:/R/library", .libPaths()))
suppressPackageStartupMessages({ library(dplyr); library(pROC); library(knitr) })

P   <- "D:/BaiduSyncdisk/research/Xu_Zhonghao"
OUT <- file.path(P, "output", "12_final_model")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

a07 <- readRDS(file.path(P, "output/10_screening_index/screening_index_memory_v3_2026-09-26.rds"))
a08 <- readRDS(file.path(P, "output/10_screening_index/A08_bootstrap_2026-09-26.rds"))
a12 <- readRDS(file.path(P, "output/11_external_validation/A12_results_2026-09-26.rds"))
a13 <- readRDS(file.path(P, "output/11_external_validation/A13_results_2026-09-26.rds"))

model_u <- a07$model_u
scenario_lab <- c(baseline_all_mean_sd="S0_baseline_all", complete_fu_mean_sd="S1_complete_fu", fixed_charls_le4="S2_fixed_le4")
model_levels <- model_u$xlevels$edu_group
coef_mat <- summary(model_u)$coefficients
ci_mat <- confint.default(model_u)
term_labels <- c(
  "(Intercept)" = "intercept",
  "age" = "age(per year of age)",
  "edu_groupPrimary school" = "education: primary school(vs no schooling)",
  "edu_groupMiddle school" = "education: middle school(vs no schooling)",
  "edu_groupHigh school or above" = "education: high school or above(vs no schooling)",
  "dep" = "depression: CES-D-10 ≥10(vs <10)",
  "adl" = "ADL: any limitation(vs no limitation)"
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

# final continuous model formula
beta <- setNames(coef_table$beta, coef_table$term)
formula_text <- sprintf(
  "logit(P(4 year memory threshold crossing)) = %.4f + %.4f×age + (%.4f)×primary school + (%.4f)×middle school + (%.4f)×high school or above + %.4f×depression + (%.4f)×ADL",
  beta[["(Intercept)"]], beta[["age"]], beta[["edu_groupPrimary school"]],
  beta[["edu_groupMiddle school"]], beta[["edu_groupHigh school or above"]],
  beta[["dep"]], beta[["adl"]]
)

# final score table
score_table <- data.frame(
  project = c("age","education","education","education","education","depression","ADL"),
  class = c("45–49 years old=0; 50–54=+2; per increase of 5 years then +2, no upper limit","no schooling(Less than primary)",
           "primary school(Primary school)","middle school(Middle school)","high school or above(High school or above)",
           "CES-D-10 ≥10","any ADL limited"),
  point_value = c("+2 / per 5 years old","0","−2","−7","−10","+3","0"),
  notes = c("45–49 years equals 0 points, no upper limit","reference group","relative to no schooling","relative to no schooling","relative to no schooling","vs CES-D-10 <10","β≈0, no points; continuous model still retains"),
  row.names = NULL
)

# risk stratification
strata_table <- a07$strata %>%
  mutate(risk_stratum = as.character(stratum),
         risk_stratum_definition = c("≤−3","−2 ～ +4","≥+5")[match(stratum, c("low","in","high"))],
         event_rate_num = as.numeric(sub("%","",event_rate))/100) %>%
  select(risk_stratum, risk_stratum_definition, n, events, event_rate, mean_score)

# internal performance
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

# external risk stratification S0
external_strata <- a12$strata %>% filter(cut_definition == "baseline_all_mean_sd") %>%
  transmute(cohort, stratum = as.character(stratum), n, events, rate, pred_mean, oe_model)

# external meta
meta_summary <- a13$cross %>% mutate(scenario = scenario_lab[cut_definition], metric_zh = ifelse(metric=="model","continuous model","integer score"))

# save CSV
write.csv(coef_table, file.path(OUT,"A14_final_model_coefficients.csv"), row.names=FALSE, fileEncoding="UTF-8")
write.csv(score_table, file.path(OUT,"A14_final_score_table.csv"), row.names=FALSE, fileEncoding="UTF-8")
write.csv(strata_table, file.path(OUT,"A14_risk_stratification.csv"), row.names=FALSE, fileEncoding="UTF-8")
write.csv(performance, file.path(OUT,"A14_performance_summary.csv"), row.names=FALSE, fileEncoding="UTF-8")
write.csv(external_strata, file.path(OUT,"A14_external_risk_stratification_S0.csv"), row.names=FALSE, fileEncoding="UTF-8")
write.csv(meta_summary, file.path(OUT,"A14_cross_cohort_meta_summary.csv"), row.names=FALSE, fileEncoding="UTF-8")

# save final model artifact
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

# format report tables
fmt <- function(x,d=3) ifelse(is.na(x),"NA",sprintf(paste0("%.",d,"f"),x))
coef_tab <- coef_table %>% mutate(
  items = term_zh,
  β = fmt(beta,4), SE = fmt(se,4), OR = fmt(OR,3),
  OR_95CI = paste0(fmt(OR_low,3),"–",fmt(OR_high,3)),
  p = ifelse(p<0.001,"<0.001",fmt(p,3))
) %>% mutate(OR = ifelse(term=="(Intercept)","—",OR),
             OR_95CI = ifelse(term=="(Intercept)","—",OR_95CI)) %>%
  select(items, β, SE, OR, OR_95CI, p)
score_tab <- score_table %>% select(project, class, point_value, notes)
strata_tab <- strata_table
fmt_ci <- function(a,lo,hi,d=3) ifelse(is.na(a), "—", paste0(fmt(a,d), " (", fmt(lo,d), "–", fmt(hi,d), ")"))
perf_tab <- performance %>% filter(scenario %in% c("internal_apparent","internal_bootstrap","external_S0")) %>%
  mutate(
    AUC_model = fmt_ci(auc_model, auc_model_low, auc_model_high),
    AUC_integer_score = fmt_ci(auc_score, auc_score_low, auc_score_high),
    calibration_slope = fmt(cal_slope,2), calibration_intercept = fmt(cal_intercept,2), Brier = fmt(brier,3)
  ) %>% select(cohort=cohort, scenario=scenario, full_risk_set_n=n_total, full_risk_set_events=events_total, model_n=n_model, events=events_model, AUC_model, AUC_integer_score, calibration_slope, calibration_intercept, Brier)
sens_tab <- performance %>% filter(scenario %in% c("external_S1_complete_fu","external_S2_fixed_le4")) %>%
  mutate(
    AUC_model = fmt_ci(auc_model, auc_model_low, auc_model_high),
    AUC_integer_score = fmt_ci(auc_score, auc_score_low, auc_score_high),
    calibration_slope = fmt(cal_slope,2), calibration_intercept = fmt(cal_intercept,2), Brier = fmt(brier,3)
  ) %>% select(cohort=cohort, scenario=scenario, full_risk_set_n=n_total, full_risk_set_events=events_total, model_n=n_model, events=events_model, AUC_model, AUC_integer_score, calibration_slope, calibration_intercept, Brier)
meta_tab <- meta_summary %>% mutate(
  combine_AUC = paste0(fmt(pooled_auc), " (", fmt(pooled_low), "–", fmt(pooled_high), ")"),
  I2 = sprintf("%.1f%%", I2)
) %>% select(cutoff=scenario, metrics=metric_zh, combine_AUC, I2, k)

md <- c(
"# A14: freeze final model and score table(v1.0, 2026-09-26)",
"",
"> **freeze statement**: from 2026-09-26 onward, final model, score table, risk strata and formula frozen.",
"> external cohorts used for validation only, **not refitted, variables not reselected, cutoffs unchanged**; any subsequent changes must be labeled post hoc and revalidated.",
"",
"## one, final continuous model",
"",
paste0("**model formula**: ", formula_text),
"",
kable(coef_tab, format="pipe", row.names=FALSE),
"",
"> reference categories: education=no schooling(Less than primary); depression=CES-D-10 <10; ADL=no limitation.",
"> outcome: CHARLS 2011 year memory score >4 participants, 2015 year memory score ≤4(4 year threshold crossing).",
"> ADL coefficient close to 0, 95%CI across 1, therefore in final simplified score table ADL score 0 points; retain in continuous model ADL to preserve prespecified full model.",
"",
"## two, final simplified score table",
"",
kable(score_tab, format="pipe", row.names=FALSE),
"",
"**total score** = age points + education points + depression points; ADL no points.",
"> age points have no upper limit(per increase of 5 years old +2).total score is not fixed −10～+20; observed development sample range −10～+20, external SHARE observed maximum +23.
> when applied to external cohorts, depression cutoff as prespecified harmonization: HRS/ELSA CES-D-8 ≥4; SHARE EURO-D ≥4(see details in A11 v2.2 mapping table).",
"",
"## three, risk stratification",
"",
kable(strata_tab, format="pipe", row.names=FALSE),
"",
"> strata derived from CHARLS tertiles of development sample integer score: low ≤−3; in −2～+4; high ≥+5.",
"> observed development sample event rate: low 5.6%, in 15.4%, high 27.3%; external cohort S0 event rates in A14_external_risk_stratification_S0.csv.",
"",
"## four, internal and external performance summary",
"",
kable(perf_tab, format="pipe", row.names=FALSE),
"",
"> internal Bootstrap adjusted AUC 0.698(0.679–0.720), calibration slope 0.98, Brier 0.1184. in the table“model n/events”is continuous model complete sample; “full risk set n/events”is integer score AUC sample used.HRS full risk set 14,051/1,842; ELSA 6,725/554; SHARE 26,879/2,930.",
"> external cohort S0 continuous model AUC: HRS 0.696, ELSA 0.763, SHARE 0.751; substantial between cohort calibration differences(HRS underestimation, SHARE overestimation, ELSA slope above one).",
"",
"## five, formal sensitivity S1/S2",
"",
kable(sens_tab, format="pipe", row.names=FALSE),
"",
"> S1=complete-fu mean−1SD; S2=fixed CHARLS ≤4.S2 AUC higher but event rate lower, depends on absolute comparability of memory scores across cohorts; primary analysis remains S0.",
"",
"## six, cross cohort random effects summary(A13)",
"",
kable(meta_tab, format="pipe", row.names=FALSE),
"",
"> cross cohort I² very high(85%–95%), pooled values are descriptive only; manuscript should emphasize cohort results, calibration differences and heterogeneity sources.",
"",
"## seven, final intended use and boundaries",
"",
"**permitted uses**: ",
"- 45 preoperative related cognitive vulnerability in middle aged and older adults aged at least**screening/triage**; ",
"- identify people who may need further MoCA/MMSE cognitive assessment; ",
"- developed using accessible data, screening index with internal and cross cohort external validation.",
"",
"**prohibited uses**: ",
"- not for dementia/MCI diagnostic tool; ",
"- does not predict postoperative delirium, PND, ICU, hospitalization or death; ",
"- do not claim validation in surgical cohorts; ",
"- not for treatment or perioperative decision making.",
"",
"## eight, limitations",
"",
"1. outcome is 4 year memory score threshold crossing, including baseline position, measurement error and regression to mean components, not equivalent to clinical dementia/MCI incidence; ",
"2. external calibration differs by cohort, HRS underestimates risk, SHARE overestimates risk, ELSA calibration slope 1.82; ",
"3. high cross cohort heterogeneity(I² 85%–95%), and k=3; ",
"4. ELSA age 90 age top coding; ADL missing HRS 79 participants, SHARE 5 participants; ",
"5. death/attrition/participants missing memory scores excluded, competing risks not addressed; ",
"6. not yet validated in surgical patients; ",
"7. S0/S1 uses within each cohort mean−1SD, and CHARLS fixed ≤4 not the same absolute outcome; S2 fixed ≤4 assumes absolute memory score comparability across cohorts.three complementary result sets, not interchangeable.",
"",
"## nine, reproducibility and files",
"",
"- model artifact: `A14_final_model_2026-09-26.rds`",
"- coefficients: `A14_final_model_coefficients.csv`",
"- score table: `A14_final_score_table.csv`",
"- risk stratification: `A14_risk_stratification.csv`",
"- performance summary: `A14_performance_summary.csv`, `A14_external_risk_stratification_S0.csv`, `A14_cross_cohort_meta_summary.csv`",
"- Session info: `A14_sessionInfo_2026-09-26.txt`",
"- freeze verification: `A14_artifact_sha256_2026-09-26.txt`",
"- scripts: `code/27_finalize_model_A14.R`",
"- independent verification: `A14_independent_verification_2026_09_26.txt`(from verify_a14 generate)",
"",
"## ten, conclusion",
"",
"> freeze final model and simplified score table, version v1.0(2026-09-26).",
"> final tool is**preoperative related cognitive vulnerability screening/triage tool**, not a postoperative outcome predictor."
)
writeLines(md, file.path(OUT,"A14_final_model_specification_2026_09_26.md"), useBytes=TRUE)

# one page scorecard
scorecard <- c(
"# preoperative related cognitive vulnerability scorecard(v1.0, 2026-09-26)",
"",
"> purpose: 45 preoperative cognitive assessment triage in age eligible middle aged and older adults; not diagnostic, does not predict postoperative delirium/PND.",
"",
"## scoring steps",
"",
"### 1. age points",
"| age | point value |",
"|---|---|",
"| 45–49 years old | 0 |",
"| 50–54 years old | +2 |",
"| 55–59 years old | +4 |",
"| 60–64 years old | +6 |",
"| 65–69 years old | +8 |",
"| 70–74 years old | +10 |",
"| 75–79 years old | +12 |",
"| 80–84 years old | +14 |",
"| 85–89 years old | +16 |",
"| 90–94 years old | +18 |",
"| per additional 5 years old | then +2 |",
"",
"### 2. education points",
"| education | point value |",
"|---|---|",
"| no schooling | 0 |",
"| primary school | −2 |",
"| middle school | −7 |",
"| high school or above | −10 |",
"",
"### 3. depression points",
"| CES-D-10 | point value |",
"|---|---|",
"| <10 | 0 |",
"| ≥10 | +3 |",
"",
"### 4. ADL",
"ADL no points(β≈0); continuous model still retains ADL.",
"",
"### 5. total score and risk strata",
"| risk stratum | total score | CHARLS observed event rate |",
"|---|---|---|",
"| low risk | ≤−3 | 5.6% |",
"| intermediate risk | −2 ～ +4 | 15.4% |",
"| high risk | ≥+5 | 27.3% |",
"",
"**example**: 70 years old, primary school, CES-D-10 ≥10 → 10 + (−2) + 3 = **11 points → high risk**.",
"",
"## important reminder",
"- high risk indicates need for further cognitive assessment only, does not establish diagnosis; ",
"- external cohort calibration differs, interpret individual probabilities cautiously; ",
"- not for predicting postoperative delirium, PND or perioperative outcomes; ",
"- not yet validated in a surgical cohort."
)
writeLines(scorecard, file.path(OUT,"A14_scorecard_2026_09_26.md"), useBytes=TRUE)

cat("A14 finalized.\n")