# ============================================================================
# 28_A14_final_freeze_checklist.R —— A14 final frozen inventory(one page version)
# input: A14 final model artifact and A14/A11–A14 audit log
# output: A14_final_frozen_inventory_2026_09_26.md, A14_final_frozen_inventory_key_numbers.csv
# ============================================================================
invisible(try(Sys.setlocale("LC_ALL","zh_CN.UTF-8"), silent=TRUE))
.libPaths(c("D:/R/library", .libPaths()))
suppressPackageStartupMessages({ library(dplyr); library(knitr) })

P <- "D:/BaiduSyncdisk/research/Xu_Zhonghao"; OUT <- file.path(P,"output","12_final_model")
a14 <- readRDS(file.path(OUT,"A14_final_model_2026-09-26.rds"))
perf <- read.csv(file.path(OUT,"A14_performance_summary.csv"), check.names=FALSE, fileEncoding="UTF-8")
meta <- read.csv(file.path(OUT,"A14_cross_cohort_meta_summary.csv"), check.names=FALSE, fileEncoding="UTF-8")
strata <- read.csv(file.path(OUT,"A14_risk_stratification.csv"), check.names=FALSE, fileEncoding="UTF-8")

status <- function(file) {
  if (!file.exists(file)) return("MISSING")
  if (any(grepl("RESULT: PASS", readLines(file, warn=FALSE)))) "PASS" else "CHECK"
}
s11 <- status(file.path(P,"output/10_screening_index/A11_v2_independent_verification_2026_09_26.txt"))
s12 <- status(file.path(P,"output/11_external_validation/A12_independent_verification_2026_09_26.txt"))
s13 <- status(file.path(P,"output/11_external_validation/A13_independent_verification_2026_09_26.txt"))
s14 <- status(file.path(OUT,"A14_independent_verification_2026_09_26.txt"))

fmt_ci <- function(a,lo,hi,d=3) sprintf(paste0("%.",d,"f (%.",d,"f–%.",d,"f)"),a,lo,hi)
perf_fmt <- perf %>%
  filter(scenario %in% c("internal_apparent","internal_bootstrap","external_S0")) %>%
  mutate(
    continuous_model_AUC = ifelse(is.na(auc_model),"—",fmt_ci(auc_model,auc_model_low,auc_model_high)),
    integer_score_AUC = ifelse(is.na(auc_score),"—",fmt_ci(auc_score,auc_score_low,auc_score_high)),
    calibration = paste0(ifelse(is.na(cal_slope),"NA",sprintf("%.2f",cal_slope))," / ",
                 ifelse(is.na(cal_intercept),"NA",sprintf("%.2f",cal_intercept))),
    Brier = ifelse(is.na(brier),"NA",sprintf("%.3f",brier))
  ) %>% select(cohort=cohort, scenario=scenario, continuous_model_AUC, integer_score_AUC, calibration, Brier)

meta_fmt <- meta %>% mutate(combine_AUC = fmt_ci(pooled_auc,pooled_low,pooled_high), I2=sprintf("%.1f%%",I2)) %>%
  select(cutoff=scenario, metrics=metric_zh, combine_AUC, I2)

key <- data.frame(
  project = c("final model version","CHARLS population at risk/events","HRS population at risk/events","ELSA population at risk/events","SHARE population at risk/events",
           "internal apparent continuous model AUC","internal apparent integer score AUC","internal Bootstrap adjustment AUC",
           "HRS S0 continuous model AUC","ELSA S0 continuous model AUC","SHARE S0 continuous model AUC",
           "HRS S0 integer score AUC","ELSA S0 integer score AUC","SHARE S0 integer score AUC",
           "S0 cross cohort pooling AUC(continuous model)","S0 cross cohort pooling I²"),
   value = c("v1.0(2026-09-26)","4,441 / 658(14.8%)","14,051 / 1,842(13.1%)","6,725 / 554(8.2%)","26,879 / 2,930(10.9%)",
             "0.702(0.681–0.722)","0.700(0.679–0.721)","0.698(0.679–0.720)",
             "0.696(0.682–0.709)","0.763(0.742–0.784)","0.751(0.742–0.761)",
             "0.692(0.678–0.705)","0.759(0.737–0.780)","0.752(0.742–0.761)",
             "0.737(0.690–0.779)","94.7%"),
   source = c("A14 artifact","A07/A14","A11/A12/A14","A11/A12/A14","A11/A12/A14",
              "A07/A14","A07/A14","A08/A14","A12/A14","A12/A14","A12/A14","A12/A14","A12/A14","A12/A14","A13/A14","A13/A14")
)

md <- c(
"# A14 final frozen inventory(v1.0, 2026-09-26)",
"",
paste0("> **status**: A11 = ",s11,"; A12 = ",s12,"; A13 = ",s13,"; A14 = ",s14,"."),
"> **SHA256**: `A14_artifact_sha256_2026-09-26.txt`(2026-09-26 itemwise verification, all matched).",
"> **freeze rules**: no further variable changes, coefficients, score table, strata and cutoffs; any change requires a new version and revalidation.",
"",
"## one, final continuous model",
"```",
a14$formula,
"```",
"- education reference: no schooling; depression reference: CES-D-10 <10; ADL reference: no limitation.",
"- outcome: CHARLS 2011 memory score >4 participants, 2015 memory score ≤4(4 year threshold crossing).",
"- ADL β≈0, 95%CI across 1; simplified score ADL=0, retain in continuous model ADL.",
"",
"## two, final simplified score table",
"",
"| project | point value |",
"|---|---|",
"| age | 45–49=0; per increase of 5 years old **+2**(no upper limit) |",
"| education | no schooling **0**; primary school **−2**; middle school **−7**; high school or above **−10** |",
"| depression | CES-D-10 ≥10 **+3**; <10 **0** |",
"| ADL | **0**(no points) |",
"",
"**total score = age points + education points + depression points.** total score has no fixed theoretical upper bound; observed in development sample −10～+20, SHARE highest +23.",
"> external cohort harmonization: HRS/ELSA CES-D-8 ≥4; SHARE EURO-D ≥4.",
"",
"## three, risk stratification",
"",
"| risk stratum | total score | CHARLS observed event rate |",
"|---|---|---:|",
"| low risk | ≤−3 | 5.6% |",
"| intermediate risk | −2～+4 | 15.4% |",
"| high risk | ≥+5 | 27.3% |",
"",
"external cohort S0 high risk/low risk: HRS 4.9×, ELSA 10.5×, SHARE 7.3×.",
"",
"## four, key performance",
"",
kable(perf_fmt, format="pipe", row.names=FALSE),
"",
"## five, cross cohort random effects(A13)",
"",
kable(meta_fmt, format="pipe", row.names=FALSE),
"",
"## six, final numbers required for manuscript",
"",
kable(key, format="pipe", row.names=FALSE),
"",
"## seven, permitted and prohibited uses",
"",
"**permitted**: 45 preoperative related cognitive vulnerability in age eligible middle aged and older adults**screening/triage**; recommend further MoCA/MMSE cognitive assessments.",
"",
"**prohibited**: dementia/MCI diagnostics; postoperative delirium/PND/ICU/hospitalization/death prediction; claim of validation in surgical cohorts; treatment or perioperative decision tool.",
"",
"## eight, limitations that must be reported",
"",
"1. outcome is 4 year memory score threshold crossing, includes regression to mean and measurement error, not equivalent to clinical dementia/MCI incidence; ",
"2. substantial external calibration differences: HRS underestimation(intercept +0.33), SHARE overestimation(intercept −0.31), ELSA slope 1.82; ",
"3. cross cohort heterogeneity I² 85%–95%, and k=3; SHARE between countries I² S0/S1≈50%, S2≈77%; ",
"4. ELSA age 90 age top coding; ADL missing HRS 79 participants, SHARE 5 participants; ",
"5. death/attrition/participants missing memory scores excluded, competing risks not addressed; ",
"6. not yet validated in surgical patients; ",
"7. S0/S1 and S2 different threshold meanings, three complementary result sets.",
"",
"## nine, verification files",
"",
"- `A14_artifact_sha256_2026-09-26.txt`",
"- `A14_independent_verification_2026_09_26.txt`",
"- `A13_independent_verification_2026_09_26.txt`",
"- `A12_independent_verification_2026_09_26.txt`",
"- `A11_v2_independent_verification_2026_09_26.txt`"
)
writeLines(md, file.path(OUT,"A14_final_frozen_inventory_2026_09_26.md"), useBytes=TRUE)
write.csv(key, file.path(OUT,"A14_final_frozen_inventory_key_numbers.csv"), row.names=FALSE, fileEncoding="UTF-8")
cat("A14 freeze checklist written.\n")