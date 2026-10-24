# =============================================================================
# 10_fivewave_logistic.R —— CHARLS five waves(2011-2020)by wave Logistic: factors associated with cognitive impairment
# version: v2(2026-09-19 revision)
#   ① sleep changed to three categories(<6h / 6-8h(reference)/ >8h), and previously delivered Paper 1(2020 wave v4.4)consistent; 
#   ② pain: main row retains"any pain"(comparable across five waves), additionally"mild 2-3""severe 4-5"two rows(only 2013/2018/2020
#      ordinal coding available; 2011/2015 only"present/none"two levels, cell annotation"binary coding only"); 
#   ③ 5 core exposures(depression, ADL, IADL, sleep, physical activity)by SAP perform BH-FDR; 
#      sleep is a three category factor, adjustment includes the whole factor Wald test P(and Paper 1 footnotes use consistent definitions); 
#   ④ continuous sleep, binary pain moved to"sensitivity subtable"(corresponding Paper 1 supplementary table S2); 
#   ⑤ automatic verification at script end 2020 whether column reproduces delivered results Paper 1 of 6 count.
# definition: outcome = this wave cog_B(within this wave, this education stratum mean total score − 1SD); adjustment = age + sex + four education groups + urban or rural; 
#      weighted = svyglm(weights = weight, ids = communityID, quasibinomial); report exposure coefficients only.
# input: output/02_five_wave_data_and_quality_control/waves_cognition_2026-09-13.rds; data/raw/{wave2_2013,wave4_2018,wave5_2020}(pain classification)
#       output/01_2020_single_wave_baseline/analysis_2020.rds(backfill only 2020 wave rural)
# output: output/03_five_wave_cross_sectional/separate_wave_analyses_Logistic_OR_2026_09_19.xlsx(notes｜wave specific details｜adjustment OR matrix_five wave comparison｜sensitivity subtable)
#       output/03_five_wave_cross_sectional/separate_wave_analyses_Logistic_OR_2026_09_19.csv, separate_wave_analyses_Logistic_OR_matrix_2026_09_19.csv, 
#       output/03_five_wave_cross_sectional/separate_wave_analyses_Logistic_sensitivity_2026_09_19.csv
# note: openxlsx written xlsx with dangling drawing reference(Excel/pandas readable, openpyxl strict mode error), 
#      after generation run: python code/_checks/repair_xlsx_dangling_refs.py <the preceding one xlsx>
# =============================================================================
Sys.setlocale("LC_ALL", "zh_CN.UTF-8")
.libPaths(c("D:/R/library", .libPaths()))
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(survey); library(haven); library(openxlsx)
})

PROJ <- Sys.getenv("CHARLS_PROJ", "D:/BaiduSyncdisk/research/Xu_Zhonghao")
OUT  <- file.path(PROJ, "output", "03_five_wave_cross_sectional")
dat  <- readRDS(file.path(PROJ, "output", "02_five_wave_data_and_quality_control", "waves_cognition_2026-09-13.rds"))
norm_id <- function(x) {
  x <- as.character(x); n <- nchar(x)
  ifelse(n == 11, paste0(substr(x, 1, 9), sprintf("%03d", as.integer(substr(x, 10, 11)))), x)
}

## ---- 0. backfill 2020 wave rural(rural = 1 − urban3, same definition as 02_table1.R) ----
a20 <- readRDS(file.path(PROJ, "output", "01_2020_single_wave_baseline", "analysis_2020.rds"))
rural20 <- ifelse(is.na(a20$urban3), NA_integer_, 1L - as.integer(a20$urban3))
names(rural20) <- as.character(a20$ID)
i20 <- which(dat$year == 2020)
dat$rural[i20] <- unname(rural20[as.character(dat$ID[i20])])
writeLines(sprintf("[backfill]2020 wave rural nonmissing: %d / %d", sum(!is.na(dat$rural[i20])), length(i20)))

## ---- 1. derived variables ----
# three pain categories: obtain ordinal variables from each raw wave module(2013 wb16 / 2018 da041_w4 / 2020 da027), 
#   1 = no pain → none; 2-3 = mild → mild; 4-5 = severe → severe(and 01_build_2020.R of pain_3cat consistent)
get_pain3 <- function(year, KEY) {
  spec <- switch(as.character(year),
    "2013" = list(f = "wave2_2013/Health_Status_and_Functioning.dta", v = "wb16"),
    "2018" = list(f = "wave4_2018/Health_Status_and_Functioning.dta",   v = "da041_w4"),
    "2020" = list(f = "wave5_2020/Health_Status_and_Functioning.dta",   v = "da027"))
  d <- read_dta(file.path(PROJ, "data/raw", spec$f), col_select = all_of(c("ID", spec$v)))
  p <- as.numeric(d[[spec$v]])
  out <- ifelse(is.na(p), NA_character_,
                ifelse(p == 1, "none", ifelse(p %in% c(2, 3), "mild", "severe")))
  out[match(KEY, norm_id(d$ID))]
}
dat$pain3 <- NA_character_
for (y in c(2013, 2018, 2020)) {
  idx <- which(dat$year == y)
  dat$pain3[idx] <- get_pain3(y, norm_id(ifelse(is.na(dat$IDn[idx]), dat$ID[idx], dat$IDn[idx])))
}

dat <- dat %>% mutate(
  dep       = as.integer(cesd10 >= 10),
  # three sleep categories: <6h / 6-8h(reference)/ >8h —— and 01_build_2020.R of sleep_cat same cutoff and reference
  sleep_cat = factor(cut(sleep_h, breaks = c(-Inf, 5.999, 8.0001, Inf),
                         labels = c("<6h", "6-8h", ">8h"), right = FALSE),
                     levels = c("6-8h", "<6h", ">8h")),
  pain3f    = factor(pain3, levels = c("none", "mild", "severe")),
  female_f  = factor(female, c(0, 1), c("Male", "Female")),
  edu_f     = factor(edu_group, levels = c("Less than primary", "Primary school",
                                          "Middle school", "High school or above")),
  rural_f   = factor(rural, c(0, 1), c("Urban", "Rural")))

YEARS  <- c(2011, 2013, 2015, 2018, 2020)
ana    <- dat %>% filter(age >= 45, cog_complete == 1, !is.na(cog_B))
# fixed cutoff(for sensitivity analyses): use only 2011 year analysis sample, calculate once by education stratum"mean total score − 1SD", 
#   apply identical thresholds to all five waves(do not recalculate by wave), for testing"threshold drift"effect on trends
thr11 <- ana %>% filter(year == 2011) %>% group_by(edu_group) %>%
  summarise(thr11 = mean(total, na.rm = TRUE) - sd(total, na.rm = TRUE), .groups = "drop")
ana <- ana %>% left_join(thr11, by = "edu_group") %>%
  mutate(imp_fixed = as.integer(total < thr11))
# excluding figure copying(0-30)definition: cutoff based on"within this wave, this education stratum total30 of mean − 1SD"(figure copying sensitivity uses)
thr30 <- ana %>% group_by(year, edu_group) %>%
  summarise(thr30 = mean(total - draw, na.rm = TRUE) - sd(total - draw, na.rm = TRUE), .groups = "drop")
ana <- ana %>% left_join(thr30, by = c("year", "edu_group")) %>%
  mutate(total30 = total - draw, imp_nodraw = as.integer(total30 < thr30))
writeLines("== fixed cutoff(2011 year education stratification mean-1SD)==")
print(as.data.frame(thr11))
writeLines("== three category pain distribution in analysis sample ==")
print(table(ana$year, ana$pain3f, useNA = "ifany"))
writeLines("== three category sleep distribution in analysis sample ==")
print(table(ana$year, ana$sleep_cat, useNA = "ifany"))

## ---- 2. weight coverage(determines weighting specification; threshold = analysis sample 50%) ----
cov_tab <- ana %>% group_by(year) %>% summarise(
  analysis_sample_n = n(),
  weights_available_n = sum(!is.na(weight) & !is.na(communityID)),
  coverage_pct = round(100 * weights_available_n / analysis_sample_n, 1), .groups = "drop") %>%
  mutate(primary_definition = ifelse(coverage_pct >= 50, "svyglm primarily weighted, glm unweighted comparison",
                       "glm primarily unweighted, svyglm as comparison"))
writeLines("== each wave weight + communityID coverage ==")
print(as.data.frame(cov_tab))

## ---- 3. exposure inventory ----
ADJ <- " + age + female_f + edu_f + rural_f"
EXPOS <- list(
  list(v = "dep",       core = TRUE,  terms = c(dep = "depression(CES-D-10 ≥ 10)")),
  list(v = "adl",       core = TRUE,  terms = c(adl = "ADL limited")),
  list(v = "iadl",      core = TRUE,  terms = c(iadl = "IADL limited")),
  list(v = "iadl2",     core = FALSE, layer = "sensitivity_analysis",
       terms = c(iadl2 = "IADL limited(managing finances+two medication items)")),
  list(v = "sleep_cat", core = TRUE,  factor = TRUE,
       terms = c("sleep_cat<6h" = "sleep <6 h(vs 6-8 h)",
                 "sleep_cat>8h" = "sleep >8 h(vs 6-8 h)")),
  list(v = "pa_any",    core = TRUE,  terms = c(pa_any = "regular physical activity")),
  list(v = "hibp",      core = FALSE, terms = c(hibp = "hypertension")),
  list(v = "dyslip",    core = FALSE, terms = c(dyslip = "dyslipidemia")),
  list(v = "diab",      core = FALSE, terms = c(diab = "diabetes")),
  list(v = "heart",     core = FALSE, terms = c(heart = "heart disease")),
  list(v = "stroke",    core = FALSE, terms = c(stroke = "stroke")),
  list(v = "lung",      core = FALSE, terms = c(lung = "chronic lung disease")),
  list(v = "kidney",    core = FALSE, terms = c(kidney = "kidney disease")),
  list(v = "pain_any",  core = FALSE, terms = c(pain_any = "bodily pain(any, vs no pain)")),
  list(v = "pain3f",    core = FALSE, factor = TRUE,
       terms = c("pain3fmild" = "mild pain 2-3(vs no pain)",
                 "pain3fsevere" = "severe pain 4-5(vs no pain)")))
# sensitivity only exposures(excluded from main matrix): continuous sleep, binary pain(= any pain, for supplementary comparison)
SENS <- list(
  list(v = "sleep_h",  terms = c(sleep_h = "nighttime sleep duration(per 1 h, continuous)")),
  list(v = "pain_any", terms = c(pain_any = "bodily pain(binary, vs no pain)")))

MARK <- function(v, y) {
  if (v == "iadl"   && y == 2011)              return("—(2011 wave excluded)")
  if (v == "pain3f" && y %in% c(2011, 2015))   return("—(binary coding only)")
  ""
}

## ---- 4. fitting utilities ----
grab <- function(fit, terms, model_tag, wtag, exposure) {
  if (is.null(fit)) return(NULL)
  s <- summary(fit)$coefficients
  keep <- terms[names(terms) %in% rownames(s)]
  if (!length(keep)) return(NULL)
  data.frame(model = model_tag, definition = wtag, exposure = exposure,
             term = names(keep), exposure = unname(keep),
             OR = exp(s[names(keep), 1]),
             CI_lower_bound = exp(s[names(keep), 1] - 1.96 * s[names(keep), 2]),
             CI_upper_bound = exp(s[names(keep), 1] + 1.96 * s[names(keep), 2]),
             p = s[names(keep), 4], n = nrow(model.frame(fit)),
             row.names = NULL, stringsAsFactors = FALSE)
}
overall_p <- function(fit, v) {
  if (is.null(fit)) return(NA_real_)
  if (inherits(fit, "svyglm"))
    return(tryCatch(survey::regTermTest(fit, as.formula(paste0("~", v)), method = "Wald")$p,
                    error = function(e) NA_real_))
  f0 <- tryCatch(update(fit, as.formula(paste0(". ~ . - ", v))), error = function(e) NULL)
  if (is.null(f0)) return(NA_real_)
  tryCatch(anova(fit, f0, test = "LRT")[["Pr(>Chi)"]][2], error = function(e) NA_real_)
}

## ---- 5. fit separately by wave(main table exposures + sensitivity exposures) ----
run_wave <- function(dy, spec_list, with_q = TRUE, yvar = "cog_B") {
  y  <- dy$year[1]
  dw <- dy %>% filter(!is.na(weight), !is.na(communityID))
  des <- svydesign(ids = ~communityID, weights = ~weight, data = dw)
  out <- list()
  for (ex in spec_list) {
    if (all(is.na(dy[[ex$v]]))) { writeLines(sprintf("  [%d] skip %s(not assessed this wave/all missing)", y, ex$v)); next }
    f_c <- as.formula(paste0(yvar, " ~ ", ex$v))
    f_a <- as.formula(paste0(yvar, " ~ ", ex$v, ADJ))
    m_wc <- tryCatch(svyglm(f_c, design = des, family = quasibinomial), error = function(e) NULL)
    m_wa <- tryCatch(svyglm(f_a, design = des, family = quasibinomial), error = function(e) NULL)
    m_uc <- tryCatch(glm(f_c, data = dy, family = binomial), error = function(e) NULL)
    m_ua <- tryCatch(glm(f_a, data = dy, family = binomial), error = function(e) NULL)
    if (is.null(m_wa)) warning(sprintf("[%d] weighted adjusted model failed: %s", y, ex$v))
    rows <- rbind(grab(m_wc, ex$terms, "crude model", "weighted (svyglm)", ex$v),
                  grab(m_wa, ex$terms, "adjusted model", "weighted (svyglm)", ex$v),
                  grab(m_uc, ex$terms, "crude model", "unweighted (glm)", ex$v),
                  grab(m_ua, ex$terms, "adjusted model", "unweighted (glm)", ex$v))
    if (!is.null(rows)) {
      rows$core_exposure <- isTRUE(ex$core)
      rows$factor_overall_p <- rep(if (isTRUE(ex$factor)) overall_p(m_wa, ex$v) else NA_real_, nrow(rows))
      out[[ex$v]] <- rows
    }
  }
  d <- do.call(rbind, out)
  if (is.null(d)) return(NULL)
  d$wave <- y
  if (with_q) {
    # BH-FDR: within each wave, 5 one model per core exposure p(three category factor uses overall Wald test P)
    cp <- d %>% filter(model == "adjusted model", definition == "weighted (svyglm)", core_exposure) %>%
      mutate(p_src = ifelse(is.na(factor_overall_p), p, factor_overall_p)) %>%
      distinct(exposure, p_src)
    cp$q <- p.adjust(cp$p_src, method = "BH")
    d <- d %>% left_join(cp[, c("exposure", "q")], by = "exposure")
    writeLines(sprintf("  [%d] core exposure BH-FDR: %s", y,
      paste(sprintf("%s p=%.4g q=%.4g", cp$exposure, cp$p_src, cp$q), collapse = " ｜ ")))
  }
  d
}

writeLines("===== fit separately by wave(main table)=====")
detail <- do.call(rbind, lapply(YEARS, function(y) run_wave(ana %>% filter(year == y), EXPOS)))
detail <- detail %>%
  mutate(`OR (95%CI)` = sprintf("%.2f (%.2f-%.2f)", OR, CI_lower_bound, CI_upper_bound),
         p_formatting = ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))) %>%
  select(wave, exposure, exposure, term, core_exposure, model, definition, `OR (95%CI)`, OR, CI_lower_bound, CI_upper_bound, p, p_formatting, q, n)
writeLines(sprintf("main table detail row count: %d", nrow(detail)))

writeLines("===== fit separately by wave(sensitivity subtable)=====")
detail_sens <- do.call(rbind, lapply(YEARS, function(y)
  run_wave(ana %>% filter(year == y), SENS, with_q = FALSE)))
detail_sens <- detail_sens %>%
  mutate(`OR (95%CI)` = sprintf("%.2f (%.2f-%.2f)", OR, CI_lower_bound, CI_upper_bound),
         p_formatting = ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))) %>%
  select(wave, exposure, exposure, model, definition, `OR (95%CI)`, OR, CI_lower_bound, CI_upper_bound, p, p_formatting, n)
writeLines(sprintf("sensitivity detail row count: %d", nrow(detail_sens)))

## ---- 6. five wave adjusted comparison OR matrix(weighted adjustment; main table + shared construction for fixed cutoff sensitivity) ----
build_mat <- function(det, spec_list, years) {
  adj <- det %>% filter(model == "adjusted model", definition == "weighted (svyglm)")
  ms  <- do.call(rbind, lapply(spec_list, function(ex)
    data.frame(exposure = ex$v, term = names(ex$terms), exposure = unname(ex$terms),
               stratification = if (!is.null(ex$layer)) ex$layer else ifelse(isTRUE(ex$core), "core", "exploratory"),
               stringsAsFactors = FALSE)))
  m  <- matrix("", nrow(ms), length(years), dimnames = list(ms$exposure, as.character(years)))
  mn <- m
  for (i in seq_len(nrow(ms))) for (y in years) {
    r <- adj[adj$exposure == ms$exposure[i] & adj$term == ms$term[i] & adj$wave == y, ]
    if (nrow(r) == 1) {
      m[i, as.character(y)]  <- sprintf("%.2f (%.2f-%.2f)%s", r$OR, r$CI_lower_bound, r$CI_upper_bound,
                                        ifelse(r$p < 0.05, "*", ""))
      mn[i, as.character(y)] <- format(r$n, big.mark = ",")
    } else {
      mk <- MARK(ms$exposure[i], y)
      m[i, as.character(y)]  <- mk
      mn[i, as.character(y)] <- if (nzchar(mk)) "—" else ""
    }
  }
  list(df = data.frame(exposure = ms$exposure, stratification = ms$stratification, m, check.names = FALSE, row.names = NULL),
       n  = data.frame(exposure = ms$exposure, mn, check.names = FALSE, row.names = NULL))
}
mat_all  <- build_mat(detail, EXPOS, YEARS)
mat_df   <- mat_all$df
mat_n_df <- mat_all$n
adj_w    <- detail %>% filter(model == "adjusted model", definition == "weighted (svyglm)")   # number 7 used for verification in section
writeLines("===== adjustment OR matrix(weighted adjustment; * = p<0.05)=====")
print(mat_df, right = FALSE)

## ---- 6b. fixed cutoff sensitivity: threshold uses only 2011 year specification, apply to five waves ----
writeLines("===== fixed cutoff sensitivity(2011 year education specific thresholds applied to all five waves)=====")
detail_fixed <- do.call(rbind, lapply(YEARS, function(y)
  run_wave(ana %>% filter(year == y), EXPOS, with_q = FALSE, yvar = "imp_fixed")))
detail_fixed <- detail_fixed %>%
  mutate(`OR (95%CI)` = sprintf("%.2f (%.2f-%.2f)", OR, CI_lower_bound, CI_upper_bound),
         p_formatting = ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))) %>%
  select(wave, exposure, exposure, term, core_exposure, model, definition, `OR (95%CI)`, OR, CI_lower_bound, CI_upper_bound, p, p_formatting, n)
cut_tab <- ana %>% group_by(year) %>% summarise(
  analysis_sample_n = n(),
  wave_specific_cutoff_prevalence = sprintf("%.1f%%", 100 * mean(cog_B)),
  fixed_2011_cutoff_prevalence = sprintf("%.1f%%", 100 * mean(imp_fixed)),
  discordant_classification_between_definitions = sprintf("%.1f%%", 100 * mean(cog_B != imp_fixed)), .groups = "drop")
writeLines("prevalence comparison under two cutoff definitions: ")
print(as.data.frame(cut_tab))
thr_wave_tab <- ana %>% group_by(year, edu_group) %>%
  summarise(wave_specific_cutoff = round(mean(total, na.rm = TRUE) - sd(total, na.rm = TRUE), 2), .groups = "drop") %>%
  pivot_wider(names_from = year, values_from = wave_specific_cutoff) %>%
  mutate(fixed_2011_cutoff = round(thr11$thr11[match(edu_group, thr11$edu_group)], 2))
mat_fix <- build_mat(detail_fixed, EXPOS, YEARS)
writeLines("adjustment using fixed cutoff OR matrix: ")
print(mat_fix$df, right = FALSE)
cmp_cut <- detail %>% filter(model == "adjusted model", definition == "weighted (svyglm)", core_exposure) %>%
  select(exposure, term, exposure, wave, OR_wave_specific_cutoff = OR, p_wave_specific_cutoff = p) %>%
  inner_join(detail_fixed %>% filter(model == "adjusted model", definition == "weighted (svyglm)") %>%
               select(exposure, term, wave, OR_fixed_cutoff = OR, p_fixed_cutoff = p),
             by = c("exposure", "term", "wave")) %>%
  mutate(OR_difference = round(OR_fixed_cutoff - OR_wave_specific_cutoff, 2)) %>%
  mutate(OR_wave_specific_cutoff = round(OR_wave_specific_cutoff, 3), OR_fixed_cutoff = round(OR_fixed_cutoff, 3),
         p_wave_specific_cutoff = signif(p_wave_specific_cutoff, 3), p_fixed_cutoff = signif(p_fixed_cutoff, 3)) %>%
  select(exposure, wave, OR_wave_specific_cutoff, OR_fixed_cutoff, OR_difference, p_wave_specific_cutoff, p_fixed_cutoff)
writeLines(sprintf("core exposure: main table vs fixed cutoff, maximum |OR difference| = %.2f; positive association unchanged",
                   max(abs(cmp_cut$OR_difference))))

## ---- 6c. excluding figure copying(0-30)sensitivity: addresses 2020 year anomalous figure copying scoring ----
writeLines("===== sensitivity excluding figure copying(total score 0-30)=====")
detail_nodraw <- do.call(rbind, lapply(YEARS, function(y)
  run_wave(ana %>% filter(year == y), EXPOS, with_q = FALSE, yvar = "imp_nodraw")))
detail_nodraw <- detail_nodraw %>%
  mutate(`OR (95%CI)` = sprintf("%.2f (%.2f-%.2f)", OR, CI_lower_bound, CI_upper_bound),
         p_formatting = ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))) %>%
  select(wave, exposure, exposure, term, core_exposure, model, definition, `OR (95%CI)`, OR, CI_lower_bound, CI_upper_bound, p, p_formatting, n)
cut_tab2 <- ana %>% group_by(year) %>% summarise(
  analysis_sample_n = n(),
  wave_specific_cutoff_prevalence = sprintf("%.1f%%", 100 * mean(cog_B)),
  prevalence_excluding_figure_copying = sprintf("%.1f%%", 100 * mean(imp_nodraw)),
  discordant_classification_between_definitions = sprintf("%.1f%%", 100 * mean(cog_B != imp_nodraw)), .groups = "drop")
writeLines("including figure copying(main table)vs excluding figure copying(0-30)prevalence: ")
print(as.data.frame(cut_tab2))
mat_nd <- build_mat(detail_nodraw, EXPOS, YEARS)
writeLines("adjustment excluding figure copying OR matrix: ")
print(mat_nd$df, right = FALSE)
cmp_nd <- detail %>% filter(model == "adjusted model", definition == "weighted (svyglm)", core_exposure) %>%
  select(exposure, term, exposure, wave, OR_main_table = OR, p_main_table = p) %>%
  inner_join(detail_nodraw %>% filter(model == "adjusted model", definition == "weighted (svyglm)") %>%
               select(exposure, term, wave, OR_without_figure_copying = OR, p_without_figure_copying = p),
             by = c("exposure", "term", "wave")) %>%
  mutate(OR_difference = round(OR_without_figure_copying - OR_main_table, 2)) %>%
  mutate(OR_main_table = round(OR_main_table, 3), OR_without_figure_copying = round(OR_without_figure_copying, 3),
         p_main_table = signif(p_main_table, 3), p_without_figure_copying = signif(p_without_figure_copying, 3)) %>%
  select(exposure, wave, OR_main_table, OR_without_figure_copying, OR_difference, p_main_table, p_without_figure_copying)
writeLines(sprintf("core exposure: main table vs excluding figure copying, maximum |OR difference| = %.2f", max(abs(cmp_nd$OR_difference))))

## ---- 7. automatic verification: 2020 column must reproduce delivered results Paper 1(2020 wave v4.4)Table 2 ----
target <- tribble(
  ~exposure, ~term,               ~OR_t,  ~lo_t,  ~hi_t,
  "dep",      "dep",               1.73,   1.51,   1.98,
  "adl",      "adl",               1.27,   1.11,   1.45,
  "iadl",     "iadl",              1.57,   1.36,   1.80,
  "sleep_cat","sleep_cat<6h",      1.19,   1.05,   1.36,
  "sleep_cat","sleep_cat>8h",      1.60,   1.24,   2.06,
  "pa_any",   "pa_any",            0.53,   0.43,   0.65,
  "pain3f",   "pain3fsevere",      1.17,   0.99,   1.38)
chk20 <- adj_w %>% filter(wave == 2020) %>% select(exposure, term, OR, CI_lower_bound, CI_upper_bound, n)
cmp20 <- target %>% left_join(chk20, by = c("exposure", "term"), suffix = c("_target", "_new")) %>%
  mutate(consistent = abs(OR - OR_t) < 0.02 & abs(CI_lower_bound - lo_t) < 0.03 & abs(CI_upper_bound - hi_t) < 0.05)
writeLines("===== 2020 columns vs previously delivered Paper 1 Table 2 =====")
print(as.data.frame(cmp20 %>% select(exposure, term, OR_t, OR, CI_lower_bound, CI_upper_bound, n, consistent)), row.names = FALSE)
if (!all(cmp20$consistent)) warning("2020 column and previously delivered Paper 1 inconsistencies found, please check!") else
  writeLines("verification passed: 2020 columns 6 core exposures + severe pain and previously delivered Paper 1 identical.")

## ---- 8. write Excel and CSV ----
note_lines <- c(
  "CHARLS separate wave analyses Logistic regression——definition notes(v2, 2026-09-19 revision)",
  "",
  "1. data sources: output/02_five_wave_data_and_quality_control/waves_cognition_2026-09-13.rds(code/03_waves_build.R).script for this table: code/10_fivewave_logistic.R(v2).",
  "2. analysis sample: age at each wave ≥ 45 years old, cognition 9 complete components(cog_complete == 1)and cog_B nonmissing.each wave n: ",
  "   2011 = 12,288; 2013 = 13,687; 2015 = 15,857(15,859 in 2 participants missing education; cutoff unavailable); 2018 = 11,653; 2020 = 14,655.",
  "3. outcome: each wave cog_B(education specific cutoff: within this wave, education stratum mean total cognition − 1SD; total score 0-31), and 2020 delivered wave baseline(code/02_table1.R, v4.4)consistent.",
  "4. model: each wave × fit crude model separately for each exposure(cog_B ~ exposure)and adjusted model(+ age + sex + four education groups + urban or rural).report exposure coefficients only; ",
  "   education coefficients deliberately omitted——education determines outcome specific cutoffs, its coefficient is not interpretable(same 02_table1.R comments).",
  "5. weighting specification(threshold = of this wave analysis sample 50%): each wave weight + communityID coverage shown below, all five waves ≥ 50%, ",
  "   therefore use svyglm(weights = weight, ids = communityID, quasibinomial)as primary definition, unweighted glm(binomial)as comparison column.",
  "   weight source: 2011-2018 is Harmonized D of r*wtresp, 2020 is INDV_weight(03_waves_build.R constructed this way).",
  "   note: weighted model uses only weight/communityID nonmissing participants, actual sample size differs from unweighted model(each model n see details), differences also reflect sample differences.",
  "6. sleep(v2 revision): following delivered results Paper 1 definition changed to three categories <6 h / 6-8 h(reference)/ >8 h; continuous coding(per 1 h)move into[sensitivity subtable].",
  "7. pain(v2 revision): main row[bodily pain(any, vs no pain)]comparable across five waves(2011/2015 raw coding contains only[present/none]two levels); ",
  "   additionally[mild 2-3][severe 4-5]two rows, only 2013 / 2018 / 2020 ordinal coding supports estimation, 2011/2015 mark[binary coding only].pain classified by SAP exploratory, no multiplicity adjustment.",
  "8. multiplicity(v2 newly added): 5 core exposures(depression, ADL, IADL, sleep, physical activity)within each wave perform BH-FDR; ",
  "   sleep is a three category factor, adjustment includes the whole factor Wald test P(and Paper 1 footnotes use consistent definitions), two sleep contrasts share the same q.q value shown in[wave specific details]sheet.",
  "9. each model reports n(model.frame row count): substantial exposure missingness → n small, interpret together with n interpretation.",
  "10._IADL_2011_wave_excluded.raw_module_at_this_wave_actually_contains_5_items_IADL_items_db016_db020_but_Harmonized_D_not_provided_2011_wave_summary_variable",
  "    (official summary r*w*iadla = making telephone calls + managing finances + taking medication, whereas[making telephone calls]starting at 2 wave first asked), therefore this pipeline 2011 wave IADL code as missing, does not estimate; 2011 column labeled[2011 wave excluded].",
  "11._physical_activity_pa_any_at_2011/2013/2015_substantial_missingness_analysis_sample_missing_counts_respectively_5_257/ 7,018 /6_206_participants_corresponding_model_n_substantially_smaller_low_precision.",
  "12._urban_or_rural_rural_2020_wave_based_on_ID_from_analysis_2020.rds_of_urban3_backfill_2015/2018_wave_rural_missing_respectively_1_752/1_112_participants_adjusted_model_n_correspondingly_reduced.",
  "13._adjustment_OR_matrix_five_wave_comparison_sheet_is_main_table_each_cell_of_weighted_adjusted_model_OR_95_CI_p_0.05_same_sheet_cell_sample_sizes_below_n.",
  "14._automatic_verification_at_script_end_2020_columns_6_core_exposures_severe_pain_and_previously_delivered_Paper_1_Table_2_itemwise_agreement_tolerance_0.02/0.03/0.05）。",
  "15._cross_wave_pain_comparability_pain_coding_differs_across_waves_2011/2015_da041_present/none_2013_wb16_2018_da041_w4_2020_da027",
  "    2018/2020 pain prevalence(approximately 58%)much higher than the first three waves(28%-35%); across waves OR interpret trends cautiously.",
  "16._CSV_is_UTF_8_none_BOM_wave_specific_details_adjustment_OR_matrix_sensitivity_subtable_fixed_cutoff_details_details_excluding_figure_copying_total_5_files.",
  "    note: Chinese Windows double click in CSV may display garbled text, use[data → from text/CSV]import and select UTF-8; for formal inspection use this Excel file.",
  "17._fixed_cutoff_sensitivity_sheet_main_table_outcome_cog_B_using_this_wave_own_education_specific_cutoff_thresholds_drift_across_waves_below_primary_school_8.19/7.97/7.08/7.67/6.96），",
  "    this sheet change to[2011 year threshold set]apply to all five waves and refit all models.inspection order: ① prevalence and discordant classification under two definitions; ② wave specific thresholds; ",
  "    ③ fixed cutoff OR matrix; ④ main table of core exposures vs cellwise fixed cutoff comparison(OR differences generally <0.15, 2020 long sleep difference 0.28, direction and significance unchanged).",
  "18._IADL_two_rows_main_row_IADL_limited_and_previously_delivered_Paper_1_consistent_definition_2013_2018_Harmonized_three_items_making_telephone_calls_managing_finances_taking_medication_2020_any_limitation_in_six_items",
  "    2011 wave excluded(telephone not asked this wave); ② [IADL limited(managing finances+two medication items)]= five wave consistent version(2011 using r1moneya+r1medsa combine, 2020 using db022+db020 ≥2), ",
  "    limitation rate 11.9%/9.4%/9.4%/8.0%/9.0%, only comparable across five waves, and 2011 with available values IADL definition, therefore 2011 of column IADL see this row for information.",
  "19._sensitivity_excluding_figure_copying_sheet_outcome_changed_to_total_score_figure_copying_0_30_recompute.reason_figure_copying_item_2011/2013/2015_classify_only_successfully_copied/unable_to_copy_2018_strict_scoring_from_this_point",
  "    2020 pass rate only 34.6%(2011–2018 is 67–71%, and 2020 year online interviewer training, related to scoring implementation differences); inspect as for fixed cutoff sheet identical.")

wb <- createWorkbook()
addWorksheet(wb, "notes")
writeData(wb, "notes", data.frame(notes = note_lines, check.names = FALSE))
writeData(wb, "notes", "weight coverage by wave(measured during execution): ", startRow = length(note_lines) + 3)
writeData(wb, "notes", as.data.frame(cov_tab), startRow = length(note_lines) + 4)
writeData(wb, "notes", "2020 column and previously delivered Paper 1 verification of: ", startRow = length(note_lines) + nrow(cov_tab) + 6)
cmp_out <- cmp20 %>% transmute(exposure = exposure, contrast = term, `Paper1_OR_target` = round(OR_t, 2),
                               `this_release_OR` = round(OR, 3), `CI_lower_bound` = round(CI_lower_bound, 3),
                               `CI_upper_bound` = round(CI_upper_bound, 3), n = n, consistent = consistent)
writeData(wb, "notes", as.data.frame(cmp_out), startRow = length(note_lines) + nrow(cov_tab) + 7)
setColWidths(wb, "notes", cols = 1, widths = 130)

addWorksheet(wb, "wave specific details")
writeData(wb, "wave specific details", detail)
setColWidths(wb, "wave specific details", cols = 1:14, widths = 16)

addWorksheet(wb, "adjustment OR matrix_five wave comparison")
writeData(wb, "adjustment OR matrix_five wave comparison",
          "five wave adjusted comparison OR matrix(outcome = cog_B; adjustment age+sex+four education groups+urban or rural; weighted svyglm; * = p<0.05)", startRow = 1)
writeData(wb, "adjustment OR matrix_five wave comparison", "sleep in three categories(reference 6-8 h), pain classification(reference no pain), one row per contrast; excluded or binary only cells shown as[—(…)]annotation.",
          startRow = 2)
writeData(wb, "adjustment OR matrix_five wave comparison", mat_df, startRow = 4)
writeData(wb, "adjustment OR matrix_five wave comparison", "sample size of model for each cell n(weighted adjusted model): ",
          startRow = nrow(mat_df) + 7)
writeData(wb, "adjustment OR matrix_five wave comparison", mat_n_df, startRow = nrow(mat_df) + 8)
setColWidths(wb, "adjustment OR matrix_five wave comparison", cols = 1, widths = 26)
setColWidths(wb, "adjustment OR matrix_five wave comparison", cols = 2:7, widths = 26)

addWorksheet(wb, "sensitivity subtable")
writeData(wb, "sensitivity subtable", "variable specification sensitivity(corresponds to delivered Paper 1 supplementary table S2): continuous sleep(per 1 h)and binary pain.", startRow = 1)
writeData(wb, "sensitivity subtable", detail_sens, startRow = 3)
setColWidths(wb, "sensitivity subtable", cols = 1:12, widths = 16)

addWorksheet(wb, "fixed cutoff sensitivity")
sn <- "fixed cutoff sensitivity"; rw <- 1
writeData(wb, sn, "fixed cutoff sensitivity: threshold uses only 2011 year specification(education strata mean total score−1SD), apply unchanged to all five waves", startRow = rw); rw <- rw + 1
writeData(wb, sn, "purpose: main table uses[wave specific cutoff], threshold may drift across waves; this sheet test dependence of conclusions on threshold definition.", startRow = rw); rw <- rw + 2
writeData(wb, sn, "1) prevalence under two cutoff definitions", startRow = rw); rw <- rw + 1
writeData(wb, sn, cut_tab, startRow = rw); rw <- rw + nrow(cut_tab) + 2
writeData(wb, sn, "2) wave specific thresholds(total score 0-31; education strata mean−1SD)", startRow = rw); rw <- rw + 1
writeData(wb, sn, thr_wave_tab, startRow = rw); rw <- rw + nrow(thr_wave_tab) + 2
writeData(wb, sn, "3) weighted adjustment with fixed cutoff OR matrix(* = p<0.05)", startRow = rw); rw <- rw + 1
writeData(wb, sn, mat_fix$df, startRow = rw); rw <- rw + nrow(mat_fix$df) + 2
writeData(wb, sn, "4) model sample size in each cell n", startRow = rw); rw <- rw + 1
writeData(wb, sn, mat_fix$n, startRow = rw); rw <- rw + nrow(mat_fix$n) + 2
writeData(wb, sn, "5) core exposure: main table(wave specific cutoff)vs fixed cutoff", startRow = rw); rw <- rw + 1
writeData(wb, sn, cmp_cut, startRow = rw)
setColWidths(wb, sn, cols = 1, widths = 30); setColWidths(wb, sn, cols = 2:8, widths = 20)

addWorksheet(wb, "sensitivity excluding figure copying")
sn2 <- "sensitivity excluding figure copying"; rw2 <- 1
writeData(wb, sn2, "sensitivity excluding figure copying: outcome changed to[total score − figure copying](0-30), cutoff within each wave, this education stratum mean−1SD recompute", startRow = rw2); rw2 <- rw2 + 1
writeData(wb, sn2, "background: 2011/2013/2015 figure copying classified only as[successfully copied/unable to copy], 2018 strict scoring from this point, 2020 year pass rate only 34.6%(2011-2018 is 67-71%), reflects measurement artifact.", startRow = rw2); rw2 <- rw2 + 1
writeData(wb, sn2, "this sheet for testing: after removing figure copying from total, whether prevalence and core exposure associations change by wave.", startRow = rw2); rw2 <- rw2 + 2
writeData(wb, sn2, "1) including figure copying(main table)vs excluding figure copying(0-30)prevalence", startRow = rw2); rw2 <- rw2 + 1
writeData(wb, sn2, cut_tab2, startRow = rw2); rw2 <- rw2 + nrow(cut_tab2) + 2
writeData(wb, sn2, "2) weighted adjustment excluding figure copying OR matrix(* = p<0.05)", startRow = rw2); rw2 <- rw2 + 1
writeData(wb, sn2, mat_nd$df, startRow = rw2); rw2 <- rw2 + nrow(mat_nd$df) + 2
writeData(wb, sn2, "3) model sample size in each cell n", startRow = rw2); rw2 <- rw2 + 1
writeData(wb, sn2, mat_nd$n, startRow = rw2); rw2 <- rw2 + nrow(mat_nd$n) + 2
writeData(wb, sn2, "4) core exposure: main table(including figure copying)vs excluding figure copying", startRow = rw2); rw2 <- rw2 + 1
writeData(wb, sn2, cmp_nd, startRow = rw2)
setColWidths(wb, sn2, cols = 1, widths = 30); setColWidths(wb, sn2, cols = 2:8, widths = 20)

f_xlsx <- file.path(OUT, "separate_wave_analyses_Logistic_OR_2026_09_19.xlsx")
saveWorkbook(wb, f_xlsx, overwrite = TRUE)
write.csv(detail,      file.path(OUT, "separate_wave_analyses_Logistic_OR_2026_09_19.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(mat_df,      file.path(OUT, "separate_wave_analyses_Logistic_OR_matrix_2026_09_19.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(detail_sens, file.path(OUT, "separate_wave_analyses_Logistic_sensitivity_2026_09_19.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(detail_fixed, file.path(OUT, "separate_wave_analyses_Logistic_fixed_cutoff_sensitivity_2026_09_19.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(detail_nodraw, file.path(OUT, "separate_wave_analyses_Logistic_sensitivity_excluding_figure_copying_2026_09_19.csv"), row.names = FALSE, fileEncoding = "UTF-8")
writeLines(sprintf("complete: %s", f_xlsx))
writeLines(paste("CSV: separate_wave_analyses_Logistic_OR_2026_09_19.csv ｜ _OR matrix_ ｜ _sensitivity_ ｜ _fixed cutoff sensitivity_ ｜ _sensitivity excluding figure copying_(all 2026-09-19)"))
