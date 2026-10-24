# ============================================================================
# 29_A15_sensitivity_analyses.R —— A15: weights, IPW, 2018 measurement differences, cutoff, missingness analysis
# basis: SAP v3 §10 sensitivity inventory, §14 T5/T6, §15–§17(A07–A14 frozen model)
# principle: final model unchanged, score table, primary cutoff and primary external validation definition; all results are sensitivity analyses only.
# input: 
#   output/02_five_wave_data_and_quality_control/waves_cognition_2026-09-13.rds
#   output/10_screening_index/screening_index_memory_v3_2026-09-26.rds
#   output/11_external_validation/A12_results_2026-09-26.rds
#   external_data/{HRS,ELSA,SHARE}/analysis/*_sensitivity_risk_sets_2026-09-26.rds
# output: output/13_A15_sensitivity_analysis/
# generate: 2026-09-26
# ============================================================================
invisible(try(Sys.setlocale("LC_ALL", "zh_CN.UTF-8"), silent = TRUE))
.libPaths(c("D:/R/library", .libPaths()))
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(survey); library(pROC)
  library(mice); library(knitr)
})
options(survey.lonely.psu = "adjust", stringsAsFactors = FALSE)
set.seed(20260926)

P <- "D:/BaiduSyncdisk/research/Xu_Zhonghao"
OUT <- file.path(P, "output", "13_A15_sensitivity_analysis")
FIG <- file.path(P, "output", "09_figures_and_tables")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# 0. general helper functions
# ---------------------------------------------------------------------------
fmt <- function(x, d = 3) ifelse(is.na(x), "NA", sprintf(paste0("%.", d, "f"), x))
fmt_ci <- function(est, lo, hi, d = 3) ifelse(is.na(est), "NA",
  sprintf(paste0("%.", d, "f (%.", d, "f–%.", d, "f)"), est, lo, hi))
fmt_p <- function(p) ifelse(is.na(p), "NA", ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))

auc_ci_unweighted <- function(y, p) {
  ok <- !is.na(y) & !is.na(p)
  y <- y[ok]; p <- p[ok]
  if (length(y) < 10 || length(unique(y)) < 2) return(c(auc = NA, low = NA, high = NA))
  r <- pROC::roc(y, p, direction = "<", quiet = TRUE)
  ci <- as.numeric(pROC::ci.auc(r))
  c(auc = as.numeric(pROC::auc(r)), low = ci[1], high = ci[3])
}

auc_ci_weighted <- function(y, p, w, cluster, B = 500, seed = 20260926) {
  ok <- !is.na(y) & !is.na(p) & !is.na(w) & !is.na(cluster)
  y <- y[ok]; p <- p[ok]; w <- w[ok]; cluster <- as.character(cluster[ok])
  if (length(y) < 10 || length(unique(y)) < 2) return(c(auc = NA, low = NA, high = NA))
  point <- tryCatch(as.numeric(pROC::auc(pROC::roc(y, p, weights = w, direction = "<", quiet = TRUE))),
                    error = function(e) NA_real_)
  ids <- unique(cluster)
  idx_list <- split(seq_along(cluster), cluster)
  set.seed(seed)
  boot <- replicate(B, {
    samp <- sample(ids, length(ids), replace = TRUE)
    idx <- unlist(idx_list[samp], use.names = FALSE)
    if (length(unique(y[idx])) < 2) return(NA_real_)
    tryCatch(as.numeric(pROC::auc(pROC::roc(y[idx], p[idx], weights = w[idx],
                                            direction = "<", quiet = TRUE))),
             error = function(e) NA_real_)
  })
  ci <- suppressWarnings(quantile(boot, c(0.025, 0.975), na.rm = TRUE))
  c(auc = point, low = unname(ci[1]), high = unname(ci[2]))
}

calibration_glm <- function(y, lp) {
  ok <- !is.na(y) & !is.na(lp)
  y <- y[ok]; lp <- lp[ok]
  fs <- glm(y ~ lp, family = binomial())
  cs <- coef(fs)[["lp"]]; cis <- confint.default(fs)["lp", ]
  fi <- glm(y ~ offset(lp), family = binomial())
  ci <- coef(fi)[["(Intercept)"]]; cii <- confint.default(fi)["(Intercept)", ]
  c(slope = unname(cs), slope_low = unname(cis[1]), slope_high = unname(cis[2]),
    intercept = unname(ci), intercept_low = unname(cii[1]), intercept_high = unname(cii[2]))
}

calibration_svy <- function(y, lp, w, cluster) {
  ok <- !is.na(y) & !is.na(lp) & !is.na(w) & !is.na(cluster)
  d <- data.frame(y = y[ok], lp = lp[ok], w = w[ok], cluster = as.character(cluster[ok]))
  des <- svydesign(ids = ~cluster, weights = ~w, data = d)
  fs <- svyglm(y ~ lp, design = des, family = quasibinomial())
  cs <- coef(fs)[["lp"]]; cis <- confint(fs)["lp", ]
  fi <- svyglm(y ~ offset(lp), design = des, family = quasibinomial())
  ci <- coef(fi)[["(Intercept)"]]; cii <- confint(fi)["(Intercept)", ]
  c(slope = unname(cs), slope_low = unname(cis[1]), slope_high = unname(cis[2]),
    intercept = unname(ci), intercept_low = unname(cii[1]), intercept_high = unname(cii[2]))
}

model_metrics <- function(section, scenario, model_type, y, p, score = NULL,
                          w = NULL, cluster = NULL, n_note = NA_integer_, B = 500) {
  ok <- !is.na(y) & !is.na(p)
  ys <- y[ok]; ps <- p[ok]
  if (!is.null(score)) { sc <- score[ok] } else { sc <- NULL }
  if (!is.null(w)) {
    ws <- w[ok]; cl <- cluster[ok]
    a <- auc_ci_weighted(ys, ps, ws, cl, B = B)
    cal <- calibration_svy(ys, qlogis(pmin(pmax(ps, 1e-6), 1 - 1e-6)), ws, cl)
    brier <- weighted.mean((ps - ys)^2, ws, na.rm = TRUE)
  } else {
    a <- auc_ci_unweighted(ys, ps)
    cal <- calibration_glm(ys, qlogis(pmin(pmax(ps, 1e-6), 1 - 1e-6)))
    brier <- mean((ps - ys)^2, na.rm = TRUE)
  }
  a_score <- if (!is.null(sc)) auc_ci_unweighted(ys, sc) else c(auc = NA, low = NA, high = NA)
  data.frame(
    section = section, scenario = scenario, model_type = model_type,
    n_model = length(ys), events_model = sum(ys, na.rm = TRUE),
    event_rate = if (!is.null(w)) weighted.mean(ys, ws, na.rm = TRUE) else mean(ys, na.rm = TRUE),
    event_rate_type = if (!is.null(w)) "weighted" else "unweighted",
    auc_model = a[["auc"]], auc_low = a[["low"]], auc_high = a[["high"]],
    auc_score = a_score[["auc"]], auc_score_low = a_score[["low"]], auc_score_high = a_score[["high"]],
    cal_slope = cal[["slope"]], cal_slope_low = cal[["slope_low"]], cal_slope_high = cal[["slope_high"]],
    cal_intercept = cal[["intercept"]], cal_intercept_low = cal[["intercept_low"]],
    cal_intercept_high = cal[["intercept_high"]], brier = brier,
    n_note = n_note, stringsAsFactors = FALSE
  )
}

coef_table <- function(fit, section, scenario, model_type, n, events) {
  co <- summary(fit)$coefficients
  est <- co[, 1]; se <- co[, 2]; pval <- co[, ncol(co)]
  ci <- tryCatch(if (inherits(fit, "svyglm")) confint(fit) else confint.default(fit),
                 error = function(e) cbind(est - 1.96 * se, est + 1.96 * se))
  data.frame(
    section = section, scenario = scenario, model_type = model_type,
    term = rownames(co), beta = est, se = se, p = pval,
    OR = exp(est), OR_low = exp(ci[, 1]), OR_high = exp(ci[, 2]),
    n = n, events = events, row.names = NULL, check.names = FALSE
  )
}

normalize_weight <- function(w) w / mean(w, na.rm = TRUE)

safe_write_csv <- function(x, file, ...) {
  for (i in seq_len(10)) {
    ok <- tryCatch({ utils::write.csv(x, file, ...); TRUE }, error = function(e) FALSE)
    if (ok) return(invisible(TRUE))
    Sys.sleep(0.5)
  }
  stop("cannot write file(retry 10 times): ", file)
}

safe_write_lines <- function(text, con, ...) {
  for (i in seq_len(10)) {
    ok <- tryCatch({ base::writeLines(text, con, ...); TRUE }, error = function(e) FALSE)
    if (ok) return(invisible(TRUE))
    Sys.sleep(0.5)
  }
  stop("cannot write file(retry 10 times): ", con)
}

# fixed formula of frozen model(A14 §17)
fixed_score <- function(d) {
  2 * pmax(floor((d$age - 45) / 5), 0) +
    ifelse(d$edu_group == "Primary school", -2,
      ifelse(d$edu_group == "Middle school", -7,
        ifelse(d$edu_group == "High school or above", -10, 0))) +
    3 * d$dep
}

fixed_lp_manual <- function(d, b) {
  e <- as.character(d$edu_group)
  b[["(Intercept)"]] + b[["age"]] * d$age +
    b[["edu_groupPrimary school"]] * (e == "Primary school") +
    b[["edu_groupMiddle school"]] * (e == "Middle school") +
    b[["edu_groupHigh school or above"]] * (e == "High school or above") +
    b[["dep"]] * d$dep + b[["adl"]] * d$adl
}

# ---------------------------------------------------------------------------
# 1. build CHARLS internal sensitivity data
# ---------------------------------------------------------------------------
w <- readRDS(file.path(P, "output/02_five_wave_data_and_quality_control/waves_cognition_2026-09-13.rds"))
w$KEY <- ifelse(is.na(w$IDn), w$ID, w$IDn)
w$ms10 <- as.numeric(w$orient) + as.numeric(w$ser7)
stopifnot(!anyDuplicated(w[, c("KEY", "year")]))

base <- w %>% filter(year == 2011, !is.na(age), age >= 45) %>%
  transmute(
    KEY, age, female, edu_group = factor(edu_group,
      levels = c("Less than primary", "Primary school", "Middle school", "High school or above")),
    rural, communityID, weight, cesd10, dep = as.integer(cesd10 >= 10),
    adl, iadl2, sleep_h, pain_any, hibp, diab, stroke,
    mem2011 = imrc + dlrc, ms10_2011 = ms10
  )
m15 <- w %>% filter(year == 2015) %>%
  transmute(KEY, mem2015 = imrc + dlrc, imrc2015 = imrc, dlrc2015 = dlrc)
m18 <- w %>% filter(year == 2018) %>%
  transmute(KEY, age2018 = age, mem2018 = imrc + dlrc, imrc2018 = imrc, dlrc2018 = dlrc)
n_ms10 <- w %>% filter(!is.na(ms10)) %>% count(KEY, name = "n_ms10")
base <- base %>% left_join(m15, by = "KEY") %>% left_join(m18, by = "KEY") %>%
  left_join(n_ms10, by = "KEY") %>% mutate(fixed5 = as.integer(!is.na(n_ms10) & n_ms10 == 5))

# primary cutoff: 2011 fixed cohort mean−1SD = 4.77 → floor = 4(A14 frozen)
pool <- base %>% filter(!is.na(mem2011), mem2011 > 4)
main <- pool %>% filter(fixed5 == 1, !is.na(mem2015), complete.cases(.[, c("age", "edu_group", "dep", "adl")])) %>%
  transmute(KEY, age, edu_group, dep, adl, mem2011, mem2015, event = as.integer(mem2015 <= 4),
            weight, communityID, fixed5, female, rural, cesd10, sleep_h, pain_any, hibp, diab, stroke,
            ms10_2011, age2018, mem2018, imrc2018, dlrc2018)
broad <- pool %>% filter(!is.na(mem2015), complete.cases(.[, c("age", "edu_group", "dep", "adl")])) %>%
  transmute(KEY, age, edu_group, dep, adl, mem2011, mem2015, event = as.integer(mem2015 <= 4),
            weight, communityID, fixed5, female, rural, cesd10, sleep_h, pain_any, hibp, diab, stroke,
            ms10_2011, age2018, mem2018, imrc2018, dlrc2018)
main$score <- fixed_score(main); broad$score <- fixed_score(broad)

stopifnot(nrow(main) == 4441, sum(main$event) == 658)
stopifnot(nrow(broad) == 8296, sum(broad$event) == 1889)

a07 <- readRDS(file.path(P, "output/10_screening_index/screening_index_memory_v3_2026-09-26.rds"))
stopifnot(all(main$KEY %in% a07$data$KEY), all(a07$data$KEY %in% main$KEY))
stopifnot(sum(main$score != a07$data$score[match(main$KEY, a07$data$KEY)]) == 0)

f <- event ~ age + edu_group + dep + adl
model_main <- a07$model_u
b_main <- coef(model_main)
lp_main_main <- fixed_lp_manual(main, b_main)
p_main_main <- plogis(lp_main_main)
lp_main_broad <- fixed_lp_manual(broad, b_main)
p_main_broad <- plogis(lp_main_broad)

# flow table
flow <- data.frame(
  step = c("2011 age in year≥45 years old", "2011 year memory score measurable", "2011 year memory score>4(unimpaired at baseline)",
           "2015 year memory score measurable", "complete four variable data", "five waves ms10 fixed cohort"),
  n = c(nrow(base), sum(!is.na(base$mem2011)), nrow(pool),
        sum(!is.na(pool$mem2015)), nrow(broad), nrow(main)),
  note = c("CHARLS 2011 year", "immediate recall(imrc)and delayed recall(dlrc)both nonmissing", "prespecified primary cutoff ≤4",
           "4 year follow up", "age, education, depression, ADL", "A14 primary sample for frozen model(complete cases)")
)
safe_write_csv(flow, file.path(OUT, "A15_sample_flow_2026_09_26.csv"), row.names = FALSE, fileEncoding = "UTF-8")

# ---------------------------------------------------------------------------
# 2. weight sensitivity(2011 baseline weights + communityID is PSU)
# ---------------------------------------------------------------------------
des_main <- svydesign(ids = ~communityID, weights = ~weight, data = main)
fit_w <- svyglm(f, design = des_main, family = quasibinomial())
weights_diag <- data.frame(
  scenario = c("unweighted", "weighted"),
  n = nrow(main), communities = length(unique(main$communityID)),
  weight_median = c(NA, median(main$weight, na.rm = TRUE)),
  weight_cv = c(NA, sd(main$weight, na.rm = TRUE) / mean(main$weight, na.rm = TRUE)),
  effective_n = c(nrow(main), sum(main$weight, na.rm = TRUE)^2 / sum(main$weight^2, na.rm = TRUE)),
  event_rate_weighted = c(mean(main$event), weighted.mean(main$event, main$weight))
)
safe_write_csv(weights_diag, file.path(OUT, "A15_weights_design_diagnostics.csv"), row.names = FALSE, fileEncoding = "UTF-8")

coef_weight <- bind_rows(
  coef_table(model_main, "weights", "primary model(unweighted)", "unweighted", nrow(main), sum(main$event)),
  coef_table(fit_w, "weights", "2011 weights+PSU(refit)", "survey_refit", nrow(main), sum(main$event))
)
safe_write_csv(coef_weight, file.path(OUT, "A15_weights_model_comparison.csv"), row.names = FALSE, fileEncoding = "UTF-8")

met_weight <- bind_rows(
  model_metrics("weights", "primary model(unweighted, frozen)", "frozen", main$event, p_main_main,
                score = main$score, n_note = nrow(main), B = 500),
  model_metrics("weights", "2011 weights+PSU(frozen model)", "frozen", main$event, p_main_main,
                score = main$score, w = main$weight, cluster = main$communityID, n_note = nrow(main), B = 500),
  model_metrics("weights", "2011 weights+PSU(refit)", "survey_refit", main$event,
                as.numeric(predict(fit_w, type = "response")), score = main$score,
                w = main$weight, cluster = main$communityID, n_note = nrow(main), B = 500)
)
safe_write_csv(met_weight, file.path(OUT, "A15_weights_performance_comparison.csv"), row.names = FALSE, fileEncoding = "UTF-8")

# ---------------------------------------------------------------------------
# 3. IPW: inverse probability weighting for selection into final analysis
# ---------------------------------------------------------------------------
sel_vars <- c("age", "female", "edu_group", "rural", "mem2011", "cesd10", "adl")
pool_ipw <- pool %>% filter(complete.cases(.[, sel_vars]))
pool_ipw$selected <- as.integer(
  pool_ipw$fixed5 == 1 & !is.na(pool_ipw$mem2015) &
  complete.cases(pool_ipw[, c("age", "edu_group", "dep", "adl")])
)
ps_formula <- as.formula(paste("selected ~", paste(sel_vars, collapse = " + ")))
ps_fit <- glm(ps_formula, data = pool_ipw, family = binomial())
pool_ipw$p_sel <- pmin(pmax(predict(ps_fit, type = "response"), 0.02), 0.98)
pool_ipw$ipw <- mean(pool_ipw$selected) / pool_ipw$p_sel
main_ipw <- pool_ipw %>% filter(selected == 1) %>%
  transmute(KEY, age, edu_group, dep, adl, mem2011, mem2015, event = as.integer(mem2015 <= 4),
            weight, communityID, ipw, p_sel, female, rural, cesd10, sleep_h, pain_any, hibp, diab, stroke)
main_ipw$score <- fixed_score(main_ipw)
main_ipw$ipw_combined <- normalize_weight(main_ipw$weight * main_ipw$ipw)
stopifnot(nrow(main_ipw) == nrow(main), sum(main_ipw$event) == sum(main$event))

ipw_diag <- data.frame(
  selection_pool_n = nrow(pool_ipw), selected_n = sum(pool_ipw$selected),
  selection_fraction = mean(pool_ipw$selected),
  p_min = min(pool_ipw$p_sel), p_median = median(pool_ipw$p_sel), p_max = max(pool_ipw$p_sel),
  ipw_min = min(main_ipw$ipw), ipw_median = median(main_ipw$ipw), ipw_max = max(main_ipw$ipw),
  effective_n = sum(main_ipw$ipw)^2 / sum(main_ipw$ipw^2)
)
safe_write_csv(ipw_diag, file.path(OUT, "A15_IPW_diagnostics.csv"), row.names = FALSE, fileEncoding = "UTF-8")

ipw_selection_coef <- coef_table(ps_fit, "IPW", "selection into final sample(logistic)", "selection_model",
                                   nrow(pool_ipw), sum(pool_ipw$selected))
safe_write_csv(ipw_selection_coef, file.path(OUT, "A15_IPW_selection_model.csv"),
               row.names = FALSE, fileEncoding = "UTF-8")

smd_weighted <- function(x, w, ref) {
  ok <- !is.na(x) & !is.na(w); x <- x[ok]; w <- w[ok]
  ref <- ref[!is.na(ref)]
  if (length(x) < 2 || length(ref) < 2) return(NA_real_)
  mx <- sum(w * x) / sum(w); vx <- sum(w * (x - mx)^2) / sum(w)
  s <- sqrt((vx + var(ref)) / 2)
  if (s == 0) return(0)
  (mx - mean(ref)) / s
}
balance_vars <- c("age", "female", "edu_group", "rural", "mem2011", "cesd10", "adl",
                  "sleep_h", "pain_any", "hibp", "diab", "stroke")
balance_ipw <- bind_rows(lapply(balance_vars, function(v) {
  x <- if (is.factor(main_ipw[[v]])) as.numeric(main_ipw[[v]]) else main_ipw[[v]]
  r <- if (is.factor(pool_ipw[[v]])) as.numeric(pool_ipw[[v]]) else pool_ipw[[v]]
  data.frame(variable = v, selected_mean = mean(x, na.rm = TRUE), pool_mean = mean(r, na.rm = TRUE),
             smd_before = smd_weighted(x, rep(1, length(x)), r),
             smd_after = smd_weighted(x, main_ipw$ipw, r))
}))
safe_write_csv(balance_ipw, file.path(OUT, "A15_IPW_balance.csv"), row.names = FALSE, fileEncoding = "UTF-8")

des_ipw <- svydesign(ids = ~communityID, weights = ~ipw, data = main_ipw)
des_ipw_comb <- svydesign(ids = ~communityID, weights = ~ipw_combined, data = main_ipw)
fit_ipw <- svyglm(f, design = des_ipw, family = quasibinomial())
fit_ipw_comb <- svyglm(f, design = des_ipw_comb, family = quasibinomial())
coef_ipw <- bind_rows(
  coef_table(model_main, "IPW", "primary model(unweighted)", "unweighted", nrow(main_ipw), sum(main_ipw$event)),
  coef_table(fit_ipw, "IPW", "IPW(selection into final sample)", "ipw_refit", nrow(main_ipw), sum(main_ipw$event)),
  coef_table(fit_ipw_comb, "IPW", "weights×IPW(refit)", "combined_refit", nrow(main_ipw), sum(main_ipw$event))
)
safe_write_csv(coef_ipw, file.path(OUT, "A15_IPW_model_comparison.csv"), row.names = FALSE, fileEncoding = "UTF-8")

met_ipw <- bind_rows(
  model_metrics("IPW", "primary model(unweighted, frozen)", "frozen", main_ipw$event,
                plogis(fixed_lp_manual(main_ipw, b_main)), score = main_ipw$score,
                n_note = nrow(main_ipw), B = 500),
  model_metrics("IPW", "IPW(frozen model)", "frozen", main_ipw$event,
                plogis(fixed_lp_manual(main_ipw, b_main)), score = main_ipw$score,
                w = main_ipw$ipw, cluster = main_ipw$communityID, n_note = nrow(main_ipw), B = 500),
  model_metrics("IPW", "IPW(refit)", "ipw_refit", main_ipw$event,
                as.numeric(predict(fit_ipw, type = "response")), score = main_ipw$score,
                w = main_ipw$ipw, cluster = main_ipw$communityID, n_note = nrow(main_ipw), B = 500),
  model_metrics("IPW", "weights×IPW(refit)", "combined_refit", main_ipw$event,
                as.numeric(predict(fit_ipw_comb, type = "response")), score = main_ipw$score,
                w = main_ipw$ipw_combined, cluster = main_ipw$communityID, n_note = nrow(main_ipw), B = 500)
)
safe_write_csv(met_ipw, file.path(OUT, "A15_IPW_performance_comparison.csv"), row.names = FALSE, fileEncoding = "UTF-8")

# ---------------------------------------------------------------------------
# 4. broader sample/fixed cohort selection sensitivity
# ---------------------------------------------------------------------------
nonfixed <- broad %>% filter(fixed5 == 0)
smd_plain <- function(x1, x2) {
  x1 <- x1[!is.na(x1)]; x2 <- x2[!is.na(x2)]
  if (length(x1) < 2 || length(x2) < 2) return(NA_real_)
  s <- sqrt((var(x1) + var(x2)) / 2)
  if (s == 0) return(0)
  (mean(x1) - mean(x2)) / s
}
selection_vars <- c("age", "female", "edu_group", "rural", "mem2011", "cesd10", "adl", "event")
selection_balance <- bind_rows(lapply(selection_vars, function(v) {
  x <- main[[v]]; y <- nonfixed[[v]]
  if (is.factor(x)) x <- as.numeric(x)
  if (is.factor(y)) y <- as.numeric(y)
  data.frame(variable = v, main_mean = mean(x, na.rm = TRUE), nonfixed_mean = mean(y, na.rm = TRUE),
             smd = smd_plain(x, y))
}))
safe_write_csv(selection_balance, file.path(OUT, "A15_missing_fixed_cohort_versus_broader_sample.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

fit_broad <- glm(f, data = broad, family = binomial())
coef_broad <- bind_rows(
  coef_table(model_main, "missing/selection", "primary model(fixed cohort)", "unweighted", nrow(main), sum(main$event)),
  coef_table(fit_broad, "missing/selection", "broader sample(five wave participation not required ms10 complete)", "broad_refit",
             nrow(broad), sum(broad$event))
)
safe_write_csv(coef_broad, file.path(OUT, "A15_missing_broader_sample_model_comparison.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

met_broad <- bind_rows(
  model_metrics("missing/selection", "fixed cohort(primary model, frozen)", "frozen", main$event, p_main_main,
                score = main$score, n_note = nrow(main), B = 500),
  model_metrics("missing/selection", "broader sample(primary model, frozen)", "frozen", broad$event, p_main_broad,
                score = broad$score, n_note = nrow(broad), B = 500),
  model_metrics("missing/selection", "broader sample(refit)", "broad_refit", broad$event,
                as.numeric(predict(fit_broad, type = "response")), score = broad$score,
                n_note = nrow(broad), B = 500)
)
safe_write_csv(met_broad, file.path(OUT, "A15_missing_broader_sample_performance_comparison.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

# ---------------------------------------------------------------------------
# 5. multiple imputation(MICE, m=20)
# ---------------------------------------------------------------------------
# F30: read each original vector directly; a same-name summarise() mean would
# overwrite that vector before subsequent quantiles. These are descriptive
# between-imputation percentiles, not pooled sampling confidence intervals.
summarise_mi_performance <- function(perf) {
  stopifnot(nrow(perf) > 1L, length(unique(perf$n_model)) == 1L)
  result <- data.frame(n_model = perf$n_model[1], events_model = mean(perf$events_model),
    event_rate = mean(perf$event_rate), event_rate_type = "unweighted")
  for (metric in c("auc_model", "auc_score", "cal_slope", "cal_intercept")) {
    values <- perf[[metric]]
    stopifnot(all(is.finite(values)))
    prefix <- if (metric == "auc_model") "auc" else metric
    result[[metric]] <- mean(values)
    result[[paste0(prefix, "_low")]] <- unname(quantile(values, .025, type = 7))
    result[[paste0(prefix, "_high")]] <- unname(quantile(values, .975, type = 7))
  }
  result$brier <- mean(perf$brier)
  result
}
mi_vars <- c("age", "female", "edu_group", "rural", "mem2011", "cesd10", "dep",
             "adl", "sleep_h", "pain_any", "hibp", "diab", "stroke", "ms10_2011", "fixed5")
set.seed(20260926)

# 5.1 missing covariates: participants with observed outcomes n=8,371
d_cov_mi <- pool %>% filter(!is.na(mem2015)) %>%
  transmute(across(all_of(mi_vars)), event = as.integer(mem2015 <= 4))
imp_cov <- mice(d_cov_mi, m = 20, maxit = 10, method = "pmm",
                printFlag = FALSE, seed = 20260926)
fits_cov <- with(imp_cov, glm(event ~ age + edu_group + dep + adl, family = binomial()))
pool_cov <- summary(pool(fits_cov), conf.int = TRUE)
coef_cov_mi <- data.frame(
  section = "missing", scenario = "broader sample MICE m=20(missing covariates; participants with observed outcomes)", model_type = "mice_cov",
  term = pool_cov$term, beta = pool_cov$estimate, se = pool_cov$std.error,
  p = pool_cov$p.value, OR = exp(pool_cov$estimate),
  OR_low = exp(pool_cov$`2.5 %`), OR_high = exp(pool_cov$`97.5 %`),
  n = nrow(d_cov_mi), events = sum(d_cov_mi$event), row.names = NULL
)
safe_write_csv(coef_cov_mi, file.path(OUT, "A15_missing_MICE_covariates_coefficients.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

perf_cov_mi <- bind_rows(lapply(1:20, function(i) {
  d <- complete(imp_cov, i)
  mm <- glm(event ~ age + edu_group + dep + adl, data = d, family = binomial())
  model_metrics("missing", paste0("MICE covariates imputation", i), "mice_cov", d$event,
                predict(mm, type = "response"), score = fixed_score(d), n_note = nrow(d), B = 0)
}))
perf_cov_sum <- summarise_mi_performance(perf_cov_mi) %>%
  mutate(section = "missing", scenario = "broader sample MICE m=20(missing covariates; participants with observed outcomes)",
             model_type = "mice_cov", interval_type = "between imputations 2.5–97.5%range(nonsampling CI)",
             n_note = nrow(d_cov_mi), .before = 1)
safe_write_csv(perf_cov_sum, file.path(OUT, "A15_missing_MICE_covariates_performance_summary.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

# 5.2 missing outcome: baseline unimpaired participants n=10,422, imputation 2015 memory score
d_out_mi <- pool %>% select(all_of(mi_vars), mem2015)
imp_out <- mice(d_out_mi, m = 20, maxit = 10, method = "pmm",
                printFlag = FALSE, seed = 20260926)
fits_out <- with(imp_out, glm(I(mem2015 <= 4) ~ age + edu_group + dep + adl, family = binomial()))
pool_out <- summary(pool(fits_out), conf.int = TRUE)
coef_out_mi <- data.frame(
  section = "missing", scenario = "broader sample MICE m=20(missing outcome, MAR; including fixed cohort auxiliary variables)", model_type = "mice_outcome",
  term = pool_out$term, beta = pool_out$estimate, se = pool_out$std.error,
  p = pool_out$p.value, OR = exp(pool_out$estimate),
  OR_low = exp(pool_out$`2.5 %`), OR_high = exp(pool_out$`97.5 %`),
  n = nrow(d_out_mi), events = NA_integer_, row.names = NULL
)
safe_write_csv(coef_out_mi, file.path(OUT, "A15_missing_MICE_outcome_coefficients.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

perf_out_mi <- bind_rows(lapply(1:20, function(i) {
  d <- complete(imp_out, i); d$event <- as.integer(d$mem2015 <= 4)
  mm <- glm(event ~ age + edu_group + dep + adl, data = d, family = binomial())
  model_metrics("missing", paste0("MICE outcome imputation", i), "mice_outcome", d$event,
                predict(mm, type = "response"), score = fixed_score(d), n_note = nrow(d), B = 0)
}))
perf_out_sum <- summarise_mi_performance(perf_out_mi) %>%
  mutate(section = "missing", scenario = "broader sample MICE m=20(missing outcome, MAR; including fixed cohort auxiliary variables)",
             model_type = "mice_outcome", interval_type = "between imputations 2.5–97.5%range(nonsampling CI)",
             n_note = nrow(d_out_mi), .before = 1)
safe_write_csv(perf_out_sum, file.path(OUT, "A15_missing_MICE_outcome_performance_summary.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

# ---------------------------------------------------------------------------
# 6. 2018 measurement differences: T09 prespecified restriction ≥60 years old / restricted to complete memory items
# ---------------------------------------------------------------------------
main18 <- main %>% mutate(mem18_complete = !is.na(imrc2018) & !is.na(dlrc2018))
scenarios_2018 <- list(
  "entire primary analysis sample" = main18,
  "2018 aged at least 60 years old" = main18 %>% filter(age2018 >= 60),
  "2018 complete memory items" = main18 %>% filter(mem18_complete),
  "2018 aged at least 60 years old with complete memory items" = main18 %>% filter(age2018 >= 60, mem18_complete)
)
coef_2018 <- bind_rows(lapply(names(scenarios_2018), function(nm) {
  d <- scenarios_2018[[nm]]
  mm <- glm(f, data = d, family = binomial())
  coef_table(mm, "2018 year measurement differences", nm, "refit_subset", nrow(d), sum(d$event))
}))
safe_write_csv(coef_2018, file.path(OUT, "A15_2018_measurement_differences_model_comparison.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

met_2018 <- bind_rows(lapply(names(scenarios_2018), function(nm) {
  d <- scenarios_2018[[nm]]
  mm <- glm(f, data = d, family = binomial())
  frozen <- plogis(fixed_lp_manual(d, b_main))
  bind_rows(
    model_metrics("2018 year measurement differences", paste0(nm, "(frozen model)"), "frozen", d$event, frozen,
                  score = d$score, n_note = nrow(d), B = 0),
    model_metrics("2018 year measurement differences", paste0(nm, "(refit)"), "refit_subset", d$event,
                  predict(mm, type = "response"), score = d$score, n_note = nrow(d), B = 0)
  )
}))
safe_write_csv(met_2018, file.path(OUT, "A15_2018_measurement_differences_performance_comparison.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

# ---------------------------------------------------------------------------
# 7. cutoff and outcome definition sensitivity
# ---------------------------------------------------------------------------
cut_rows <- list(); cut_met <- list()
for (t in 3:6) {
  d <- main %>% filter(mem2011 > t) %>% mutate(event_t = as.integer(mem2015 <= t))
  mm <- glm(event_t ~ age + edu_group + dep + adl, data = d, family = binomial())
  cut_rows[[length(cut_rows) + 1]] <- coef_table(mm, "cutoff", paste0("internal cutoff ≤", t), "refit_cut",
                                                 nrow(d), sum(d$event_t))
  cut_met[[length(cut_met) + 1]] <- model_metrics("cutoff", paste0("internal cutoff ≤", t, "(refit)"),
                                                   "refit_cut", d$event_t,
                                                   predict(mm, type = "response"),
                                                   score = fixed_score(d), n_note = nrow(d), B = 0)
}
coef_cut <- bind_rows(cut_rows)
met_cut <- bind_rows(cut_met)
cut2015 <- floor(mean(main$mem2015, na.rm = TRUE) - sd(main$mem2015, na.rm = TRUE))
ref_cut <- base %>% filter(fixed5 == 1, !is.na(mem2011))
cut2011_raw <- mean(ref_cut$mem2011, na.rm = TRUE) - sd(ref_cut$mem2011, na.rm = TRUE)
cut2011 <- floor(cut2011_raw)
cut_note <- data.frame(
  definition = c("2011 mean−1SD(primary cutoff)", "2015 mean−1SD(reestimate in same risk set)", "exact continuous cutoff mean−1SD(not rounded)"),
  raw_cut = c(cut2011_raw, mean(main$mem2015) - sd(main$mem2015), cut2011_raw),
  integer_cut = c(cut2011, cut2015, 4),
  note = c("A14 frozen primary cutoff ≤4(entire fixed cohort mean−1SD)", "2015 reestimate raw=4.797; by SAP floor rule is ≤4(if rounded to ≤5, see ≤5 sensitivity)", "memory scores are integers, exact cutoff equivalent to ≤4")
)
safe_write_csv(coef_cut, file.path(OUT, "A15_cutoff_internal_model_comparison.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
safe_write_csv(met_cut, file.path(OUT, "A15_cutoff_internal_performance_comparison.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
safe_write_csv(cut_note, file.path(OUT, "A15_cutoff_definition_summary.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

# continuous change and baseline adjustment
main$delta <- main$mem2015 - main$mem2011
main$decline2 <- as.integer(main$delta <= -2)
lm_delta <- lm(delta ~ age + edu_group + dep + adl, data = main)
lm_adj <- lm(mem2015 ~ mem2011 + age + edu_group + dep + adl, data = main)
lm_tab <- function(fit, label) {
  co <- summary(fit)$coefficients
  ci <- confint.default(fit)
  data.frame(model = label, term = rownames(co), beta = co[, 1], se = co[, 2],
             t = co[, 3], p = co[, 4], beta_low = ci[, 1], beta_high = ci[, 2],
             row.names = NULL)
}
outcome_cont <- bind_rows(lm_tab(lm_delta, "unadjusted continuous change delta"),
                          lm_tab(lm_adj, "after baseline adjustment mem2015"))
fit_decline2 <- glm(decline2 ~ age + edu_group + dep + adl, data = main, family = binomial())
fit_adj_event <- glm(event ~ mem2011 + age + edu_group + dep + adl, data = main, family = binomial())
outcome_alt <- bind_rows(
  coef_table(model_main, "outcome_definition", "primary outcome: 4 year threshold crossing ≤4", "unweighted", nrow(main), sum(main$event)),
  coef_table(fit_decline2, "outcome_definition", "alternative: memory decline ≥2 points", "alt_binary", nrow(main), sum(main$decline2)),
  coef_table(fit_adj_event, "outcome_definition", "threshold crossing after baseline adjustment(including mem2011)", "baseline_adjusted",
             nrow(main), sum(main$event))
)
met_outcome_alt <- bind_rows(
  model_metrics("outcome_definition", "primary outcome: 4 year threshold crossing ≤4", "unweighted", main$event,
                predict(model_main, newdata = main, type = "response"), score = main$score,
                n_note = nrow(main), B = 0),
  model_metrics("outcome_definition", "alternative: memory decline ≥2 points", "alt_binary", main$decline2,
                predict(fit_decline2, type = "response"), score = main$score,
                n_note = nrow(main), B = 0),
  model_metrics("outcome_definition", "threshold crossing after baseline adjustment(including mem2011)", "baseline_adjusted", main$event,
                predict(fit_adj_event, type = "response"), score = main$score,
                n_note = nrow(main), B = 0)
)
safe_write_csv(outcome_cont, file.path(OUT, "A15_outcome_definition_continuous_change_model.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
safe_write_csv(outcome_alt, file.path(OUT, "A15_outcome_definition_model_comparison.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
safe_write_csv(met_outcome_alt, file.path(OUT, "A15_outcome_definition_performance_comparison.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

# external cutoff/mapping sensitivity summary(organize directly A12 formal results, do not recompute external primary model)
a12 <- readRDS(file.path(P, "output/11_external_validation/A12_results_2026-09-26.rds"))
ext_cut <- a12$perf %>%
  transmute(cohort, cut_definition, n_total, events_total, n_model, events_model,
            auc_model, auc_model_low, auc_model_high, auc_score_full,
            cal_slope, cal_intercept, brier)
safe_write_csv(ext_cut, file.path(OUT, "A15_cutoff_external_summary.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
ext_extra <- bind_rows(
  a12$extra_dep %>% transmute(cohort, sensitivity = "alternative depression cutoff", detail = dep_rule,
                              n, events, auc_model, auc_model_low, auc_model_high,
                              auc_score, cal_slope, cal_intercept, brier),
  a12$extra_score_alt %>% transmute(cohort, sensitivity = "score_alt(positional education coding)", detail = NA_character_,
                                    n, events, auc_model = NA_real_, auc_model_low = NA_real_,
                                    auc_model_high = NA_real_, auc_score = auc_score_alt,
                                    cal_slope = NA_real_, cal_intercept = NA_real_, brier = NA_real_),
  a12$extra_elsa_age %>% transmute(cohort, sensitivity = "ELSA exclude age 90 top coding", detail = NA_character_,
                                   n, events, auc_model, auc_model_low, auc_model_high,
                                   auc_score, cal_slope, cal_intercept, brier)
)
safe_write_csv(ext_extra, file.path(OUT, "A15_cutoff_additional_external_sensitivity.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

# ---------------------------------------------------------------------------
# 8. external cohort ADL missing data sensitivity(integer scores unaffected; continuous model ADL=0 vs complete cases)
# ---------------------------------------------------------------------------
adl_ext <- bind_rows(lapply(c("HRS", "ELSA", "SHARE"), function(co) {
  x <- readRDS(file.path(P, "external_data", co, "analysis",
                         paste0(co, "_sensitivity_risk_sets_2026-09-26.rds")))[[paste0(co, "__baseline_all_mean_sd")]]
  x_obs <- x; x_zero <- x; x_zero$adl[is.na(x_zero$adl)] <- 0
  p_obs <- plogis(fixed_lp_manual(x_obs, b_main))
  p_zero <- plogis(fixed_lp_manual(x_zero, b_main))
  a_obs <- auc_ci_unweighted(x$event, p_obs); a_zero <- auc_ci_unweighted(x$event, p_zero)
  b_obs <- mean((p_obs - x$event)^2, na.rm = TRUE); b_zero <- mean((p_zero - x$event)^2, na.rm = TRUE)
  data.frame(cohort = co, n_total = nrow(x), events_total = sum(x$event),
             adl_missing = sum(is.na(x$adl)),
             n_model_complete = sum(!is.na(p_obs)), n_model_zero = sum(!is.na(p_zero)),
             auc_observed = a_obs[["auc"]], auc_observed_low = a_obs[["low"]], auc_observed_high = a_obs[["high"]],
             auc_adl0 = a_zero[["auc"]], auc_adl0_low = a_zero[["low"]], auc_adl0_high = a_zero[["high"]],
             delta_auc = a_zero[["auc"]] - a_obs[["auc"]],
             brier_observed = b_obs, brier_adl0 = b_zero, row.names = NULL)
}))
safe_write_csv(adl_ext, file.path(OUT, "A15_missing_external_ADL_sensitivity.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

# ---------------------------------------------------------------------------
# 9. summary table and stability assessment
# ---------------------------------------------------------------------------
summary_internal <- bind_rows(met_weight, met_ipw, met_broad, perf_cov_sum, perf_out_sum,
                              met_2018, met_cut, met_outcome_alt) %>%
  mutate(domain = section, .before = 1) %>%
  select(domain, scenario, model_type, n_model, events_model, event_rate, event_rate_type,
         auc_model, auc_low, auc_high, auc_score, cal_slope, cal_slope_low, cal_slope_high,
         cal_intercept, cal_intercept_low, cal_intercept_high, brier, interval_type, n_note)
safe_write_csv(summary_internal, file.path(OUT, "A15_internal_sensitivity_summary.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

summary_external <- bind_rows(
  ext_cut %>% mutate(domain = "external cutoff", scenario = paste0(cohort, " / ", cut_definition),
                     model_type = "fixed_external", .before = 1) %>%
    transmute(domain, scenario, model_type, n_model, events_model, event_rate = events_model / n_model,
              auc_model, auc_low = auc_model_low, auc_high = auc_model_high, auc_score = auc_score_full,
              cal_slope, cal_slope_low = NA_real_, cal_slope_high = NA_real_,
              cal_intercept, cal_intercept_low = NA_real_, cal_intercept_high = NA_real_,
              brier, n_note = NA_integer_),
  ext_extra %>% mutate(domain = "external supplement", scenario = paste0(cohort, " / ", sensitivity),
                       model_type = "fixed_external", .before = 1) %>%
    transmute(domain, scenario, model_type, n_model = n, events_model = events, event_rate = events / n,
              auc_model, auc_low = auc_model_low, auc_high = auc_model_high, auc_score,
              cal_slope, cal_slope_low = NA_real_, cal_slope_high = NA_real_,
              cal_intercept, cal_intercept_low = NA_real_, cal_intercept_high = NA_real_,
              brier, n_note = NA_integer_)
)
safe_write_csv(summary_external, file.path(OUT, "A15_external_sensitivity_summary.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

key_direction <- function(tab, scenario, term) {
  z <- tab %>% filter(.data$scenario == !!scenario, .data$term == !!term)
  if (nrow(z) == 0) return(NA_character_)
  ifelse(is.na(z$OR[1]), NA_character_,
         ifelse(z$OR[1] > 1 & z$OR_low[1] > 1, "OR>1 and CI>1",
                ifelse(z$OR[1] < 1 & z$OR_high[1] < 1, "OR<1 and CI<1",
                       ifelse(z$OR[1] > 1, "OR>1(CI including 1)", "OR<1(CI including 1)"))))
}
coef_all <- bind_rows(coef_weight, coef_ipw, coef_broad, coef_2018, coef_cut, outcome_alt, coef_cov_mi, coef_out_mi)
key_terms <- c("age", "edu_groupHigh school or above", "dep", "adl")
key_or <- bind_rows(lapply(unique(coef_all$scenario), function(sc) {
  bind_rows(lapply(key_terms, function(tm) {
    z <- coef_all %>% filter(.data$scenario == sc, .data$term == tm)
    if (nrow(z) == 0) return(NULL)
    data.frame(scenario = sc, model_type = z$model_type[1], term = tm,
               OR = z$OR[1], OR_low = z$OR_low[1], OR_high = z$OR_high[1],
               direction = key_direction(coef_all, sc, tm))
  }))
}))
safe_write_csv(key_or, file.path(OUT, "A15_key_OR_directional_stability.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

stability <- data.frame(
  domain = c("weights", "IPW", "2018 year measurement differences", "cutoff", "missing/selection", "overall"),
  conclusion = c(
    "weighted refit AUC 0.701, unweighted primary model AUC 0.702; age, education, depression and ADL association direction consistent with primary model.",
    "IPW refitted model AUC 0.702; selection differences after IPW maximum after |SMD|=0.042; key OR direction consistent with primary model.",
    "primary outcome 2011→2015 do not use 2018 year data; restricted to 2018 aged at least 60 age restriction or complete memory items OR consistent direction, AUC 0.677–0.698(slight decline).",
    "internal cutoff ≤3–≤6 of AUC 0.699–0.711; external S0/S1/S2 and alternative cutoff AUC all≥0.69.",
    "broader sample n=8,296, event rate 22.8%, refit AUC 0.729; MICE covariates/outcome AUC 0.728/0.734; stable direction, but absolute risk and calibration vary by population.",
    "main conclusions directionally stable under prespecified sensitivities; absolute risk, calibration intercept and AUC affected by analysis population/and cutoff definition; interpret using A14 frozen model(v1.0)as reference."
  ),
  status = c("stable", "stable", "stable(does not directly affect primary outcome)", "stable", "stable direction; interpret absolute risk and calibration cautiously", "passed")
)
safe_write_csv(stability, file.path(OUT, "A15_assessment_of_conclusion_stability.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

# ---------------------------------------------------------------------------
# 10. figure: internal sensitivity AUC and key OR direction
# ---------------------------------------------------------------------------
plot_tab <- summary_internal %>%
  filter(!grepl("imputation", scenario), !grepl("MICE", scenario), !is.na(auc_model)) %>%
  distinct(domain, scenario, .keep_all = TRUE) %>%
  mutate(label = paste0(domain, " | ", scenario))
plot_tab <- plot_tab[order(plot_tab$auc_model), ]
# F29: include every confidence interval, including the decline-outcome analysis.
auc_xlim <- range(c(plot_tab$auc_low, plot_tab$auc_high), na.rm = TRUE) + c(-0.02, 0.02)
png(file.path(FIG, "FigureS8_A15_sensitivity_2026-09-26.png"), width = 22, height = 12, units = "in", res = 300)
par(mfrow = c(1, 2), mar = c(4.5, 25, 3, 1.5), family = "sans")
yy <- seq_len(nrow(plot_tab))
plot(plot_tab$auc_model, yy, xlim = auc_xlim, ylim = c(0.5, nrow(plot_tab) + 0.5),
     pch = 16, col = "#1F4E78", xlab = "AUC (continuous model)", ylab = "",
     main = "A15 internal sensitivity: discrimination", yaxt = "n", cex = 0.85)
segments(plot_tab$auc_low, yy, plot_tab$auc_high, yy, col = "#9DC3E6", lwd = 2)
points(plot_tab$auc_model, yy, pch = 16, col = "#1F4E78", cex = 0.9)
abline(v = summary_internal$auc_model[summary_internal$scenario == "primary model(unweighted, frozen)"][1], lty = 2, col = "#C00000")
axis(2, at = yy, labels = plot_tab$label, las = 1, cex.axis = 0.65)
grid(col = "grey90")
plot_or <- key_or %>% filter(term != "adl", !is.na(OR), !grepl("imputation", scenario))
plot_or$term_plot <- ifelse(plot_or$term == "age", "Age (per 5 years)",
                     ifelse(plot_or$term == "dep", "Depression", "Education high"))
plot_or$OR_plot <- ifelse(plot_or$term == "age", plot_or$OR^5, plot_or$OR)
plot_or$low_plot <- ifelse(plot_or$term == "age", plot_or$OR_low^5, plot_or$OR_low)
plot_or$high_plot <- ifelse(plot_or$term == "age", plot_or$OR_high^5, plot_or$OR_high)
plot_or <- plot_or[order(plot_or$term_plot, plot_or$OR_plot), ]
y2 <- seq_len(nrow(plot_or))
cols <- c("Age (per 5 years)" = "#1F4E78", "Depression" = "#C00000", "Education high" = "#2E7D32")
plot(plot_or$OR_plot, y2, log = "x", xlim = range(c(plot_or$low_plot, plot_or$high_plot, 1), na.rm = TRUE),
     ylim = c(0.5, nrow(plot_or) + 0.5), pch = 16, col = cols[plot_or$term_plot],
     xlab = "OR (log scale)", ylab = "", main = "A15: direction stability of key predictors", yaxt = "n", cex = 0.85)
segments(plot_or$low_plot, y2, plot_or$high_plot, y2, col = cols[plot_or$term_plot], lwd = 2)
points(plot_or$OR_plot, y2, pch = 16, col = cols[plot_or$term_plot], cex = 0.9)
abline(v = 1, lty = 2, col = "grey40")
axis(2, at = y2, labels = paste0(plot_or$term_plot, " | ", plot_or$scenario), las = 1, cex.axis = 0.55)
grid(col = "grey90")
legend("topright", legend = names(cols), col = cols, pch = 16, bty = "n", cex = 0.8)
dev.off()
tryCatch({
  tiff(file.path(FIG, "FigureS8_A15_sensitivity_2026-09-26.tif"), width = 22, height = 12, units = "in", res = 600, compression = "lzw")
  par(mfrow = c(1, 2), mar = c(4.5, 25, 3, 1.5), family = "sans")
  yy <- seq_len(nrow(plot_tab)); plot(plot_tab$auc_model, yy, xlim = auc_xlim, ylim = c(0.5, nrow(plot_tab) + 0.5), pch = 16, col = "#1F4E78", xlab = "AUC (continuous model)", ylab = "", main = "A15 internal sensitivity: discrimination", yaxt = "n", cex = 0.85); segments(plot_tab$auc_low, yy, plot_tab$auc_high, yy, col = "#9DC3E6", lwd = 2); points(plot_tab$auc_model, yy, pch = 16, col = "#1F4E78", cex = 0.9); abline(v = summary_internal$auc_model[summary_internal$scenario == "primary model(unweighted, frozen)"][1], lty = 2, col = "#C00000"); axis(2, at = yy, labels = plot_tab$label, las = 1, cex.axis = 0.65); grid(col = "grey90")
  plot(plot_or$OR_plot, y2, log = "x", xlim = range(c(plot_or$low_plot, plot_or$high_plot, 1), na.rm = TRUE), ylim = c(0.5, nrow(plot_or) + 0.5), pch = 16, col = cols[plot_or$term_plot], xlab = "OR (log scale)", ylab = "", main = "A15: direction stability of key predictors", yaxt = "n", cex = 0.85); segments(plot_or$low_plot, y2, plot_or$high_plot, y2, col = cols[plot_or$term_plot], lwd = 2); points(plot_or$OR_plot, y2, pch = 16, col = cols[plot_or$term_plot], cex = 0.9); abline(v = 1, lty = 2, col = "grey40"); axis(2, at = y2, labels = paste0(plot_or$term_plot, " | ", plot_or$scenario), las = 1, cex.axis = 0.55); grid(col = "grey90"); legend("topright", legend = names(cols), col = cols, pch = 16, bty = "n", cex = 0.8)
  dev.off()
}, error = function(e) { if (length(dev.list())) dev.off() })

# ---------------------------------------------------------------------------
# 11. save RDS / sessionInfo
# ---------------------------------------------------------------------------
saveRDS(list(main = main, broad = broad, main_ipw = main_ipw, pool_ipw = pool_ipw,
             model_main = model_main, fit_w = fit_w, fit_ipw = fit_ipw, fit_ipw_comb = fit_ipw_comb,
             fit_broad = fit_broad, imp_cov = imp_cov, imp_out = imp_out,
             flow = flow, weights_diag = weights_diag, coef_weight = coef_weight,
             met_weight = met_weight, ipw_diag = ipw_diag, ipw_selection_coef = ipw_selection_coef, balance_ipw = balance_ipw,
             coef_ipw = coef_ipw, met_ipw = met_ipw, selection_balance = selection_balance,
             coef_broad = coef_broad, met_broad = met_broad, coef_cov_mi = coef_cov_mi,
             perf_cov_sum = perf_cov_sum, coef_out_mi = coef_out_mi, perf_out_sum = perf_out_sum,
             coef_2018 = coef_2018, met_2018 = met_2018, coef_cut = coef_cut,
             met_cut = met_cut, cut_note = cut_note, outcome_cont = outcome_cont,
             outcome_alt = outcome_alt, met_outcome_alt = met_outcome_alt,
             ext_cut = ext_cut, ext_extra = ext_extra, adl_ext = adl_ext,
             summary_internal = summary_internal, summary_external = summary_external,
             key_or = key_or, stability = stability),
        file.path(OUT, "A15_results_2026-09-26.rds"))
safe_write_lines(capture.output(sessionInfo()), file.path(OUT, "A15_sessionInfo_2026-09-26.txt"), useBytes = TRUE)

cat("A15 analyses completed.\n")

# ---------------------------------------------------------------------------
# 12. A15 report
# ---------------------------------------------------------------------------
kable_md <- function(df) paste(capture.output(kable(df, format = "pipe", row.names = FALSE)), collapse = "\n")
get_metric <- function(tab, scenario, col) {
  z <- tab %>% filter(.data$scenario == !!scenario)
  if (nrow(z) == 0) return(NA_real_)
  as.numeric(z[[col]][1])
}
get_or_triplet <- function(sc, term, pow = 1) {
  z <- key_or %>% filter(.data$scenario == !!sc, .data$term == !!term)
  if (nrow(z) == 0) return(rep(NA_real_, 3))
  c(z$OR[1]^pow, z$OR_low[1]^pow, z$OR_high[1]^pow)
}
main_auc <- get_metric(summary_internal, "primary model(unweighted, frozen)", "auc_model")
w_auc <- get_metric(summary_internal, "2011 weights+PSU(frozen model)", "auc_model")
w_refit_auc <- get_metric(summary_internal, "2011 weights+PSU(refit)", "auc_model")
ipw_auc <- get_metric(summary_internal, "IPW(frozen model)", "auc_model")
ipw_refit_auc <- get_metric(summary_internal, "IPW(refit)", "auc_model")
broad_frozen_auc <- get_metric(summary_internal, "broader sample(primary model, frozen)", "auc_model")
broad_refit_auc <- get_metric(summary_internal, "broader sample(refit)", "auc_model")
cut_auc_range <- range(met_cut$auc_model, na.rm = TRUE)
ext_s0 <- ext_cut %>% filter(cut_definition == "baseline_all_mean_sd") %>% arrange(cohort)
ext_s1 <- ext_cut %>% filter(cut_definition == "complete_fu_mean_sd") %>% arrange(cohort)
ext_s2 <- ext_cut %>% filter(cut_definition == "fixed_charls_le4") %>% arrange(cohort)
max_smd_after <- max(abs(balance_ipw$smd_after), na.rm = TRUE)
max_smd_before <- max(abs(balance_ipw$smd_before), na.rm = TRUE)

or_show <- bind_rows(lapply(c("primary model(unweighted)", "2011 weights+PSU(refit)",
                              "IPW(selection into final sample)", "broader sample(five wave participation not required ms10 complete)",
                              "broader sample MICE m=20(missing outcome, MAR; including fixed cohort auxiliary variables)",
                              "threshold crossing after baseline adjustment(including mem2011)"), function(sc) {
  a <- get_or_triplet(sc, "age", 5); d <- get_or_triplet(sc, "dep", 1); e <- get_or_triplet(sc, "edu_groupHigh school or above", 1)
  data.frame(model = sc,
    age5 = fmt_ci(a[1], a[2], a[3]), depression = fmt_ci(d[1], d[2], d[3]),
    education_high = fmt_ci(e[1], e[2], e[3]))
}))

weight_show <- met_weight %>%
  transmute(scenario = scenario, n = format(n_model, big.mark = ","), events = format(events_model, big.mark = ","),
            `event_rate_see_note` = sprintf("%.1f%%", 100*event_rate),
            AUC = fmt_ci(auc_model, auc_low, auc_high), calibration_slope = fmt(cal_slope, 2),
            calibration_intercept = fmt(cal_intercept, 2), `Brier_score` = fmt(brier, 4))
ipw_show <- met_ipw %>%
  transmute(scenario = scenario, n = format(n_model, big.mark = ","), events = format(events_model, big.mark = ","),
            `event_rate_see_note` = sprintf("%.1f%%", 100*event_rate),
            AUC = fmt_ci(auc_model, auc_low, auc_high), calibration_slope = fmt(cal_slope, 2),
            calibration_intercept = fmt(cal_intercept, 2), `Brier_score` = fmt(brier, 4))
t2018_show <- met_2018 %>%
  transmute(scenario = scenario, n = format(n_model, big.mark = ","), events = format(events_model, big.mark = ","),
            event_rate = sprintf("%.1f%%", 100*event_rate),
            AUC = fmt_ci(auc_model, auc_low, auc_high), calibration_slope = fmt(cal_slope, 2),
            calibration_intercept = fmt(cal_intercept, 2), `Brier_score` = fmt(brier, 4))
cut_show <- met_cut %>%
  transmute(cutoff = scenario, n = format(n_model, big.mark = ","), events = format(events_model, big.mark = ","),
            event_rate = sprintf("%.1f%%", 100*event_rate),
            AUC = fmt_ci(auc_model, auc_low, auc_high), calibration_slope = fmt(cal_slope, 2),
            calibration_intercept = fmt(cal_intercept, 2), `Brier_score` = fmt(brier, 4))
missing_show <- bind_rows(
  met_broad %>% transmute(scenario = scenario, n = format(n_model, big.mark = ","), events = format(round(events_model), big.mark = ","),
                          event_rate = sprintf("%.1f%%", 100*event_rate),
                          AUC = fmt_ci(auc_model, auc_low, auc_high), calibration_slope = fmt(cal_slope, 2),
                          calibration_intercept = fmt(cal_intercept, 2), `Brier_score` = fmt(brier, 4)),
  perf_cov_sum %>% transmute(scenario = scenario, n = format(n_model, big.mark = ","), events = format(round(events_model), big.mark = ","),
                             event_rate = sprintf("%.1f%%", 100*event_rate),
                             AUC = paste0(fmt(auc_model, 6), "(between imputations ", fmt(auc_low, 6), "–", fmt(auc_high, 6), ")"),
                             calibration_slope = fmt(cal_slope, 2), calibration_intercept = fmt(cal_intercept, 2), `Brier_score` = fmt(brier, 4)),
  perf_out_sum %>% transmute(scenario = scenario, n = format(n_model, big.mark = ","), events = paste0(format(events_model, big.mark = ","), "(mean across imputations)"),
                             event_rate = sprintf("%.1f%%", 100*event_rate),
                             AUC = paste0(fmt(auc_model, 6), "(between imputations ", fmt(auc_low, 6), "–", fmt(auc_high, 6), ")"),
                             calibration_slope = fmt(cal_slope, 2), calibration_intercept = fmt(cal_intercept, 2), `Brier_score` = fmt(brier, 4))
)
ext_show <- summary_external %>%
  filter(grepl("baseline_all_mean_sd|complete_fu_mean_sd|fixed_charls_le4", scenario)) %>%
  transmute(scenario = scenario, n = n_model, events = events_model, event_rate = sprintf("%.1f%%", 100*event_rate),
            AUC = fmt_ci(auc_model, auc_low, auc_high), calibration_slope = fmt(cal_slope, 2),
            calibration_intercept = fmt(cal_intercept, 2), Brier = fmt(brier, 4))

report <- c(
"# A15 sensitivity analysis results(2026-09-26)", "",
"> task: SAP v3 §10 and execution board A15——weights(survey design), inverse probability weighting(IPW), 2018 year measurement differences, cutoff, missingness analysis.",
"> boundaries: unchanged A14 frozen model(v1.0), score table, primary cutoff and primary external validation definition; all results in this report are sensitivity analyses.primary analysis remains A14 of CHARLS fixed cohort complete case analysis, other samples or imputation results do not replace primary analysis.",
"> conclusion: under prespecified sensitivity definitions, **main conclusions directionally stable**; absolute risk, calibration intercept and AUC(area under curve)varies by analysis population or cutoff.", "",
"## one, analysis definitions and samples", "",
kable_md(flow %>% mutate(n = format(n, big.mark = ",")) %>% rename(steps = step, notes = note)), "",
"broader sample definition: requires only 2011 year memory score>4, 2015 measurable memory in year and age, education, depression, ADL complete four variable data, not required 2013/2018/2020 year ms10 complete.",
"terminology: in this report“primary model”all refer to A14 frozen model(v1.0)in original 4,441 participant complete case primary sample performance; “refit”all are sensitivity analyses, do not replace primary model.",
paste0("primary analysis sample n=", format(nrow(main), big.mark = ","), ", event count=", format(sum(main$event), big.mark = ","), "(", sprintf("%.1f%%", 100*mean(main$event)), ")."),
paste0("broader sample n=", format(nrow(broad), big.mark = ","), ", event count=", format(sum(broad$event), big.mark = ","), "(", sprintf("%.1f%%", 100*mean(broad$event)), ")."), "",
"## two, weight sensitivity(2011 year baseline individual weights + communityID is PSU(primary sampling unit))", "",
kable_md(weight_show), "",
paste0("survey weighted event rate is ", sprintf("%.1f%%", 100*weighted.mean(main$event, main$weight)), "; weighted frozen model AUC=", fmt(w_auc), ", unweighted primary model AUC=", fmt(main_auc), "."),
"age, education, depression and ADL of OR(odds ratio)direction consistent with primary model; weight design in `A15_weights_design_diagnostics.csv`.",
"note: calibration slope in weighted refit row/intercept is apparent calibration of weighted model on same data; frozen model row evaluates A14 calibration of fixed formula in weighted population.",
"in the table“event count”is the unweighted event count; weighted row reports survey weighted event rate(658/4,441 = 14.8%, after weighting 15.1%).", "",
"## three, IPW(inverse probability weighting): probability weighting for selection into final sample", "",
paste0("IPW target pool n=", format(ipw_diag$selection_pool_n, big.mark = ","), "(2011 age in year≥45 years old, baseline memory score>4 and complete selection model variables), final analysis sample n=", format(ipw_diag$selected_n, big.mark = ","), "; "),
paste0("range of predicted selection probabilities ", fmt(ipw_diag$p_min, 3), "–", fmt(ipw_diag$p_max, 3), "; stabilized IPW effective sample size of ", format(round(ipw_diag$effective_n), big.mark = ","), "."),
"selection model includes age, sex, education, urban or rural, 2011 memory score, CES-D-10 and ADL; coefficients and OR see `A15_IPW_selection_model.csv`.",
paste0("IPW maximum before |SMD|(standardized mean difference)=", fmt(max_smd_before, 3), ", IPW maximum after |SMD|=", fmt(max_smd_after, 3), "."),
kable_md(ipw_show), "",
paste0("IPW frozen model AUC=", fmt(ipw_auc), "; IPW refitted model AUC=", fmt(ipw_refit_auc), "."),
"IPW age after, education, depression values OR direction consistent with primary model, ADL association remains close to zero; unmeasured selection and death cannot be addressed by IPW adjustment.",
"IPW target pool includes only complete selection model variables 10,326 participants; ADL participants missing selection variables excluded IPW.",
"in the table“event count”is the unweighted event count; IPW event rate in row is IPW weighted event rate(unweighted 14.8%, after weighting 16.6%).", "",
"in subsequent alternative outcome table“IPW(selection into final sample)”and in this section IPW refitted model is the same model.", "",
"## four, 2018 year measurement differences", "",
"primary outcome is 2011→2015 memory score threshold crossing, **2018 year measurement differences do not enter primary outcome directly**.",
"this section follows T09 prespecified definition, restricted respectively to 2018 aged at least 60 years old, 2018 complete immediate and delayed recall in year, and satisfying both criteria; refit same four variable model(age, education, depression, ADL).",
kable_md(t2018_show), "",
"under three restrictions OR direction consistent with primary model; AUC 0.677–0.698, slightly lower than primary sample, but retains moderate discrimination; T09 reported 2015→2018 common item change under restriction ≥60 direction consistent in age eligible or memory complete participants(see T09 report).",
"notes: T09 reported n=3,761(≥60 years old)is ms10 five wave fixed cohort; this table n=2,872(≥60 years old)is final memory outcome risk sample(unimpaired at baseline, 2015 memory available, complete four variable data), sample definitions differ, sample sizes should not be directly compared.",
"therefore cannot treat 2018 module differences as alternative explanation of primary conclusion, but also must not state“complete measurement equivalence demonstrated”.", "",
"## five, cutoff and outcome definition sensitivity", "",
kable_md(cut_show), "",
"in the table“internal cutoff ≤t”denotes baseline memory score>t, follow up memory score≤t, both baseline and follow up thresholds change.",
paste0("internal cutoff ≤3–≤6 of AUC range is ", fmt(cut_auc_range[1]), "–", fmt(cut_auc_range[2]), "."),
paste0("2011 primary cutoff raw=", fmt(cut_note$raw_cut[1], 3), ", as prespecified floor(round down)rule is ≤", cut_note$integer_cut[1], "; 2015 reestimate within same risk set raw=", fmt(cut_note$raw_cut[2], 3), ", by floor rule remains ≤", cut_note$integer_cut[2], "."),
"supplementary: if ordinary rounding used, 2015 cutoff is ≤5; A15 full results reported ≤5 cutoff sensitivity(AUC 0.706), primary conclusion unchanged.",
"external cutoff definition: S0=full baseline sample in each cohort mean−1SD(mean−1 standard deviations); S1=complete-fu mean−1SD; S2=fixed CHARLS ≤4.also check alternative depression cutoff, score_alt(positional education coding)and ELSA age 90 top coding, see `A15_cutoff_external_summary.csv` and `A15_cutoff_additional_external_sensitivity.csv`.",
"", "### alternative outcome definition", "",
kable_md(or_show %>% rename(model = model, `age_per_5_years` = age5, depression = depression, `education_high_school_or_above` = education_high)), "",
"unadjusted continuous change(Δ=2015 memory score−2011 memory score)affected by baseline position and regression to mean: depression, ADL and Δ association direction may differ from primary outcome OR opposite direction.",
"in baseline adjusted change and threshold crossing models(both include 2011 memory score)in, age, education, depression direction agrees with primary outcome, ADL association remains close to zero.",
"therefore A15 unadjusted continuous change should not replace primary outcome; primary outcome remains 4 year threshold crossing.",
"note: baseline adjusted threshold crossing model additionally includes 2011 memory score, its AUC not comparable with primary model AUC direct comparison; model only assesses regression to mean and directional sensitivity.", "",
"## six, missingness and selection analysis", "",
kable_md(missing_show), "",
"in the table“fixed cohort(primary model, frozen)”means A14 frozen model in original 4,441 participant primary sample; “broader sample(primary model, frozen)”means applying A14 frozen model to 8,296 participant broader sample; “broader sample(refit)”means refitting same four variable model in broader sample.",
"MICE parentheses in this row show 20 across imputations 2.5%–97.5% between imputation range, not a sampling distribution or confidence interval; parentheses in other unweighted rows show DeLong 95% CI, weighted/IPW row represents cluster bootstrap 95% CI.",
"MICE continuous model refitted and evaluated in each imputed dataset, represents apparent performance; calibration slope approximately 1, intercept approximately 0 does not represent external calibration.integer score applies the original fixed rules.event count in missing outcome analysis is 20 mean across imputations, not an observed integer event count.",
paste0("broader sample(five wave participation not required ms10 complete)frozen model AUC=", fmt(broad_frozen_auc), "; broader sample refit AUC=", fmt(broad_refit_auc), "."),
paste0("MICE(multiple imputation, m=20)missing covariates AUC=", fmt(perf_cov_sum$auc_model), "; missing outcome AUC=", fmt(perf_out_sum$auc_model), "."),
"event rate increases after relaxing fixed cohort restriction, calibration intercept changes, fixed cohort selection affects absolute risk; but age, education, depression direction stable.",
"external ADL missing data sensitivity: HRS and SHARE in, complete case analysis and ADL=0 compared with imputation AUC differences are small; integer score unaffected by ADL missingness.",
"external ADL missing Brier scores not directly comparable across different samples, emphasize AUC differences primarily.",
"outcome imputation includes fixed cohort indicator as auxiliary variable; result is MAR(missing at random)sensitivity under assumption, does not imply ignorable missingness mechanism.",
"MICE imputation variables include age, sex, education, urban or rural, 2011 memory score, CES-D, ADL, sleep, pain, hypertension, diabetes, stroke, 2011 ms10 and fixed cohort indicator; analysis model remains four variables.",
"FigureS8 left panel excludes MICE rows(interval represents between imputation range), MICE see results table in this section.", "",
"### SAP §10 remaining items", "",
"| SAP items | A15 handling |", "|---|---|",
"| S6 excludes disputed variables(pain/sleep/PA) | final A14 model already excludes these variables; primary analysis already satisfies, no separate model. |",
"| S7 outcome ③ confirmed at two consecutive waves | memory scores only 2011/2015 two waves, cannot construct two consecutive wave confirmation; continuous change and decline already used≥2 points as alternative. |",
"| S8 ADL/IADL count version | ADL binary term retained in continuous model but scores 0 points; IADL already in A01/T1 exclude from primary score, do not reopen primary model. |",
"| S9 SHARE alternative outcomes | SHARE external validation already uses recall scores, see A12/A15 external summary. |",
"| S10 baseline adjustment | baseline adjusted change and threshold crossing models completed. |", "",
"## seven, overall stability assessment", "",
kable_md(stability %>% rename(domain = domain, conclusion = conclusion, assessment = status)), "",
"**overall conclusion**: A15 no evidence sufficient to overturn A14 directional stability of final model.",
"weights, IPW, 2018 year measurement difference restrictions, internal cutoffs and MICE in missingness analysis, age, education, depression and ADL association directions remain consistent; external cohort S0/S1/S2 discrimination stable under alternative cutoffs; ",
"but broader sample and outcome imputation show event rate, calibration intercept and absolute risk changes; manuscript should state: tool is for screening/triage, not for fixed absolute risk interpretation.",
"limitations: IPW cannot adjust unmeasured selection and death; outcome multiple imputation assumes MAR(missing at random)as assumption, competing risks unaddressed; results for sensitivity interpretation only, do not replace primary analysis.", "",
"## eight, output file", "",
"- `A15_sample_flow_2026_09_26.csv`",
"- `A15_weights_design_diagnostics.csv`, `A15_weights_model_comparison.csv`, `A15_weights_performance_comparison.csv`",
"- `A15_IPW_diagnostics.csv`, `A15_IPW_selection_model.csv`, `A15_IPW_balance.csv`, `A15_IPW_model_comparison.csv`, `A15_IPW_performance_comparison.csv`",
"- `A15_2018_measurement_differences_model_comparison.csv`, `A15_2018_measurement_differences_performance_comparison.csv`",
"- `A15_cutoff_internal_model_comparison.csv`, `A15_cutoff_internal_performance_comparison.csv`, `A15_cutoff_definition_summary.csv`, `A15_cutoff_external_summary.csv`, `A15_cutoff_additional_external_sensitivity.csv`",
"- `A15_outcome_definition_continuous_change_model.csv`, `A15_outcome_definition_model_comparison.csv`, `A15_outcome_definition_performance_comparison.csv`",
"- `A15_missing_fixed_cohort_versus_broader_sample.csv`, `A15_missing_broader_sample_model_comparison.csv`, `A15_missing_broader_sample_performance_comparison.csv`",
"- `A15_missing_MICE_covariates_coefficients.csv`, `A15_missing_MICE_covariates_performance_summary.csv`, `A15_missing_MICE_outcome_coefficients.csv`, `A15_missing_MICE_outcome_performance_summary.csv`",
"- `A15_missing_external_ADL_sensitivity.csv`",
"- `A15_internal_sensitivity_summary.csv`, `A15_external_sensitivity_summary.csv`, `A15_key_OR_directional_stability.csv`, `A15_assessment_of_conclusion_stability.csv`",
"- `FigureS8_A15_sensitivity_2026-09-26.png/.tif`",
"- `A15_results_2026-09-26.rds`, `A15_sessionInfo_2026-09-26.txt`",
""
)
safe_write_lines(report, file.path(OUT, "A15_sensitivity_analysis_report_2026_09_26.md"), useBytes = TRUE)
cat("A15 report written.\n")
