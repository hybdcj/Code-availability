# ============================================================================
# 14_attrition_selection.R —— T07: attrition, fixed cohort selection bias and IPW sensitivity analysis
# version: v1.0(2026-09-20)
# input: output/02_five_wave_data_and_quality_control/waves_cognition_2026-09-13.rds
# output: 
#   output/04_longitudinal_cohort/attrition_pattern_2026-09-20.csv
#   output/04_longitudinal_cohort/attrition_balance_2026-09-20.csv
#   output/04_longitudinal_cohort/selection_model_2026-09-20.csv
#   output/04_longitudinal_cohort/ipw_diagnostics_2026-09-20.csv
#   output/04_longitudinal_cohort/ipw_longitudinal_models_2026-09-20.csv
#   output/04_longitudinal_cohort/ipw_cumulative_incidence_2026-09-20.csv
#   output/04_longitudinal_cohort/attrition_selection_2026-09-20.md
#   output/09_figures_and_tables/FigureS_attrition_selection_2026-09-20.png/.tif
#
# important definition: 
#   target population = 2011 age in year ≥45 years old and baseline ms10(orientation+serial subtraction of 7)scorable participants.
#   meets criterion and all five waves ms10 scorable participants constitute age eligible fixed cohort n=5,672.
#   old longitudinal analysis previously used 5,884 participants, including 212 participants aged<45 or participants with missing age; 
#   2026-09-20 already changed T04-T06 harmonized to age eligible cohort definition.
# ============================================================================

invisible(try(Sys.setlocale("LC_ALL", "zh_CN.UTF-8"), silent = TRUE))
if (dir.exists("D:/R/library")) .libPaths(c("D:/R/library", .libPaths()))
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(survey); library(survival) })
options(survey.lonely.psu = "adjust")

P <- Sys.getenv("CHARLS_PROJ", "D:/BaiduSyncdisk/research/Xu_Zhonghao")
OUT <- file.path(P, "output", "04_longitudinal_cohort")
FIG_DIR <- file.path(P, "output", "09_figures_and_tables")
if (!dir.exists(OUT)) dir.create(OUT, recursive = TRUE)
if (!dir.exists(FIG_DIR)) dir.create(FIG_DIR, recursive = TRUE)

w <- readRDS(file.path(P, "output", "02_five_wave_data_and_quality_control", "waves_cognition_2026-09-13.rds"))
w$KEY <- ifelse(is.na(w$IDn), w$ID, w$IDn)
w$ms10 <- as.numeric(w$orient) + as.numeric(w$ser7)
stopifnot(!anyDuplicated(w[, c("KEY", "year")]))

# number of waves with scorable cognition.
nwaves <- w %>% filter(!is.na(ms10)) %>% count(KEY, name = "n_waves")
b2011 <- w %>% filter(year == 2011)
eligible_all <- b2011 %>% filter(!is.na(age), age >= 45)
eligible_ms10 <- eligible_all %>% filter(!is.na(ms10)) %>% left_join(nwaves, by = "KEY") %>%
  mutate(fixed = as.integer(n_waves == 5))
eligible_all <- eligible_all %>% left_join(nwaves, by = "KEY") %>% mutate(fixed = as.integer(!is.na(n_waves) & n_waves == 5))

stopifnot(nrow(eligible_all) == 16831)
stopifnot(nrow(eligible_ms10) == 13710, sum(eligible_ms10$fixed) == 5672)

# ------------------------- 1. attrition patterns and retention by wave -------------------------
pattern_tab <- eligible_ms10 %>% count(n_waves, name = "n") %>%
  mutate(percent = 100 * n / sum(n)) %>%
  arrange(n_waves)
colnames(pattern_tab) <- c("waves_with_ms10", "n", "percent")

ids <- eligible_ms10$KEY
retention <- w %>% filter(KEY %in% ids) %>%
  group_by(year) %>%
  summarise(interviewed_n = n(), ms10_n = sum(!is.na(ms10)), .groups = "drop") %>%
  mutate(retained_percent = 100 * ms10_n / nrow(eligible_ms10))
write.csv(pattern_tab, file.path(OUT, "attrition_pattern_2026-09-20.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(retention, file.path(OUT, "attrition_retention_2026-09-20.csv"), row.names = FALSE, fileEncoding = "UTF-8")

# ------------------------- 2. baseline comparison of included and excluded participants -------------------------
balance_vars <- c("age", "female", "rural", "ms10", "cesd10", "adl", "sleep_h", "pain_any", "hibp", "diab")
balance_labels <- c(age = "age(years old)", female = "female", rural = "rural", ms10 = "baseline ms10",
                    cesd10 = "CES-D-10 score", adl = "ADL limited", sleep_h = "sleep duration(hours)",
                    pain_any = "pain", hibp = "hypertension", diab = "diabetes")

smd_cont <- function(x1, x2) {
  x1 <- x1[!is.na(x1)]; x2 <- x2[!is.na(x2)]
  if (length(x1) < 2 || length(x2) < 2) return(NA_real_)
  s <- sqrt((var(x1) + var(x2)) / 2)
  if (s == 0) return(0)
  (mean(x1) - mean(x2)) / s
}
wmean <- function(x, w) sum(x * w, na.rm = TRUE) / sum(w[!is.na(x)], na.rm = TRUE)
wvar <- function(x, w) {
  keep <- !is.na(x) & !is.na(w)
  x <- x[keep]; w <- w[keep]
  mu <- sum(w * x) / sum(w)
  sum(w * (x - mu)^2) / sum(w)
}
smd_weighted <- function(x1, w1, x2) {
  keep1 <- !is.na(x1) & !is.na(w1); keep2 <- !is.na(x2)
  x1 <- x1[keep1]; w1 <- w1[keep1]; x2 <- x2[keep2]
  if (length(x1) < 2 || length(x2) < 2) return(NA_real_)
  s <- sqrt((wvar(x1, w1) + var(x2)) / 2)
  if (s == 0) return(0)
  (wmean(x1, w1) - mean(x2)) / s
}

fixed <- eligible_ms10 %>% filter(fixed == 1)
nonfixed <- eligible_ms10 %>% filter(fixed == 0)
balance <- do.call(rbind, lapply(balance_vars, function(v) {
  x1 <- fixed[[v]]; x2 <- nonfixed[[v]]
  data.frame(
    variable = v, label = balance_labels[[v]],
    fixed_n = sum(!is.na(x1)), fixed_mean = mean(x1, na.rm = TRUE), fixed_sd = sd(x1, na.rm = TRUE),
    nonfixed_n = sum(!is.na(x2)), nonfixed_mean = mean(x2, na.rm = TRUE), nonfixed_sd = sd(x2, na.rm = TRUE),
    smd = smd_cont(x1, x2),
    stringsAsFactors = FALSE
  )
}))
write.csv(balance, file.path(OUT, "attrition_balance_2026-09-20.csv"), row.names = FALSE, fileEncoding = "UTF-8")

# ------------------------- 3. selection probability model -------------------------
pred_vars <- c("age", "female", "edu_group", "rural", "ms10", "cesd10", "adl", "sleep_h", "pain_any", "hibp", "diab")
ps_data <- eligible_ms10 %>% filter(complete.cases(.[, pred_vars])) %>% mutate(fixed_int = as.integer(fixed))
ps_formula <- as.formula(paste("fixed_int ~", paste(c("age", "female", "edu_group", "rural", "ms10", "cesd10", "adl", "sleep_h", "pain_any", "hibp", "diab"), collapse = " + ")))
ps_fit <- glm(ps_formula, data = ps_data, family = binomial())
ps_data$p_hat <- pmin(pmax(predict(ps_fit, type = "response"), 0.02), 0.98)
ps_data$ipw_stab <- mean(ps_data$fixed_int) / ps_data$p_hat
ps_fixed <- ps_data %>% filter(fixed_int == 1)

selection_coef <- summary(ps_fit)$coefficients
selection_tab <- data.frame(
  term = rownames(selection_coef), estimate = selection_coef[, 1], se = selection_coef[, 2],
  p = selection_coef[, ncol(selection_coef)],
  OR = exp(selection_coef[, 1]),
  OR_lower = exp(selection_coef[, 1] - 1.96 * selection_coef[, 2]),
  OR_upper = exp(selection_coef[, 1] + 1.96 * selection_coef[, 2]),
  stringsAsFactors = FALSE
)
write.csv(selection_tab, file.path(OUT, "selection_model_2026-09-20.csv"), row.names = FALSE, fileEncoding = "UTF-8")

# IPW diagnostics; conventional effective sample size sum(w)^2 / sum(w^2).
ipw_weight <- ps_fixed$ipw_stab
ipw_diag <- data.frame(
  selection_model_n = nrow(ps_data), fixed_n = nrow(ps_fixed), selection_fraction = mean(ps_data$fixed_int),
  p_min = min(ps_data$p_hat), p_max = max(ps_data$p_hat), p_median = median(ps_data$p_hat),
  ipw_min = min(ipw_weight), ipw_median = median(ipw_weight), ipw_max = max(ipw_weight),
  ipw_mean = mean(ipw_weight), effective_n = sum(ipw_weight)^2 / sum(ipw_weight^2),
  stringsAsFactors = FALSE
)
write.csv(ipw_diag, file.path(OUT, "ipw_diagnostics_2026-09-20.csv"), row.names = FALSE, fileEncoding = "UTF-8")

# using IPW reweight fixed cohort, compare against baseline complete case pool.
target <- ps_data %>% filter(fixed_int == 0 | fixed_int == 1)
balance_ipw <- do.call(rbind, lapply(balance_vars, function(v) {
  # target mean is the full eligible complete-case pool; IPW reweights fixed cases to it.
  x_pool <- target[[v]]
  x_fixed <- ps_fixed[[v]]
  w_ipw <- ps_fixed$ipw_stab
  data.frame(
    variable = v, label = balance_labels[[v]],
    fixed_unweighted = mean(x_fixed, na.rm = TRUE),
    pool_unweighted = mean(x_pool, na.rm = TRUE),
    fixed_ipw_weighted = wmean(x_fixed, w_ipw),
    smd_before = smd_cont(x_fixed, x_pool),
    smd_after = smd_weighted(x_fixed, w_ipw, x_pool),
    stringsAsFactors = FALSE
  )
}))
balance_ipw$smd_abs_before <- abs(balance_ipw$smd_before)
balance_ipw$smd_abs_after <- abs(balance_ipw$smd_after)
write.csv(balance_ipw, file.path(OUT, "attrition_balance_ipw_2026-09-20.csv"), row.names = FALSE, fileEncoding = "UTF-8")

# ------------------------- 4. construct outcomes for age eligible fixed cohort -------------------------
fixed_ids <- eligible_ms10$KEY[eligible_ms10$fixed == 1]
# continuous cognition score matrix.
score_wide <- w %>% filter(KEY %in% fixed_ids) %>%
  select(KEY, year, ms10) %>%
  pivot_wider(names_from = year, values_from = ms10, names_prefix = "ms") %>%
  arrange(KEY)
# 2011 year education specific fixed cutoff.
base11 <- w %>% filter(KEY %in% fixed_ids, year == 2011) %>% select(KEY, edu_group)
thr <- w %>% filter(KEY %in% fixed_ids, year == 2011) %>%
  group_by(edu_group) %>% summarise(anchor = mean(ms10, na.rm = TRUE) - sd(ms10, na.rm = TRUE), .groups = "drop")
# impairment status matrix.
status_wide <- w %>% filter(KEY %in% fixed_ids) %>%
  select(KEY, year, ms10, edu_group) %>%
  left_join(thr, by = "edu_group") %>%
  mutate(impaired = as.integer(ms10 < anchor)) %>%
  select(KEY, year, impaired) %>%
  pivot_wider(names_from = year, values_from = impaired, names_prefix = "im") %>%
  arrange(KEY)
wide <- score_wide %>% left_join(status_wide, by = "KEY")

# baseline exposures and ten year cognitive change.
base_cov <- ps_data %>% select(KEY, ipw_stab, p_hat, fixed_int)
fixed_base <- w %>% filter(KEY %in% fixed_ids, year == 2011) %>%
  transmute(KEY, age, female, edu_group, rural,
            dep = as.integer(cesd10 >= 10), adl, sleep_h, pain = pain_any, pa = pa_any,
            hibp, diab, weight, communityID)
fixed_base <- fixed_base %>% left_join(base_cov, by = "KEY") %>% left_join(
  wide %>% transmute(
    KEY,
    delta = ms2020 - ms2011,
    first_index = case_when(
      im2013 == 1 ~ 1L, im2015 == 1 ~ 2L, im2018 == 1 ~ 3L, im2020 == 1 ~ 4L,
      TRUE ~ NA_integer_
    ),
    confirmed_index = case_when(
      im2013 == 1 & im2015 == 1 ~ 1L,
      im2015 == 1 & im2018 == 1 ~ 2L,
      im2018 == 1 & im2020 == 1 ~ 3L,
      TRUE ~ NA_integer_
    ),
    endpoint_2020 = im2020,
    ms10_2011 = ms2011
  ), by = "KEY"
) %>% mutate(
  time = ifelse(is.na(first_index), 9, c(2, 4, 7, 9)[first_index]),
  event = as.integer(!is.na(first_index)),
  ipw_combined = weight * ipw_stab
)
# standardize within fixed cohort IPW combined weights; constant weight scaling does not change regression estimates.
fixed_base$ipw_combined <- fixed_base$ipw_combined / mean(fixed_base$ipw_combined, na.rm = TRUE)

# ------------------------- 5. IPW longitudinal model -------------------------
exposures <- c("dep", "adl", "sleep_h", "pain", "pa", "hibp", "diab")
labels <- c(dep = "depression(CES-D-10 ≥10)", adl = "ADL limited", sleep_h = "sleep duration(per 1 hours)",
            pain = "pain", pa = "regular physical activity", hibp = "hypertension", diab = "diabetes")
extract_est <- function(fit, term, type = c("ratio", "beta")) {
  type <- match.arg(type); b <- coef(fit)[term]; se <- sqrt(vcov(fit)[term, term])
  if (type == "ratio") c(exp(b), exp(b - 1.96 * se), exp(b + 1.96 * se))
  else c(b, b - 1.96 * se, b + 1.96 * se)
}
ipw_rows <- list()
for (e in exposures) {
  vars_inc <- c(e, "age", "female", "edu_group", "endpoint_2020", "time", "event", "ipw_combined", "communityID")
  vars_change <- c(e, "age", "female", "edu_group", "delta", "ipw_combined", "communityID")
  d_inc <- fixed_base[complete.cases(fixed_base[, vars_inc]), ]
  d_change <- fixed_base[complete.cases(fixed_base[, vars_change]), ]
  des_inc <- svydesign(ids = ~communityID, weights = ~ipw_combined, data = d_inc)
  des_change <- svydesign(ids = ~communityID, weights = ~ipw_combined, data = d_change)
  rhs <- paste(e, "+ age + female + edu_group")
  f1 <- as.formula(paste("endpoint_2020 ~", rhs))
  f2 <- as.formula(paste("Surv(time, event) ~", rhs))
  f3 <- as.formula(paste("delta ~", rhs))
  wl <- svyglm(f1, design = des_inc, family = quasibinomial())
  wc <- svycoxph(f2, design = des_inc)
  wch <- svyglm(f3, design = des_change)
  u1 <- glm(f1, data = d_inc, family = binomial())
  u2 <- coxph(f2, data = d_inc)
  u3 <- lm(f3, data = d_change)
  specs <- list(
    list(model = "endpoint_2020_logistic", label = "2020 endpoint assessment", measure = "OR",
         w = extract_est(wl, e, "ratio"), u = extract_est(u1, e, "ratio"),
         n = nrow(d_inc), events = sum(d_inc$endpoint_2020), communities = length(unique(d_inc$communityID))),
    list(model = "first_event_cox", label = "first impairment Cox", measure = "HR",
         w = extract_est(wc, e, "ratio"), u = extract_est(u2, e, "ratio"),
         n = nrow(d_inc), events = sum(d_inc$event), communities = length(unique(d_inc$communityID))),
    list(model = "ten_year_change_lm", label = "ten year change(fixed cohort complete cases)", measure = "beta",
         w = extract_est(wch, e, "beta"), u = extract_est(u3, e, "beta"),
         n = nrow(d_change), events = nrow(d_change), communities = length(unique(d_change$communityID)))
  )
  for (z in specs) {
    ipw_rows[[length(ipw_rows) + 1L]] <- data.frame(
      model = z$model, model_label = z$label, exposure = e, exposure_label = labels[[e]],
      effect_measure = z$measure, n = z$n, events = z$events, communities = z$communities,
      ipw_est = z$w[1], ipw_lower = z$w[2], ipw_upper = z$w[3],
      unweighted_est = z$u[1], unweighted_lower = z$u[2], unweighted_upper = z$u[3],
      stringsAsFactors = FALSE
    )
  }
}
ipw_models <- bind_rows(ipw_rows)
write.csv(ipw_models, file.path(OUT, "ipw_longitudinal_models_2026-09-20.csv"), row.names = FALSE, fileEncoding = "UTF-8")

# ------------------------- 6. IPW cumulative incidence -------------------------
time_points <- c(0,2,4,7,9); first_limits <- c(0L,1L,2L,3L,4L); conf_limits <- c(0L,1L,2L,3L,3L)
ipw_inc <- list()
for (k in seq_along(time_points)) {
  d <- fixed_base
  d$e_first <- as.integer(!is.na(d$first_index) & d$first_index <= first_limits[k])
  d$e_conf <- as.integer(!is.na(d$confirmed_index) & d$confirmed_index <= conf_limits[k])
  des <- svydesign(ids = ~communityID, weights = ~ipw_combined, data = d[!is.na(d$ipw_combined), ])
  m1 <- svymean(~e_first, des, na.rm = TRUE); m2 <- svymean(~e_conf, des, na.rm = TRUE)
  ipw_inc[[k]] <- data.frame(
    time_years = rep(time_points[k], 2), outcome = c("first_episode", "confirmed_first_episode"),
    ipw_cumulative_incidence = c(as.numeric(m1), as.numeric(m2)),
    unweighted_cumulative_incidence = c(mean(d$e_first, na.rm = TRUE), mean(d$e_conf, na.rm = TRUE)),
    stringsAsFactors = FALSE
  )
}
ipw_inc <- bind_rows(ipw_inc)
write.csv(ipw_inc, file.path(OUT, "ipw_cumulative_incidence_2026-09-20.csv"), row.names = FALSE, fileEncoding = "UTF-8")

# ------------------------- 7. figures and tables: retention rate + SMD + IPW forest plot -------------------------
draw_attrition_figure <- function() {
  par(mfrow = c(1, 3), mar = c(5, 7, 3, 1), family = "sans")
  # Retention.
  plot(retention$year, retention$retained_percent, type = "b", pch = 16, lwd = 2, col = "#1F4E78",
       ylim = c(0, 100), xaxt = "n", ylab = "ms10 available (%)", xlab = "Year", main = "Follow-up availability")
  axis(1, at = retention$year); grid(col = "grey88")
  # SMD before/after.
  z <- balance_ipw[order(abs(balance_ipw$smd_before)), ]
  yy <- seq_len(nrow(z)); xlim <- range(c(z$smd_before, z$smd_after, -0.1, 0.1), na.rm = TRUE)
  plot(NA, xlim = xlim, ylim = c(0.5, nrow(z) + 0.5), yaxt = "n", xlab = "Standardized mean difference", ylab = "", main = "Covariate balance")
  abline(v = c(-0.1, 0, 0.1), lty = c(2,1,2), col = c("grey60","grey40","grey60"))
  segments(z$smd_before, yy - 0.10, z$smd_after, yy + 0.10, col = "grey80")
  points(z$smd_before, yy - 0.10, pch = 17, col = "#C00000", cex = 0.9)
  points(z$smd_after, yy + 0.10, pch = 16, col = "#1F4E78", cex = 0.9)
  axis(2, at = yy, labels = z$label, las = 1, cex.axis = 0.75)
  legend("topright", legend = c("Before IPW", "After IPW"), col = c("#C00000","#1F4E78"), pch = c(17,16), bty = "n", cex = 0.75)
  # IPW vs unweighted Cox.
  z2 <- ipw_models[ipw_models$model == "first_event_cox", ]
  z2 <- z2[order(z2$ipw_est), ]; yy <- seq_len(nrow(z2))
  xlim2 <- range(c(z2$ipw_lower, z2$ipw_upper, z2$unweighted_lower, z2$unweighted_upper, 1), na.rm = TRUE)
  plot(NA, xlim = xlim2, ylim = c(0.5, nrow(z2)+0.5), yaxt = "n", xlab = "Hazard ratio", ylab = "", main = "IPW sensitivity: first event")
  abline(v = 1, lty = 3, col = "grey60")
  segments(z2$unweighted_lower, yy - 0.1, z2$unweighted_upper, yy - 0.1, col = "#C00000", lwd = 2)
  points(z2$unweighted_est, yy - 0.1, pch = 17, col = "#C00000")
  segments(z2$ipw_lower, yy + 0.1, z2$ipw_upper, yy + 0.1, col = "#1F4E78", lwd = 2)
  points(z2$ipw_est, yy + 0.1, pch = 16, col = "#1F4E78")
  axis(2, at = yy, labels = z2$exposure_label, las = 1, cex.axis = 0.75)
  legend("topright", legend = c("Unweighted", "IPW (selection-adjusted)"), col = c("#C00000","#1F4E78"), pch = c(17,16), bty = "n", cex = 0.72)
}
png(file.path(FIG_DIR, "FigureS_attrition_selection_2026-09-20.png"), width = 15, height = 6, units = "in", res = 300)
draw_attrition_figure(); dev.off()
tryCatch({ tiff(file.path(FIG_DIR, "FigureS_attrition_selection_2026-09-20.tif"), width = 15, height = 6, units = "in", res = 600, compression = "lzw"); draw_attrition_figure(); dev.off() }, error = function(e) { if (length(dev.list())) dev.off() })

# ------------------------- 8. report -------------------------
fmt_pct <- function(x) sprintf("%.1f%%", 100*x)
report <- c(
  "# T07: attrition, fixed cohort selection bias and IPW sensitivity analysis(2026-09-20)",
  "",
  "> target population: 2011 age in year ≥45 years old and baseline ms10 scorable participants.",
  sprintf("> analyzable at baseline n=%d; fixed cohort meeting age eligibility n=%d; old 5,884 participant cohort includes 212 participants aged<45 or participants with missing age.", nrow(eligible_ms10), nrow(fixed_base)),
  "",
  "## one, attrition patterns",
  "",
  "| number of scorable waves out of five | participant count | proportion |",
  "|---:|---:|---:|",
  paste0("| ", pattern_tab[[1]], " | ", pattern_tab[[2]], " | ", sprintf("%.1f%%", pattern_tab[[3]]), " |"),
  "",
  "| year | actual interviewed count | ms10 scorable | relative to 2011 retain |",
  "|---:|---:|---:|---:|",
  paste0("| ", retention$year, " | ", retention$interviewed_n, " | ", retention$ms10_n, " | ", fmt_pct(retention$retained_percent/100), " |"),
  "",
  "## two, baseline differences between included and excluded participants",
  "",
  "| metrics | fixed cohort | excluded | SMD |",
  "|---|---:|---:|---:|",
  paste0("| ", balance$label, " | ", sprintf("%.2f", balance$fixed_mean), " | ", sprintf("%.2f", balance$nonfixed_mean), " | ", sprintf("%.2f", balance$smd), " |"),
  "",
  "> SMD absolute value closer to 0 indicates better balance; |SMD|≥0.10 indicates meaningful baseline differences.",
  "",
  "## three, selection probability and IPW",
  "",
  sprintf("- selection model complete cases: %d participants; of whom in fixed cohort %d participants.", ipw_diag$selection_model_n, ipw_diag$fixed_n),
  sprintf("- truncated predicted selection probability: median %.2f, range %.2f–%.2f.", ipw_diag$p_median, ipw_diag$p_min, ipw_diag$p_max),
  sprintf("- stabilized IPW: median %.2f, maximum %.2f, effective sample size %.0f.", ipw_diag$ipw_median, ipw_diag$ipw_max, ipw_diag$effective_n),
  "",
  "| metrics | unweighted fixed cohort | unweighted full pool | IPW reweighted | SMD before | SMD after |",
  "|---|---:|---:|---:|---:|---:|",
  paste0("| ", balance_ipw$label, " | ", sprintf("%.2f", balance_ipw$fixed_unweighted), " | ", sprintf("%.2f", balance_ipw$pool_unweighted), " | ", sprintf("%.2f", balance_ipw$fixed_ipw_weighted), " | ", sprintf("%.2f", balance_ipw$smd_before), " | ", sprintf("%.2f", balance_ipw$smd_after), " |"),
  "",
  "## four, IPW longitudinal model",
  "",
  "| model | exposure | n | IPW estimate | unweighted |",
  "|---|---|---:|---:|---:|",
  paste0("| ", ipw_models$model_label, " | ", ipw_models$exposure_label, " | ", ipw_models$n, " | ",
         sprintf("%.2f (%.2f–%.2f)", ipw_models$ipw_est, ipw_models$ipw_lower, ipw_models$ipw_upper), " | ",
         sprintf("%.2f (%.2f–%.2f)", ipw_models$unweighted_est, ipw_models$unweighted_lower, ipw_models$unweighted_upper), " |"),
  "",
  "## five, conclusion",
  "",
  sprintf("- age eligible fixed cohort as fraction of analyzable baseline participants %.1f%%.", 100 * nrow(fixed_base) / nrow(eligible_ms10)),
  sprintf("- IPW after SMD mean absolute value from %.2f changes to %.2f.", mean(abs(balance_ipw$smd_before)), mean(abs(balance_ipw$smd_after))),
  "- fixed cohort is substantially younger, healthier, higher baseline cognition, selection bias is present.",
  "- IPW adjusts observed baseline factors only, cannot adjust unmeasured selection, death and missing cognitive measurement.",
  "- T04–T06 rerun with identical age eligible definition; in this report 5,672 participants consistent with current longitudinal primary analysis.",
  "",
  "## six, output file",
  "",
  "- `attrition_pattern_2026-09-20.csv`",
  "- `attrition_retention_2026-09-20.csv`",
  "- `attrition_balance_2026-09-20.csv`",
  "- `attrition_balance_ipw_2026-09-20.csv`",
  "- `selection_model_2026-09-20.csv`",
  "- `ipw_diagnostics_2026-09-20.csv`",
  "- `ipw_longitudinal_models_2026-09-20.csv`",
  "- `ipw_cumulative_incidence_2026-09-20.csv`",
  "- `output/09_figures_and_tables/FigureS_attrition_selection_2026-09-20.png/.tif`",
  "",
  sprintf("*generation time: %s*", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
)
writeLines(report, file.path(OUT, "attrition_selection_2026-09-20.md"), useBytes = TRUE)

cat("\n========== T07 attrition and selection bias analysis complete ==========\n")
cat("baseline≥45 years old and ms10 analyzable: ", nrow(eligible_ms10), "; fixed cohort meeting age eligibility: ", nrow(fixed_base), "\n", sep = "")
cat("selection model complete cases: ", nrow(ps_data), "; fixed cohort: ", nrow(ps_fixed), "\n", sep = "")
cat("IPW mean|SMD|: ", sprintf("%.2f -> %.2f", mean(abs(balance_ipw$smd_before)), mean(abs(balance_ipw$smd_after))), "\n", sep = "")
cat("main figure: ", file.path(FIG_DIR, "FigureS_attrition_selection_2026-09-20.png"), "\n", sep = "")
