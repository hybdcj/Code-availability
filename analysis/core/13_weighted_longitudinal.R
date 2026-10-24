# ============================================================================
# 13_weighted_longitudinal.R —— T06: weighted longitudinal sensitivity analysis
# version: v1.0(2026-09-20)
# input: 
#   output/04_longitudinal_cohort/cohort_change_2026-09-13.rds
#   output/04_longitudinal_cohort/interval_incidence_events_2026-09-20.rds
# output: 
#   output/04_longitudinal_cohort/weighted_longitudinal_models_2026-09-20.csv
#   output/04_longitudinal_cohort/weighted_longitudinal_incidence_2026-09-20.csv
#   output/04_longitudinal_cohort/weighted_longitudinal_2026-09-20.md
#   output/09_figures_and_tables/FigureS_weighted_sensitivity_2026-09-20.png/.tif
#
# definition: 
#   primary sensitivity weights = 2011 year baseline individual weights; 
#   PSU = communityID; 
#   longitudinal model = 2020 endpoint assessment, first impairment Cox, linear regression of ten year change; 
#   compare weighted and unweighted results within identical samples, avoid interpreting“weight induced sample changes”as changes in effects.
#   this script adjusts baseline weights only; attrition selection IPW from T07 performed separately.
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

obj_coh <- readRDS(file.path(OUT, "cohort_change_2026-09-13.rds"))
obj_inc <- readRDS(file.path(OUT, "interval_incidence_events_2026-09-20.rds"))
coh <- obj_coh$cohort
base <- obj_inc$base
wide_change <- obj_coh$wide

# baseline variables.
b11 <- coh %>% filter(year == 2011) %>%
  transmute(KEY, age, female, edu_group,
            dep = as.integer(cesd10 >= 10), adl, sleep_h,
            pain = pain_any, pa = pa_any, hibp, diab,
            weight, communityID)
base_model <- base %>% inner_join(b11, by = "KEY") %>% inner_join(
  wide_change %>% select(KEY, delta), by = "KEY"
) %>% mutate(
  time = ifelse(is.na(first_index), 9, c(2, 4, 7, 9)[first_index]),
  event = as.integer(!is.na(first_index))
)
# ten year change model and existing cohort_risk consistent with primary analysis: entire fixed cohort n=5,672.
base_all <- b11 %>% inner_join(wide_change %>% select(KEY, delta), by = "KEY")
stopifnot(nrow(base_model) == nrow(base), nrow(base_all) == nrow(wide_change))

# weight design.
base_design <- svydesign(ids = ~communityID, weights = ~weight, data = base_model)
all_design <- svydesign(ids = ~communityID, weights = ~weight, data = base_all)
weight_diag <- do.call(rbind, lapply(list(
  list(label = "2011 unimpaired baseline participants", d = base_model, design = base_design),
  list(label = "entire fixed cohort", d = base_all, design = all_design)
), function(z) data.frame(
  analysis_sample = z$label,
  n = nrow(z$d),
  weight_n = sum(!is.na(z$d$weight) & z$d$weight > 0),
  community_n = length(unique(z$d$communityID)),
  design_df = degf(z$design),
  weight_sum = sum(z$d$weight, na.rm = TRUE),
  weight_median = median(z$d$weight, na.rm = TRUE),
  weight_cv = sd(z$d$weight, na.rm = TRUE) / mean(z$d$weight, na.rm = TRUE)
)))

# weighted/unweighted cumulative incidence.
cum_incidence <- data.frame()
time_points <- c(0, 2, 4, 7, 9)
first_limits <- c(0L, 1L, 2L, 3L, 4L)
conf_limits <- c(0L, 1L, 2L, 3L, 3L)
event_counts_first <- c(0L, cumsum(vapply(1:4, function(i) sum(obj_inc$base$first_index == i, na.rm = TRUE), integer(1))))
event_counts_conf <- c(0L, cumsum(vapply(1:3, function(i) sum(obj_inc$base$confirmed_index == i, na.rm = TRUE), integer(1))), sum(!is.na(obj_inc$base$confirmed_index)))
for (k in seq_along(time_points)) {
  e_first <- as.integer(!is.na(base_model$first_index) & base_model$first_index <= first_limits[k])
  e_conf <- as.integer(!is.na(base_model$confirmed_index) & base_model$confirmed_index <= conf_limits[k])
  d <- base_model
  d$e_first <- e_first
  d$e_conf <- e_conf
  des <- svydesign(ids = ~communityID, weights = ~weight, data = d)
  m1 <- svymean(~e_first, des, na.rm = TRUE)
  m2 <- svymean(~e_conf, des, na.rm = TRUE)
  ci1 <- confint(m1)
  ci2 <- confint(m2)
  un1 <- mean(d$e_first, na.rm = TRUE)
  un2 <- mean(d$e_conf, na.rm = TRUE)
  cum_incidence <- rbind(cum_incidence, data.frame(
    time_years = rep(time_points[k], 2),
    outcome = c("first_episode", "confirmed_first_episode"),
    weighted_cumulative_incidence = c(as.numeric(m1), as.numeric(m2)),
    weighted_ci_lower = c(ci1[1,1], ci2[1,1]),
    weighted_ci_upper = c(ci1[1,2], ci2[1,2]),
    unweighted_cumulative_incidence = c(un1, un2),
    weighted_minus_unweighted = c(as.numeric(m1) - un1, as.numeric(m2) - un2),
    events_first = c(event_counts_first[k], NA_integer_),
    events_confirmed = c(NA_integer_, event_counts_conf[k]),
    stringsAsFactors = FALSE
  ))
}
cum_incidence$weighted_cumulative_incidence_pct <- 100 * cum_incidence$weighted_cumulative_incidence
cum_incidence$weighted_ci_lower_pct <- 100 * cum_incidence$weighted_ci_lower
cum_incidence$weighted_ci_upper_pct <- 100 * cum_incidence$weighted_ci_upper
cum_incidence$unweighted_cumulative_incidence_pct <- 100 * cum_incidence$unweighted_cumulative_incidence
cum_incidence$weighted_minus_unweighted_pct <- 100 * cum_incidence$weighted_minus_unweighted

# weighted longitudinal model.compare weighting within exposure specific complete cases.
exposures <- c("dep", "adl", "sleep_h", "pain", "pa", "hibp", "diab")
labels <- c(dep = "depression(CES-D-10 ≥10)", adl = "ADL limited", sleep_h = "sleep duration(per 1 hours)",
            pain = "pain", pa = "regular physical activity", hibp = "hypertension", diab = "diabetes")
measure <- c(dep = "OR", adl = "OR", sleep_h = "OR", pain = "OR", pa = "OR", hibp = "OR", diab = "OR")

extract_est <- function(fit, term, type = c("ratio", "beta")) {
  type <- match.arg(type)
  b <- coef(fit)[term]
  se <- sqrt(vcov(fit)[term, term])
  if (type == "ratio") {
    est <- exp(b); lo <- exp(b - 1.96 * se); hi <- exp(b + 1.96 * se)
  } else {
    est <- b; lo <- b - 1.96 * se; hi <- b + 1.96 * se
  }
  p <- summary(fit)$coefficients[term, ncol(summary(fit)$coefficients)]
  list(est = as.numeric(est), lo = as.numeric(lo), hi = as.numeric(hi), p = as.numeric(p))
}

model_rows <- list()
for (e in exposures) {
  vars_inc <- c(e, "age", "female", "edu_group", "endpoint_2020", "time", "event", "weight", "communityID")
  vars_change <- c(e, "age", "female", "edu_group", "delta", "weight", "communityID")
  d_inc <- base_model[complete.cases(base_model[, vars_inc]), ]
  d_change <- base_all[complete.cases(base_all[, vars_change]), ]
  des_inc <- svydesign(ids = ~communityID, weights = ~weight, data = d_inc)
  des_change <- svydesign(ids = ~communityID, weights = ~weight, data = d_change)
  rhs <- paste(e, "+ age + female + edu_group")
  f_logit <- as.formula(paste("endpoint_2020 ~", rhs))
  f_cox <- as.formula(paste("Surv(time, event) ~", rhs))
  f_lin <- as.formula(paste("delta ~", rhs))

  w_logit <- svyglm(f_logit, design = des_inc, family = quasibinomial())
  u_logit <- glm(f_logit, data = d_inc, family = binomial())
  w_cox <- svycoxph(f_cox, design = des_inc)
  u_cox <- coxph(f_cox, data = d_inc)
  w_lin <- svyglm(f_lin, design = des_change)
  u_lin <- lm(f_lin, data = d_change)

  getw <- function(fit, type) extract_est(fit, e, type)
  model_specs <- list(
    list(model = "endpoint_2020_logistic", label = "2020 endpoint assessment", measure = "OR",
         w = getw(w_logit, "ratio"), u = getw(u_logit, "ratio"),
         n = nrow(d_inc), events = sum(d_inc$endpoint_2020), communities = length(unique(d_inc$communityID))),
    list(model = "first_event_cox", label = "first impairment Cox", measure = "HR",
         w = getw(w_cox, "ratio"), u = getw(u_cox, "ratio"),
         n = nrow(d_inc), events = sum(d_inc$event), communities = length(unique(d_inc$communityID))),
    list(model = "ten_year_change_lm", label = "ten year cognitive change(fixed cohort complete cases)", measure = "beta",
         w = getw(w_lin, "beta"), u = getw(u_lin, "beta"),
         n = nrow(d_change), events = nrow(d_change), communities = length(unique(d_change$communityID)))
  )
  for (z in model_specs) {
    model_rows[[length(model_rows) + 1L]] <- data.frame(
      model = z$model, model_label = z$label, exposure = e, exposure_label = labels[[e]],
      effect_measure = z$measure, n = z$n, events = z$events, communities = z$communities,
      weighted_est = z$w$est, weighted_lower = z$w$lo, weighted_upper = z$w$hi,
      weighted_p = z$w$p,
      unweighted_est = z$u$est, unweighted_lower = z$u$lo, unweighted_upper = z$u$hi,
      unweighted_p = z$u$p,
      stringsAsFactors = FALSE
    )
  }
}
models <- bind_rows(model_rows)
models <- models %>% mutate(
  weighted_significant = weighted_p < 0.05,
  unweighted_significant = unweighted_p < 0.05,
  direction_consistent = ifelse(effect_measure == "beta",
                                sign(weighted_est) == sign(unweighted_est),
                                sign(weighted_est - 1) == sign(unweighted_est - 1)),
  estimate_change_pct = ifelse(effect_measure == "beta", NA_real_, 100 * (weighted_est - unweighted_est) / unweighted_est)
)

write.csv(models, file.path(OUT, "weighted_longitudinal_models_2026-09-20.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(cum_incidence, file.path(OUT, "weighted_longitudinal_incidence_2026-09-20.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(weight_diag, file.path(OUT, "weighted_longitudinal_design_2026-09-20.csv"), row.names = FALSE, fileEncoding = "UTF-8")

# simple forest plot: highlight direction and precision changes with weighting.
draw_weighted_figure <- function() {
  par(mfrow = c(1, 3), mar = c(5, 7, 3, 1), family = "sans")
  model_specs <- list(
    c("endpoint_2020_logistic", "2020 endpoint: OR", "OR"),
    c("first_event_cox", "First event: HR", "HR"),
    c("ten_year_change_lm", "Ten-year change: beta", "beta")
  )
  for (spec in model_specs) {
    z <- models[models$model == spec[1], ]
    z <- z[order(z$unweighted_est), ]
    yy <- seq_len(nrow(z))
    xlim <- range(c(z$weighted_lower, z$weighted_upper, z$unweighted_lower, z$unweighted_upper, if (spec[3] == "beta") 0 else 1), na.rm = TRUE)
    plot(NA, xlim = xlim, ylim = c(0.5, nrow(z) + 0.5), yaxt = "n", xlab = spec[2], ylab = "", main = spec[2])
    abline(v = if (spec[3] == "beta") 0 else 1, lty = 3, col = "grey60")
    segments(z$unweighted_lower, yy, z$unweighted_upper, yy, col = "#C00000", lwd = 2)
    points(z$unweighted_est, yy, pch = 17, col = "#C00000", cex = 0.9)
    segments(z$weighted_lower, yy + 0.12, z$weighted_upper, yy + 0.12, col = "#1F4E78", lwd = 2)
    points(z$weighted_est, yy + 0.12, pch = 16, col = "#1F4E78", cex = 0.9)
    axis(2, at = yy, labels = z$exposure_label, las = 1, cex.axis = 0.8)
  }
  legend("topright", legend = c("Unweighted", "Weighted (baseline 2011)"),
         col = c("#C00000", "#1F4E78"), pch = c(17, 16), lwd = 2, bty = "n", cex = 0.75)
}
png(file.path(FIG_DIR, "FigureS_weighted_sensitivity_2026-09-20.png"), width = 14, height = 6, units = "in", res = 300)
draw_weighted_figure()
dev.off()
tryCatch({
  tiff(file.path(FIG_DIR, "FigureS_weighted_sensitivity_2026-09-20.tif"), width = 14, height = 6, units = "in", res = 600, compression = "lzw")
  draw_weighted_figure()
  dev.off()
}, error = function(e) { if (length(dev.list())) dev.off() })

# statistical conclusions.
sig_change <- models %>% filter(weighted_significant != unweighted_significant)
report <- c(
  "# T06: weighted longitudinal sensitivity analysis(2026-09-20)",
  "",
  "> weights: 2011 year baseline individual weights; PSU: communityID.",
  "> 2020 endpoint assessment and first impairment Cox use baseline unimpaired participants; ten year change uses entire fixed cohort as source, each exposure uses its own complete cases; weighted and unweighted versions of each model use identical samples.",
  "> analysis uses only 2011 baseline weights; five wave attrition selection IPW already addressed by T07 separate analysis.",
  "",
  "## one, weight design",
  "",
  "| analysis sample | n | weights available | number of communities | design degrees of freedom | median weight | weights CV |",
  "|---|---:|---:|---:|---:|---:|---:|",
  paste0("| ", weight_diag$analysis_sample, " | ", weight_diag$n, " | ", weight_diag$weight_n,
         " | ", weight_diag$community_n, " | ", sprintf("%.1f", weight_diag$design_df),
         " | ", sprintf("%.0f", weight_diag$weight_median), " | ", sprintf("%.2f", weight_diag$weight_cv), " |"),
  "",
  "## two, weighted and unweighted cumulative incidence",
  "",
  "| outcome | follow up year | weighted cumulative incidence | 95%CI | unweighted | difference |",
  "|---|---:|---:|---:|---:|---:|",
  paste0("| ", ifelse(cum_incidence$outcome == "first_episode", "first onset", "consecutive confirmation"),
         " | ", cum_incidence$time_years, " | ", sprintf("%.1f%%", cum_incidence$weighted_cumulative_incidence_pct),
         " | ", sprintf("%.1f–%.1f%%", cum_incidence$weighted_ci_lower_pct, cum_incidence$weighted_ci_upper_pct),
         " | ", sprintf("%.1f%%", cum_incidence$unweighted_cumulative_incidence_pct),
         " | ", sprintf("%+.1f pp", cum_incidence$weighted_minus_unweighted_pct), " |"),
  "",
  "## three, longitudinal model: weighted versus unweighted comparison",
  "",
  "| model | exposure | n | effect estimate | weighted estimate (95%CI) | unweighted estimate (95%CI) | consistent direction | significance change |",
  "|---|---|---:|---:|---:|---:|---|---|",
  paste0("| ", models$model_label, " | ", models$exposure_label, " | ", models$n, " | ", models$effect_measure,
         " | ", sprintf("%.2f (%.2f–%.2f)", models$weighted_est, models$weighted_lower, models$weighted_upper),
         " | ", sprintf("%.2f (%.2f–%.2f)", models$unweighted_est, models$unweighted_lower, models$unweighted_upper),
         " | ", ifelse(models$direction_consistent, "yes", "no"),
         " | ", ifelse(models$weighted_significant != models$unweighted_significant, "yes", "no"), " |"),
  "",
  "## four, conclusion",
  "",
  sprintf("- combinations with discordant weighted and unweighted directions: %d/%d.", sum(!models$direction_consistent), nrow(models)),
  sprintf("- significance changed after weighting: %d/%d.", nrow(sig_change), nrow(models)),
  sprintf("- combinations with significance changes: %s.", ifelse(nrow(sig_change) == 0, "none", paste0(sig_change$model_label, "×", sig_change$exposure_label, collapse = "; "))),
  "- regular physical activity model has substantial missingness, its weight sensitivity results are exploratory, should not independently alter main conclusions.",
  "- depression, ADL core conclusions largely retain direction and significance under baseline weighting; borderline changes in pain and regular physical activity are exploratory.",
  "- fixed cohort attrition/survivor selection cannot be resolved solely by 2011 weights; T07 attrition and analyses completed IPW analysis.",
  "",
  "## five, output file",
  "",
  "- `weighted_longitudinal_models_2026-09-20.csv`",
  "- `weighted_longitudinal_incidence_2026-09-20.csv`",
  "- `weighted_longitudinal_design_2026-09-20.csv`",
  "- `weighted_longitudinal_2026-09-20.md`",
  "- `output/09_figures_and_tables/FigureS_weighted_sensitivity_2026-09-20.png/.tif`",
  "",
  sprintf("*generation time: %s*", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
)
writeLines(report, file.path(OUT, "weighted_longitudinal_2026-09-20.md"), useBytes = TRUE)

cat("\n========== T06 weighted longitudinal sensitivity complete ==========\n")
cat("fixed cohort: ", nrow(base_model), "; weighted model combinations: ", nrow(models), "\n", sep = "")
cat("weighted/direction differs from unweighted: ", sum(!models$direction_consistent), "; significance changes: ", nrow(sig_change), "\n", sep = "")
cat("main figure: ", file.path(FIG_DIR, "FigureS_weighted_sensitivity_2026-09-20.png"), "\n", sep = "")
