# ============================================================================
# 12_cumulative_incidence.R —— T05: fixed cohort cumulative incidence curves and numbers at risk
# version: v1.0(2026-09-20)
# input: output/04_longitudinal_cohort/interval_incidence_events_2026-09-20.rds
# output: 
#   output/04_longitudinal_cohort/cumulative_incidence_2026-09-20.csv
#   output/04_longitudinal_cohort/cumulative_incidence_risk_table_2026-09-20.csv
#   output/04_longitudinal_cohort/cumulative_incidence_km_sensitivity_2026-09-20.csv
#   output/04_longitudinal_cohort/cumulative_incidence_2026-09-20.md
#   output/09_figures_and_tables/Figure4_cumulative_incidence_2026-09-20.png
#   output/09_figures_and_tables/Figure4_cumulative_incidence_2026-09-20.tif
#
# notes: 
#   primary curve: at scheduled follow up times 2, 4, 7, 9 year discrete cumulative incidence; 
#   sensitivity curve: impairment confirmed at two consecutive waves; 
#   KM: unknown event time, use interval midpoint approximation 1, 3, 5.5, 8 year, sensitivity analysis only.
# ============================================================================

invisible(try(Sys.setlocale("LC_ALL", "zh_CN.UTF-8"), silent = TRUE))
if (dir.exists("D:/R/library")) .libPaths(c("D:/R/library", .libPaths()))
suppressPackageStartupMessages(library(survival))

P <- Sys.getenv("CHARLS_PROJ", "D:/BaiduSyncdisk/research/Xu_Zhonghao")
OUT <- file.path(P, "output", "04_longitudinal_cohort")
FIG_DIR <- file.path(P, "output", "09_figures_and_tables")
if (!dir.exists(OUT)) dir.create(OUT, recursive = TRUE)
if (!dir.exists(FIG_DIR)) dir.create(FIG_DIR, recursive = TRUE)

obj <- readRDS(file.path(OUT, "interval_incidence_events_2026-09-20.rds"))
base <- obj$base
n_base <- nrow(base)
stopifnot(n_base > 0, all(c("first_index", "confirmed_index", "endpoint_2020") %in% names(base)))

times <- c(0, 2, 4, 7, 9)
interval_len <- c(2, 2, 3, 2)
midpoints <- c(1, 3, 5.5, 8)

binomial_ci <- function(events, n) {
  if (n <= 0) return(c(NA_real_, NA_real_))
  as.numeric(binom.test(events, n)$conf.int)
}
fmt_pct <- function(x) ifelse(is.na(x), "—", sprintf("%.1f%%", 100 * x))

# primary curve: first onset.
events_first <- vapply(1:4, function(i) sum(base$first_index == i, na.rm = TRUE), integer(1))
cum_first <- c(0L, cumsum(events_first))
ci_first <- t(vapply(cum_first, function(x) binomial_ci(x, n_base), numeric(2)))
risk_first <- n_base - cum_first
# t=9 number at risk is participants without first event before final follow up; thereafter study ends for all.
risk_first[5] <- n_base - cum_first[5]

# sensitivity curve: confirmed at two consecutive waves.2018-2020 no next wave, cannot confirm new events.
events_conf <- vapply(1:3, function(i) sum(base$confirmed_index == i, na.rm = TRUE), integer(1))
cum_conf <- c(0L, cumsum(events_conf), sum(events_conf))
ci_conf <- t(vapply(cum_conf, function(x) binomial_ci(x, n_base), numeric(2)))
risk_conf <- c(n_base, n_base - cumsum(events_conf), NA_integer_)
risk_conf_status <- c("estimable", "estimable", "estimable", "estimable", "unconfirmed")

curve_tab <- data.frame(
  time_years = rep(times, 2),
  outcome = rep(c("first_episode", "confirmed_first_episode"), each = length(times)),
  outcome_label = rep(c("first impairment within nine years", "impairment confirmed at two consecutive waves"), each = length(times)),
  cumulative_events = c(cum_first, cum_conf),
  cumulative_incidence = c(cum_first / n_base, cum_conf / n_base),
  ci_lower = c(ci_first[, 1], ci_conf[, 1]),
  ci_upper = c(ci_first[, 2], ci_conf[, 2]),
  event_free_at_time = c(risk_first, risk_conf),
  risk_status = c(rep("estimable", 5), risk_conf_status),
  stringsAsFactors = FALSE
)
curve_tab$cumulative_incidence_pct <- 100 * curve_tab$cumulative_incidence
curve_tab$ci_lower_pct <- 100 * curve_tab$ci_lower
curve_tab$ci_upper_pct <- 100 * curve_tab$ci_upper

# number at risk table: retain only two rows needed for manuscript figure.
risk_tab <- data.frame(
  time_years = times,
  at_risk_first_episode = risk_first,
  cumulative_events_first_episode = cum_first,
  at_risk_confirmed = risk_conf,
  cumulative_events_confirmed = cum_conf,
  confirmed_risk_status = risk_conf_status,
  stringsAsFactors = FALSE
)

# KM sensitivity: place event time at interval midpoint.
km_time <- ifelse(is.na(base$first_index), 9, midpoints[base$first_index])
km_event <- as.integer(!is.na(base$first_index))
km_fit <- survfit(Surv(km_time, km_event) ~ 1)
km_s <- summary(km_fit, times = c(2, 4, 7, 9), extend = TRUE)
km_tab <- data.frame(
  time_years = km_s$time,
  n_risk = km_s$n.risk,
  n_event = km_s$n.event,
  cumulative_incidence_km = 1 - km_s$surv,
  ci_lower = 1 - km_s$upper,
  ci_upper = 1 - km_s$lower,
  stringsAsFactors = FALSE
)
km_tab$cumulative_incidence_pct <- 100 * km_tab$cumulative_incidence_km
km_tab$ci_lower_pct <- 100 * km_tab$ci_lower
km_tab$ci_upper_pct <- 100 * km_tab$ci_upper

# output data.
write.csv(curve_tab, file.path(OUT, "cumulative_incidence_2026-09-20.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(risk_tab, file.path(OUT, "cumulative_incidence_risk_table_2026-09-20.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(km_tab, file.path(OUT, "cumulative_incidence_km_sensitivity_2026-09-20.csv"), row.names = FALSE, fileEncoding = "UTF-8")
saveRDS(list(curve = curve_tab, risk = risk_tab, km = km_tab, fit = km_fit),
        file.path(OUT, "cumulative_incidence_2026-09-20.rds"))

# main figure: English labels, for subsequent submission.risk table place in x below axis.
draw_figure <- function() {
  layout(matrix(c(1, 2), ncol = 1), heights = c(4.1, 1.7))
  # Upper panel: cumulative incidence curves.
  par(mar = c(2.8, 8.0, 3.2, 1.2), family = "sans")
  ylim <- c(0, 50)
  plot(times, 100 * cum_first / n_base, type = "n", xaxt = "n", yaxt = "n",
       xlim = c(0, 9.4), ylim = ylim,
       xlab = "", ylab = "Cumulative incidence (%)",
       main = "Cumulative Incidence of Cognitive Impairment (Baseline Normal)")
  axis(1, at = times, labels = FALSE, tcl = -0.25)
  axis(2, at = seq(0, 50, 10), las = 1)
  grid(col = "grey88", lty = 3)
  # First episode: 95% CI and line.
  segments(times, 100 * ci_first[, 1], times, 100 * ci_first[, 2], col = "#1F4E78", lwd = 1.5)
  lines(times, 100 * cum_first / n_base, type = "s", lwd = 2, col = "#1F4E78")
  points(times, 100 * cum_first / n_base, pch = 16, cex = 0.9, col = "#1F4E78")
  # Confirmed episode: last time is carry-forward, not a new confirmable event.
  ctimes <- c(0, 2, 4, 7)
  segments(ctimes, 100 * ci_conf[1:4, 1], ctimes, 100 * ci_conf[1:4, 2], col = "#C00000", lwd = 1.5)
  lines(ctimes, 100 * cum_conf[1:4] / n_base, type = "s", lwd = 2, col = "#C00000")
  points(ctimes, 100 * cum_conf[1:4] / n_base, pch = 17, cex = 0.9, col = "#C00000")
  points(9, 100 * cum_conf[5] / n_base, pch = 1, cex = 1.35, col = "#C00000", lwd = 2)
  legend("topleft",
         legend = c("First impairment episode (any wave)",
                    "Confirmed in two consecutive waves",
                    "Not confirmable at 9 years"),
         col = c("#1F4E78", "#C00000", "#C00000"),
         pch = c(16, 17, 1), lwd = c(2, 2, NA), bty = "n", cex = 0.82)
  # Lower panel: risk set and cumulative events table.
  par(mar = c(1.0, 8.0, 0.2, 1.2), family = "sans")
  plot(NA, xlim = c(0, 9.4), ylim = c(0, 5.6), axes = FALSE, xlab = "", ylab = "")
  mtext("No. at risk", side = 2, line = 5.2, at = 5.35, las = 1, adj = 0, font = 2, cex = 0.82)
  mtext("First episode", side = 2, line = 5.2, at = 4.15, las = 1, adj = 0, col = "#1F4E78", font = 2, cex = 0.80)
  mtext("Confirmed", side = 2, line = 5.2, at = 3.00, las = 1, adj = 0, col = "#C00000", font = 2, cex = 0.80)
  mtext("First-event count", side = 2, line = 5.2, at = 1.85, las = 1, adj = 0, font = 2, cex = 0.80)
  text(times, rep(4.15, length(times)), risk_first, cex = 0.78, col = "#1F4E78")
  conf_risk_display <- c(as.character(risk_conf[1:4]), "NE")
  text(times, rep(3.00, length(times)), conf_risk_display, cex = 0.78, col = "#C00000")
  text(times, rep(1.85, length(times)), cum_first, cex = 0.78, col = "black")
  text(times, rep(0.25, length(times)), times, cex = 0.82, font = 2)
  mtext("Follow-up (years)", side = 1, line = 1.85, cex = 0.82)
}

png(file.path(FIG_DIR, "Figure4_cumulative_incidence_2026-09-20.png"), width = 10, height = 7.2, units = "in", res = 300)
draw_figure()
dev.off()

tiff_ok <- TRUE
tryCatch({
  tiff(file.path(FIG_DIR, "Figure4_cumulative_incidence_2026-09-20.tif"), width = 10, height = 7.2, units = "in", res = 600, compression = "lzw")
  draw_figure()
  dev.off()
}, error = function(e) {
  tiff_ok <<- FALSE
  if (length(dev.list())) dev.off()
})

# Chinese report.
fmt_int <- function(x) ifelse(is.na(x), "—", format(as.integer(x), big.mark = ",", scientific = FALSE, trim = TRUE))
report <- c(
  "# T05: fixed cohort cumulative incidence curves(2026-09-20)",
  "",
  "> input: T04 of `interval_incidence_events_2026-09-20.rds`",
  sprintf("> primary outcome: 2011 unimpaired baseline participants(n=%s)first occurrence `impaired_fixed=1`", format(n_base, big.mark = ",")),
  "> **note: ** this curve represents 2011 year≥45 age eligible and all five waves ms10 conditional cumulative incidence in fully scorable fixed cohort, cannot extrapolate to entire CHARLS population incidence.",
  "> primary curve at scheduled follow up times 2, 4, 7, 9 year display of discrete cumulative incidence; number at risk means event free participants at that time.KM use interval midpoint event times for consistency check.",
  "",
  "## one, primary curve data",
  "",
  "| year | follow up year | cumulative events | cumulative incidence | 95%CI | number still event free at time point |",
  "|---|---:|---:|---:|---:|---:|",
  paste0("| ", c("2011","2013","2015","2018","2020"), " | ", times, " | ", cum_first, " | ",
         fmt_pct(cum_first / n_base), " | ",
         ifelse(is.na(ci_first[,1]), "—", sprintf("%.1f–%.1f%%", 100*ci_first[,1], 100*ci_first[,2])),
         " | ", fmt_int(risk_first), " |"),
  "",
  "## two, cumulative incidence confirmed at two consecutive waves(sensitivity)",
  "",
  "| year | follow up year | cumulative confirmed events | cumulative incidence | 95%CI | number still unconfirmed at time point |",
  "|---|---:|---:|---:|---:|---:|",
  paste0("| ", c("2011","2013","2015","2018","2020"), " | ", times, " | ", cum_conf, " | ",
         fmt_pct(cum_conf / n_base), " | ",
         ifelse(is.na(ci_conf[,1]), "—", sprintf("%.1f–%.1f%%", 100*ci_conf[,1], 100*ci_conf[,2])),
         " | ", ifelse(is.na(risk_conf), "unconfirmed", fmt_int(risk_conf)), " |"),
  "",
  "> 2018–2020 interval cannot confirm new two wave events, because requires 2022 next wave data in year; time point marked as“unconfirmed”, is not 0 events.",
  "> confirmed event time is first impaired wave of consecutive pair, confirmation depends on subsequent wave.",
  "",
  "## three, KM consistency check(fixed cohort has no censoring)",
  "",
  "KM use interval midpoints 1, 3, 5.5, 8 years as event times.because this fixed cohort within 9 years has no interim censoring, KM estimate agrees with discrete cumulative incidence counts; table checks computational consistency, not an independent sensitivity finding.",
  "",
  "| follow up year | number at risk | KM cumulative incidence | 95%CI |",
  "|---:|---:|---:|---:|",
  paste0("| ", km_tab$time_years, " | ", km_tab$n_risk, " | ", fmt_pct(km_tab$cumulative_incidence_km), " | ",
         sprintf("%.1f–%.1f%%", km_tab$ci_lower_pct, km_tab$ci_upper_pct), " |"),
  "",
  "## four, conclusion",
  "",
  sprintf("- cumulative incidence of first onset: 2 year %s, 4 year %s, 7 year %s, 9 year %s.",
          fmt_pct(cum_first[2]/n_base), fmt_pct(cum_first[3]/n_base),
          fmt_pct(cum_first[4]/n_base), fmt_pct(cum_first[5]/n_base)),
  sprintf("- cumulative incidence confirmed at two consecutive waves: 2 year %s, 4 year %s, 7 year %s; 9 year cannot add confirmed events.",
          fmt_pct(cum_conf[2]/n_base), fmt_pct(cum_conf[3]/n_base), fmt_pct(cum_conf[4]/n_base)),
  "- main figure suitable for Figure 4 cumulative incidence panel, but must also show numbers at risk.",
  sprintf("- single threshold crossing proportion is high, consecutively confirmed curve better demonstrates robustness; do not describe first onset in manuscript as %s persistent cognitive impairment rate.", fmt_pct(cum_first[5] / n_base)),
  "",
  "## five, output file",
  "",
  "- `cumulative_incidence_2026-09-20.csv`",
  "- `cumulative_incidence_risk_table_2026-09-20.csv`",
  "- `cumulative_incidence_km_sensitivity_2026-09-20.csv`",
  "- `cumulative_incidence_2026-09-20.rds`",
  "- `output/09_figures_and_tables/Figure4_cumulative_incidence_2026-09-20.png`",
  "- `output/09_figures_and_tables/Figure4_cumulative_incidence_2026-09-20.tif`",
  "",
  sprintf("*generation time: %s*", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
)
writeLines(report, file.path(OUT, "cumulative_incidence_2026-09-20.md"), useBytes = TRUE)

cat("\n========== T05 cumulative incidence curves complete ==========\n")
cat("cumulative incidence of first onset: ", paste(fmt_pct(cum_first / n_base), collapse = " | "), "\n", sep = "")
cat("consecutively confirmed cumulative incidence: ", paste(fmt_pct(cum_conf / n_base), collapse = " | "), "\n", sep = "")
cat("main figure: ", file.path(FIG_DIR, "Figure4_cumulative_incidence_2026-09-20.png"), "\n", sep = "")
cat("TIFF: ", ifelse(tiff_ok, file.path(FIG_DIR, "Figure4_cumulative_incidence_2026-09-20.tif"), "generation failed"), "\n", sep = "")
