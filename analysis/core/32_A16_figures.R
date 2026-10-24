# ============================================================================
# 32_A16_figures.R —— A16: Figure 1–5 Chinese and English submission versions(PNG 300 dpi + TIF 600 dpi)
# input: A07/A08/A09/A12/A13/A14/A15 frozen results
# output: output/14_A16_figures_and_tables/{Chinese,English}/Figure*_ZH/EN.png/.tif
# generate: 2026-09-27
# ============================================================================
invisible(try(Sys.setlocale("LC_ALL", "zh_CN.UTF-8"), silent = TRUE))
.libPaths(c("D:/R/library", .libPaths()))
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(pROC); library(readr) })
P <- "D:/BaiduSyncdisk/research/Xu_Zhonghao"
OUT <- file.path(P, "output", "14_A16_figures_and_tables")
ZH <- file.path(OUT, "Chinese"); EN <- file.path(OUT, "English"); DATA <- file.path(OUT, "data")
for (d in c(OUT, ZH, EN, DATA)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

save_both <- function(draw_fun, basename, lang) {
  dir <- if (lang == "ZH") ZH else EN
  for (ext in c("png", "tif")) {
    path <- file.path(dir, paste0(basename, ".", ext))
    ok <- FALSE
    for (i in seq_len(10)) {
      ok <- tryCatch({
        if (ext == "png") {
          png(path, width = 12, height = 7, units = "in", res = 300)
        } else {
          tiff(path, width = 12, height = 7, units = "in", res = 600, compression = "lzw")
        }
        draw_fun(lang)
        dev.off()
        TRUE
      }, error = function(e) {
        if (length(dev.list())) dev.off()
        FALSE
      })
      if (ok) break
      Sys.sleep(1)
    }
    if (!ok) stop("cannot write figure(retry 10 times): ", path)
  }
}
safe_write_csv <- function(x, file, ...) {
  for (i in seq_len(10)) {
    ok <- tryCatch({ utils::write.csv(x, file, ...); TRUE }, error = function(e) FALSE)
    if (ok) return(invisible(TRUE))
    Sys.sleep(1)
  }
  stop("cannot write file(retry 10 times): ", file)
}
parse_leading_num <- function(x) suppressWarnings(as.numeric(sub(" .*$", "", trimws(x))))
parse_est <- function(x) {
  x <- sub(", p=.*$", "", trimws(x))
  m <- regmatches(x, regexec("^([-0-9.]+) \\(([-0-9.]+)-([-0-9.]+)\\)$", x))[[1]]
  if (length(m) == 4) return(as.numeric(m[2:4]))
  c(NA_real_, NA_real_, NA_real_)
}

# ---------------------------------------------------------------------------
# Figure 1: study flowchart
# ---------------------------------------------------------------------------
draw_flow <- function(lang) {
  zh <- lang == "ZH"
  par(mar = c(1, 1, 1, 1), xpd = NA, family = "sans")
  plot(NA, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE, xlab = "", ylab = "")
  box <- function(x, y, w, h, label, col = "white", border = "#2C3E50", lwd = 1.5, cex = 0.85) {
    rect(x - w/2, y - h/2, x + w/2, y + h/2, col = col, border = border, lwd = lwd)
    text(x, y, label, cex = cex, font = 1)
  }
  ar <- function(x1, y1, x2, y2, label = NULL, col = "#2C3E50") {
    arrows(x1, y1, x2, y2, length = 0.08, lwd = 1.5, col = col)
    if (!is.null(label)) text((x1+x2)/2 + 0.01, (y1+y2)/2, label, cex = 0.72, col = col, srt = 0)
  }
  title_zh <- "CHARLS development and external validation flow"
  title_en <- "Development and external validation flow"
  text(0.5, 0.98, if (zh) title_zh else title_en, cex = 1.3, font = 2)
  # development column
  text(0.25, 0.925, if (zh) "CHARLS development cohort" else "CHARLS development cohort", cex = 1.0, font = 2)
  boxes <- list(
    c(0.25, 0.845, if (zh) "2011 age in year≥45 years old\nn=16,831" else "2011 age ≥45 years\nn=16,831"),
    c(0.25, 0.725, if (zh) "2011 year memory score measurable\nn=13,444" else "2011 memory measured\nn=13,444"),
    c(0.25, 0.605, if (zh) "unimpaired at baseline(memory score>4)\nn=10,422" else "Baseline unimpaired (memory >4)\nn=10,422"),
    c(0.25, 0.485, if (zh) "2015 year memory score measurable\nn=8,371" else "2015 memory measured\nn=8,371"),
    c(0.25, 0.365, if (zh) "complete four variable data\nn=8,296" else "Four variables complete\nn=8,296"),
    c(0.25, 0.225, if (zh) "final development sample(five waves ms10 fixed cohort)\nn=4,441; events=658(14.8%)" else "Final development sample (five-wave ms10 fixed cohort)\nn=4,441; events=658 (14.8%)")
  )
  for (i in seq_along(boxes)) {
    b <- boxes[[i]]
    box(as.numeric(b[1]), as.numeric(b[2]), 0.42, 0.095,
        b[3], col = if (i == length(boxes)) "#D6EAF8" else "white", cex = 0.82)
    if (i < length(boxes)) ar(as.numeric(b[1]), as.numeric(b[2]) - 0.048,
                              as.numeric(boxes[[i+1]][1]), as.numeric(boxes[[i+1]][2]) + 0.048)
  }
  # external validation
  text(0.76, 0.925, if (zh) "external validation cohorts(S0)" else "External validation cohorts (S0)", cex = 1.0, font = 2)
  ext <- list(
    c(0.76, 0.82, if (zh) "HRS w10→w12\nn=14,051; events=1,842(13.1%)" else "HRS w10→w12\nn=14,051; events=1,842 (13.1%)"),
    c(0.76, 0.68, if (zh) "ELSA w4→w6\nn=6,725; events=554(8.2%)" else "ELSA w4→w6\nn=6,725; events=554 (8.2%)"),
    c(0.76, 0.54, if (zh) "SHARE w4→w6\nn=26,879; events=2,930(10.9%)" else "SHARE w4→w6\nn=26,879; events=2,930 (10.9%)")
  )
  for (b in ext) box(as.numeric(b[1]), as.numeric(b[2]), 0.42, 0.105, b[3], col = "#FDEBD0", cex = 0.82)
  ar(0.47, 0.225, 0.55, 0.70, if (zh) "frozen A14 model\n without refitting" else "Frozen A14 model\nno refitting")
  box(0.76, 0.36, 0.42, 0.13,
      if (zh) "S0: full baseline sample in each cohort mean−1SD\nS1/S2: sensitivity analysis\n external AUC 0.692–0.763" else
        "S0: baseline-all mean−1SD\nS1/S2: sensitivity analyses\nExternal AUC 0.692–0.763",
      col = "#D5F5E3", cex = 0.80)
  text(0.25, 0.135, if (zh) "fixed five waves ms10 excluded by completeness restriction 3,855 participants" else "3,855 excluded by five-wave ms10 completeness", cex = 0.72, col = "#7F8C8D")
  text(0.5, 0.055,
       if (zh) "note: CHARLS outcome is 2011→2015 memory score threshold crossing(≤4); external cohorts use 4 year follow up across cohort specific thresholds.ms10(orientation + serial subtraction of 7)used only to define fixed five wave cohort." else
         "Note: CHARLS outcome is 2011→2015 memory cross-threshold (≤4); external cohorts use 4-year follow-up and cohort-specific thresholds. ms10 (orientation + serial 7s) was used only to define the five-wave fixed cohort.",
       cex = 0.66, adj = 0.5)
}
save_both(draw_flow, "Figure1_study_flow_ZH", "ZH")
save_both(draw_flow, "Figure1_study_flow_EN", "EN")

# ---------------------------------------------------------------------------
# Figure 2: candidate factor stability
# ---------------------------------------------------------------------------
orm <- read.csv(file.path(P, "output/03_five_wave_cross_sectional/separate_wave_analyses_Logistic_OR_matrix_2026_09_19.csv"),
                check.names = FALSE, fileEncoding = "UTF-8")
risk <- read.csv(file.path(P, "output/04_longitudinal_cohort/cohort_risk_2026-09-13.csv"),
                 check.names = FALSE, fileEncoding = "UTF-8")
iadl2 <- read.csv(file.path(P, "output/04_longitudinal_cohort/iadl2_longitudinal_models_2026-09-26.csv"),
                  check.names = FALSE, fileEncoding = "UTF-8")
or_stats <- function(pattern) {
  idx <- grep(pattern, orm[[1]])
  if (length(idx) == 0) return(c(NA, NA, NA))
  vals <- unlist(orm[idx[1], 3:7], use.names = FALSE)
  nums <- parse_leading_num(vals); nums <- nums[!is.na(nums)]
  if (length(nums) == 0) return(c(NA, NA, NA))
  c(min(nums), median(nums), max(nums))
}
risk_est <- function(rowname) {
  idx <- which(risk[[1]] == rowname)
  if (length(idx) == 0) return(c(NA, NA, NA))
  parse_est(risk[idx[1], 2])
}
risk_hr <- function(rowname) {
  idx <- which(risk[[1]] == rowname)
  if (length(idx) == 0) return(c(NA, NA, NA))
  parse_est(risk[idx[1], 3])
}
iadl_est <- function(model) {
  x <- iadl2$unweighted[iadl2$exposure == "iadl2" & iadl2$model == model]
  if (length(x) == 0) return(c(NA, NA, NA))
  parse_est(x[1])
}
plot_candidate_stability <- function(lang) {
  zh <- lang == "ZH"
  par(mfrow = c(1, 2), mar = c(4.5, 6, 3, 2), family = "sans")
  rows <- list(
    list(name_zh = "depression", name_en = "Depression", role = "Core", or = or_stats("depression"), lr = risk_est("depression(CES-D-10 ≥10)"), hr = risk_hr("depression(CES-D-10 ≥10)")),
    list(name_zh = "ADL limited", name_en = "ADL limitation", role = "Core", or = or_stats("ADL limited"), lr = risk_est("ADL limited"), hr = risk_hr("ADL limited")),
    list(name_zh = "IADL limited", name_en = "IADL limitation", role = "Supplementary", or = or_stats("IADL limited"), lr = iadl_est("logistic"), hr = iadl_est("cox")),
    list(name_zh = "sleep <6 h", name_en = "Sleep <6 h", role = "Incremental", or = or_stats("sleep <6"), lr = risk_est("sleep duration(per 1 hours)"), hr = risk_hr("sleep duration(per 1 hours)")),
    list(name_zh = "sleep >8 h", name_en = "Sleep >8 h", role = "Incremental", or = or_stats("sleep >8"), lr = c(NA, NA, NA), hr = c(NA, NA, NA)),
    list(name_zh = "regular physical activity", name_en = "Physical activity", role = "Incremental", or = or_stats("regular physical activity"), lr = risk_est("regular physical activity"), hr = risk_hr("regular physical activity")),
    list(name_zh = "hypertension", name_en = "Hypertension", role = "Background", or = or_stats("hypertension"), lr = risk_est("hypertension"), hr = risk_hr("hypertension")),
    list(name_zh = "diabetes", name_en = "Diabetes", role = "Incremental", or = or_stats("diabetes"), lr = risk_est("diabetes"), hr = risk_hr("diabetes")),
    list(name_zh = "stroke", name_en = "Stroke", role = "Incremental", or = or_stats("stroke"), lr = c(NA, NA, NA), hr = c(NA, NA, NA)),
    list(name_zh = "pain", name_en = "Pain", role = "Excluded", or = or_stats("bodily pain"), lr = risk_est("pain"), hr = risk_hr("pain"))
  )
  cols <- c(Core = "#1F4E78", Supplementary = "#2E7D32", Incremental = "#E67E22", Background = "#7D3C98", Excluded = "#C0392B")
  # Panel A: cross-sectional OR range
  yy <- seq_along(rows)
  plot(NA, xlim = c(0.45, 2.6), ylim = c(0.5, length(rows) + 0.5), log = "x",
       yaxt = "n", xlab = if (zh) "OR(log scale)" else "OR (log scale)", ylab = "",
       main = if (zh) "(a) five wave cross sectional OR range" else "(a) Cross-sectional OR range across five waves", bty = "n")
  abline(v = 1, lty = 2, col = "grey60")
  for (i in seq_along(rows)) {
    z <- rows[[i]]; if (all(is.na(z$or))) next
    segments(z$or[1], i, z$or[3], i, col = cols[z$role], lwd = 3)
    points(z$or[2], i, pch = 18, cex = 1.3, col = cols[z$role])
  }
  axis(2, at = yy, labels = if (zh) vapply(rows, `[[`, "", "name_zh") else vapply(rows, `[[`, "", "name_en"),
       las = 1, cex.axis = 0.75)
  legend("bottomright", if (zh) c("core", "supplementary", "incremental", "background", "exclude") else names(cols), col = cols, lwd = 3, pch = 18, bty = "n", cex = 0.7)
  # Panel B: longitudinal OR/HR forest (only factors with longitudinal estimates)
  rows_b <- rows[!vapply(rows, function(z) all(is.na(z$lr)) && all(is.na(z$hr)), logical(1))]
  lab <- if (zh) vapply(rows_b, `[[`, "", "name_zh") else vapply(rows_b, `[[`, "", "name_en")
  y2 <- rev(seq_along(rows_b))
  plot(NA, xlim = c(0.45, 2.6), ylim = c(0.5, length(rows_b) + 0.5), log = "x",
       yaxt = "n", xlab = if (zh) "OR / HR(log scale)" else "OR / HR (log scale)", ylab = "",
       main = if (zh) "(b) longitudinal effect estimate" else "(b) Longitudinal effect estimates", bty = "n")
  abline(v = 1, lty = 2, col = "grey60")
  for (i in seq_along(rows_b)) {
    z <- rows_b[[i]]; y <- y2[i]
    if (!all(is.na(z$lr))) {
      segments(z$lr[2], y + 0.12, z$lr[3], y + 0.12, col = "#1F4E78", lwd = 2)
      points(z$lr[1], y + 0.12, pch = 19, col = "#1F4E78", cex = 1.1)
    }
    if (!all(is.na(z$hr))) {
      segments(z$hr[2], y - 0.12, z$hr[3], y - 0.12, col = "#C0392B", lwd = 2)
      points(z$hr[1], y - 0.12, pch = 17, col = "#C0392B", cex = 1.1)
    }
  }
  axis(2, at = y2, labels = lab, las = 1, cex.axis = 0.75)
  legend("bottomright", c(if (zh) "endpoint assessment OR" else "Endpoint OR", if (zh) "first impairment HR" else "First impairment HR"),
         col = c("#1F4E78", "#C0392B"), pch = c(19, 17), lwd = 2, bty = "n", cex = 0.75)
}
save_both(plot_candidate_stability, "Figure2_candidate_stability_ZH", "ZH")
save_both(plot_candidate_stability, "Figure2_candidate_stability_EN", "EN")

safe_write_csv(data.frame(pattern = c("depression", "ADL limited", "IADL limited", "sleep <6", "sleep >8", "regular physical activity", "hypertension", "diabetes", "stroke", "bodily pain")),
               file.path(DATA, "Figure2_candidate_patterns.csv"), row.names = FALSE, fileEncoding = "UTF-8")

# ---------------------------------------------------------------------------
# Figure 3: risk stratification
# ---------------------------------------------------------------------------
a14_strata <- read.csv(file.path(P, "output/12_final_model/A14_risk_stratification.csv"),
                       check.names = FALSE, fileEncoding = "UTF-8")
a09_strata <- read.csv(file.path(P, "output/10_screening_index/A09_strata_observed_vs_predicted_2026-09-26.csv"),
                       check.names = FALSE, fileEncoding = "UTF-8")
draw_risk_strata <- function(lang) {
  zh <- lang == "ZH"
  par(mar = c(5, 5, 3, 2), family = "sans", cex.lab = 1.2)
  bp <- barplot(a09_strata$observed_rate, names.arg = if (zh) c("low risk", "intermediate risk", "high risk") else c("Low", "Intermediate", "High"),
                col = c("#C6DBEF", "#6BAED6", "#2171B5"), ylim = c(0, max(a09_strata$observed_rate) * 1.35),
                ylab = if (zh) "4 year incident cognitive impairment rate" else "4-year incident cognitive impairment",
                xlab = if (zh) "screening index risk strata" else "Screening index risk stratum",
                main = if (zh) "screening index risk stratification(CHARLS development cohort)" else "Screening index risk stratification (CHARLS development cohort)")
  text(bp, a09_strata$observed_rate, sprintf("%.1f%%\n(n=%s)", 100*a09_strata$observed_rate, format(a09_strata$n, big.mark=",")), pos = 3, cex = 0.95)
  points(bp, a09_strata$predicted_rate, pch = 18, cex = 1.6, col = "#C0392B")
  legend("topleft", if (zh) c("observed_event_rate", "model_predicted_rate") else c("Observed event rate", "Model-predicted risk"),
         pch = c(15, 18), col = c("#2171B5", "#C0392B"), bty = "n", cex = 0.9)
}
save_both(draw_risk_strata, "Figure3_risk_stratification_ZH", "ZH")
save_both(draw_risk_strata, "Figure3_risk_stratification_EN", "EN")

# ---------------------------------------------------------------------------
# Figure 4: internal validation ROC / calibration / DCA
# ---------------------------------------------------------------------------
a07 <- readRDS(file.path(P, "output/10_screening_index/screening_index_memory_v3_2026-09-26.rds"))
a09_cal <- read.csv(file.path(P, "output/10_screening_index/A09_calibration_deciles_2026-09-26.csv"),
                    check.names = FALSE, fileEncoding = "UTF-8")
a09_dca <- read.csv(file.path(P, "output/10_screening_index/A09_dca_2026-09-26.csv"),
                    check.names = FALSE, fileEncoding = "UTF-8")
d <- a07$data; m <- a07$model_u
d$p_hat <- predict(m, type = "response")
r <- roc(d$incident_2015, d$p_hat, direction = "<", quiet = TRUE)
auc_v <- as.numeric(auc(r)); ci_v <- as.numeric(ci.auc(r))
cal_slope <- glm(incident_2015 ~ qlogis(pmin(pmax(p_hat, 1e-6), 1 - 1e-6)), data = d, family = binomial())$coefficients[[2]]
cal_int <- glm(incident_2015 ~ offset(qlogis(pmin(pmax(p_hat, 1e-6), 1 - 1e-6))), data = d, family = binomial())$coefficients[[1]]
draw_internal <- function(lang) {
  zh <- lang == "ZH"
  par(mfrow = c(1, 3), mar = c(5, 5, 3, 2), family = "sans", cex.lab = 1.0, cex.axis = 0.9)
  plot(1 - r$specificities, r$sensitivities, type = "l", lwd = 2.5, col = "#2C7FB8",
       xlab = if (zh) "1 − specificity" else "1 − Specificity", ylab = if (zh) "sensitivity" else "Sensitivity",
       main = if (zh) "(a) ROC curve" else "(a) ROC curve", xlim = c(0, 1), ylim = c(0, 1))
  abline(0, 1, lty = 3, col = "grey60")
  text(0.60, 0.18, sprintf("AUC = %.3f\n(95%% CI %.3f–%.3f)", auc_v, ci_v[1], ci_v[3]), cex = 0.9)
  plot(a09_cal$mean_prediction, a09_cal$observed_rate, pch = 19, cex = 1.2, col = "#31A354",
       xlim = c(0, max(a09_cal$observed_CI_high) * 1.05), ylim = c(0, max(a09_cal$observed_CI_high) * 1.05),
       xlab = if (zh) "predicted probability" else "Predicted probability", ylab = if (zh) "observed_event_rate" else "Observed event rate",
       main = if (zh) "(b) calibration plot(deciles)" else "(b) Calibration plot (deciles)")
  segments(a09_cal$mean_prediction, a09_cal$observed_CI_low, a09_cal$mean_prediction, a09_cal$observed_CI_high, col = "#31A354")
  abline(0, 1, lty = 3, col = "grey60")
  lines(lowess(a09_cal$mean_prediction, a09_cal$observed_rate), col = "#756BB1", lwd = 1.5)
  text(max(a09_cal$mean_prediction) * 0.42, max(a09_cal$observed_CI_high) * 0.95,
       sprintf(if (zh) "calibration slope = %.2f\n intercept = %.2f" else "Calibration slope = %.2f\nIntercept = %.2f", cal_slope, cal_int), cex = 0.9, adj = 0)
  plot(a09_dca$threshold, a09_dca$net_benefit_model, type = "l", lwd = 2.5, col = "#2C7FB8",
       ylim = c(min(a09_dca$net_benefit_model, -0.02), max(a09_dca$net_benefit_model, mean(d$incident_2015)) * 1.05),
       xlab = if (zh) "threshold probability" else "Threshold probability", ylab = if (zh) "net_benefit" else "Net benefit",
       main = if (zh) "(c) decision curve analysis" else "(c) Decision curve analysis")
  lines(a09_dca$threshold, a09_dca$net_benefit_intervene_in_all, lty = 2, col = "#E6550D", lwd = 2)
  abline(h = 0, lty = 3, col = "grey60")
  legend("topright", if (zh) c("screening index", "intervene in all", "no intervention") else c("Screening index", "Treat all", "Treat none"),
         lty = c(1, 2, 3), lwd = c(2.5, 2, 1), col = c("#2C7FB8", "#E6550D", "grey60"), bty = "n", cex = 0.8)
}
save_both(draw_internal, "Figure4_internal_validation_ZH", "ZH")
save_both(draw_internal, "Figure4_internal_validation_EN", "EN")

# ---------------------------------------------------------------------------
# Figure 5: external validation ROC / calibration / AUC forest
# ---------------------------------------------------------------------------
a12 <- readRDS(file.path(P, "output/11_external_validation/A12_results_2026-09-26.rds"))
model_u <- a07$model_u; model_levels <- model_u$xlevels$edu_group
cohorts <- c("HRS", "ELSA", "SHARE")
cohort_cols <- c(HRS = "#1F77B4", ELSA = "#2CA02C", SHARE = "#D62728")
risk_sets <- list()
for (co in cohorts) {
  x <- readRDS(file.path(P, "external_data", co, "analysis",
                         paste0(co, "_sensitivity_risk_sets_2026-09-26.rds")))[[paste0(co, "__baseline_all_mean_sd")]]
  x$edu_group <- factor(as.character(x$edu_group), levels = model_levels)
  x$p_hat <- as.numeric(predict(model_u, newdata = data.frame(age = x$age, edu_group = x$edu_group,
                                                              dep = x$dep, adl = x$adl), type = "response"))
  risk_sets[[co]] <- x
}
draw_external <- function(lang) {
  zh <- lang == "ZH"
  par(mfrow = c(1, 3), mar = c(4.8, 4.8, 3, 1.5), family = "sans", cex.lab = 1.0, cex.axis = 0.85)
  # (a) ROC
  plot(NA, xlim = c(0, 1), ylim = c(0, 1), xlab = if (zh) "1 − specificity" else "1 − Specificity",
       ylab = if (zh) "sensitivity" else "Sensitivity", main = if (zh) "(a) external cohort ROC(S0)" else "(a) External ROC (S0 model)", bty = "n")
  abline(0, 1, lty = 2, col = "grey70")
  for (co in cohorts) {
    d <- risk_sets[[co]]; d <- d[!is.na(d$p_hat), ]
    r <- roc(d$event, d$p_hat, direction = "<", quiet = TRUE)
    lines(1 - r$specificities, r$sensitivities, col = cohort_cols[co], lwd = 2)
  }
  a <- a12$perf[a12$perf$cut_definition == "baseline_all_mean_sd", ]
  a <- a[match(cohorts, a$cohort), ]
  legend("bottomright", sprintf("%s AUC %.3f (%.3f–%.3f)", a$cohort, a$auc_model, a$auc_model_low, a$auc_model_high),
         col = cohort_cols[a$cohort], lwd = 2, bty = "n", cex = 0.70)
  # (b) Calibration
  dec <- a12$deciles[a12$deciles$cut_definition == "baseline_all_mean_sd", ]
  ymax <- max(c(dec$observed_high, dec$pred_mean), na.rm = TRUE) * 1.05
  plot(NA, xlim = c(0, ymax), ylim = c(0, ymax), xlab = if (zh) "predicted_risk" else "Predicted risk",
       ylab = if (zh) "observed risk" else "Observed risk", main = if (zh) "(b) calibration plot(S0)" else "(b) Calibration by decile (S0 model)", bty = "n")
  abline(0, 1, lty = 2, col = "grey70")
  for (co in cohorts) {
    d <- dec[dec$cohort == co, ]
    arrows(d$pred_mean, d$observed_low, d$pred_mean, d$observed_high, angle = 90, code = 3, length = 0.03, col = cohort_cols[co])
    lines(d$pred_mean, d$observed, type = "b", pch = 19, col = cohort_cols[co], lwd = 2)
  }
  legend("bottomright", cohorts, col = cohort_cols[cohorts], lwd = 2, pch = 19, bty = "n", cex = 0.75)
  # (c) AUC forest
  m <- a12$meta[a12$meta$cut_definition == "baseline_all_mean_sd" & a12$meta$metric == "model", ]
  labels <- c(cohorts, if (zh) "random effects pooling" else "Pooled RE")
  est <- c(a$auc_model, m$pooled_auc); lo <- c(a$auc_model_low, m$pooled_low); hi <- c(a$auc_model_high, m$pooled_high)
  y <- rev(seq_along(labels))
  plot(NA, xlim = c(0.62, 0.86), ylim = c(0.5, length(labels) + 0.5), yaxt = "n",
       xlab = if (zh) "AUC(95% CI)" else "AUC (95% CI)", ylab = "", main = if (zh) "(c) AUC forest plot(S0)" else "(c) AUC forest (S0 model)", bty = "n")
  abline(v = 0.5, lty = 3, col = "grey70")
  segments(lo, y, hi, y, col = c(cohort_cols[a$cohort], "black"), lwd = 2)
  points(est, y, pch = c(19, 19, 19, 18), col = c(cohort_cols[a$cohort], "black"), cex = c(1, 1, 1, 1.4))
  axis(2, at = y, labels = labels, las = 1, cex.axis = 0.8)
}
save_both(draw_external, "Figure5_external_validation_ZH", "ZH")
save_both(draw_external, "Figure5_external_validation_EN", "EN")

a_s0_csv <- a12$perf[a12$perf$cut_definition == "baseline_all_mean_sd", ]
a_s0_csv <- a_s0_csv[match(cohorts, a_s0_csv$cohort), ]
safe_write_csv(data.frame(cohort = cohorts, n = a_s0_csv$n_total, events = a_s0_csv$events_total,
                          auc = a_s0_csv$auc_model),
               file.path(DATA, "Figure5_external_summary.csv"), row.names = FALSE, fileEncoding = "UTF-8")
cat("A16 figures written.\n")
