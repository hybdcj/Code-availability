# ============================================================================
# 25_a12_figures_report.R —— A12 figures tables and results report
# input: output/11_external_validation/A12_results_2026-09-26.rds
# output: Figure5, FigureS1, FigureS2(PNG 300 dpi + TIF 600 dpi)
#       A12_external_validation_results_2026_09_26.md
# ============================================================================
invisible(try(Sys.setlocale("LC_ALL", "zh_CN.UTF-8"), silent = TRUE))
.libPaths(c("D:/R/library", .libPaths()))
suppressPackageStartupMessages({ library(dplyr); library(pROC); library(knitr) })

P   <- "D:/BaiduSyncdisk/research/Xu_Zhonghao"
OUT <- file.path(P, "output", "11_external_validation")
FIG <- file.path(P, "output", "09_figures_and_tables")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)

res <- readRDS(file.path(OUT, "A12_results_2026-09-26.rds"))
perf <- res$perf; strata <- res$strata; deciles <- res$deciles; thresholds <- res$thresholds
dca <- res$dca; share_country <- res$share_country; share_boot <- res$share_boot; meta <- res$meta
extra_dep <- res$extra_dep; extra_score_alt <- res$extra_score_alt; extra_elsa_age <- res$extra_elsa_age

a07 <- readRDS(file.path(P, "output/10_screening_index/screening_index_memory_v3_2026-09-26.rds"))
model_u <- a07$model_u; model_levels <- model_u$xlevels$edu_group
cohorts <- c("HRS","ELSA","SHARE")
cohort_cols <- c(HRS = "#1f77b4", ELSA = "#2ca02c", SHARE = "#d62728")
scenario_lab <- c(baseline_all_mean_sd = "S0_baseline_all", complete_fu_mean_sd = "S1_complete_fu", fixed_charls_le4 = "S2_fixed_le4")

add_pred <- function(d, model) {
  d <- as.data.frame(d)
  d$edu_group <- factor(as.character(d$edu_group), levels = model_levels)
  d$p_hat <- as.numeric(predict(model, newdata = data.frame(age = d$age, edu_group = d$edu_group,
                                                            dep = d$dep, adl = d$adl), type = "response"))
  d$lp <- qlogis(pmin(pmax(d$p_hat, 1e-6), 1 - 1e-6))
  d
}
risk_sets <- list()
for (co in cohorts) {
  x <- readRDS(file.path(P, "external_data", co, "analysis", paste0(co, "_sensitivity_risk_sets_2026-09-26.rds")))
  risk_sets[[co]] <- add_pred(x[[paste(co, "baseline_all_mean_sd", sep = "__")]], model_u)
}

# ---------- plotting functions ----------
draw_main <- function() {
  par(mfrow = c(1,3), mar = c(4.2,4.2,3,1), oma = c(0,0,0,0))
  # (a) ROC
  plot(NA, xlim = c(0,1), ylim = c(0,1), xlab = "1 - Specificity", ylab = "Sensitivity",
       main = "(a) External ROC (S0 model)", bty = "n")
  abline(0,1,lty = 2,col = "grey70")
  for (co in cohorts) {
    d <- risk_sets[[co]]; d <- d[!is.na(d$p_hat), ]
    r <- roc(d$event, d$p_hat, direction = "<", quiet = TRUE)
    lines(1 - r$specificities, r$sensitivities, col = cohort_cols[co], lwd = 2)
  }
  a <- perf %>% filter(cut_definition == "baseline_all_mean_sd")
  leg <- sprintf("%s AUC %.3f (%.3f-%.3f)", a$cohort, a$auc_model, a$auc_model_low, a$auc_model_high)
  legend("bottomright", legend = leg, col = cohort_cols[a$cohort], lwd = 2, bty = "n", cex = 0.72)
  # (b) Calibration deciles
  dec <- deciles %>% filter(cut_definition == "baseline_all_mean_sd")
  ymax <- max(c(dec$observed_high, dec$pred_mean), na.rm = TRUE) * 1.05
  plot(NA, xlim = c(0, ymax), ylim = c(0, ymax), xlab = "Predicted risk", ylab = "Observed risk",
       main = "(b) Calibration by decile (S0 model)", bty = "n")
  abline(0,1,lty = 2,col = "grey70")
  for (co in cohorts) {
    d <- dec[dec$cohort == co, ]
    arrows(d$pred_mean, d$observed_low, d$pred_mean, d$observed_high, angle = 90, code = 3,
           length = 0.03, col = cohort_cols[co])
    lines(d$pred_mean, d$observed, type = "b", pch = 19, col = cohort_cols[co], lwd = 2)
  }
  legend("bottomright", legend = cohorts, col = cohort_cols[cohorts], lwd = 2, pch = 19, bty = "n", cex = 0.8)
  # (c) Forest plot
  a <- perf %>% filter(cut_definition == "baseline_all_mean_sd") %>% arrange(match(cohort, cohorts))
  m <- meta %>% filter(cut_definition == "baseline_all_mean_sd", metric == "model")
  labels <- c(a$cohort, "Pooled RE")
  est <- c(a$auc_model, m$pooled_auc); lo <- c(a$auc_model_low, m$pooled_low); hi <- c(a$auc_model_high, m$pooled_high)
  y <- rev(seq_along(labels))
  plot(NA, xlim = c(0.62, 0.86), ylim = c(0.5, length(labels)+0.5), yaxt = "n",
       xlab = "AUC (95% CI)", ylab = "", main = "(c) AUC forest (S0 model)", bty = "n")
  abline(v = 0.5, lty = 3, col = "grey70")
  segments(lo, y, hi, y, col = c(cohort_cols[a$cohort], "black"), lwd = 2)
  points(est, y, pch = c(19,19,19,18), col = c(cohort_cols[a$cohort], "black"), cex = c(1,1,1,1.4))
  axis(2, at = y, labels = labels, las = 1, cex.axis = 0.85)
}
draw_sens <- function() {
  par(mfrow = c(1,3), mar = c(4,3.2,3,1))
  for (cd in c("baseline_all_mean_sd","complete_fu_mean_sd","fixed_charls_le4")) {
    a <- perf %>% filter(cut_definition == cd) %>% arrange(match(cohort, cohorts))
    m <- meta %>% filter(cut_definition == cd, metric == "model")
    est <- c(a$auc_model, m$pooled_auc); lo <- c(a$auc_model_low, m$pooled_low); hi <- c(a$auc_model_high, m$pooled_high)
    labels <- c(a$cohort, "Pooled")
    y <- rev(seq_along(labels))
    plot(NA, xlim = c(0.62, 0.88), ylim = c(0.5, length(labels)+0.5), yaxt = "n",
         xlab = "AUC", ylab = "", main = scenario_lab[[cd]], bty = "n")
    abline(v = 0.5, lty = 3, col = "grey70")
    segments(lo, y, hi, y, col = c(cohort_cols[a$cohort], "black"), lwd = 2)
    points(est, y, pch = 19, col = c(cohort_cols[a$cohort], "black"), cex = c(1,1,1,1.3))
    axis(2, at = y, labels = labels, las = 1, cex.axis = 0.8)
  }
}
draw_share_country <- function() {
  d <- share_country %>% arrange(auc_model)
  y <- seq_len(nrow(d))
  par(mar = c(4,7,3,1))
  plot(NA, xlim = c(0.6, 0.9), ylim = c(0.5, nrow(d)+0.5), yaxt = "n",
       xlab = "AUC (95% CI)", ylab = "", main = "SHARE country-level AUC (S0 model)", bty = "n")
  abline(v = 0.5, lty = 3, col = "grey70")
  segments(d$auc_model_low, y, d$auc_model_high, y, col = "#d62728", lwd = 2)
  points(d$auc_model, y, pch = 19, col = "#d62728", cex = 0.9)
  axis(2, at = y, labels = d$country_label, las = 1, cex.axis = 0.78)
}
draw_cal_sens <- function() {
  par(mfrow = c(1,3), mar = c(4,3.5,3,1))
  for (cd in c("baseline_all_mean_sd","complete_fu_mean_sd","fixed_charls_le4")) {
    dec <- deciles %>% filter(cut_definition == cd)
    ymax <- max(c(dec$observed_high, dec$pred_mean), na.rm = TRUE) * 1.05
    plot(NA, xlim = c(0, ymax), ylim = c(0, ymax), xlab = "Predicted risk", ylab = "Observed risk",
         main = scenario_lab[[cd]], bty = "n")
    abline(0,1,lty = 2,col = "grey70")
    for (co in cohorts) {
      d <- dec[dec$cohort == co, ]
      arrows(d$pred_mean, d$observed_low, d$pred_mean, d$observed_high, angle = 90, code = 3,
             length = 0.03, col = cohort_cols[co])
      lines(d$pred_mean, d$observed, type = "b", pch = 19, col = cohort_cols[co], lwd = 1.7)
    }
    if (cd == "baseline_all_mean_sd") legend("bottomright", legend = cohorts, col = cohort_cols[cohorts], lwd = 2, pch = 19, bty = "n", cex = 0.75)
  }
}
draw_dca <- function() {
  par(mfrow = c(1,3), mar = c(4,4,3,1))
  for (co in cohorts) {
    d <- dca %>% filter(cohort == co, cut_definition == "baseline_all_mean_sd")
    ylim <- range(d$net_benefit, na.rm = TRUE)
    plot(NA, xlim = c(0.05,0.40), ylim = ylim, xlab = "Threshold probability", ylab = "Net benefit",
         main = paste0(co, " (S0)"), bty = "n")
    abline(h = 0, lty = 3, col = "grey70")
    for (st in c("model","all")) {
      dd <- d[d$strategy == st, ]
      lines(dd$threshold, dd$net_benefit, col = ifelse(st == "model", cohort_cols[co], "grey40"),
            lwd = 2, lty = ifelse(st == "model", 1, 2))
    }
    legend("topright", legend = c("Model","Treat all"), col = c(cohort_cols[co],"grey40"),
           lwd = 2, lty = c(1,2), bty = "n", cex = 0.75)
  }
}
save_both <- function(base, width, height, draw) {
  png(paste0(base, ".png"), width = width, height = height, units = "in", res = 300)
  draw(); dev.off()
  tiff(paste0(base, ".tif"), width = width, height = height, units = "in", res = 600, compression = "lzw")
  draw(); dev.off()
}
save_both(file.path(FIG, "Figure5_external_validation_2026-09-26"), 13, 4.5, draw_main)
save_both(file.path(FIG, "FigureS1_external_sensitivity_forest_2026-09-26"), 12, 4.2, draw_sens)
save_both(file.path(FIG, "FigureS2_SHARE_country_forest_2026-09-26"), 8, 5.5, draw_share_country)
save_both(file.path(FIG, "FigureS3_external_calibration_S0_S1_S2_2026-09-26"), 12, 4.2, draw_cal_sens)
save_both(file.path(FIG, "FigureS4_external_DCA_S0_2026-09-26"), 12, 4.2, draw_dca)

# ---------- table formatting ----------
fmt3 <- function(x) ifelse(is.na(x), "NA", sprintf("%.3f", x))
fmt2 <- function(x) ifelse(is.na(x), "NA", sprintf("%.2f", x))
fmt_pct <- function(x) ifelse(is.na(x), "NA", sprintf("%.1f%%", 100*x))
rng <- function(a,b) paste0(fmt3(a), " (", fmt3(b), ")", ifelse(is.na(a) | is.na(b), "", ""))
rng2 <- function(a,b) paste0(fmt2(a), " (", fmt2(b), ")")

perf_s0 <- perf %>% filter(cut_definition == "baseline_all_mean_sd") %>% arrange(match(cohort, cohorts))
tab_main <- data.frame(
  cohort = perf_s0$cohort,
  model_n = perf_s0$n_model,
  events = perf_s0$events_model,
  continuous_model_AUC = paste0(fmt3(perf_s0$auc_model), " (", fmt3(perf_s0$auc_model_low), "–", fmt3(perf_s0$auc_model_high), ")"),
  integer_score_AUC = paste0(fmt3(perf_s0$auc_score_full), " (", fmt3(perf_s0$auc_score_full_low), "–", fmt3(perf_s0$auc_score_full_high), ")"),
  calibration_slope = paste0(fmt2(perf_s0$cal_slope), " (", fmt2(perf_s0$cal_slope_low), " to ", fmt2(perf_s0$cal_slope_high), ")"),
  calibration_intercept = paste0(fmt2(perf_s0$cal_intercept), " (", fmt2(perf_s0$cal_intercept_low), " to ", fmt2(perf_s0$cal_intercept_high), ")"),
  Brier = fmt3(perf_s0$brier), standardized_Brier = fmt3(perf_s0$scaled_brier)
)
tab_strata <- strata %>% filter(cut_definition == "baseline_all_mean_sd") %>%
  mutate(stratification = as.character(stratum),
         observed_event_proportion = fmt_pct(rate), predicted_risk = fmt_pct(pred_mean), OE = fmt2(oe_model)) %>%
  select(cohort = cohort, stratification, n, events = events, observed_event_proportion, predicted_risk, OE)
tab_thr <- thresholds %>% filter(cut_definition == "baseline_all_mean_sd") %>%
  mutate(threshold = paste0(threshold*100, "%"), sensitivity = fmt_pct(sensitivity), specificity = fmt_pct(specificity),
         PPV = fmt_pct(ppv), NPV = fmt_pct(npv), net_benefit = fmt3(nb_model), net_benefit_of_intervening_in_all = fmt3(nb_all)) %>%
  mutate(positive_rate = fmt_pct(pos_rate)) %>% select(cohort = cohort, threshold, positive_count = n_positive, positive_rate, sensitivity, specificity, PPV, NPV, net_benefit, net_benefit_of_intervening_in_all)
tab_meta <- meta %>% mutate(cut = scenario_lab[cut_definition], metric = ifelse(metric == "model", "continuous model", "integer score"),
                            pooled = paste0(fmt3(pooled_auc), " (", fmt3(pooled_low), "–", fmt3(pooled_high), ")"),
                            I2 = sprintf("%.1f%%", I2)) %>%
  select(cutoff = cut, metrics = metric, combine_AUC = pooled, I2, k)
tab_sens <- perf %>% filter(cut_definition != "baseline_all_mean_sd") %>%
  mutate(cutoff = scenario_lab[cut_definition],
         continuous_model_AUC = paste0(fmt3(auc_model), " (", fmt3(auc_model_low), "–", fmt3(auc_model_high), ")"),
         integer_score_AUC = paste0(fmt3(auc_score_full), " (", fmt3(auc_score_full_low), "–", fmt3(auc_score_full_high), ")"),
         calibration_slope = paste0(fmt2(cal_slope), " (", fmt2(cal_slope_low), " to ", fmt2(cal_slope_high), ")"), calibration_intercept = paste0(fmt2(cal_intercept), " (", fmt2(cal_intercept_low), " to ", fmt2(cal_intercept_high), ")"), Brier = fmt3(brier)) %>%
  select(cohort = cohort, cutoff, n = n_model, events = events_model, continuous_model_AUC, integer_score_AUC, calibration_slope, calibration_intercept, Brier)
tab_share <- share_country %>% arrange(desc(n)) %>%
  mutate(AUC = paste0(fmt3(auc_model), " (", fmt3(auc_model_low), "–", fmt3(auc_model_high), ")"),
         event_rate = fmt_pct(event_rate), calibration_slope = fmt2(cal_slope), calibration_intercept = fmt2(cal_intercept)) %>%
  select(country = country_label, n, events = events, event_rate, AUC, calibration_slope, calibration_intercept)
tab_extra <- bind_rows(
  extra_dep %>% transmute(cohort = cohort, sensitivity_analysis = "depression cutoff", rules = dep_rule, n, events = events,
                          continuous_model_AUC = paste0(fmt3(auc_model), " (", fmt3(auc_model_low), "–", fmt3(auc_model_high), ")"),
                          integer_score_AUC = paste0(fmt3(auc_score), " (", fmt3(auc_score_low), "–", fmt3(auc_score_high), ")")),
  extra_score_alt %>% transmute(cohort = cohort, sensitivity_analysis = "score_alt", rules = "positional education mapping", n, events = events,
                                continuous_model_AUC = "—", integer_score_AUC = paste0(fmt3(auc_score_alt), " (", fmt3(auc_low), "–", fmt3(auc_high), ")")),
  extra_elsa_age %>% transmute(cohort = cohort, sensitivity_analysis = "ELSA exclude age=90", rules = "age top coding sensitivity", n, events = events,
                               continuous_model_AUC = paste0(fmt3(auc_model), " (", fmt3(auc_model_low), "–", fmt3(auc_model_high), ")"),
                               integer_score_AUC = paste0(fmt3(auc_score), " (", fmt3(auc_score_low), "–", fmt3(auc_score_high), ")"))
)
tab_dca <- dca %>% filter(cut_definition == "baseline_all_mean_sd", strategy %in% c("model","all"),
                          round(threshold, 2) %in% c(0.10,0.15,0.20,0.25,0.30)) %>%
  mutate(threshold = paste0(threshold*100, "%")) %>%
  select(cohort = cohort, threshold, strategy = strategy, net_benefit = net_benefit) %>%
  mutate(net_benefit = fmt3(net_benefit))

# ---------- report ----------
m_s0_model <- meta %>% filter(cut_definition == "baseline_all_mean_sd", metric == "model")
m_s1_model <- meta %>% filter(cut_definition == "complete_fu_mean_sd", metric == "model")
m_s2_model <- meta %>% filter(cut_definition == "fixed_charls_le4", metric == "model")
meta_text <- sprintf(
  "> random effects use DerSimonian–Laird method, at logit(AUC) pool on scale.S0 continuous model pooling AUC = %.3f(%.3f–%.3f), but **I² = %.1f%%**; S1 combine AUC = %.3f(%.3f–%.3f), I² = %.1f%%; S2 combine AUC = %.3f(%.3f–%.3f), I² = %.1f%%.very high heterogeneity, pooled values are descriptive only, manuscript should emphasize cohort results and heterogeneity sources.",
  m_s0_model$pooled_auc, m_s0_model$pooled_low, m_s0_model$pooled_high, m_s0_model$I2,
  m_s1_model$pooled_auc, m_s1_model$pooled_low, m_s1_model$pooled_high, m_s1_model$I2,
  m_s2_model$pooled_auc, m_s2_model$pooled_low, m_s2_model$pooled_high, m_s2_model$I2)
sb_model <- share_boot[share_boot$cut_definition == "baseline_all_mean_sd" & share_boot$metric == "auc_model", ]
sb_score <- share_boot[share_boot$cut_definition == "baseline_all_mean_sd" & share_boot$metric == "auc_score", ]
share_boot_text <- sprintf(
  "> by country cluster bootstrap(%d times): continuous model AUC **%.3f(%.3f–%.3f)**, integer score AUC **%.3f(%.3f–%.3f)**.point estimate and conventional AUC close; cluster bootstrap CI wider, reflects country clustering uncertainty.emphasize country level differences in manuscript.",
  sb_model$B, sb_model$orig_auc, sb_model$boot_low, sb_model$boot_high,
  sb_score$orig_auc, sb_score$boot_low, sb_score$boot_high)
share_s0_fig <- risk_sets[["SHARE"]]
n_missing_adl_share <- sum(is.na(share_s0_fig$p_hat))
events_missing_adl_share <- sum(share_s0_fig$event[is.na(share_s0_fig$p_hat)])
share_auc_text <- paste0(
  sprintf("> S0 SHARE substantial country heterogeneity: country AUC range approximately %.3f–%.3f, event rate range %.1f%%–%.1f%%.", 
          min(share_country$auc_model), max(share_country$auc_model),
          100*min(share_country$event_rate), 100*max(share_country$event_rate)),
  "country level calibration slope/intercept based on smaller samples, exploratory result, not a main conclusion.
",
  sprintf("> country table uses continuous model complete cases n=%s; all SHARE risk set n=%s(%d participants ADL missing, %d events).",
          format(sum(!is.na(share_s0_fig$p_hat)), big.mark = ","), format(nrow(share_s0_fig), big.mark = ","),
          n_missing_adl_share, events_missing_adl_share))
md <- c(
"# A12: formal external validation results(2026-09-26)",
"",
"> primary analysis: **S0 = all eligible baseline participants in each cohort mean−1SD**; formal sensitivity: **S1 = complete-fu mean−1SD**, **S2 = fixed CHARLS ≤4**.",
"> model: CHARLS fixed model(age + education + depression + ADL), **without refitting**.",
"> cohort: HRS(w10→w12), ELSA(w4→w6), SHARE(w4→w6).",
"",
"## methods overview",
"",
"- use CHARLS fixed model(age + education + depression + ADL), in external cohorts**without refitting, do not reselect variables**; ",
"- primary analysis S0 = all eligible baseline participants in each cohort mean−1SD; S1 = complete-fu mean−1SD; S2 = fixed CHARLS ≤4; ",
"- AUC using DeLong calculated using method 95%CI; calibration slope = glm(event ~ logit(p)) in logit(p) coefficients of; calibration intercept = glm(event ~ offset(logit(p))) intercept of(calibration-in-the-large); Brier and standardized Brier calculated in model complete cases; ",
"- cross cohort pooling uses DerSimonian–Laird random effects model(logit AUC scale), SE from each cohort 95%CI back calculate; SHARE pooling uses country cluster bootstrap CI(fixed seed 20260926), HRS/ELSA use DeLong CI; ",
"- risk stratification uses CHARLS fixed cutoff low ≤−3｜in −2～4｜high ≥5.",
"",
"## one, primary analysis(S0): external performance of continuous model",
"",
kable(tab_main, format = "pipe", row.names = FALSE),
"",
"> continuous model analysis uses ADL complete cases: HRS 13,972(missing 79 participants, including events 9), ELSA 6,725(missing 0), SHARE 26,874(missing 5, events 0).integer score AUC use the full population at risk.
> among model complete cases, integer score AUC respectively HRS 0.693, ELSA 0.759, SHARE 0.752, almost identical to continuous model; integer score in table above AUC use the full population at risk.main table SHARE of CI is DeLong CI, country clustering not considered; by country cluster bootstrap of CI see section six.",
"> **main findings**: continuous model AUC is HRS 0.696, ELSA 0.763, SHARE 0.751; integer score AUC is 0.692, 0.759, 0.752.similar point estimates(approximately 0.70–0.76), but high cross cohort statistical heterogeneity(I² 85%–95%, see section four).",
"> **substantial calibration heterogeneity**: HRS calibration slope 1.23, intercept +0.33(systematic risk underestimation); ELSA slope 1.82, intercept −0.19; SHARE slope 1.09, intercept −0.31(systematic risk overestimation).",
"",
"## two, S0 risk stratification(fixed score cutoffs)",
"",
kable(tab_strata, format = "pipe", row.names = FALSE),
"",
"> stratification cutoffs are CHARLS fixed cutoff low ≤−3｜in −2～4｜high ≥5.observed event rates use full population at risk; predicted risk and O/E calculated in model complete cases.",
"> HRS high risk event rate 36.3%, ELSA 31.7%, SHARE 27.1%; high risk/high to low risk ratios approximately **4.9×, 10.5×, 7.3×**.",
"",
"## three, S0 threshold performance(predicted probability)",
"",
kable(tab_thr, format = "pipe", row.names = FALSE),
"",
"> threshold performance based on model complete cases; corresponding DCA see FigureS4.",
"",
"## four, cross cohort random effects summary",
"",
kable(tab_meta, format = "pipe", row.names = FALSE),
"",
meta_text,
"",
"## five, formal sensitivity analyses(S1/S2)",
"",
kable(tab_sens, format = "pipe", row.names = FALSE),
"",
"> S1/S2 of table n and events are**continuous model complete sample**; corresponding full risk set event count: S2 HRS 945, ELSA 324, SHARE 2,256, S1 ELSA 706, HRS/SHARE and S0 identical.
> ELSA only cohort materially affected by cutoff definition: S0 ≤6 vs S1 ≤7; S2 fixed ≤4 increases in all three cohorts AUC, but event rates substantially decrease(full risk set HRS 6.1%, ELSA 4.5%, SHARE 7.9%), depends on absolute comparability of memory scores across cohorts.report all three result sets.",
"",
"## six, SHARE country level and clustering",
"",
kable(tab_share, format = "pipe", row.names = FALSE),
"",
share_auc_text,
share_boot_text,
"",
"## seven, additional sensitivity analyses",
"",
kable(tab_extra, format = "pipe", row.names = FALSE),
"",
"> depression cutoff and ELSA age90 rows use model complete cases; score_alt rows use full risk set.depression cutoff changed to CES-D-8 ≥3 / EURO-D ≥5 after, AUC largely unchanged; positional education mapping score_alt of AUC is HRS 0.707, ELSA 0.755, SHARE 0.756; ELSA exclude top coded ages 90 after AUC slight decline(model 0.757, integer score 0.753).",
"",
"## eight, decision curve(S0, selected net benefit results)",
"",
kable(tab_dca, format = "pipe", row.names = FALSE),
"",
"> DCA based on model complete cases; at 10%–30% within threshold range, S0 continuous model net benefit exceeds“intervene in all”strategy; approximately 35% then gradually approaches 0.",
"",
"## nine, interpretation and limitations",
"",
"1. **intended role**: tool assesses preoperative related cognitive vulnerability**screening/triage**tool, does not predict postoperative delirium, PND, ICU or hospitalization outcomes.",
"2. **outcome definition**: external cohort outcome is 4 year follow up“threshold crossing”, including baseline position, measurement error and regression to mean components, not equivalent to clinical dementia/MCI incidence.",
"3. **discrimination**: continuous model AUC approximately 0.70–0.76, moderate; cannot claim high accuracy prediction.",
"4. **calibration**: HRS systematic risk underestimation, SHARE systematic risk overestimation, ELSA slope substantially >1; report calibration differences accurately, must not report only AUC.",
"5. **heterogeneity**: cross cohort pooling I² 85%–95%, indicates substantial cohort differences; SHARE country heterogeneity also substantial.",
"6. **measurement**: ELSA age at 90 age top coding; cross cohort memory word lists, language and administration differ; S2 fixed ≤4 assumes absolute comparability, requires caution.",
"7. **missing**: continuous model analysis excludes ADL missing HRS 79 participants, SHARE 5 participants; integer score analysis unaffected.",
"8. **clinical validation**: not yet validated in a surgical cohort, cannot claim improved perioperative decisions or clinical outcomes.",
"9. **follow up selection**: analysis restricted to 4 participants with measurable memory after years; death, attrition or missing memory excluded, possible selection bias, death not treated as a competing risk.",
"10._meaning_of_thresholds_S0/S1_uses_within_each_cohort_mean_1SD_and_CHARLS_fixed_4_not_the_same_absolute_outcome_S2_fixed_4_addresses_same_absolute_threshold_issue_but_assumes_absolute_memory_score_comparability_across_cohorts.three_complementary_result_sets_not_interchangeable.",
"",
"## ten, output file",
"",
"- `A12_primary_results_three_cohorts_S0_S1_S2.csv`",
"- `A12_risk_stratification_S0_S1_S2.csv`",
"- `A12_calibration_deciles_S0_S1_S2.csv`",
"- `A12_threshold_performance_S0_S1_S2.csv`",
"- `A12_DCA_S0_S1_S2.csv`",
"- `A12_SHARE_country_level.csv`, `A12_SHARE_country_clustering_bootstrap.csv`",
"- `A12_cross_cohort_random_effects.csv`",
"- `A12_sensitivity_depression_cutoff.csv`, `A12_sensitivity_score_alt.csv`, `A12_sensitivity_ELSA_age90.csv`",
"- `A12_sessionInfo_2026-09-26.txt`",
"- `A12_independent_verification_2026_09_26.txt`",
"- `Figure5_external_validation_2026-09-26.png/.tif`",
"- `FigureS1_external_sensitivity_forest_2026-09-26.png/.tif`",
"- `FigureS2_SHARE_country_forest_2026-09-26.png/.tif`",
"- `FigureS3_external_calibration_S0_S1_S2_2026-09-26.png/.tif`",
"- `FigureS4_external_DCA_S0_2026-09-26.png/.tif`"
)
writeLines(md, file.path(OUT, "A12_external_validation_results_2026_09_26.md"), useBytes = TRUE)
cat("Figures and report written.\n")