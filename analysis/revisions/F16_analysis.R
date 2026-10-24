# ============================================================================
# F16 analysis: external cohort(HRS / ELSA / SHARE)broader population definition
# and A11/A12 only difference: baseline memory availability no longer required, baseline unimpaired status no longer required
# cutoff definition and A11 v2 primary definition(baseline_all_mean_sd)consistent, for direct comparison
# ============================================================================
invisible(try(Sys.setlocale("LC_ALL", "zh_CN.UTF-8"), silent = TRUE))
suppressPackageStartupMessages({ library(data.table); library(pROC) })
if (dir.exists("D:/R/library")) .libPaths(c("D:/R/library", .libPaths()))
TMP <- "C:/y03tmp/tmp"; OUT <- "C:/y03tmp/out"
con <- file(file.path(OUT, "F16_report.txt"), open = "wt", encoding = "UTF-8")
say <- function(...) { s <- paste0(..., collapse = ""); cat(s, "\n"); cat(s, "\n", file = con) }
pt <- function(x) say(paste(capture.output(print(x)), collapse = "\n"))

B <- c(intercept = -4.19524508778947, age = 0.0492039250378111,
       eduPri = -0.290105771000328, eduMid = -0.994711800202235,
       eduHigh = -1.43588388035287, dep = 0.412166484954086, adl = -0.0110119344848014)
CO <- c("HRS", "ELSA", "SHARE")
LB <- c(U = "U_full_population", N = "N_unimpaired_at_baseline", I = "I_impaired_at_baseline", M = "M_missing_baseline_memory")

prep <- function(co) {
  d <- as.data.table(readRDS(file.path(TMP, paste0("F16_base_", co, ".rds"))))
  d <- d[!is.na(age) & age >= 45]
  d[, score := 2L * pmax(floor((age - 45) / 5), 0) + edu_pts + 3L * dep]
  d[, lp := B[["intercept"]] + B[["age"]] * age + B[["dep"]] * dep + B[["adl"]] * adl +
      B[["eduPri"]] * as.integer(edu_pts == -2) +
      B[["eduMid"]] * as.integer(edu_pts == -7) +
      B[["eduHigh"]] * as.integer(edu_pts == -10)]
  d[, phat := plogis(lp)]
  cs <- d[!is.na(mem_base)]
  d[, cut := mean(cs$mem_base) - sd(cs$mem_base)]
  d[, cut_int := floor(cut)]
  d[, event := as.integer(mem_fu <= cut_int)]
  d[, low0 := as.integer(mem_base <= cut_int)]
  d[, cohort := co]
  d[]
}
D <- lapply(CO, prep); names(D) <- CO

say("cutoff verification(expected HRS 6.244/6, ELSA 6.77/6, SHARE 5.28/5): ")
pt(rbindlist(lapply(CO, function(co) data.table(
  cohort = co, original_cutoff = round(unique(D[[co]]$cut), 3), integer_threshold = unique(D[[co]]$cut_int),
  mem_base_available = sum(!is.na(D[[co]]$mem_base)), complete_age_education_and_depression = nrow(D[[co]])))))

say("\n========== 1. pipeline verification(reproduction A11 primary risk set)==========")
say("expected: HRS 14051 / events 1842; ELSA 6725 / 554; SHARE 26879 / 2930")
ver <- rbindlist(lapply(CO, function(co) {
  d <- D[[co]][!is.na(mem_base) & mem_base > cut_int & !is.na(mem_fu)]
  data.table(cohort = co, risk_set_n = nrow(d), events = sum(d$event),
             AUC_scorecard = round(as.numeric(auc(roc(d$event, d$score, quiet = TRUE, direction = "<"))), 4))
}))
pt(ver); fwrite(ver, file.path(OUT, "F16_verify_A11.csv"), bom = TRUE)

fl <- rbindlist(lapply(CO, function(co) {
  d <- D[[co]]; fu <- !is.na(d$mem_fu)
  data.table(cohort = co, complete_age_education_and_depression = nrow(d), follow_up_memory_scorable = sum(fu),
             U_full_population = sum(fu),
             U_missing_baseline_memory = sum(fu & is.na(d$mem_base)),
             N_unimpaired_at_baseline = sum(fu & !is.na(d$mem_base) & d$mem_base > d$cut_int),
             I_impaired_at_baseline = sum(fu & !is.na(d$mem_base) & d$mem_base <= d$cut_int))
}))
say("\n========== 2. sample flow ==========")
pt(fl); fwrite(fl, file.path(OUT, "F16_flow.csv"), bom = TRUE)

aucrow <- function(d, v, co, arm, vlab) {
  ok <- !is.na(d[[v]]) & !is.na(d$event)
  if (sum(ok) < 30 || length(unique(d$event[ok])) < 2) return(NULL)
  r <- roc(d$event[ok], d[[v]][ok], quiet = TRUE, direction = "<")
  ci <- as.numeric(ci.auc(r))
  data.table(cohort = co, definition = arm, predictors = vlab, n = sum(ok), events = sum(d$event[ok]),
             event_rate = round(100 * mean(d$event[ok]), 2),
             AUC = round(ci[2], 4), AUC_low = round(ci[1], 4), AUC_high = round(ci[3], 4))
}
mk_sets <- function(d) {
  list(U = d,
       N = d[!is.na(mem_base) & mem_base > cut_int],
       I = d[!is.na(mem_base) & mem_base <= cut_int],
       M = d[is.na(mem_base)])
}
aucT <- rbindlist(lapply(CO, function(co) {
  st <- mk_sets(D[[co]][!is.na(mem_fu)])
  rbindlist(lapply(names(st), function(s) {
    if (nrow(st[[s]]) == 0) return(NULL)
    rbind(aucrow(st[[s]], "score", co, LB[[s]], "three item scorecard"),
          aucrow(st[[s]], "phat", co, LB[[s]], "continuous model"))
  }))
}))
say("\n========== 3. each definition AUC ==========")
pt(aucT); fwrite(aucT, file.path(OUT, "F16_auc.csv"), bom = TRUE)

dec <- rbindlist(lapply(CO, function(co) {
  d <- D[[co]][!is.na(mem_fu) & !is.na(mem_base)]
  r_low <- roc(d$event, d$low0, quiet = TRUE, direction = "<")
  r_scr <- roc(d$event, d$score, quiet = TRUE, direction = "<")
  fit <- glm(event ~ score + low0, data = d, family = binomial())
  r_cmb <- roc(d$event, predict(fit, type = "response"), quiet = TRUE, direction = "<")
  data.table(cohort = co, n = nrow(d), events = sum(d$event),
             baseline_low_scorers_as_fraction_of_events = round(100 * mean(d$low0[d$event == 1]), 1),
             AUC_baseline_low_score_only = round(as.numeric(auc(r_low)), 4),
             AUC_scorecard_only = round(as.numeric(auc(r_scr)), 4),
             AUC_combine = round(as.numeric(auc(r_cmb)), 4),
             P_combined_approach_outperforms_scorecard = signif(roc.test(r_cmb, r_scr, paired = TRUE)$p.value, 3))
}))
say("\n========== 4. decomposition: predictive value of baseline low score ==========")
pt(dec); fwrite(dec, file.path(OUT, "F16_decomposition.csv"), bom = TRUE)

scr <- rbindlist(lapply(CO, function(co) {
  d <- D[[co]][!is.na(mem_fu)]
  s <- d$score; ev <- d$event; n <- length(ev)
  out <- rbindlist(lapply(c(0.10, 0.20, 0.30), function(fr) {
    k <- floor(fr * n); th <- sort(s, decreasing = TRUE)[k]; pos <- s >= th
    sub <- d[pos]; yld <- sum(sub$event == 1 & sub$low0 == 0, na.rm = TRUE)
    data.table(cohort = co, target_screening_fraction = paste0(fr * 100, "%"), threshold = th, n = length(ev), events = sum(ev), cases_detected = sum(pos), actual_screening_fraction = round(100 * mean(pos), 2), TP = sum(pos & ev == 1), FP = sum(pos & ev == 0), TN = sum(!pos & ev == 0), FN = sum(!pos & ev == 1),
               sensitivity = round(100 * sum(ev[pos]) / sum(ev), 1),
               specificity = round(100 * sum(!pos & ev == 0) / sum(ev == 0), 1),
               positive_predictive_value = round(100 * mean(ev[pos]), 1),
               negative_predictive_value = round(100 * mean(ev[!pos] == 0), 1),
               fraction_already_low_scoring_at_baseline_among_detected_cases = round(100 * mean(sub$low0 == 1, na.rm = TRUE), 1),
               true_incident_events = yld, per_detected_1_incident_case_number_needed_to_screen = if (yld > 0) round(n / yld, 1) else NA_real_)
  }))
  c2 <- coords(roc(ev, s, quiet = TRUE, direction = "<"), "best", best.method = "youden",
               ret = c("threshold", "sensitivity", "specificity"))
  c2 <- c2[which.max(c2$specificity), , drop = FALSE]
  pos <- s >= c2$threshold; sub <- d[pos]
  yld <- sum(sub$event == 1 & sub$low0 == 0, na.rm = TRUE)
  rbind(out, data.table(cohort = co, target_screening_fraction = "optimal Youden cutoff", threshold = c2$threshold, n = length(ev), events = sum(ev), cases_detected = sum(pos), actual_screening_fraction = round(100 * mean(pos), 2), TP = sum(pos & ev == 1), FP = sum(pos & ev == 0), TN = sum(!pos & ev == 0), FN = sum(!pos & ev == 1),
                        sensitivity = round(100 * c2$sensitivity, 1), specificity = round(100 * c2$specificity, 1),
                        positive_predictive_value = round(100 * mean(ev[pos]), 1),
                        negative_predictive_value = round(100 * mean(ev[!pos] == 0), 1),
                        fraction_already_low_scoring_at_baseline_among_detected_cases = round(100 * mean(sub$low0 == 1, na.rm = TRUE), 1),
                        true_incident_events = yld, per_detected_1_incident_case_number_needed_to_screen = if (yld > 0) round(n / yld, 1) else NA_real_))
}))
say("\n========== 5. broader population screening simulation ==========")
pt(scr); fwrite(scr, file.path(OUT, "F16_screening.csv"), bom = TRUE)

strat_auc <- function(s, co, kind, lab) {
  if (nrow(s) < 30 || length(unique(s$event)) < 2) return(NULL)
  data.table(cohort = co, stratification = kind, strata = lab, n = nrow(s),
             event_rate = round(100 * mean(s$event), 1),
             AUC = round(as.numeric(auc(roc(s$event, s$score, quiet = TRUE, direction = "<"))), 4))
}
st <- rbindlist(lapply(CO, function(co) {
  d <- D[[co]][!is.na(mem_fu)]
  d[, ag := fcase(age < 60, "45-59", age < 70, "60-69", age < 80, "70-79", default = "80+")]
  out <- rbindlist(lapply(c("45-59", "60-69", "70-79", "80+"),
                          function(a) strat_auc(d[ag == a], co, "age_group", a)))
  if (!all(is.na(d$sex))) {
    sx <- d[!is.na(sex)]
    sx[, sexl := paste0("values", as.character(as.integer(sex)))]
    out <- rbind(out, rbindlist(lapply(unique(sx$sexl),
                                       function(v) strat_auc(sx[sexl == v], co, "sex", v))))
    say("   [", co, "] raw sex value counts: ",
        paste(names(table(as.integer(sx$sex))), as.integer(table(as.integer(sx$sex))),
              sep = "=", collapse = "  "))
  }
  eg <- d[!is.na(mem_base)]
  eg[, elab := fcase(edu_pts == 0, "no schooling", edu_pts == -2, "primary school", edu_pts == -7, "middle school",
                     edu_pts == -10, "high school or above", default = "other")]
  out <- rbind(out, rbindlist(lapply(unique(eg$elab),
                                     function(v) strat_auc(eg[elab == v], co, "education", v))))
  out
}))
say("\n========== 6. stratification stability ==========")
pt(st); fwrite(st, file.path(OUT, "F16_strata.csv"), bom = TRUE)

cal <- rbindlist(lapply(CO, function(co) {
  d0 <- D[[co]][!is.na(mem_fu) & !is.na(phat)]
  rbindlist(lapply(c("U", "N"), function(s) {
    d <- if (s == "U") d0 else d0[!is.na(mem_base) & mem_base > cut_int]
    fit <- glm(event ~ lp, data = d, family = binomial())
    data.table(cohort = co, definition = LB[[s]], n = nrow(d),
               expected_event_rate = round(100 * mean(d$phat), 1), actual_event_rate = round(100 * mean(d$event), 1),
               calibration_slope = round(coef(fit)[2], 3), calibration_intercept = round(coef(fit)[1], 3),
               Brier = round(mean((d$phat - d$event)^2), 4))
  }))
}))
say("\n========== 7. calibration ==========")
pt(cal); fwrite(cal, file.path(OUT, "F16_calibration.csv"), bom = TRUE)

sh <- D[["SHARE"]][!is.na(mem_fu) & !is.na(country_label)]
cn <- sh[, .(n = .N, event_rate = round(100 * mean(event), 1),
             broader_population_AUC = round(as.numeric(auc(roc(event, score, quiet = TRUE, direction = "<"))), 3),
             baseline_low_score_fraction = round(100 * mean(low0 == 1, na.rm = TRUE), 1)),
         by = country_label][order(-n)]
cn <- cn[n >= 200]
say("\n========== 8. SHARE by country(n>=200, broader population definition)==========")
pt(cn); fwrite(cn, file.path(OUT, "F16_share_country.csv"), bom = TRUE)

close(con)
cat("\n complete ->", file.path(OUT, "F16_report.txt"), "\n")
