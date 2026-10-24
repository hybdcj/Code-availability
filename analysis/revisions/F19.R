# ============================================================================
# F19: external broader population definition × external S2 sensitivity threshold(fixed CHARLS <=4)
# purpose: additional external fixed threshold sensitivity and post hoc broader population analysis; does not replace S0 primary analysis
# verification: N definition(baseline >4)should reproduce A11 S2 —— HRS 15521/945/0.7475, 
#       ELSA 7203/324/0.8025, SHARE 28391/2256/0.7769
# ============================================================================
invisible(try(Sys.setlocale("LC_ALL", "zh_CN.UTF-8"), silent = TRUE))
suppressPackageStartupMessages({ library(data.table); library(pROC) })
if (dir.exists("D:/R/library")) .libPaths(c("D:/R/library", .libPaths()))
TMP <- "C:/y03tmp/tmp"; OUT <- "C:/y03tmp/out"
con <- file(file.path(OUT, "F19_report.txt"), open = "wt", encoding = "UTF-8")
say <- function(...) { s <- paste0(..., collapse = ""); cat(s, "\n"); cat(s, "\n", file = con) }
pt <- function(x) say(paste(capture.output(print(x)), collapse = "\n"))

B <- c(intercept = -4.19524508778947, age = 0.0492039250378111,
       eduPri = -0.290105771000328, eduMid = -0.994711800202235,
       eduHigh = -1.43588388035287, dep = 0.412166484954086, adl = -0.0110119344848014)
CO <- c("HRS", "ELSA", "SHARE")
K <- 4L

prep <- function(co) {
  d <- as.data.table(readRDS(file.path(TMP, paste0("F16_base_", co, ".rds"))))
  d <- d[!is.na(age) & age >= 45]
  d[, score := 2L * pmax(floor((age - 45) / 5), 0) + edu_pts + 3L * dep]
  d[, lp := B[["intercept"]] + B[["age"]] * age + B[["dep"]] * dep + B[["adl"]] * adl +
      B[["eduPri"]] * as.integer(edu_pts == -2) +
      B[["eduMid"]] * as.integer(edu_pts == -7) +
      B[["eduHigh"]] * as.integer(edu_pts == -10)]
  d[, phat := plogis(lp)]
  d[, event := as.integer(mem_fu <= K)]
  d[, low0 := as.integer(mem_base <= K)]
  d[]
}
D <- lapply(CO, prep); names(D) <- CO

say("threshold: fixed CHARLS <= ", K, "(external S2 sensitivity definition; S0 remains primary analysis)")
say("\n========== 1. pipeline verification(reproduction A11 S2 of N definition)==========")
say("expected: HRS 15521 / 945 / AUC 0.7475; ELSA 7203 / 324 / 0.8025; SHARE 28391 / 2256 / 0.7769")
ver <- rbindlist(lapply(CO, function(co) {
  d <- D[[co]][!is.na(mem_fu) & !is.na(mem_base) & mem_base > K]
  data.table(cohort = co, n = nrow(d), events = sum(d$event),
             AUC_scorecard = round(as.numeric(auc(roc(d$event, d$score, quiet = TRUE, direction = "<"))), 4))
}))
pt(ver); fwrite(ver, file.path(OUT, "F19_verify_S2.csv"), bom = TRUE)

say("\n========== 2. each definition AUC(fixed <=4)==========")
aucrow <- function(d, v, co, arm, vlab) {
  ok <- !is.na(d[[v]]) & !is.na(d$event)
  if (sum(ok) < 30 || length(unique(d$event[ok])) < 2) return(NULL)
  r <- roc(d$event[ok], d[[v]][ok], quiet = TRUE, direction = "<")
  ci <- as.numeric(ci.auc(r))
  data.table(cohort = co, definition = arm, predictors = vlab, n = sum(ok), events = sum(d$event[ok]),
             event_rate = round(100 * mean(d$event[ok]), 2),
             AUC = round(ci[2], 4), AUC_low = round(ci[1], 4), AUC_high = round(ci[3], 4))
}
mk <- function(d) list(U = d, N = d[!is.na(mem_base) & mem_base > K],
                       I = d[!is.na(mem_base) & mem_base <= K], M = d[is.na(mem_base)])
LB <- c(U = "U_full_population", N = "N_unimpaired_at_baseline", I = "I_impaired_at_baseline", M = "M_missing_baseline_memory")
aucT <- rbindlist(lapply(CO, function(co) {
  st <- mk(D[[co]][!is.na(mem_fu)])
  rbindlist(lapply(names(st), function(s) {
    if (nrow(st[[s]]) == 0) return(NULL)
    rbind(aucrow(st[[s]], "score", co, LB[[s]], "three item scorecard"),
          aucrow(st[[s]], "phat", co, LB[[s]], "continuous model"))
  }))
}))
pt(aucT); fwrite(aucT, file.path(OUT, "F19_auc.csv"), bom = TRUE)

say("\n========== 3. decomposition ==========")
dec <- rbindlist(lapply(CO, function(co) {
  d <- D[[co]][!is.na(mem_fu) & !is.na(mem_base)]
  r_low <- roc(d$event, d$low0, quiet = TRUE, direction = "<")
  r_scr <- roc(d$event, d$score, quiet = TRUE, direction = "<")
  fit <- glm(event ~ score + low0, data = d, family = binomial())
  r_cmb <- roc(d$event, predict(fit, type = "response"), quiet = TRUE, direction = "<")
  data.table(cohort = co, n = nrow(d), events = sum(d$event),
             baseline_low_scorers_as_fraction_of_events = round(100 * mean(d$low0[d$event == 1]), 1),
             event_rate_among_baseline_low_scorers = round(100 * mean(d[low0 == 1]$event), 1),
             event_rate_among_baseline_unimpaired_participants = round(100 * mean(d[low0 == 0]$event), 1),
             AUC_baseline_low_score_only = round(as.numeric(auc(r_low)), 4),
             AUC_scorecard_only = round(as.numeric(auc(r_scr)), 4),
             AUC_combine = round(as.numeric(auc(r_cmb)), 4),
             P_combined_approach_outperforms_scorecard = signif(roc.test(r_cmb, r_scr, paired = TRUE)$p.value, 3))
}))
pt(dec); fwrite(dec, file.path(OUT, "F19_decomposition.csv"), bom = TRUE)

say("\n========== 4. broader population screening simulation ==========")
scr <- rbindlist(lapply(CO, function(co) {
  d <- D[[co]][!is.na(mem_fu)]
  s <- d$score; ev <- d$event; n <- length(ev)
  rbindlist(lapply(c(0.10, 0.20, 0.30), function(fr) {
    k <- floor(fr * n); th <- sort(s, decreasing = TRUE)[k]; pos <- s >= th
    sub <- d[pos]; yld <- sum(sub$event == 1 & sub$low0 == 0, na.rm = TRUE)
    data.table(cohort = co, target_screening_fraction = paste0(fr * 100, "%"), threshold = th, n = length(ev), events = sum(ev), cases_detected = sum(pos), actual_screening_fraction = round(100 * mean(pos), 2), TP = sum(pos & ev == 1), FP = sum(pos & ev == 0), TN = sum(!pos & ev == 0), FN = sum(!pos & ev == 1),
               sensitivity = round(100 * sum(ev[pos]) / sum(ev), 1),
               specificity = round(100 * sum(!pos & ev == 0) / sum(ev == 0), 1),
               negative_predictive_value = round(100 * mean(ev[!pos] == 0), 1), positive_predictive_value = round(100 * mean(ev[pos]), 1),
               fraction_already_low_scoring_at_baseline_among_detected_cases = round(100 * mean(sub$low0 == 1, na.rm = TRUE), 1),
               true_incident_events = yld, per_detected_1_incident_case_number_needed_to_screen = if (yld > 0) round(n / yld, 1) else NA_real_)
  }))
}))
pt(scr); fwrite(scr, file.path(OUT, "F19_screening.csv"), bom = TRUE)

say("\n========== 5. summary linking to main manuscript definition ==========")
pt(rbindlist(lapply(CO, function(co) {
  d0 <- D[[co]][!is.na(mem_fu)]; d1 <- d0[!is.na(mem_base) & mem_base > K]
  data.table(cohort = co,
             N_unimpaired_at_baseline_n = nrow(d1), N_events = sum(d1$event),
             N_AUC = round(as.numeric(auc(roc(d1$event, d1$score, quiet = TRUE, direction = "<"))), 4),
             U_full_population_n = nrow(d0), U_events = sum(d0$event),
             U_AUC = round(as.numeric(auc(roc(d0$event, d0$score, quiet = TRUE, direction = "<"))), 4))
})))

close(con)
cat("\n complete ->", file.path(OUT, "F19_report.txt"), "\n")
