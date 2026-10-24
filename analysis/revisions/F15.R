# F32(2026-10-16): Retain the historical common complete sample calculation below for traceability; CHARLS The integer score now uses its own sample from the same directory F32_correct_population_eligibility_and_denominators.R.Do not merely remove the lp restriction and run the old script, because the downstream continuous model, denominators and headings require separate handling.
# ============================================================================
# F15  broad-population application of the frozen 3-item memory risk score
# frozen model A14(2011 -> 2015 development): three item scorecard + continuous model
# outcome: follow up wave memory score < education strata mean-1SD cutoff(cutoff from 2011 wave definition)
# definition: U full population(no baseline memory restriction)| N unimpaired at baseline(current)| I impaired at baseline
# window: W1 = 2011 -> 2015(development cohort)| W2 = 2015 -> 2018(independent sample)
# ============================================================================
invisible(try(Sys.setlocale("LC_ALL", "zh_CN.UTF-8"), silent = TRUE))
suppressPackageStartupMessages({ library(data.table); library(pROC) })
if (dir.exists("D:/R/library")) .libPaths(c("D:/R/library", .libPaths()))
TMP <- "C:/y03tmp/tmp"; OUT <- "C:/y03tmp/out"
con <- file(file.path(OUT, "F15_report.txt"), open = "wt", encoding = "UTF-8")
say <- function(...) { s <- paste0(..., collapse = ""); cat(s, "\n"); cat(s, "\n", file = con) }
pt  <- function(x) say(paste(capture.output(print(x)), collapse = "\n"))

B <- c(intercept = -4.19524508778947, age = 0.0492039250378111,
       eduPri = -0.290105771000328, eduMid = -0.994711800202235,
       eduHigh = -1.43588388035287, dep = 0.412166484954086, adl = -0.0110119344848014)
LV <- c("Less than primary", "Primary school", "Middle school", "High school or above")
G3 <- c("A_none", "B_primid", "C_high")

## ---------------- 1. data ---------------- 
w <- as.data.table(readRDS(file.path(TMP, "waves_cognition.rds")))
g <- w[year %in% c(2011, 2015, 2018), .(ID = IDn, year = as.numeric(year),
        mem = as.numeric(imrc) + as.numeric(dlrc),
        age = as.numeric(age), edu = edu_group, cesd = as.numeric(cesd10),
        adl = as.numeric(adl), female = as.numeric(female), wt = as.numeric(weight))]
say("available years in five wave data: ", paste(sort(unique(w$year)), collapse = ", "))
say("years used here: ", paste(sort(unique(g$year)), collapse = ", "),
    "| row count ", nrow(g))

## ---------------- 2. education specific cutoff(2011 definition, same F14)---------------- 
d11 <- g[year == 2011 & !is.na(mem) & !is.na(edu) & age >= 50]
d11[, g3 := fcase(edu == LV[1], G3[1], edu %in% LV[2:3], G3[2], edu == LV[4], G3[3])]
cut3 <- d11[!is.na(g3), .(n = .N, mean = mean(mem), sd = sd(mem)), by = g3][order(g3)]
cut3[, cutoff := mean - sd]
say("\n========== 1. education specific cutoff(2011 wave, mean-1SD)==========")
pt(cut3[, .(strata = g3, n, mean = round(mean, 2), SD = round(sd, 2), cutoff = round(cutoff, 2))])
fwrite(cut3, file.path(OUT, "F15_cutoffs.csv"), bom = TRUE)

## ---------------- 3. construct window dataset ---------------- 
build <- function(y0, y1) {
  a <- g[year == y0, .(ID, mem0 = mem, age0 = age, edu0 = edu, cesd0 = cesd,
                       adl0 = adl, wt0 = wt, female0 = female)]
  b <- g[year == y1, .(ID, mem1 = mem)]
  d <- merge(a, b, by = "ID", all = FALSE)
  d[, g3 := fcase(edu0 == LV[1], G3[1], edu0 %in% LV[2:3], G3[2], edu0 == LV[4], G3[3])]
  d <- merge(d, cut3[, .(g3, cutoff)], by = "g3", all.x = TRUE)
  d[, edu_f := factor(edu0, levels = LV)]
  d[, score := 2L * pmax(floor((age0 - 45) / 5), 0) +
        (-2L) * as.integer(edu_f == LV[2]) + (-7L) * as.integer(edu_f == LV[3]) +
        (-10L) * as.integer(edu_f == LV[4]) + 3L * as.integer(cesd0 >= 10)]
  d[, lp := B[["intercept"]] + B[["age"]] * age0 + B[["dep"]] * as.integer(cesd0 >= 10) +
        B[["adl"]] * adl0 + B[["eduPri"]] * as.integer(edu_f == LV[2]) +
        B[["eduMid"]] * as.integer(edu_f == LV[3]) + B[["eduHigh"]] * as.integer(edu_f == LV[4])]
  d[, phat := plogis(lp)]
  d[, event := as.integer(mem1 < cutoff)]
  d[, low0 := as.integer(mem0 < cutoff)]
  d[]
}
w1 <- build(2011, 2015)
w2 <- build(2015, 2018)

# same inclusion criteria as y04 national weighted distribution: age >= 45(original file includes spouses/child rows, ages may be as low as 10 years old)
core <- function(d) !is.na(d$score) & !is.na(d$lp) & !is.na(d$event) & !is.na(d$age0) & d$age0 >= 45
getc <- function(nm) { d <- get(nm, envir = .GlobalEnv); d[core(d)] }

## ---------------- 4. sample flow ---------------- 
flow <- data.table()
for (nm in c("w1", "w2")) {
  d <- get(nm)
  lab <- ifelse(nm == "w1", "2011 -> 2015(development cohort)", "2015 -> 2018(independent sample)")
  ok <- core(d)
  flow <- rbind(flow, data.table(
    window = lab,
    two_waves_ID_matchable = nrow(d),
    follow_up_memory_scorable = sum(!is.na(d$mem1)),
    follow_up_scorable_and_cutoff_available = sum(!is.na(d$mem1) & !is.na(d$cutoff)),
    three_item_score_calculable = sum(ok),
    because_younger_than_45_years_excluded = sum(!is.na(d$score) & !is.na(d$lp) & !is.na(d$event) &
                             !is.na(d$age0) & d$age0 < 45),
    excluded_for_missing_age = sum(!is.na(d$score) & !is.na(d$lp) & !is.na(d$event) & is.na(d$age0)),
    baseline_memory_scorable = sum(!is.na(d$mem0)),
    U_full_population = sum(ok),
    U_including_missing_baseline_memory = sum(ok & is.na(d$mem0)),
    N_unimpaired_at_baseline = sum(ok & d$low0 == 0, na.rm = TRUE),
    I_impaired_at_baseline = sum(ok & d$low0 == 1, na.rm = TRUE)))
}
say("\n========== 2. sample flow ==========")
pt(flow)
fwrite(flow, file.path(OUT, "F15_flow.csv"), bom = TRUE)

## ---------------- 5. AUC main table ---------------- 
aucrow <- function(d, v, win, arm, vlab) {
  ok <- !is.na(d[[v]]) & !is.na(d$event)
  r <- roc(d$event[ok], d[[v]][ok], quiet = TRUE, direction = "<")
  ci <- as.numeric(ci.auc(r))
  ev <- sum(d$event[ok] == 1); n <- sum(ok)
  data.table(window = win, definition = arm, predictors = vlab, n = n, events = ev,
             event_rate = round(100 * ev / n, 2),
             AUC = round(ci[2], 4), AUC_low = round(ci[1], 4), AUC_high = round(ci[3], 4))
}
armset <- function(d, win, out) {
  ok <- core(d)
  U <- d[ok]
  N <- d[ok & !is.na(low0) & low0 == 0]
  I <- d[ok & !is.na(low0) & low0 == 1]
  Bm <- d[ok & !is.na(low0)]
  Mu <- d[ok & is.na(low0)]
  sets <- list(U_full_population = U, N_unimpaired_at_baseline = N, I_impaired_at_baseline = I,
               M_missing_baseline_memory = Mu, B_baseline_memory_available = Bm)
  for (s in names(sets)) {
    if (nrow(sets[[s]]) == 0) next
    out <- rbind(out, aucrow(sets[[s]], "score", win, s, "three item scorecard"))
    out <- rbind(out, aucrow(sets[[s]], "phat", win, s, "continuous model"))
  }
  out
}
W1L <- "2011 -> 2015(development cohort)"; W2L <- "2015 -> 2018(independent sample)"
aucT <- data.table()
aucT <- armset(w1, W1L, aucT)
aucT <- armset(w2, W2L, aucT)
say("\n========== 3. each definition AUC(frozen model, four predictors)==========")
pt(aucT)
fwrite(aucT, file.path(OUT, "F15_auc.csv"), bom = TRUE)

## ---------------- 6. decomposition: predictive value of baseline low score alone ---------------- 
say("\n========== 4. decomposition: baseline low score indicator vs scorecard ==========")
decT <- data.table()
for (nm in c("w1", "w2")) {
  d <- getc(nm)[!is.na(low0)]
  lab <- ifelse(nm == "w1", W1L, W2L)
  r_low <- roc(d$event, d$low0, quiet = TRUE, direction = "<")
  r_scr <- roc(d$event, d$score, quiet = TRUE, direction = "<")
  fit <- glm(event ~ score + low0, data = d, family = binomial())
  pp <- predict(fit, type = "response")
  r_cmb <- roc(d$event, pp, quiet = TRUE, direction = "<")
  dt1 <- roc.test(r_scr, r_low, paired = TRUE)
  dt2 <- roc.test(r_cmb, r_scr, paired = TRUE)
  decT <- rbind(decT, data.table(
    window = lab, n = nrow(d), events = sum(d$event),
    fraction_of_events_already_low_scoring_at_baseline = round(100 * mean(d$low0[d$event == 1]), 1),
    AUC_baseline_low_score_only = round(as.numeric(auc(r_low)), 4),
    AUC_scorecard_only = round(as.numeric(auc(r_scr)), 4),
    AUC_scorecard_plus_baseline_low_score = round(as.numeric(auc(r_cmb)), 4),
    P_baseline_indicator_outperforms_scorecard = signif(dt1$p.value, 3),
    P_combined_approach_outperforms_scorecard = signif(dt2$p.value, 3)))
}
pt(decT)
fwrite(decT, file.path(OUT, "F15_decomposition.csv"), bom = TRUE)

## ---------------- 7. full population risk stratification ---------------- 
say("\n========== 5. full population risk stratification(score takes integer values)==========")
strT <- data.table()
for (nm in c("w1", "w2")) {
  d <- copy(getc(nm))
  lab <- ifelse(nm == "w1", W1L, W2L)
  d[, strt := cut(score, breaks = quantile(score, seq(0, 1, .25)), include.lowest = TRUE,
                  labels = c("Q1 lowest", "Q2", "Q3", "Q4 highest"))]
  tab <- d[, .(n = .N, event_rate = round(100 * mean(event), 1),
               including_baseline_low_scorers = round(100 * mean(low0 == 1, na.rm = TRUE), 1)),
           by = strt][order(strt)]
  say("-- ", lab, " --"); pt(tab)
  strT <- rbind(strT, cbind(window = lab, tab))
}
fwrite(strT, file.path(OUT, "F15_strata.csv"), bom = TRUE)

## ---------------- 8. full population screening simulation ---------------- 
say("\n========== 6. full population screening simulation(select score cutoff by fraction)==========")
scrT <- data.table()
for (nm in c("w1", "w2")) {
  d <- copy(getc(nm))
  lab <- ifelse(nm == "w1", W1L, W2L)
  s <- d$score; ev <- d$event; n <- length(ev)
  for (fr in c(0.05, 0.10, 0.20, 0.30)) {
    k <- floor(fr * n); th <- sort(s, decreasing = TRUE)[k]; pos <- s >= th
    sens <- sum(ev[pos]) / sum(ev); spec <- sum(!pos & ev == 0) / sum(ev == 0)
    ppv <- ifelse(sum(pos) > 0, mean(ev[pos]), NA_real_)
    npv <- mean(ev[!pos] == 0)
    sub <- d[pos]
    yld_new <- sum(sub$event == 1 & sub$low0 == 0, na.rm = TRUE)
    scrT <- rbind(scrT, data.table(
      window = lab, target_screening_fraction = paste0(fr * 100, "%"), threshold = th, n = length(ev), events = sum(ev), cases_detected = sum(pos), actual_screening_fraction = round(100 * mean(pos), 2), TP = sum(pos & ev == 1), FP = sum(pos & ev == 0), TN = sum(!pos & ev == 0), FN = sum(!pos & ev == 1),
      sensitivity = round(100 * sens, 1), specificity = round(100 * spec, 1),
      positive_predictive_value = round(100 * ppv, 1), negative_predictive_value = round(100 * npv, 1),
      fraction_already_low_scoring_at_baseline_among_detected_cases = round(100 * mean(sub$low0 == 1, na.rm = TRUE), 1),
      true_incident_events = yld_new, per_detected_1_incident_case_number_needed_to_screen = if (yld_new > 0) round(n / yld_new, 1) else NA_real_))
  }
  co <- coords(roc(ev, s, quiet = TRUE, direction = "<"), "best", best.method = "youden",
               ret = c("threshold", "sensitivity", "specificity"))
  co <- co[which.max(co$specificity), , drop = FALSE]
  pos <- s >= co$threshold
  sub <- d[pos]
  yld_new <- sum(sub$event == 1 & sub$low0 == 0, na.rm = TRUE)
  scrT <- rbind(scrT, data.table(
    window = lab, target_screening_fraction = "optimal Youden cutoff", threshold = co$threshold, n = length(ev), events = sum(ev), cases_detected = sum(pos), actual_screening_fraction = round(100 * mean(pos), 2), TP = sum(pos & ev == 1), FP = sum(pos & ev == 0), TN = sum(!pos & ev == 0), FN = sum(!pos & ev == 1),
    sensitivity = round(100 * co$sensitivity, 1), specificity = round(100 * co$specificity, 1),
    positive_predictive_value = round(100 * mean(ev[pos]), 1), negative_predictive_value = round(100 * mean(ev[!pos] == 0), 1),
    fraction_already_low_scoring_at_baseline_among_detected_cases = round(100 * mean(sub$low0 == 1, na.rm = TRUE), 1),
    true_incident_events = yld_new, per_detected_1_incident_case_number_needed_to_screen = if (yld_new > 0) round(n / yld_new, 1) else NA_real_))
}
pt(scrT)
fwrite(scrT, file.path(OUT, "F15_screening.csv"), bom = TRUE)

## ---------------- 9. stability across age groups(full population)---------------- 
say("\n========== 7. stability across age groups(full population definition)==========")
ageT <- data.table()
for (nm in c("w1", "w2")) {
  d <- copy(getc(nm))
  lab <- ifelse(nm == "w1", W1L, W2L)
  d[, ag := fcase(age0 < 60, "45-59", age0 < 70, "60-69", age0 < 80, "70-79",
                  default = "80+")]
  for (a in c("45-59", "60-69", "70-79", "80+")) {
    s <- d[ag == a]
    if (nrow(s) < 30 || length(unique(s$event)) < 2) next
    r <- roc(s$event, s$score, quiet = TRUE, direction = "<")
    ageT <- rbind(ageT, data.table(window = lab, age_group = a, n = nrow(s),
                                   events = sum(s$event),
                                   event_rate = round(100 * mean(s$event), 1),
                                   AUC_three_item_scorecard = round(as.numeric(auc(r)), 4)))
  }
}
pt(ageT)
fwrite(ageT, file.path(OUT, "F15_agegroup.csv"), bom = TRUE)

## ---------------- 10. calibration(full population vs unimpaired)---------------- 
say("\n========== 8. calibration ==========")
calT <- data.table()
for (nm in c("w1", "w2")) {
  d0 <- copy(getc(nm))
  lab <- ifelse(nm == "w1", W1L, W2L)
  for (s in c("U_full_population", "N_unimpaired_at_baseline")) {
    d <- if (s == "U_full_population") d0 else d0[!is.na(low0) & low0 == 0]
    fit <- glm(event ~ lp, data = d, family = binomial())
    calT <- rbind(calT, data.table(window = lab, definition = s, n = nrow(d),
                                   expected_event_rate = round(100 * mean(d$phat), 1),
                                   actual_event_rate = round(100 * mean(d$event), 1),
                                   calibration_slope = round(coef(fit)[2], 3),
                                   calibration_intercept = round(coef(fit)[1], 3),
                                   Brier = round(mean((d$phat - d$event)^2), 4)))
  }
}
pt(calT)
fwrite(calT, file.path(OUT, "F15_calibration.csv"), bom = TRUE)

## ---------------- 11. weighted sensitivity(full population)---------------- 
say("\n========== 9. weighted sensitivity(baseline weights, full population definition)==========")
wt_auc <- function(p, y, wgt) {
  cs <- y == 1; ct <- y == 0
  num <- 0
  for (i in which(cs)) num <- num + wgt[i] * sum(wgt[ct] * ((p[i] > p[ct]) + 0.5 * (p[i] == p[ct])))
  num / (sum(wgt[cs]) * sum(wgt[ct]))
}
wtT <- data.table()
for (nm in c("w1", "w2")) {
  d <- getc(nm)[!is.na(wt0)]
  lab <- ifelse(nm == "w1", W1L, W2L)
  wtT <- rbind(wtT, data.table(window = lab, n = nrow(d),
                               weighted_AUC_scorecard = round(wt_auc(d$score, d$event, d$wt0), 4),
                               unweighted_AUC_scorecard =
                                 round(as.numeric(auc(roc(d$event, d$score, quiet = TRUE, direction = "<"))), 4)))
}
pt(wtT)
fwrite(wtT, file.path(OUT, "F15_weighted.csv"), bom = TRUE)

## ---------------- 12. figure ---------------- 
## ---------------- 12b. stratification stability(within full population definition)---------------- 
say("\n========== 10. within stratum stability(full population definition)==========")
say("(notes: age and education are score components, within stratum AUC necessarily lower, used here to check directional consistency)")
strataT <- data.table()
for (nm in c("w1", "w2")) {
  d <- copy(getc(nm))
  lab <- ifelse(nm == "w1", W1L, W2L)
  d[, sex := fcase(female0 == 1, "female", female0 == 0, "male", default = NA_character_)]
  d[, edulab := fcase(g3 == G3[1], "no schooling", g3 == G3[2], "primary school/middle school",
                      g3 == G3[3], "high school or above", default = NA_character_)]
  for (sv in c("sex", "edulab")) {
    for (lv in na.omit(unique(d[[sv]]))) {
      s <- d[get(sv) == lv]
      if (nrow(s) < 30 || length(unique(s$event)) < 2) next
      r <- roc(s$event, s$score, quiet = TRUE, direction = "<")
      strataT <- rbind(strataT, data.table(
        window = lab, stratification_variable = ifelse(sv == "sex", "sex", "education"),
        strata = lv, n = nrow(s), event_rate = round(100 * mean(s$event), 1),
        AUC_three_item_scorecard = round(as.numeric(auc(r)), 4)))
    }
  }
  s <- d[age0 < 50]
  if (nrow(s) >= 30 && length(unique(s$event)) >= 2) {
    r <- roc(s$event, s$score, quiet = TRUE, direction = "<")
    strataT <- rbind(strataT, data.table(window = lab, stratification_variable = "minimum age",
                                         strata = "45-49", n = nrow(s),
                                          event_rate = round(100 * mean(s$event), 1),
                                          AUC_three_item_scorecard = round(as.numeric(auc(r)), 4)))
  }
  say("-- ", lab, ": age range ", min(d$age0), "-", max(d$age0),
      ", of which 45-49 years old ", sum(d$age0 < 50), " participants --")
}
pt(strataT)
fwrite(strataT, file.path(OUT, "F15_strata_within.csv"), bom = TRUE)

## ---------------- 13. figure ---------------- 
png(file.path(OUT, "F15_fig_auc.png"), width = 1800, height = 720, res = 140)
par(mfrow = c(1, 2), mar = c(5, 13, 3.5, 2))
for (nm in c("w1", "w2")) {
  lab <- ifelse(nm == "w1", "W1 2011->2015", "W2 2015->2018")
  sub <- aucT[window == ifelse(nm == "w1", W1L, W2L) & predictors == "three item scorecard"]
  sub <- sub[match(c("U_full_population", "N_unimpaired_at_baseline", "I_impaired_at_baseline", "M_missing_baseline_memory"), sub$definition)]
  bp <- barplot(rev(sub$AUC), names.arg = rev(sub$definition), horiz = TRUE, las = 1,
                xlim = c(0.4, 1), col = "steelblue", main = lab, xlab = "AUC")
  text(rev(sub$AUC), bp, labels = sprintf("%.3f", rev(sub$AUC)), pos = 4, cex = 0.9)
  abline(v = 0.5, lty = 3)
}
dev.off()

close(con)
cat("\n complete ->", file.path(OUT, "F15_report.txt"), "\n")
