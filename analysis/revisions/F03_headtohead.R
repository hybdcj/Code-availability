# ============================================================================
# F03 full version: temporal validation(A)＋ direct comparison of two weights in identical participants(B)
# population: CHARLS 2018 baseline → 2020 outcome; 2018 unimpaired participants enter risk population
# two outcome definitions: 
#   D1 absolute cutoff: 2018 memory score mean−1SD rounding(≤3), 2020 retain same absolute score
#   D2 equipercentile anchoring: take D1 cutoff in 2018 distribution percentile p, 2020 use same percentile of own distribution
# two weight sets: frozen(2011→2015, A14)vs reestimate(2018→2020 current fit)
# optimistic bias: Harrell bootstrap adjustment ＋ nested bootstrap of ΔAUC interval
# output: C:/y03tmp/out/F03_report.txt, F03_*.csv, F03_*.png
# ============================================================================
invisible(try(Sys.setlocale("LC_ALL", "zh_CN.UTF-8"), silent = TRUE))
suppressPackageStartupMessages({ library(data.table); library(pROC) })
if (dir.exists("D:/R/library")) .libPaths(c("D:/R/library", .libPaths()))
TMP <- "C:/y03tmp/tmp"; OUT <- "C:/y03tmp/out"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
con <- file(file.path(OUT, "F03_report.txt"), open = "wt", encoding = "UTF-8")
say <- function(...) { s <- paste0(..., collapse = ""); cat(s, "\n"); cat(s, "\n", file = con) }
pt  <- function(x) say(paste(capture.output(print(x)), collapse = "\n"))

BETA_OLD <- c("(Intercept)" = -4.19524508778947, age = 0.0492039250378111,
              eduP = -0.290105771000328, eduM = -0.994711800202235, eduH = -1.43588388035287,
              dep = 0.412166484954086, adl = -0.0110119344848014)
LV <- c("Less than primary", "Primary school", "Middle school", "High school or above")

## ---------------- dataset ----------------
w <- as.data.table(readRDS(file.path(TMP, "waves_cognition.rds")))
a18 <- w[year == 2018, .(ID = IDn, mem18 = as.numeric(imrc) + as.numeric(dlrc),
                         age18 = as.numeric(age), edu18 = edu_group,
                         cesd18 = as.numeric(cesd10), adl18 = as.numeric(adl),
                         wt18 = as.numeric(weight), comm18 = as.character(communityID))]
a20 <- w[year == 2020, .(ID = IDn, mem20 = as.numeric(imrc) + as.numeric(dlrc),
                         wt20 = as.numeric(weight))]
m <- merge(a18, a20, by = "ID")
m <- m[!is.na(mem18) & !is.na(mem20) & !is.na(age18) & !is.na(edu18) &
         !is.na(cesd18) & !is.na(adl18)]
m[, dep := as.integer(cesd18 >= 10)]
m[, edu_group := factor(edu18, levels = LV)]
say("analysis population n = ", nrow(m), "(of which 2018 year <45 years old: ", sum(m$age18 < 45), ")")
m <- m[age18 >= 45]
say("restrict 2018 year ≥45 after age restriction n = ", nrow(m))

## ---------------- outcome definition ----------------
cut18 <- floor(mean(m$mem18) - sd(m$mem18))
p_anchor <- mean(m$mem18 <= cut18)
cut20 <- as.numeric(quantile(m$mem20, p_anchor, type = 1))
say("\n2018 memory score mean=", round(mean(m$mem18), 2), " SD=", round(sd(m$mem18), 2),
    " → D1 cutoff = ≤", cut18, "(this cutoff in 2018 distribution percentile = ", round(100 * p_anchor, 1), "%)")
say("2020 memory score mean=", round(mean(m$mem20), 2), " SD=", round(sd(m$mem20), 2),
    " → D2 equipercentile cutoff = ≤", cut20)

defs <- list(D1 = list(k18 = cut18, k20 = cut18, lab = paste0("D1 absolute cutoff ≤", cut18)),
             D2 = list(k18 = cut18, k20 = cut20, lab = paste0("D2 equipercentile anchoring(2020 ≤", cut20, ")")))

## ---------------- helper functions ----------------
pred_old <- function(d) {
  plogis(BETA_OLD[["(Intercept)"]] + BETA_OLD[["age"]] * d$age18 + BETA_OLD[["dep"]] * d$dep +
           BETA_OLD[["adl"]] * d$adl18 +
           BETA_OLD[["eduP"]] * as.integer(d$edu_group == LV[2]) +
           BETA_OLD[["eduM"]] * as.integer(d$edu_group == LV[3]) +
           BETA_OLD[["eduH"]] * as.integer(d$edu_group == LV[4]))
}
fit_new <- function(d) glm(event ~ age18 + edu_group + dep + adl18, data = d, family = binomial())
calc <- function(d, p, tag = "") {
  r <- roc(d$event, p, quiet = TRUE, direction = "<")
  ci <- as.numeric(ci.auc(r))
  lp <- qlogis(pmin(pmax(p, 1e-9), 1 - 1e-9))
  cal <- glm(d$event ~ lp, family = binomial())
  data.table(label = tag, n = nrow(d), events = sum(d$event),
             AUC = round(ci[2], 3), AUC_low = round(ci[1], 3), AUC_high = round(ci[3], 3),
             calibration_slope = round(coef(cal)[2], 3), calibration_intercept = round(coef(cal)[1], 3),
             Brier = round(mean((p - d$event)^2), 4),
             OE_ratio = round(sum(d$event) / sum(p), 2))
}
wt_auc <- function(p, y, wgt) {
  cs <- y == 1; ct <- y == 0
  num <- 0
  for (i in which(cs)) num <- num + wgt[i] * sum(wgt[ct] * ((p[i] > p[ct]) + 0.5 * (p[i] == p[ct])))
  num / (sum(wgt[cs]) * sum(wgt[ct]))
}
nb_curve <- function(p, y, th) sapply(th, function(t) sum(p >= t & y == 1)/length(y) -
                                        sum(p >= t & y == 0)/length(y) * (t/(1-t)))

## ---------------- A: temporal validation of frozen model ----------------
say("\n========== A. frozen model in 2018→2020 temporal validation of ==========")
A <- data.table(); stA <- NULL; d1 <- NULL
for (nm in names(defs)) {
  dfn <- defs[[nm]]
  d <- copy(m[mem18 > dfn$k18])[, event := as.integer(mem20 <= dfn$k20)]
  d[, phat := pred_old(d)]
  A <- rbind(A, calc(d, d$phat, dfn$lab))
  if (nm == "D1") {
    d1 <- copy(d)
    d1[, score := 2L * pmax(floor((age18 - 45)/5), 0) +
          (-2L) * as.integer(edu_group == LV[2]) + (-7L) * as.integer(edu_group == LV[3]) +
          (-10L) * as.integer(edu_group == LV[4]) + 3L * dep]
    stA <- d1[, .(n = .N, events = sum(event), event_rate = paste0(round(100*mean(event), 1), "%"),
                  mean_score = round(mean(score), 1)),
              by = .(risk_stratum = cut(score, c(-Inf, -3.5, 4.5, Inf), labels = c("low (≤-3)", "in (-2~+4)", "high (≥+5)")))][order(risk_stratum)]
  }
}
pt(A)
say("\n frozen score risk stratification(", defs$D1$lab, "): "); pt(stA)
say("(comparison: apparent development sample performance AUC 0.702, Bootstrap adjustment 0.698; this period n=", nrow(d1),
    ", events=", sum(d1$event), ", event rate=", round(100*mean(d1$event), 2), "%)")
fwrite(A, file.path(OUT, "F03_A_temporal_validation.csv"), bom = TRUE)
fwrite(stA, file.path(OUT, "F03_A_frozen_strata.csv"), bom = TRUE)

## ---------------- B: direct comparison of two weights in identical participants ----------------
say("\n========== B. direct comparison of two weights in identical participants ==========")
Bset <- 2000; set.seed(20261014)
B <- data.table(); DCAtab <- data.table(); stAll <- NULL; dtab <- NULL
for (nm in names(defs)) {
  dfn <- defs[[nm]]
  d <- copy(m[mem18 > dfn$k18])[, event := as.integer(mem20 <= dfn$k20)]
  d[, phat_old := pred_old(d)]
  fit <- fit_new(d)
  d[, phat_new := predict(fit, type = "response")]
  say("\n---- ", dfn$lab, "(n=", nrow(d), ", events=", sum(d$event), ")----")
  say("reestimated coefficients(2018→2020): "); pt(round(coef(fit), 4))

  rowA <- calc(d, d$phat_old, "old weights(2011 fit)")
  rowB <- calc(d, d$phat_new, "new weights(2018→2020 fit, apparent)")
  r_old <- roc(d$event, d$phat_old, quiet = TRUE, direction = "<")
  r_new <- roc(d$event, d$phat_new, quiet = TRUE, direction = "<")
  pdel <- tryCatch(roc.test(r_old, r_new, paired = TRUE)$p.value, error = function(e) NA_real_)

  opt <- rep(NA_real_, Bset); dA <- rep(NA_real_, Bset)
  for (b in seq_len(Bset)) {
    idx <- sample.int(nrow(d), replace = TRUE); db <- d[idx]
    f <- tryCatch(glm(event ~ age18 + edu_group + dep + adl18, data = db, family = binomial()),
                  error = function(e) NULL)
    if (is.null(f) || length(unique(db$event)) < 2) next
    p_orig <- predict(f, newdata = d, type = "response")
    p_boot <- predict(f, type = "response")
    opt[b] <- as.numeric(auc(roc(db$event, p_boot, quiet = TRUE, direction = "<"))) -
              as.numeric(auc(roc(d$event, p_orig, quiet = TRUE, direction = "<")))
    dA[b]  <- as.numeric(auc(roc(db$event, p_boot, quiet = TRUE, direction = "<"))) -
              as.numeric(auc(roc(db$event, pred_old(db), quiet = TRUE, direction = "<")))
  }
  auc_new_corrected <- as.numeric(auc(r_new)) - mean(opt, na.rm = TRUE)
  say("\n old weights AUC=", round(rowA$AUC, 3), " (", round(rowA$AUC_low, 3), "-", round(rowA$AUC_high, 3),
      ")｜calibration slope=", rowA$calibration_slope, " intercept=", rowA$calibration_intercept, "｜Brier=", rowA$Brier)
  say("new weights AUC=", round(rowB$AUC, 3), " (", round(rowB$AUC_low, 3), "-", round(rowB$AUC_high, 3),
      ")｜calibration slope=", rowB$calibration_slope, " intercept=", rowB$calibration_intercept, "｜Brier=", rowB$Brier)
  say("new weights Harrell adjustment AUC=", round(auc_new_corrected, 3), "(optimism ", round(mean(opt, na.rm = TRUE), 4), ")")
  say("apparent ΔAUC=", round(rowB$AUC - rowA$AUC, 4), "｜paired DeLong p=", format(pdel, digits = 3))
  say("nested bootstrap ΔAUC median=", round(median(dA, na.rm = TRUE), 4),
      ", 95% interval=[", round(quantile(dA, .025, na.rm = TRUE), 4), ", ",
      round(quantile(dA, .975, na.rm = TRUE), 4), "]")
  yo <- function(r) { co <- coords(r, "best", best.method = "youden",
                                   ret = c("sensitivity", "specificity"))
                      data.table(diagnostic_sensitivity = round(co$sensitivity, 3), specificity = round(co$specificity, 3),
                                 Youden = round(co$sensitivity + co$specificity - 1, 3)) }
  say("\n Youden cutoff: old weights"); pt(yo(r_old)); say("Youden cutoff: new weights"); pt(yo(r_new))

  ## integerize new weights(retain A14 rules)
  bb <- coef(fit)
  tb <- c(age5 = 5 * bb[["age18"]], eduP = bb[["edu_groupPrimary school"]],
          eduM = bb[["edu_groupMiddle school"]], eduH = bb[["edu_groupHigh school or above"]],
          dep = bb[["dep"]], adl = bb[["adl18"]])
  pts <- round(10 * tb / max(abs(tb)))
  say("\n rescale new weights β: "); pt(round(tb, 3)); say("integer points from new weights: "); pt(pts)
  d[, score_new := as.integer(pts[["age5"]]) * pmax(floor((d$age18 - 45)/5), 0) +
        as.integer(pts[["eduP"]]) * as.integer(d$edu_group == LV[2]) +
        as.integer(pts[["eduM"]]) * as.integer(d$edu_group == LV[3]) +
        as.integer(pts[["eduH"]]) * as.integer(d$edu_group == LV[4]) +
        as.integer(pts[["dep"]]) * d$dep]
  d[, score_old := 2L * pmax(floor((d$age18 - 45)/5), 0) +
        (-2L) * as.integer(d$edu_group == LV[2]) + (-7L) * as.integer(d$edu_group == LV[3]) +
        (-10L) * as.integer(d$edu_group == LV[4]) + 3L * d$dep]
  say("\n integer score AUC: old ", round(as.numeric(auc(roc(d$event, d$score_old, quiet = TRUE, direction = "<"))), 3),
      "｜new ", round(as.numeric(auc(roc(d$event, d$score_new, quiet = TRUE, direction = "<"))), 3))
  LVL <- c("low (≤-3)", "in (-2~+4)", "high (≥+5)")
  d[, strata_old := cut(score_old, c(-Inf, -3.5, 4.5, Inf), labels = LVL)]
  d[, strata_new := cut(score_new, c(-Inf, -3.5, 4.5, Inf), labels = LVL)]
  s_old <- d[, .(n = .N, events = sum(event), event_rate = round(100*mean(event), 1)), by = .(risk_stratum = strata_old)]
  s_old[, score := "old (2011)"]
  s_new <- d[, .(n = .N, events = sum(event), event_rate = round(100*mean(event), 1)), by = .(risk_stratum = strata_new)]
  s_new[, score := "new (2018-2020)"]
  stk <- rbind(s_old[, .(score, risk_stratum, n, events, event_rate)][order(risk_stratum)],
               s_new[, .(score, risk_stratum, n, events, event_rate)][order(risk_stratum)])
  pt(stk)
  say("\n by“high risk (≥+5)”classify positive: ")
  for (sc in c("score_old", "score_new")) {
    pos <- d[[sc]] >= 5
    say("  ", sc, ": positive n=", sum(pos), "(", round(100*mean(pos), 1), "%)",
        "｜sensitivity=", round(100*sum(d$event[pos])/sum(d$event), 1), "%",
        "｜specificity=", round(100*sum(!pos & d$event == 0)/sum(d$event == 0), 1), "%",
        "｜event rate among screen positives=", round(100*mean(d$event[pos]), 1), "%")
  }
  ## model updating hierarchy: ① unchanged ② intercept only recalibration(1 parameters)③ slope＋intercept recalibration(2 parameters)④ fully reestimated
  lp_old <- qlogis(pmin(pmax(d$phat_old, 1e-9), 1 - 1e-9))
  cal_int <- glm(event ~ 1, offset = lp_old, data = d, family = binomial())
  d[, phat_int := plogis(lp_old + coef(cal_int)[1])]
  rowC <- calc(d, d$phat_int, "frozen＋intercept only recalibration")
  cal_sl <- glm(event ~ lp_old, data = d, family = binomial())
  d[, phat_slope := predict(cal_sl, type = "response")]
  rowD <- calc(d, d$phat_slope, "frozen＋slope and intercept recalibration")
  say("\n intercept only recalibration: shift=", round(coef(cal_int)[1], 3), "｜Brier=", rowC$Brier,
      "｜AUC unchanged(monotonic transformation)=", rowC$AUC)
  say("slope＋intercept recalibration: slope=", round(coef(cal_sl)[2], 3), " intercept=", round(coef(cal_sl)[1], 3),
      "｜Brier=", rowD$Brier, "｜AUC(monotonic transformation)=", rowD$AUC)

  dw <- d[!is.na(wt18)]
  w_old <- wt_auc(dw$phat_old, dw$event, dw$wt18)
  w_int <- wt_auc(dw$phat_int, dw$event, dw$wt18)
  w_slope <- wt_auc(dw$phat_slope, dw$event, dw$wt18)
  w_new <- wt_auc(dw$phat_new, dw$event, dw$wt18)
  say("\n weighted sensitivity(2018 baseline weights, available n=", nrow(dw), "): unchanged=", round(w_old, 3),
      "｜intercept only recalibration=", round(w_int, 3), "｜slope and intercept recalibration=", round(w_slope, 3),
      "｜fully reestimated=", round(w_new, 3))

  th <- seq(0.02, 0.30, 0.02)
  DCAtab <- rbind(DCAtab, data.table(definition = dfn$lab, threshold = th,
                                     frozen_unchanged = nb_curve(d$phat_old, d$event, th),
                                     frozen_intercept_only = nb_curve(d$phat_int, d$event, th),
                                     frozen_slope_and_intercept = nb_curve(d$phat_slope, d$event, th),
                                     fully_reestimated = nb_curve(d$phat_new, d$event, th)))
  stAll <- rbind(stAll, stk)
  nr <- data.table(
    definition = dfn$lab,
    weights = c("① frozen unchanged(2011 fit)", "② frozen＋intercept only recalibration", "③ frozen＋slope and intercept recalibration",
             "④ fully reestimated(2018-2020 fit, apparent)", "④ fully reestimated(Harrell adjustment)"),
    AUC = c(rowA$AUC, rowC$AUC, rowD$AUC, rowB$AUC, round(auc_new_corrected, 3)),
    calibration_slope = c(rowA$calibration_slope, rowC$calibration_slope, rowD$calibration_slope, rowB$calibration_slope, NA_real_),
    calibration_intercept = c(rowA$calibration_intercept, rowC$calibration_intercept, rowD$calibration_intercept, rowB$calibration_intercept, NA_real_),
    Brier = c(rowA$Brier, rowC$Brier, rowD$Brier, rowB$Brier, NA_real_),
    weighted_AUC = c(round(w_old, 3), round(w_int, 3), round(w_slope, 3), round(w_new, 3), NA_real_),
    dAUC = c(NA_real_, 0, round(rowD$AUC - rowA$AUC, 4), round(rowB$AUC - rowA$AUC, 4), NA_real_),
    DeLong_p = c(NA_real_, NA_real_, NA_real_, round(pdel, 4), NA_real_),
    dAUC_boot_low = c(NA_real_, NA_real_, NA_real_, round(quantile(dA, .025, na.rm = TRUE), 4), NA_real_),
    dAUC_boot_high = c(NA_real_, NA_real_, NA_real_, round(quantile(dA, .975, na.rm = TRUE), 4), NA_real_))
  dtab <- rbind(dtab, nr)
}
fwrite(dtab, file.path(OUT, "F03_B_headtohead.csv"), bom = TRUE)
fwrite(stAll, file.path(OUT, "F03_B_integer_strata.csv"), bom = TRUE)
fwrite(DCAtab, file.path(OUT, "F03_B_dca.csv"), bom = TRUE)

## ---------------- figure ----------------
d <- copy(m[mem18 > defs$D1$k18])[, event := as.integer(mem20 <= defs$D1$k20)]
d[, phat_old := pred_old(d)]
d[, phat_new := predict(fit_new(d), type = "response")]
d[, lp_old := qlogis(pmin(pmax(phat_old, 1e-9), 1 - 1e-9))]
d[, phat_int := plogis(lp_old + coef(glm(event ~ 1, offset = lp_old, family = binomial()))[1])]

png(file.path(OUT, "F03_fig_DCA.png"), width = 1100, height = 850, res = 120)
th <- seq(0, 0.4, 0.01)
plot(th, nb_curve(d$phat_old, d$event, th), type = "l", lwd = 2, col = "grey40", lty = 2,
     xlab = "Threshold probability", ylab = "Net benefit", ylim = c(-0.02, 0.10),
     main = "DCA: frozen (2011) vs refit (2018-2020) weights")
lines(th, nb_curve(d$phat_int, d$event, th), lwd = 2, col = "black")
lines(th, nb_curve(d$phat_new, d$event, th), lwd = 2, col = "red", lty = 1)
abline(h = 0, lty = 3)
legend("topright", c("Frozen 2011 weights (as published)", "Frozen + intercept recalibration",
                     "Refit 2018-2020 weights", "Treat none"),
       col = c("grey40", "black", "red", "black"), lty = c(2, 1, 1, 3), lwd = c(2, 2, 2, 1), bty = "n")
dev.off()

png(file.path(OUT, "F03_fig_calibration.png"), width = 1500, height = 520, res = 120)
par(mfrow = c(1, 3))
PVS  <- c(phat_old = "1) Frozen 2011 (as published)", phat_int = "2) Frozen + intercept recalibration",
          phat_new = "4) Refit 2018-2020")
for (pv in names(PVS)) {
  br <- unique(quantile(d[[pv]], seq(0, 1, 0.1)))
  d[, grp := cut(d[[pv]], br, include.lowest = TRUE)]
  agg <- d[, .(pred = mean(d[[pv]]), obs = mean(event), n = .N), by = grp][order(pred)]
  plot(agg$pred, agg$obs, pch = 19, cex = 1.2, xlim = c(0, 0.20), ylim = c(0, 0.20),
       xlab = "Predicted risk", ylab = "Observed risk", main = PVS[[pv]])
  abline(0, 1, lty = 3)
  text(agg$pred, agg$obs, labels = agg$n, pos = 3, cex = 0.7)
}
dev.off()

close(con)
cat("\n complete ->", file.path(OUT, "F03_report.txt"), "\n")
