# F26 LEGACY: BDSI implementation/input is a historical version; for corrected analyses run F26_bdsi_audit.R.
# original numerical results retained for traceability, cannot directly serve as original model BDSI valid comparison; Tu see also limitations in F25.
# ============================================================================
# F14: direct comparison with established models(including sensitivity, significance, sensitivity/specificity)
# three arms: ① our three item scorecard  ② BDSI(Barnes 2014, 7 items)
#       ③ Tu 2026 six factor model(including memory＋orientation＋grip strength)
# window: CHARLS 2015 → 2018; outcome: memory score falls to[education strata mean−1SD cutoff]below
# output: C:/y03tmp/out/F14_report.txt, F14_*.csv
# ============================================================================
invisible(try(Sys.setlocale("LC_ALL", "zh_CN.UTF-8"), silent = TRUE))
suppressPackageStartupMessages({ library(data.table); library(pROC); library(haven) })
if (dir.exists("D:/R/library")) .libPaths(c("D:/R/library", .libPaths()))
TMP <- "C:/y03tmp/tmp"; OUT <- "C:/y03tmp/out"
con <- file(file.path(OUT, "F14_report.txt"), open = "wt", encoding = "UTF-8")
say <- function(...) { s <- paste0(..., collapse = ""); cat(s, "\n"); cat(s, "\n", file = con) }
pt  <- function(x) say(paste(capture.output(print(x)), collapse = "\n"))

B_TU <- c(age_lt65 = -0.53, age_ge75 = 0.26, male = 0.38,
          edu_primid = 0.45, edu_none = 0.03, grip = -0.16, memory = -0.33, orient = -0.15)
LV <- c("Less than primary", "Primary school", "Middle school", "High school or above")

## ---------------- 1. cognition and covariates ----------------
w <- as.data.table(readRDS(file.path(TMP, "waves_cognition.rds")))
g <- w[year %in% c(2011, 2015, 2018), .(ID = IDn, year,
                                        mem = as.numeric(imrc) + as.numeric(dlrc),
                                        orient = as.numeric(orient), age = as.numeric(age),
                                        female = as.numeric(female), edu_group = edu_group,
                                        cesd10 = as.numeric(cesd10), adl = as.numeric(adl),
                                        iadl2 = as.numeric(iadl2),
                                        stroke = as.numeric(stroke), diab = as.numeric(diab),
                                        weight = as.numeric(weight))]
d15 <- g[year == 2015, .(ID, mem15 = mem, orient15 = orient, age15 = age, female15 = female,
                         edu15 = edu_group, cesd15 = cesd10, adl15 = adl, iadl2_15 = iadl2,
                         stroke15 = stroke, diab15 = diab, wt15 = weight)]
d18 <- g[year == 2018, .(ID, mem18 = mem)]

## ---------------- 2. physical measurements and grip strength(including sentinel value cleaning)----------------
b <- as.data.table(read_dta(file.path(TMP, "w3_biomarker.dta"),
                            col_select = any_of(c("ID", "qi002", "ql002", "qc003", "qc004", "qc005", "qc006"))))
b[, `:=`(height = fifelse(qi002 >= 120 & qi002 <= 220, qi002, NA_real_),
         weight = fifelse(ql002 >= 25 & ql002 <= 200, ql002, NA_real_))]
for (v in c("qc003", "qc004", "qc005", "qc006"))
  b[[v]] <- fifelse(b[[v]] > 0 & b[[v]] <= 100, b[[v]], NA_real_)
b[, `:=`(bmi = weight / (height / 100)^2,
         grip_max = pmax(qc003, qc004, qc005, qc006, na.rm = TRUE))]
b[, grip_alt := rowMeans(cbind(pmax(qc003, qc004, na.rm = TRUE), pmax(qc005, qc006, na.rm = TRUE)),
                         na.rm = TRUE)]
b <- b[, .(ID = as.character(ID), bmi, grip_max, grip_alt)]
say("physical measurement data: height available ", sum(!is.na(b$bmi)), " participants(BMI median ",
    round(median(b$bmi, na.rm = TRUE), 2), ")")

## ---------------- 3. education specific cutoff(2011 reconstruct, same F07)----------------
d11 <- g[year == 2011 & !is.na(mem) & !is.na(edu_group) & age >= 50]
d11[, edu3 := fcase(edu_group == LV[1], "No formal schooling",
                    edu_group %in% LV[2:3], "Primary/middle",
                    edu_group == LV[4], "High school and above")]
cut3 <- d11[, .(n = .N, cutoff = mean(mem) - sd(mem)), by = edu3]

## ---------------- 4. analysis sample ----------------
mm <- merge(merge(d15, d18, by = "ID"), b, by = "ID", all.x = TRUE)
mm[, edu3 := fcase(edu15 == LV[1], "No formal schooling",
                   edu15 %in% LV[2:3], "Primary/middle",
                   edu15 == LV[4], "High school and above")]
mm <- merge(mm, cut3[, .(edu3, cutoff)], by = "edu3", all.x = TRUE)
say("\n sample flow: 2015+2018 both memory assessments scorable ", nrow(mm[!is.na(mem15) & !is.na(mem18)]),
    " → ≥50 years old ", nrow(mm[!is.na(mem15) & !is.na(mem18) & age15 >= 50]),
    " → complete variables for three arms ", nrow(mm[!is.na(mem15) & !is.na(mem18) & age15 >= 50 & !is.na(cutoff) &
                                   !is.na(bmi) & !is.na(grip_max) & !is.na(orient15) &
                                   !is.na(cesd15) & !is.na(iadl2_15) & !is.na(stroke15) & !is.na(diab15)]))
d <- mm[!is.na(mem15) & !is.na(mem18) & age15 >= 50 & !is.na(cutoff) & !is.na(bmi) &
          !is.na(grip_max) & !is.na(orient15) & !is.na(cesd15) & !is.na(iadl2_15) &
          !is.na(stroke15) & !is.na(diab15)]
d <- d[mem15 >= cutoff]                        # 2015 normal memory
d[, event := as.integer(mem18 < cutoff)]
say("analysis sample n = ", nrow(d), "｜events = ", sum(d$event), "(", round(100 * mean(d$event), 2), "%)")

## ---------------- 5. three arm scores ----------------
d[, edu_f := factor(edu15, levels = LV)]
d[, male := as.integer(female15 == 0)]
## ① our three item scorecard
d[, score_us := 2L * pmax(floor((age15 - 45) / 5), 0) +
      (-2L) * as.integer(edu_f == LV[2]) + (-7L) * as.integer(edu_f == LV[3]) +
      (-10L) * as.integer(edu_f == LV[4]) + 3L * as.integer(cesd15 >= 10)]
## ② BDSI(Barnes 2014): age 1 points/years old (65–79), <12 year education 9, stroke 6, diabetes 3, BMI<18.5 8, needs help managing finances/taking medication 10, depression 6
d[, `:=`(b_age = pmax(0, pmin(age15, 79) - 64),
         b_edu = 9L * as.integer(edu_f != LV[4]),      # <12 year ≈ below high school
         b_stroke = 6L * as.integer(stroke15 == 1),
         b_diab = 3L * as.integer(diab15 == 1),
         b_bmi = 8L * as.integer(bmi < 18.5),
         b_func = 10L * as.integer(iadl2_15 == 1),
         b_dep = 6L * as.integer(cesd15 >= 10))]
d[, score_bdsi := b_age + b_edu + b_stroke + b_diab + b_bmi + b_func + b_dep]
## ③ Tu 2026 six factor linear predictor
d[, lp_tu := B_TU[["age_lt65"]] * as.integer(age15 < 65) + B_TU[["age_ge75"]] * as.integer(age15 >= 75) +
      B_TU[["male"]] * male + B_TU[["edu_primid"]] * as.integer(edu3 == "Primary/middle") +
      B_TU[["edu_none"]] * as.integer(edu3 == "No formal schooling") +
      B_TU[["grip"]] * grip_max + B_TU[["memory"]] * mem15 + B_TU[["orient"]] * orient15]

say("\nBDSI score distribution: mean ", round(mean(d$score_bdsi), 1), " median ", median(d$score_bdsi),
    " range ", min(d$score_bdsi), "-", max(d$score_bdsi))
say("BDSI means of component items: age ", round(mean(d$b_age), 2), "｜education ", round(mean(d$b_edu), 2),
    "｜stroke ", round(mean(d$b_stroke), 2), "｜diabetes ", round(mean(d$b_diab), 2),
    "｜BMI ", round(mean(d$b_bmi), 2), "｜function ", round(mean(d$b_func), 2),
    "｜depression ", round(mean(d$b_dep), 2))
say("our three item card: mean ", round(mean(d$score_us), 1), " range ", min(d$score_us), "-", max(d$score_us))

## ---------------- 6. primary results: AUC and pairwise significance ----------------
arms <- c(score_us = "our three item scorecard", score_bdsi = "BDSI(7 items)", lp_tu = "Tu 2026(6 factor, includes cognitive testing)")
res <- data.table()
rocs <- list()
for (v in names(arms)) {
  r <- roc(d$event, d[[v]], quiet = TRUE, direction = "<")
  rocs[[v]] <- r
  ci <- as.numeric(ci.auc(r))
  res <- rbind(res, data.table(model = arms[[v]], AUC = round(ci[2], 4),
                               low = round(ci[1], 4), high = round(ci[3], 4)))
}
say("\n========== AUC ==========")
pt(res)
fwrite(res, file.path(OUT, "F14_auc.csv"), bom = TRUE)

say("\n========== paired comparison with our three item card ==========")
cmp <- data.table()
for (v in c("score_bdsi", "lp_tu")) {
  dt <- roc.test(rocs[["score_us"]], rocs[[v]], paired = TRUE)
  set.seed(20261015); B <- 2000; dv <- rep(NA_real_, B)
  for (i in seq_len(B)) {
    idx <- sample.int(nrow(d), replace = TRUE); db <- d[idx]
    if (length(unique(db$event)) < 2) next
    dv[i] <- as.numeric(auc(roc(db$event, db$score_us, quiet = TRUE, direction = "<"))) -
             as.numeric(auc(roc(db$event, db[[v]],  quiet = TRUE, direction = "<")))
  }
  cmp <- rbind(cmp, data.table(comparison = arms[[v]],
                               this_study_AUC = round(as.numeric(auc(rocs[["score_us"]])), 4),
                               comparison_AUC = round(as.numeric(auc(rocs[[v]])), 4),
                               ΔAUC = round(as.numeric(auc(rocs[["score_us"]])) - as.numeric(auc(rocs[[v]])), 4),
                               DeLong_P = signif(dt$p.value, 3),
                               boot_low = round(quantile(dv, .025, na.rm = TRUE), 4),
                               boot_high = round(quantile(dv, .975, na.rm = TRUE), 4)))
}
pt(cmp)
fwrite(cmp, file.path(OUT, "F14_pairwise.csv"), bom = TRUE)

## ---------------- 7. sensitivity／specificity(target screening fraction(including ties, report actual fraction)+ Youden)----------------
oper <- function(v, lab) {
  s <- d[[v]]; ev <- d$event; n <- length(ev)
  out <- data.table()
  for (fr in c(0.10, 0.20, 0.30)) {
    k <- floor(fr * n); th <- sort(s, decreasing = TRUE)[k]
    pos <- s >= th
    out <- rbind(out, data.table(model = lab, definition = paste0("target screening ", fr * 100, "%"),
                                 threshold = th, n = length(ev), events = sum(ev), cases_detected = sum(pos), actual_screening_fraction = round(100 * mean(pos), 2), TP = sum(pos & ev == 1), FP = sum(pos & ev == 0), TN = sum(!pos & ev == 0), FN = sum(!pos & ev == 1), sensitivity = round(100 * sum(ev[pos]) / sum(ev), 1),
                                 specificity = round(100 * sum(!pos & ev == 0) / sum(ev == 0), 1),
                                 negative_predictive_value = round(100 * mean(ev[!pos] == 0), 1), event_rate_among_screen_positives = round(100 * mean(ev[pos]), 1)))
  }
  co <- coords(roc(ev, s, quiet = TRUE, direction = "<"), "best", best.method = "youden",
               ret = c("threshold", "sensitivity", "specificity"))
  co <- co[which.max(co$specificity), , drop = FALSE]
  pos <- s >= co$threshold
  out <- rbind(out, data.table(model = lab, definition = "optimal Youden cutoff",
                               threshold = round(co$threshold, 2), n = length(ev), events = sum(ev), cases_detected = sum(pos), actual_screening_fraction = round(100 * mean(pos), 2), TP = sum(pos & ev == 1), FP = sum(pos & ev == 0), TN = sum(!pos & ev == 0), FN = sum(!pos & ev == 1),
                               sensitivity = round(100 * co$sensitivity, 1),
                               specificity = round(100 * co$specificity, 1),
                               negative_predictive_value = round(100 * mean(ev[!pos] == 0), 1), event_rate_among_screen_positives = round(100 * mean(ev[pos]), 1)))
  out
}
op <- rbind(oper("score_us", "our three item scorecard"), oper("score_bdsi", "BDSI(7 items)"),
            oper("lp_tu", "Tu 2026(6 factor)"))
say("\n========== sensitivity／specificity ==========")
pt(op)
fwrite(op, file.path(OUT, "F14_operating.csv"), bom = TRUE)

## ---------------- 8. risk stratification(cohort specific quartiles)----------------
say("\n========== event rates by risk stratum(cohort specific quartiles)==========")
for (v in names(arms)) {
  q <- cut(d[[v]], quantile(d[[v]], seq(0, 1, .25)), include.lowest = TRUE, labels = 1:4)
  cat("", arms[[v]], ": ", paste(d[, .(r = paste0(round(100 * mean(event), 1), "%")), by = q][order(q)]$r,
                                collapse = "  →  "), "\n")
  say("  ", arms[[v]], ": ", paste(d[, .(r = paste0(round(100 * mean(event), 1), "%")), by = q][order(q)]$r,
                                   collapse = "  →  "))
}

## ---------------- 9. sensitivity analysis ----------------
say("\n========== sensitivity analysis ==========")
## (a) outcome uses this study global cutoff(≤3)
mu15 <- mean(mm$mem15, na.rm = TRUE); sd15 <- sd(mm$mem15, na.rm = TRUE); k15 <- floor(mu15 - sd15)
d2 <- copy(d)[, event2 := as.integer(mem18 <= k15)]
say("(a) outcome uses this study global cutoff ≤", k15, "(events ", sum(d2$event2), ", ",
    round(100 * mean(d2$event2), 2), "%)")
for (v in names(arms))
  say("    ", arms[[v]], " AUC = ",
      round(as.numeric(auc(roc(d2$event2, d2[[v]], quiet = TRUE, direction = "<"))), 4))
## (b) BDSI strict depression threshold CES-D≥16
d3 <- copy(d)[, score_bdsi_b := b_age + b_edu + b_stroke + b_diab + b_bmi + b_func +
                 6L * as.integer(cesd15 >= 16)]
say("(b) BDSI depression item changed to CES-D≥16: AUC = ",
    round(as.numeric(auc(roc(d3$event, d3$score_bdsi_b, quiet = TRUE, direction = "<"))), 4),
    "(original CES-D≥10 is ", round(as.numeric(auc(rocs[["score_bdsi"]])), 4), ")")
## (c) grip strength changed to mean of best hand values(effect Tu arm)
d4 <- copy(d)[, lp_tu_b := B_TU[["age_lt65"]] * as.integer(age15 < 65) + B_TU[["age_ge75"]] * as.integer(age15 >= 75) +
                 B_TU[["male"]] * male + B_TU[["edu_primid"]] * as.integer(edu3 == "Primary/middle") +
                 B_TU[["edu_none"]] * as.integer(edu3 == "No formal schooling") +
                 B_TU[["grip"]] * grip_alt + B_TU[["memory"]] * mem15 + B_TU[["orient"]] * orient15]
say("(c) Tu arm grip strength uses mean of best hand values: AUC = ",
    round(as.numeric(auc(roc(d4$event, d4$lp_tu_b, quiet = TRUE, direction = "<"))), 4))
## (d) weighted(2015 baseline weights)
dw <- d[!is.na(wt15)]
wt_auc <- function(p, y, wgt) {
  cs <- y == 1; ct <- y == 0
  num <- 0
  for (i in which(cs)) num <- num + wgt[i] * sum(wgt[ct] * ((p[i] > p[ct]) + 0.5 * (p[i] == p[ct])))
  num / (sum(wgt[cs]) * sum(wgt[ct]))
}
say("(d) weighted sensitivity(2015 baseline weights, available n = ", nrow(dw), "): ")
for (v in names(arms))
  say("    ", arms[[v]], " weighted AUC = ", round(wt_auc(dw[[v]], dw$event, dw$wt15), 4))

close(con)
cat("\n complete ->", file.path(OUT, "F14_report.txt"), "\n")
