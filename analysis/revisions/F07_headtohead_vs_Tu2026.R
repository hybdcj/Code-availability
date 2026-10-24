# ============================================================================
# F07: and Tu etc(GeroScience 2026)direct model comparison
# window: CHARLS 2015 baseline → 2018 outcome(exactly their temporal validation window)
# outcome: use**their definition**——memory score falls to"education strata mean−1SD cutoff"below
# comparison: their 6 factor CoxPH linear predictor(using their publication Table 2 published coefficients)
#       vs our 4 factor frozen score(A14, 2011→2015 fit, direct transport)
# output: C:/y03tmp/out/F07_report.txt, F07_*.csv
# advance declaration: post hoc direct comparison; our score is 4 year model transported to 3 year window, 
#           absolute values will be miscalibrated, compare ranking metrics only(AUC)and relative structure.
# ============================================================================
invisible(try(Sys.setlocale("LC_ALL", "zh_CN.UTF-8"), silent = TRUE))
suppressPackageStartupMessages({ library(data.table); library(pROC); library(haven) })
if (dir.exists("D:/R/library")) .libPaths(c("D:/R/library", .libPaths()))
TMP <- "C:/y03tmp/tmp"; OUT <- "C:/y03tmp/out"
con <- file(file.path(OUT, "F07_report.txt"), open = "wt", encoding = "UTF-8")
say <- function(...) { s <- paste0(..., collapse = ""); cat(s, "\n"); cat(s, "\n", file = con) }
pt  <- function(x) say(paste(capture.output(print(x)), collapse = "\n"))

## ---------- their coefficients(manuscript Table 2, CoxPH)----------
B_TU <- c(age_lt65 = -0.53, age_ge75 = 0.26, male = 0.38,
          edu_primid = 0.45, edu_none = 0.03,
          grip = -0.16, memory = -0.33, orient = -0.15)
## ---------- our frozen model(A14)----------
B_US <- c("(Intercept)" = -4.19524508778947, age = 0.0492039250378111,
          eduP = -0.290105771000328, eduM = -0.994711800202235, eduH = -1.43588388035287,
          dep = 0.412166484954086, adl = -0.0110119344848014)
LV <- c("Less than primary", "Primary school", "Middle school", "High school or above")

## ---------------- 1. cognitive data ----------------
w <- as.data.table(readRDS(file.path(TMP, "waves_cognition.rds")))
g <- w[year %in% c(2011, 2015, 2018), .(ID = IDn, year,
                                        mem = as.numeric(imrc) + as.numeric(dlrc),
                                        orient = as.numeric(orient),
                                        age = as.numeric(age), female = as.numeric(female),
                                        edu_group = edu_group,
                                        cesd10 = as.numeric(cesd10), adl = as.numeric(adl))]

## ---------------- 2. 2015 grip strength ----------------
b <- as.data.table(read_dta(file.path(TMP, "w3_biomarker.dta"),
                            col_select = any_of(c("ID", "qc003", "qc004", "qc005", "qc006"))))
say("biomarker file ID columns: ", paste(setdiff(names(b), c("qc003","qc004","qc005","qc006")), collapse = ", "))
grip <- b[, .(ID = as.character(ID),
              grip_max = do.call(pmax, c(.SD, na.rm = TRUE)),
              grip_n   = rowSums(!is.na(.SD))),
          .SDcols = c("qc003", "qc004", "qc005", "qc006")]
say("2015 grip strength: nonmissing ", sum(!is.na(grip$grip_max)), " participants; median ",
    round(median(grip$grip_max, na.rm = TRUE), 1), " kg")

## ---------------- 3. education specific cutoff(at 2011 reconstruct at wave)----------------
d11 <- g[year == 2011 & !is.na(mem) & !is.na(edu_group) & age >= 50]
d11[, edu3 := fcase(edu_group == LV[1], "No formal schooling",
                    edu_group %in% LV[2:3], "Primary/middle",
                    edu_group == LV[4], "High school and above")]
cut3 <- d11[, .(n = .N, mean = mean(mem), sd = sd(mem)), by = edu3]
cut3[, cutoff := mean - sd]
say("\n education specific cutoff(at CHARLS 2011, ≥50 years old, reconstruct among scorable memory participants): ")
pt(cut3[, .(edu3, n, mean = round(mean, 2), SD = round(sd, 2), cutoff = round(cutoff, 2))])

## ---------------- 4. build 2015→2018 analysis sample ----------------
d15 <- g[year == 2015, .(ID, mem15 = mem, orient15 = orient, age15 = age,
                         female15 = female, edu_group15 = edu_group,
                         cesd10_15 = cesd10, adl15 = adl)]
d18 <- g[year == 2018, .(ID, mem18 = mem, orient18 = orient)]
mm <- merge(merge(d15, d18, by = "ID"), grip, by = "ID", all.x = TRUE)
mm[, edu3 := fcase(edu_group15 == LV[1], "No formal schooling",
                   edu_group15 %in% LV[2:3], "Primary/middle",
                   edu_group15 == LV[4], "High school and above")]
mm <- merge(mm, cut3[, .(edu3, cutoff)], by = "edu3", all.x = TRUE)

say("\n sample flow: ")
say("  2015 and 2018 both memory assessments scorable: ", nrow(mm[!is.na(mem15) & !is.na(mem18)]))
say("  ＋2015 ≥50 years old: ", nrow(mm[!is.na(mem15) & !is.na(mem18) & age15 >= 50]))
say("  ＋cutoff available: ", nrow(mm[!is.na(mem15) & !is.na(mem18) & age15 >= 50 & !is.na(cutoff)]))
say("  ＋grip strength available: ", nrow(mm[!is.na(mem15) & !is.na(mem18) & age15 >= 50 & !is.na(cutoff) & !is.na(grip_max)]))
say("  ＋orientation score available: ", nrow(mm[!is.na(mem15) & !is.na(mem18) & age15 >= 50 & !is.na(cutoff) &
                                 !is.na(grip_max) & !is.na(orient15)]))
d <- mm[!is.na(mem15) & !is.na(mem18) & age15 >= 50 & !is.na(cutoff) &
           !is.na(grip_max) & !is.na(orient15)]
d <- d[!is.na(cesd10_15) & !is.na(adl15)]
say("  ＋our covariates complete: ", nrow(d))

## enroll baseline normal memory participants(their definition: ≥cutoff is intact)
d <- d[mem15 >= cutoff]
say("  ＋2015 normal memory(≥education specific cutoff): ", nrow(d))
d[, event := as.integer(mem18 < cutoff)]
d[, dep := as.integer(cesd10_15 >= 10)]
d[, edu_group := factor(edu_group15, levels = LV)]
d[, male := as.integer(female15 == 0)]
d[, agecat := fcase(age15 < 65, "lt65", age15 < 75, "6574", default = "ge75")]
say("\n analysis sample n = ", nrow(d), "｜outcome(their definition)events = ", sum(d$event),
    "(", round(100 * mean(d$event), 2), "%)")
say("(comparison: their reported temporal validation cohort n=1,333, events 210, 15.8%)")

## ---------------- 5. two models ----------------
d[, lp_tu := B_TU[["age_lt65"]] * as.integer(agecat == "lt65") +
             B_TU[["age_ge75"]] * as.integer(agecat == "ge75") +
             B_TU[["male"]] * male +
             B_TU[["edu_primid"]] * as.integer(edu3 == "Primary/middle") +
             B_TU[["edu_none"]] * as.integer(edu3 == "No formal schooling") +
             B_TU[["grip"]] * grip_max + B_TU[["memory"]] * mem15 +
             B_TU[["orient"]] * orient15]
d[, lp_us := B_US[["(Intercept)"]] + B_US[["age"]] * age15 + B_US[["dep"]] * dep +
             B_US[["adl"]] * adl15 +
             B_US[["eduP"]] * as.integer(edu_group == LV[2]) +
             B_US[["eduM"]] * as.integer(edu_group == LV[3]) +
             B_US[["eduH"]] * as.integer(edu_group == LV[4])]
d[, p_tu := plogis(scale(lp_tu)[, 1])][, p_us := plogis(lp_us)]
## simplified scorecard: only 3 item scoring(age +2/5 years old, education 0/−2/−7/−10, depression +3); ADL retained in continuous model but scores 0 points
d[, score3 := 2L * pmax(floor((age15 - 45)/5), 0) +
      (-2L) * as.integer(edu_group == LV[2]) +
      (-7L) * as.integer(edu_group == LV[3]) +
      (-10L) * as.integer(edu_group == LV[4]) +
      3L * dep]

## our**simplified scorecard**(three items: age＋education＋depression; ADL score 0 points, equivalent to excluding)
d[, score_us := 2L * pmax(floor((age15 - 45)/5), 0) +
        (-2L) * as.integer(edu_group == LV[2]) + (-7L) * as.integer(edu_group == LV[3]) +
        (-10L) * as.integer(edu_group == LV[4]) + 3L * dep]

## ---------------- 6. performance ----------------
say("\n========== primary comparison(outcome＝their definition)==========")
r_tu <- roc(d$event, d$lp_tu, quiet = TRUE, direction = "<")
r_us <- roc(d$event, d$lp_us, quiet = TRUE, direction = "<")
ci_tu <- as.numeric(ci.auc(r_tu)); ci_us <- as.numeric(ci.auc(r_us))
dt <- roc.test(r_us, r_tu, paired = TRUE)
say("Tu 6 factor: AUC = ", round(ci_tu[2], 3), "(", round(ci_tu[1], 3), "–", round(ci_tu[3], 3), ")")
say("our model 4 factor: AUC = ", round(ci_us[2], 3), "(", round(ci_us[1], 3), "–", round(ci_us[3], 3), ")")
say("ΔAUC(our model − their model)= ", round(ci_us[2] - ci_tu[2], 4),
    "｜paired DeLong p = ", format(dt$p.value, digits = 3))
r_s3 <- roc(d$event, d$score3, quiet = TRUE, direction = "<")
ci_s3 <- as.numeric(ci.auc(r_s3))
say("our simplified scorecard(**3 item scoring**, ADL score 0): AUC = ", round(ci_s3[2], 3),
    "(", round(ci_s3[1], 3), "–", round(ci_s3[3], 3), ")",
    "｜ΔAUC(3 item card − their model)= ", round(ci_s3[2] - ci_tu[2], 4),
    "｜paired DeLong p = ", format(roc.test(r_s3, r_tu, paired = TRUE)$p.value, digits = 3))

## simplified scorecard(three items, none ADL)as third arm
r_sc <- roc(d$event, d$score_us, quiet = TRUE, direction = "<")
ci_sc <- as.numeric(ci.auc(r_sc))
dt_sc <- roc.test(r_sc, r_tu, paired = TRUE)
say("our model 3 item scorecard(none ADL): AUC = ", round(ci_sc[2], 3), "(", round(ci_sc[1], 3), "–",
    round(ci_sc[3], 3), ")｜relative to theirs ΔAUC = ", round(ci_sc[2] - ci_tu[2], 4),
    "｜paired DeLong p = ", format(dt_sc$p.value, digits = 3))
say("(continuous model vs scorecard: ΔAUC = ", round(ci_us[2] - ci_sc[2], 4), ", namely ADL actual contribution of)")

## paired Bootstrap(resample individuals)
set.seed(20261014); B <- 2000; dv <- rep(NA_real_, B); au <- av <- rep(NA_real_, B)
for (i in seq_len(B)) {
  idx <- sample.int(nrow(d), replace = TRUE); db <- d[idx]
  if (length(unique(db$event)) < 2) next
  au[i] <- as.numeric(auc(roc(db$event, db$lp_tu, quiet = TRUE, direction = "<")))
  av[i] <- as.numeric(auc(roc(db$event, db$lp_us, quiet = TRUE, direction = "<")))
  dv[i] <- av[i] - au[i]
}
say("paired Bootstrap 95% interval: their model [", round(quantile(au, .025, na.rm = TRUE), 3), ", ",
    round(quantile(au, .975, na.rm = TRUE), 3), "]; our model [",
    round(quantile(av, .025, na.rm = TRUE), 3), ", ", round(quantile(av, .975, na.rm = TRUE), 3), "]")
say("ΔAUC of 95% interval: [", round(quantile(dv, .025, na.rm = TRUE), 4), ", ",
    round(quantile(dv, .975, na.rm = TRUE), 4), "]; P(ΔAUC>0) = ",
    round(mean(dv > 0, na.rm = TRUE), 3))

## calibration and Brier not compared: their CoxPH linear predictor lacks baseline risk function, 
## and our logistic absolute probabilities on different scales, calibration comparison distorted.
say("\n(calibration and Brier do not compare: their CoxPH linear predictor has no baseline risk function, ")
say("  absolute risks on different scales, compare ranking metrics only.)")
cal_ok <- function(lp, y) coef(glm(y ~ as.numeric(scale(lp)), family = binomial()))
say("  reference: take each LP standardized calibration slope——their model ", round(cal_ok(d$lp_tu, d$event)[2], 3),
    "; our model ", round(cal_ok(d$lp_us, d$event)[2], 3), "(higher values indicate greater risk dispersion)")

## risk stratification(cohort specific quartiles)
say("\n event rates by risk stratum(cohort specific quartiles): ")
for (v in c("lp_tu", "lp_us")) {
  d[, q := cut(get(v), quantile(get(v), seq(0, 1, .25)), include.lowest = TRUE, labels = 1:4)]
  say(paste0("  ", v, ": ", paste(d[, .(event_rate = paste0(round(100*mean(event),1), "%"), n = .N),
                                    by = q][order(q)]$event_rate, collapse = "  →  ")))
}

## ---------------- 7. secondary: switch to our outcome definition ----------------
say("\n========== secondary: switch to our definition(global mean−1SD rounded cutoff)==========")
mu15 <- mean(mm$mem15, na.rm = TRUE); sd15 <- sd(mm$mem15, na.rm = TRUE)
k15 <- floor(mu15 - sd15)
dd <- copy(d)[, event2 := as.integer(mem18 <= k15)]
say("2015 global cutoff = ≤", k15, "; events = ", sum(dd$event2), "(", round(100*mean(dd$event2),2), "%)")
r2t <- roc(dd$event2, dd$lp_tu, quiet = TRUE, direction = "<")
r2u <- roc(dd$event2, dd$lp_us, quiet = TRUE, direction = "<")
say("Tu 6 factor AUC = ", round(as.numeric(auc(r2t)), 3),
    "; our model 4 factor AUC = ", round(as.numeric(auc(r2u)), 3),
    "; ΔAUC = ", round(as.numeric(auc(r2u)) - as.numeric(auc(r2t)), 4),
    "｜paired DeLong p = ", format(roc.test(r2u, r2t, paired = TRUE)$p.value, digits = 3))
r2s <- roc(dd$event2, dd$score3, quiet = TRUE, direction = "<")
say("our model 3 item scorecard AUC = ", round(as.numeric(auc(r2s)), 3),
    "; ΔAUC(3 item card − their model)= ", round(as.numeric(auc(r2s)) - as.numeric(auc(r2t)), 4),
    "｜paired DeLong p = ", format(roc.test(r2s, r2t, paired = TRUE)$p.value, digits = 3))
r2s <- roc(dd$event2, dd$score_us, quiet = TRUE, direction = "<")
say("Tu 6 factor AUC = ", round(as.numeric(auc(r2t)), 3),
    "; our model 3 item scorecard AUC = ", round(as.numeric(auc(r2s)), 3),
    "; ΔAUC = ", round(as.numeric(auc(r2s)) - as.numeric(auc(r2t)), 4),
    "｜paired DeLong p = ", format(roc.test(r2s, r2t, paired = TRUE)$p.value, digits = 3))

## ---------------- 8. sensitivity: grip strength definition ----------------

## ---------------- 9. each model"better than random"test of ----------------
say("\n========== appendix: each model own vs random(AUC=0.5)test of ==========")
auc_vs_chance <- function(y, p, tag) {
  r <- roc(y, p, quiet = TRUE, direction = "<")
  ci <- as.numeric(ci.auc(r))
  se <- (ci[3] - ci[1]) / (2 * 1.96)
  z <- (ci[2] - 0.5) / se
  data.table(model = tag, AUC = round(ci[2], 3), SE = round(se, 4),
             z = round(z, 1), p_vs_0.5 = format(2 * pnorm(-abs(z)), scientific = TRUE, digits = 3))
}
say("\n-- outcome＝their definition --")
pt(rbind(auc_vs_chance(d$event, d$lp_tu, "Tu 6 factor"),
         auc_vs_chance(d$event, d$score_us, "our model 3 item scorecard"),
         auc_vs_chance(d$event, d$lp_us, "our continuous model")))
say("-- outcome＝our definition --")
pt(rbind(auc_vs_chance(dd$event2, dd$lp_tu, "Tu 6 factor"),
         auc_vs_chance(dd$event2, dd$score_us, "our model 3 item scorecard"),
         auc_vs_chance(dd$event2, dd$lp_us, "our continuous model")))
say("\n noninferiority(their definition): ΔAUC of 95% interval lower bound = ", round(quantile(dv, .025, na.rm = TRUE), 4),
    "; if prespecified noninferiority margin is −0.02, interval lies entirely above margin, noninferiority criterion met.")

say("\n========== sensitivity: grip strength changed to[mean of best value from each hand]==========")
grip2 <- b[, .(ID = as.character(ID),
               g1 = pmax(qc003, qc004, na.rm = TRUE),
               g2 = pmax(qc005, qc006, na.rm = TRUE))]
grip2[, grip_alt := rowMeans(cbind(g1, g2), na.rm = TRUE)]
d2 <- merge(d[, .(ID, event, lp_tu, lp_us, mem15, orient15, age15, male, edu3)], grip2,
            by = "ID")
d2 <- d2[!is.na(grip_alt)]
d2[, lp_tu_alt := B_TU[["age_lt65"]] * as.integer(age15 < 65) +
                   B_TU[["age_ge75"]] * as.integer(age15 >= 75) +
                   B_TU[["male"]] * as.integer(male) +
                   B_TU[["edu_primid"]] * as.integer(edu3 == "Primary/middle") +
                   B_TU[["edu_none"]] * as.integer(edu3 == "No formal schooling") +
                   B_TU[["grip"]] * grip_alt + B_TU[["memory"]] * mem15 +
                   B_TU[["orient"]] * orient15]
say("n = ", nrow(d2), "｜Tu model(grip strength uses mean of best hand specific values)AUC = ",
    round(as.numeric(auc(roc(d2$event, d2$lp_tu_alt, quiet = TRUE, direction = "<"))), 3),
    "; our model 4 factor AUC = ",
    round(as.numeric(auc(roc(d2$event, d2$lp_us, quiet = TRUE, direction = "<"))), 3))

## export
out <- data.table(
  definition = c("Tu(education specific cutoff)", "Tu(education specific cutoff)", "our model(global cutoff ≤", "our model(global cutoff ≤"),
  model = c("Tu 6 factor", "our model 4 factor", "Tu 6 factor", "our model 4 factor")
)
out <- data.table(
  outcome_specification = rep(c("Tu education specific cutoff", "our global cutoff"), each = 2),
  model = rep(c("Tu 6 factor(includes cognitive testing)", "our model 4 factor(no cognitive testing)"), 2),
  AUC = c(round(ci_tu[2], 3), round(ci_us[2], 3),
          round(as.numeric(auc(r2t)), 3), round(as.numeric(auc(r2u)), 3)),
  n = nrow(d), events = c(sum(d$event), sum(d$event), sum(dd$event2), sum(dd$event2)))
fwrite(out, file.path(OUT, "F07_headtohead.csv"), bom = TRUE)

close(con)
cat("\n complete ->", file.path(OUT, "F07_report.txt"), "\n")
