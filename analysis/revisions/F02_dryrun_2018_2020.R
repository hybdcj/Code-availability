# ============================================================================
# F02 feasibility check: whether usable CHARLS 2018→2020 develop/validation model
# input: waves_cognition_2026-09-13.rds(data separately by wave), 2020 original HSF(word list version)
# output: C:/y03tmp/out/F02_report.txt, F02_*.csv
# key points: frozen model = A14(2011→2015, 4 year memory threshold crossing)
#       outcome definition: 2018 memory score(immediate recall trial 1 + delayed, 0–20)→ 2020 ties crossing threshold
#       cutoff: following original design use 2018 own mean−1SD rounding; also report"apply unchanged 2011 cutoff ≤4"as comparison
# ============================================================================
invisible(try(Sys.setlocale("LC_ALL", "zh_CN.UTF-8"), silent = TRUE))
suppressPackageStartupMessages({ library(data.table); library(pROC); library(haven) })
if (dir.exists("D:/R/library")) .libPaths(c("D:/R/library", .libPaths()))

TMP <- "C:/y03tmp/tmp"; OUT <- "C:/y03tmp/out"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
con <- file(file.path(OUT, "F02_report.txt"), open = "wt", encoding = "UTF-8")
say <- function(...) { s <- paste0(..., collapse = ""); cat(s, "\n"); cat(s, "\n", file = con) }
pt <- function(x) say(paste(capture.output(print(x)), collapse = "\n"))

## ---------- frozen model coefficients(A14_final_model_coefficients.csv) ----------
beta <- c("(Intercept)" = -4.19524508778947,
          age = 0.0492039250378111,
          eduPrimary = -0.290105771000328,
          eduMiddle = -0.994711800202235,
          eduHigh = -1.43588388035287,
          dep = 0.412166484954086,
          adl = -0.0110119344848014)

## ---------- 1. data ----------
w <- as.data.table(readRDS(file.path(TMP, "waves_cognition.rds")))
say("five wave dataset row count: ", nrow(w), "｜year: ", paste(sort(unique(w$year)), collapse = ", "))

d18 <- w[year == 2018]; d20 <- w[year == 2020]
say("2018 row count: ", nrow(d18), "｜2020 row count: ", nrow(d20))

a18 <- d18[, .(ID = IDn, imrc18 = as.numeric(imrc), dlrc18 = as.numeric(dlrc),
               age18 = as.numeric(age), edu18 = edu_group, cesd18 = as.numeric(cesd10),
               adl18 = as.numeric(adl), wt18 = as.numeric(weight))]
a18[, mem18 := imrc18 + dlrc18]
a20 <- d20[, .(ID = IDn, imrc20 = as.numeric(imrc), dlrc20 = as.numeric(dlrc))]
a20[, mem20 := imrc20 + dlrc20]
m <- merge(a18, a20, by = "ID", all = FALSE)
say("2018 and 2020 complete at both waves(by ID): ", nrow(m), " participants")

## ---------- 2. sequential sample flow ----------
say("\n========== 2. sample flow ==========")
say("① 2018 memory score calculable: ", sum(!is.na(m$mem18)))
say("② 2020 memory score calculable: ", sum(!is.na(m$mem20)))
say("③ scorable at both waves: ", sum(!is.na(m$mem18) & !is.na(m$mem20)))
say("④ ③ and 2018 all four covariates complete: ",
    sum(!is.na(m$mem18) & !is.na(m$mem20) & !is.na(m$age18) & !is.na(m$edu18) &
          !is.na(m$cesd18) & !is.na(m$adl18)))
say("   2018 missing education: ", sum(is.na(m$edu18)), "｜CES-D missing: ", sum(is.na(m$cesd18)),
    "｜ADL missing: ", sum(is.na(m$adl18)), "｜missing age: ", sum(is.na(m$age18)))

## ---------- 3. comparison of memory distributions across two waves ----------
say("\n========== 3. memory distribution comparison(same participants)==========")
pp <- m[!is.na(mem18) & !is.na(mem20)]
say("paired participant count: ", nrow(pp))
tab <- data.table(
  metrics = c("memory score (0-20)", "immediate recall (0-10)", "delayed recall (0-10)"),
  m2018 = round(c(mean(pp$mem18), mean(pp$imrc18), mean(pp$dlrc18)), 2),
  sd2018 = round(c(sd(pp$mem18), sd(pp$imrc18), sd(pp$dlrc18)), 2),
  m2020 = round(c(mean(pp$mem20), mean(pp$imrc20), mean(pp$dlrc20)), 2),
  sd2020 = round(c(sd(pp$mem20), sd(pp$imrc20), sd(pp$dlrc20)), 2),
  diff = round(c(mean(pp$mem20 - pp$mem18), mean(pp$imrc20 - pp$imrc18),
                 mean(pp$dlrc20 - pp$dlrc18)), 2),
  diff_sd = round(c(sd(pp$mem20 - pp$mem18), sd(pp$imrc20 - pp$imrc18),
                    sd(pp$dlrc20 - pp$dlrc18)), 2))
pt(tab)
say("paired difference t test(memory score): p = ",
    format(t.test(pp$mem20, pp$mem18, paired = TRUE)$p.value, scientific = TRUE, digits = 3))
fwrite(tab, file.path(OUT, "F02_memory_distribution.csv"), bom = TRUE)

## ---------- 4. cutoffs and events + transport frozen model ----------
say("\n========== 4. cutoff, events and frozen model transport ==========")
cut18 <- floor(mean(pp$mem18) - sd(pp$mem18))
say("2018 memory score: mean=", round(mean(pp$mem18), 2), " SD=", round(sd(pp$mem18), 2),
    " → mean−1SD rounding = ≤", cut18)
say("(original development definition 2011 cutoff = ≤4)")
say("2018 impaired fraction(≤", cut18, "): ", round(100 * mean(pp$mem18 <= cut18), 1), "%")

res <- data.table()
d18_keep <- NULL
for (k in c(cut18, 4)) {
  nm <- if (k == cut18) "2018 own cutoff" else "apply unchanged 2011 cutoff"
  d <- m[!is.na(mem18) & !is.na(mem20) & !is.na(age18) & !is.na(edu18) &
           !is.na(cesd18) & !is.na(adl18)]
  d <- d[mem18 > k]
  d[, event := as.integer(mem20 <= k)]
  d[, dep := as.integer(cesd18 >= 10)]
  d[, edu_group := factor(edu18, levels = c("Less than primary", "Primary school",
                                            "Middle school", "High school or above"))]
  d[, lp := beta[["(Intercept)"]] + beta[["age"]] * age18 + beta[["dep"]] * dep +
        beta[["adl"]] * adl18 +
        beta[["eduPrimary"]] * as.integer(edu_group == "Primary school") +
        beta[["eduMiddle"]] * as.integer(edu_group == "Middle school") +
        beta[["eduHigh"]] * as.integer(edu_group == "High school or above")]
  d[, phat := plogis(lp)]
  r <- roc(d$event, d$phat, quiet = TRUE, direction = "<")
  ci <- as.numeric(ci.auc(r))
  cal <- glm(event ~ lp, data = d, family = binomial())
  brier <- mean((d$phat - d$event)^2)
  say("\n-- ", nm, "(≤", k, ")--")
  say("   population at risk n=", nrow(d), "(2018 unimpaired)｜events=", sum(d$event),
      "｜event rate=", round(100 * mean(d$event), 2), "%")
  say("   frozen model: AUC=", round(ci[2], 3), " (", round(ci[1], 3), "–", round(ci[3], 3),
      ")｜calibration slope=", round(coef(cal)[2], 3), " intercept=", round(coef(cal)[1], 3),
      "｜Brier=", round(brier, 4))
  res <- rbind(res, data.table(definition = nm, cutoff = k, n = nrow(d), events = sum(d$event),
                               event_rate = round(100 * mean(d$event), 2),
                               AUC = round(ci[2], 3), AUC_low = round(ci[1], 3), AUC_high = round(ci[3], 3),
                               calibration_slope = round(coef(cal)[2], 3), calibration_intercept = round(coef(cal)[1], 3),
                               Brier = round(brier, 4)))
  if (k == cut18) d18_keep <- copy(d)
}
fwrite(res, file.path(OUT, "F02_frozen_model_transport.csv"), bom = TRUE)

## ---------- 5. reestimate weights in same participants(preview)----------
say("\n========== 5. reestimate in identical participants: 2018 weights vs 2011 frozen weights(preview)==========")
d <- copy(d18_keep)
fit <- glm(event ~ age18 + edu_group + dep + adl18, data = d, family = binomial())
say("reestimated coefficients: "); pt(round(coef(fit), 4))
d[, phat_new := predict(fit, type = "response")]
r_old <- roc(d$event, d$phat, quiet = TRUE, direction = "<")
r_new <- roc(d$event, d$phat_new, quiet = TRUE, direction = "<")
ci_old <- as.numeric(ci.auc(r_old)); ci_new <- as.numeric(ci.auc(r_new))
say("old weights AUC=", round(ci_old[2], 3), " (", round(ci_old[1], 3), "–", round(ci_old[3], 3), ")")
say("new weights AUC=", round(ci_new[2], 3), " (", round(ci_new[1], 3), "–", round(ci_new[3], 3), ")")
say("paired DeLong: p = ", format(roc.test(r_old, r_new, paired = TRUE)$p.value, digits = 3))
best <- function(r) {
  co <- coords(r, "best", best.method = "youden", ret = c("threshold", "sensitivity", "specificity"))
  data.table(threshold = round(co$threshold, 4), diagnostic_sensitivity = round(co$sensitivity, 3), specificity = round(co$specificity, 3))
}
say("\n old weights(Youden cutoff): "); pt(best(r_old))
say("new weights(Youden cutoff): "); pt(best(r_new))

## ---------- 6. 2020 word list version A/B and memory score ----------
say("\n========== 6. 2020 word list version A/B and memory score ==========")
x5 <- as.data.table(read_dta(file.path(TMP, "w5_hsf.dta"),
                             col_select = all_of(c("ID", paste0("xwordlist_", 1:10, "_")))))
say("2020 original file ID example: ", paste(head(as.character(x5$ID), 3), collapse = ", "))
x5[, ver := fifelse(as.character(xwordlist_1_) == "\u6cb3\u6d41", "A",
                    fifelse(as.character(xwordlist_1_) == "\u5929\u7a7a", "B", "not administered/missing"))]
x5 <- x5[, .(ID = as.character(ID), ver)]
mm <- merge(m[!is.na(mem18) & !is.na(mem20)], x5, by = "ID", all.x = TRUE)
pt(mm[, .(n = .N, mem2018 = round(mean(mem18), 2), mem2020 = round(mean(mem20), 2),
          immediate_2020 = round(mean(imrc20), 2), delayed_2020 = round(mean(dlrc20), 2)), by = ver][order(ver)])

close(con)
cat("\n complete ->", file.path(OUT, "F02_report.txt"), "\n")
