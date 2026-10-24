# F32(2026-10-16): Retain the historical common complete sample calculation below for traceability; CHARLS The integer score now uses its own sample from the same directory F32_correct_population_eligibility_and_denominators.R.Do not merely remove the lp restriction and run the old script, because the downstream continuous model, denominators and headings require separate handling.
# ============================================================================
# F15c: take F15 of"full population definition"replace with**original main manuscript threshold(global ≤4)**, and first round R02
#        broader sample of(n=8,296, continuous model AUC 0.726)align, resolve incompatible definitions
# window: CHARLS 2011 -> 2015; outcome: 2015 memory score <= 4(baseline > 4)
# ============================================================================
invisible(try(Sys.setlocale("LC_ALL", "zh_CN.UTF-8"), silent = TRUE))
suppressPackageStartupMessages({ library(data.table); library(pROC) })
if (dir.exists("D:/R/library")) .libPaths(c("D:/R/library", .libPaths()))
TMP <- "C:/y03tmp/tmp"; OUT <- "C:/y03tmp/out"
con <- file(file.path(OUT, "F15c_report.txt"), open = "wt", encoding = "UTF-8")
say <- function(...) { s <- paste0(..., collapse = ""); cat(s, "\n"); cat(s, "\n", file = con) }
pt <- function(x) say(paste(capture.output(print(x)), collapse = "\n"))

B <- c(intercept = -4.19524508778947, age = 0.0492039250378111,
       eduPri = -0.290105771000328, eduMid = -0.994711800202235,
       eduHigh = -1.43588388035287, dep = 0.412166484954086, adl = -0.0110119344848014)
LV <- c("Less than primary", "Primary school", "Middle school", "High school or above")
K <- 4L   # main manuscript threshold: 2011 five wave fixed cohort mean-1SD = 4.7647 -> <= 4

w <- as.data.table(readRDS(file.path(TMP, "waves_cognition.rds")))
g <- w[year %in% c(2011, 2015), .(ID = IDn, year = as.numeric(year),
        mem = as.numeric(imrc) + as.numeric(dlrc),
        age = as.numeric(age), edu = edu_group, cesd = as.numeric(cesd10),
        adl = as.numeric(adl), female = as.numeric(female), wt = as.numeric(weight))]
a <- g[year == 2011, .(ID, mem0 = mem, age0 = age, edu0 = edu, cesd0 = cesd,
                       adl0 = adl, wt0 = wt, female0 = female)]
b <- g[year == 2015, .(ID, mem1 = mem)]
d <- merge(a, b, by = "ID")
d[, edu_f := factor(edu0, levels = LV)]
d[, score := 2L * pmax(floor((age0 - 45) / 5), 0) +
      (-2L) * as.integer(edu_f == LV[2]) + (-7L) * as.integer(edu_f == LV[3]) +
      (-10L) * as.integer(edu_f == LV[4]) + 3L * as.integer(cesd0 >= 10)]
d[, lp := B[["intercept"]] + B[["age"]] * age0 + B[["dep"]] * as.integer(cesd0 >= 10) +
      B[["adl"]] * adl0 + B[["eduPri"]] * as.integer(edu_f == LV[2]) +
      B[["eduMid"]] * as.integer(edu_f == LV[3]) + B[["eduHigh"]] * as.integer(edu_f == LV[4])]
d[, phat := plogis(lp)]
d[, event := as.integer(mem1 <= K)]
d[, low0 := as.integer(mem0 <= K)]
ok <- !is.na(d$score) & !is.na(d$lp) & !is.na(d$event) & !is.na(d$age0) & d$age0 >= 45
d <- d[ok]

say("threshold: global <= ", K, "(main manuscript definition; derived from 2011 five wave fixed cohort mean-1SD = 4.7647)")
say("sample flow: two waves ID matchable ", nrow(merge(a, b, by = "ID")),
    " -> age/score/lp/complete outcomes and >=45 years old: ", nrow(d))
say("A11/R02 record: 2011 year >=45 years old and baseline >4 participants 10,422 participants; 2015 measurable memory 8,371; ")
say("             complete four variable data 8,296 participants(= broader sample); including events 1,889(22.8%)")

arm <- function(x, lab) {
  if (nrow(x) < 30 || length(unique(x$event)) < 2) return(NULL)
  r1 <- roc(x$event, x$score, quiet = TRUE, direction = "<")
  r2 <- roc(x$event, x$phat, quiet = TRUE, direction = "<")
  c1 <- as.numeric(ci.auc(r1)); c2 <- as.numeric(ci.auc(r2))
  data.table(definition = lab, n = nrow(x), events = sum(x$event),
             event_rate = round(100 * mean(x$event), 2),
             AUC_scorecard = round(c1[2], 4), scorecard_low = round(c1[1], 4), scorecard_high = round(c1[3], 4),
             AUC_continuous = round(c2[2], 4), continuous_low = round(c2[1], 4), continuous_high = round(c2[3], 4))
}
res <- rbindlist(list(
  arm(d, "U full population(no baseline memory restriction)"),
  arm(d[!is.na(mem0) & mem0 > K], "N unimpaired at baseline(main manuscript broader sample definition)"),
  arm(d[!is.na(mem0) & mem0 <= K], "I impaired at baseline"),
  arm(d[is.na(mem0)], "M missing baseline memory")), fill = TRUE)
say("\n========== global <=4 definitions under threshold ==========")
pt(res)
fwrite(res, file.path(OUT, "F15c_auc_globalcut.csv"), bom = TRUE)

say("\n========== and first round R02 broader sample comparison ==========")
pt(data.table(
  metrics = c("n", "events", "event_rate", "continuous model AUC", "scorecard AUC"),
  R02_broader_sample = c("8,296", "1,889", "22.8%", "0.726", "not reported"),
  this_release_N_definition = c(format(nrow(d[!is.na(mem0) & mem0 > K]), big.mark = ","),
                 format(sum(d[!is.na(mem0) & mem0 > K]$event), big.mark = ","),
                 paste0(round(100 * mean(d[!is.na(mem0) & mem0 > K]$event), 1), "%"),
                 as.character(res[definition == "N unimpaired at baseline(main manuscript broader sample definition)"]$AUC_continuous),
                 as.character(res[definition == "N unimpaired at baseline(main manuscript broader sample definition)"]$AUC_scorecard))))

say("\n========== sequential relaxation(2011 -> 2015 window, global <=4 threshold)==========")
uq <- d[!is.na(mem0)]
say("① original manuscript development sample(complete five wave common items): 4,441 participants / 658 events / apparent AUC 0.702")
say("② broader sample(remove five wave requirement, baseline >4): ", nrow(uq[mem0 > K]), " participants / ",
    sum(uq[mem0 > K]$event), " events / AUC 0.726(R02 record 8,296 / 1,889)")
say("③ full population(then remove baseline >4 require, this is F15 definition of): ", nrow(d), " participants / ",
    sum(d$event), " events")
say("④ full population + no baseline memory restriction: " , nrow(d), " participants / ", sum(d$event), " events")

say("\n========== decomposition: baseline low scorers as fraction of full population events ==========")
dd <- d[!is.na(mem0)]
say("baseline among events <=4 participant fraction = ", round(100 * mean(dd$low0[dd$event == 1]), 1), "%")
say("of which: event rate among baseline low scorers = ", round(100 * mean(dd[low0 == 1]$event), 1),
    "%｜event rate among baseline unimpaired participants = ", round(100 * mean(dd[low0 == 0]$event), 1), "%")

close(con)
cat("\n complete ->", file.path(OUT, "F15c_report.txt"), "\n")
