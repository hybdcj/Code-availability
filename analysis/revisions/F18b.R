# F32(2026-10-16): Retain the historical common complete sample calculation below for traceability; CHARLS The integer score now uses its own sample from the same directory F32_correct_population_eligibility_and_denominators.R.Do not merely remove the lp restriction and run the old script, because the downstream continuous model, denominators and headings require separate handling.
# ============================================================================
# F18b: add description"step number ③ level compared with ② newly included at level 3,457 participants"
# definition: main manuscript threshold <=4, CHARLS 2011 -> 2015
# ============================================================================
invisible(try(Sys.setlocale("LC_ALL", "zh_CN.UTF-8"), silent = TRUE))
suppressPackageStartupMessages({ library(data.table) })
if (dir.exists("D:/R/library")) .libPaths(c("D:/R/library", .libPaths()))
TMP <- "C:/y03tmp/tmp"; OUT <- "C:/y03tmp/out"
con <- file(file.path(OUT, "F18b_report.txt"), open = "wt", encoding = "UTF-8")
say <- function(...) { s <- paste0(..., collapse = ""); cat(s, "\n"); cat(s, "\n", file = con) }
pt <- function(x) say(paste(capture.output(print(x)), collapse = "\n"))
LV <- c("Less than primary", "Primary school", "Middle school", "High school or above")
B <- c(intercept = -4.19524508778947, age = 0.0492039250378111,
       eduPri = -0.290105771000328, eduMid = -0.994711800202235,
       eduHigh = -1.43588388035287, dep = 0.412166484954086, adl = -0.0110119344848014)
K <- 4L

w <- as.data.table(readRDS(file.path(TMP, "waves_cognition.rds")))
g <- w[year %in% c(2011, 2015), .(ID = IDn, year = as.numeric(year),
        mem = as.numeric(imrc) + as.numeric(dlrc), age = as.numeric(age),
        edu = edu_group, cesd = as.numeric(cesd10), adl = as.numeric(adl))]
a <- g[year == 2011, .(ID, mem0 = mem, age0 = age, edu0 = edu, cesd0 = cesd, adl0 = adl)]
b <- g[year == 2015, .(ID, mem1 = mem)]
d <- merge(a, b, by = "ID")
d[, edu_f := factor(edu0, levels = LV)]
d[, score := 2L * pmax(floor((age0 - 45) / 5), 0) +
      (-2L) * as.integer(edu_f == LV[2]) + (-7L) * as.integer(edu_f == LV[3]) +
      (-10L) * as.integer(edu_f == LV[4]) + 3L * as.integer(cesd0 >= 10)]
d[, lp := B[["intercept"]] + B[["age"]] * age0 + B[["dep"]] * as.integer(cesd0 >= 10) +
      B[["adl"]] * adl0 + B[["eduPri"]] * as.integer(edu_f == LV[2]) +
      B[["eduMid"]] * as.integer(edu_f == LV[3]) + B[["eduHigh"]] * as.integer(edu_f == LV[4])]
d <- d[!is.na(score) & !is.na(lp) & !is.na(age0) & age0 >= 45 & !is.na(mem1)]
d[, event := as.integer(mem1 <= K)]

grp <- function(x, lab) data.table(
  groups = lab, n = nrow(x),
  mean_age = round(mean(x$age0), 1),
  fraction_with_no_schooling = round(100 * mean(x$edu0 == LV[1]), 1),
  baseline_memory_measurable = sum(!is.na(x$mem0)),
  mean_baseline_memory = round(mean(x$mem0, na.rm = TRUE), 2),
  events = sum(x$event), event_rate = round(100 * mean(x$event), 1),
  mean_risk_score = round(mean(x$score), 2))

T <- rbindlist(list(
  grp(d[!is.na(mem0) & mem0 > K], "① broader sample(baseline>4, R02 number ② level)"),
  grp(d[!is.na(mem0) & mem0 <= K], "② baseline low score(newly added)"),
  grp(d[is.na(mem0)], "③ baseline memory not measurable(newly added)"),
  grp(d, "total(full population, number ③ level)")))
say("using main manuscript threshold <=4, CHARLS 2011->2015.")
say("number ③ level 11,753 participants − number ② level 8,296 participants = newly added ", nrow(d) - nrow(d[!is.na(mem0) & mem0 > K]), " participants")
say("verification: ② + ③ = ", nrow(d[!is.na(mem0) & mem0 <= K]) + nrow(d[is.na(mem0)]))
pt(T)
fwrite(T, file.path(OUT, "F18b_newcomers.csv"), bom = TRUE)
close(con)
cat("\n complete\n")
