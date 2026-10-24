# F32(2026-10-16): Retain the historical common complete sample calculation below for traceability; CHARLS The integer score now uses its own sample from the same directory F32_correct_population_eligibility_and_denominators.R.Do not merely remove the lp restriction and run the old script, because the downstream continuous model, denominators and headings require separate handling.
# ============================================================================
# F20: older age(>=90 years old)participant count in**broader population definition**recalculation under(update first round R09)
# four cohorts: CHARLS 2011->2015, 2015->2018; HRS, ELSA, SHARE
# also provide"main manuscript threshold <=4"and"each cohort own mean-1SD"event counts under two definitions
# ============================================================================
invisible(try(Sys.setlocale("LC_ALL", "zh_CN.UTF-8"), silent = TRUE))
suppressPackageStartupMessages({ library(data.table) })
if (dir.exists("D:/R/library")) .libPaths(c("D:/R/library", .libPaths()))
TMP <- "C:/y03tmp/tmp"; OUT <- "C:/y03tmp/out"
con <- file(file.path(OUT, "F20_report.txt"), open = "wt", encoding = "UTF-8")
say <- function(...) { s <- paste0(..., collapse = ""); cat(s, "\n"); cat(s, "\n", file = con) }
pt <- function(x) say(paste(capture.output(print(x)), collapse = "\n"))
LV <- c("Less than primary", "Primary school", "Middle school", "High school or above")
G3 <- c("A_none", "B_primid", "C_high")
K <- 4L
B <- c(intercept = -4.19524508778947, age = 0.0492039250378111,
       eduPri = -0.290105771000328, eduMid = -0.994711800202235,
       eduHigh = -1.43588388035287, dep = 0.412166484954086, adl = -0.0110119344848014)

## ---------- CHARLS: construct wide data for two windows ----------
w <- as.data.table(readRDS(file.path(TMP, "waves_cognition.rds")))
g <- w[year %in% c(2011, 2015, 2018), .(ID = IDn, year = as.numeric(year),
        mem = as.numeric(imrc) + as.numeric(dlrc), age = as.numeric(age),
        edu = edu_group, cesd = as.numeric(cesd10), adl = as.numeric(adl))]
d11 <- g[year == 2011 & !is.na(mem) & !is.na(edu) & age >= 50]
d11[, g3 := fcase(edu == LV[1], G3[1], edu %in% LV[2:3], G3[2], edu == LV[4], G3[3])]
cut3 <- d11[!is.na(g3), .(cutoff = mean(mem) - sd(mem)), by = g3]

make_win <- function(y0, y1) {
  a <- g[year == y0, .(ID, mem0 = mem, age0 = age, edu0 = edu, cesd0 = cesd, adl0 = adl)]
  b <- g[year == y1, .(ID, mem1 = mem)]
  d <- merge(a, b, by = "ID")
  d[, g3 := fcase(edu0 == LV[1], G3[1], edu0 %in% LV[2:3], G3[2], edu0 == LV[4], G3[3])]
  d <- merge(d, cut3[, .(g3, cutoff)], by = "g3", all.x = TRUE)
  d[, edu_f := factor(edu0, levels = LV)]
  d[, score := 2L * pmax(floor((age0 - 45) / 5), 0) +
        (-2L) * as.integer(edu_f == LV[2]) + (-7L) * as.integer(edu_f == LV[3]) +
        (-10L) * as.integer(edu_f == LV[4]) + 3L * as.integer(cesd0 >= 10)]
  d[, lp := B[["intercept"]] + B[["age"]] * age0 + B[["dep"]] * as.integer(cesd0 >= 10) +
      B[["adl"]] * adl0 + B[["eduPri"]] * as.integer(edu_f == LV[2]) +
      B[["eduMid"]] * as.integer(edu_f == LV[3]) + B[["eduHigh"]] * as.integer(edu_f == LV[4])]
  d <- d[!is.na(score) & !is.na(lp) & !is.na(age0) & age0 >= 45 & !is.na(mem1) & !is.na(cutoff)]
  d[]
}
W1 <- make_win(2011, 2015); W2 <- make_win(2015, 2018)

rowU <- function(d, co, win) {
  data.table(cohort = co, window = win, definition = "U full population", n = nrow(d),
             n90 = sum(d$age0 >= 90), pct90 = round(100 * mean(d$age0 >= 90), 2),
             events_main_manuscript_threshold_4 = sum(d$mem1 <= K), events_education_strata = sum(d$mem1 < d$cutoff),
             age_range = paste0(min(d$age0), "-", max(d$age0)))
}
rowN <- function(d, co, win) {
  s <- d[!is.na(mem0) & mem0 >= cutoff]
  data.table(cohort = co, window = win, definition = "N unimpaired at baseline", n = nrow(s),
             n90 = sum(s$age0 >= 90), pct90 = round(100 * mean(s$age0 >= 90), 2),
             events_main_manuscript_threshold_4 = sum(s$mem1 <= K), events_education_strata = sum(s$mem1 < s$cutoff),
             age_range = paste0(min(s$age0), "-", max(s$age0)))
}
ch <- rbind(rowU(W1, "CHARLS", "2011->2015"), rowN(W1, "CHARLS", "2011->2015"),
            rowU(W2, "CHARLS", "2015->2018"), rowN(W2, "CHARLS", "2015->2018"))
say("========== CHARLS ==========")
pt(ch); fwrite(ch, file.path(OUT, "F20_charls.csv"), bom = TRUE)

## ---------- external cohort ----------
mk_ext <- function(co) {
  d <- as.data.table(readRDS(file.path(TMP, paste0("F16_base_", co, ".rds"))))
  d <- d[!is.na(age) & age >= 45]
  cs <- d[!is.na(mem_base)]
  d[, ci := floor(mean(cs$mem_base) - sd(cs$mem_base))]
  d[, score := 2L * pmax(floor((age - 45) / 5), 0) + edu_pts + 3L * dep]
  d[, lp := B[["intercept"]] + B[["age"]] * age + B[["dep"]] * dep + B[["adl"]] * adl +
      B[["eduPri"]] * as.integer(edu_pts == -2) + B[["eduMid"]] * as.integer(edu_pts == -7) +
      B[["eduHigh"]] * as.integer(edu_pts == -10)]
  # and F16 of U consistent definition: requires calculable score and available follow up memory only(not required ADL)
  d[!is.na(score) & !is.na(mem_fu)]
}
ext_row <- function(d, co, lab) {
  data.table(cohort = co, window = "baseline->follow up", definition = lab, n = nrow(d),
             n90 = sum(d$age >= 90), pct90 = round(100 * mean(d$age >= 90), 2),
             events_main_manuscript_threshold_4 = sum(d$mem_fu <= K), events_own_cutoff = sum(d$mem_fu <= d$ci),
             age_range = paste0(min(d$age), "-", max(d$age)))
}
ext <- rbindlist(lapply(c("HRS", "ELSA", "SHARE"), function(co) {
  d <- mk_ext(co)
  rbind(ext_row(d, co, "U full population"),
        ext_row(d[!is.na(mem_base) & mem_base > ci], co, "N unimpaired at baseline"))
}))
say("\n========== external cohort ==========")
pt(ext); fwrite(ext, file.path(OUT, "F20_external.csv"), bom = TRUE)

## ---------- older age subgroup(>=90)own events ----------
sub90 <- function(d, memvar, agevar, co) {
  s <- d[!is.na(d[[agevar]]) & d[[agevar]] >= 90 & !is.na(d[[memvar]])]
  data.table(cohort = co, n = nrow(s), events_main_manuscript_threshold_4 = sum(s[[memvar]] <= K),
             event_rate = round(100 * mean(s[[memvar]] <= K), 1),
             age_range = paste0(min(s[[agevar]]), "-", max(s[[agevar]])))
}
col1 <- sub90(W1, "mem1", "age0", "CHARLS 2011->2015")
col2 <- sub90(W2, "mem1", "age0", "CHARLS 2015->2018")
col3 <- rbindlist(lapply(c("HRS", "ELSA", "SHARE"), function(co) {
  sub90(mk_ext(co), "mem_fu", "age", co)
}))
hh <- rbind(col1, col2, col3)
say("\n========== older age subgroup(>=90 years old)own ==========")
pt(hh); fwrite(hh, file.path(OUT, "F20_90plus_subgroup.csv"), bom = TRUE)

close(con)
cat("\n complete ->", file.path(OUT, "F20_report.txt"), "\n")
