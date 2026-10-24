# ============================================================================
# F02b: memory change comparisons between adjacent waves + 2020 results using own cutoff
# purpose: determine 2018→2020 of +0.92 upward score shift reflects"word list/practice effect"or true change
# output: C:/y03tmp/out/F02b_report.txt
# ============================================================================
invisible(try(Sys.setlocale("LC_ALL", "zh_CN.UTF-8"), silent = TRUE))
suppressPackageStartupMessages({ library(data.table); library(pROC) })
if (dir.exists("D:/R/library")) .libPaths(c("D:/R/library", .libPaths()))
TMP <- "C:/y03tmp/tmp"; OUT <- "C:/y03tmp/out"
con <- file(file.path(OUT, "F02b_report.txt"), open = "wt", encoding = "UTF-8")
say <- function(...) { s <- paste0(..., collapse = ""); cat(s, "\n"); cat(s, "\n", file = con) }
pt <- function(x) say(paste(capture.output(print(x)), collapse = "\n"))

w <- as.data.table(readRDS(file.path(TMP, "waves_cognition.rds")))
g <- w[, .(ID = IDn, year, mem = as.numeric(imrc) + as.numeric(dlrc))]

## ---------- 1. each adjacent wave pair: memory change within same participants ----------
say("========== 1. adjacent wave memory change(same participants)==========")
prs <- list(c(2011, 2013), c(2013, 2015), c(2015, 2018), c(2018, 2020))
out <- data.table()
for (p in prs) {
  a <- g[year == p[1], .(ID, mem_a = mem)]
  b <- g[year == p[2], .(ID, mem_b = mem)]
  d <- merge(a, b, by = "ID")
  d <- d[!is.na(mem_a) & !is.na(mem_b)]
  out <- rbind(out, data.table(
    waves = paste0(p[1], "→", p[2]), interval_in_years = p[2] - p[1], n = nrow(d),
    preceding_wave_mean = round(mean(d$mem_a), 2), subsequent_wave_mean = round(mean(d$mem_b), 2),
    change = round(mean(d$mem_b - d$mem_a), 2), change_SD = round(sd(d$mem_b - d$mem_a), 2),
    preceding_wave_SD = round(sd(d$mem_a), 2), subsequent_wave_SD = round(sd(d$mem_b), 2)))
}
pt(out)
say("\n(note: 2018, 2020 all are 3 learning trials; 2011/2015 is 1 trials, 2013 is 3 times)")

## ---------- 2. memory means by wave/SD and mean−1SD cutoff ----------
say("\n========== 2. wave memory distributions and mean−1SD cutoff ==========")
d18 <- g[year == 2018]; d20 <- g[year == 2020]
tab2 <- data.table(
  waves = c(2011, 2013, 2015, 2018, 2020),
  scorable_n = sapply(list(2011, 2013, 2015, 2018, 2020), function(y) g[year == y & !is.na(mem), .N]),
  mean = sapply(list(2011, 2013, 2015, 2018, 2020), function(y) round(g[year == y, mean(mem, na.rm = TRUE)], 2)),
  SD = sapply(list(2011, 2013, 2015, 2018, 2020), function(y) round(g[year == y, sd(mem, na.rm = TRUE)], 2)))
tab2[, mean_minus_1SD := round(mean - SD, 2)]
tab2[, rounded_cutoff := floor(mean - SD)]
pt(tab2)

## ---------- 3. 2018→2020: results under three cutoffs ----------
say("\n========== 3. 2018→2020 risk populations and events under three cutoffs ==========")
beta <- c("(Intercept)" = -4.19524508778947, age = 0.0492039250378111,
          eduP = -0.290105771000328, eduM = -0.994711800202235, eduH = -1.43588388035287,
          dep = 0.412166484954086, adl = -0.0110119344848014)
a18 <- w[year == 2018, .(ID = IDn, mem18 = as.numeric(imrc) + as.numeric(dlrc),
                         age18 = as.numeric(age), edu18 = edu_group,
                         cesd18 = as.numeric(cesd10), adl18 = as.numeric(adl))]
a20 <- w[year == 2020, .(ID = IDn, mem20 = as.numeric(imrc) + as.numeric(dlrc))]
m <- merge(a18, a20, by = "ID")
m <- m[!is.na(mem18) & !is.na(mem20) & !is.na(age18) & !is.na(edu18) &
         !is.na(cesd18) & !is.na(adl18)]
m[, dep := as.integer(cesd18 >= 10)]
m[, edu_group := factor(edu18, levels = c("Less than primary", "Primary school",
                                          "Middle school", "High school or above"))]
m[, lp := beta[["(Intercept)"]] + beta[["age"]] * age18 + beta[["dep"]] * dep +
      beta[["adl"]] * adl18 +
      beta[["eduP"]] * as.integer(edu_group == "Primary school") +
      beta[["eduM"]] * as.integer(edu_group == "Middle school") +
      beta[["eduH"]] * as.integer(edu_group == "High school or above")]
m[, phat := plogis(lp)]

res <- data.table()
for (spec in list(list("2018 own mean−1SD", floor(mean(m$mem18) - sd(m$mem18))),
                  list("apply unchanged 2011 cutoff", 4L),
                  list("2020 own mean−1SD", floor(mean(m$mem20) - sd(m$mem20))))) {
  k <- spec[[2]]
  d <- m[mem18 > k]
  d[, event := as.integer(mem20 <= k)]
  r <- roc(d$event, d$phat, quiet = TRUE, direction = "<")
  ci <- as.numeric(ci.auc(r))
  cal <- glm(event ~ lp, data = d, family = binomial())
  res <- rbind(res, data.table(cutoff_definition = spec[[1]], cutoff = k, population_at_risk = nrow(d),
                               events = sum(d$event), event_rate = round(100 * mean(d$event), 2),
                               AUC = round(ci[2], 3), calibration_slope = round(coef(cal)[2], 3),
                               calibration_intercept = round(coef(cal)[1], 3)))
}
pt(res)
say("\n2020 counts by score category(same participants ", nrow(m), " participants): ")
pt(m[, .(n = .N, fraction = round(100 * .N / nrow(m), 1)), by = .(category = cut(mem20, c(-Inf, 3, 4, 5, Inf)))])

close(con)
cat("complete ->", file.path(OUT, "F02b_report.txt"), "\n")
