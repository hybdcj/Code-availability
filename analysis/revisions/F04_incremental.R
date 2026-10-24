# ============================================================================
# F04 plan A: same participants(2018→2020)incremental variable comparison
# purpose: originally classified as"incremental/exploration/exclude"candidate variables, at"memory score threshold crossing"this outcome, 
#       in the past decade population, whether they improve the existing four item model
# prespecified inventory: M0 baseline; M1 sleep; M2 physical activity; M3 IADL; M4 chronic disease burden; 
#       M5 sex+urban or rural; M6 stroke+diabetes; M7 add all four domains; M8 continuous depression variable; M9 spline for age
# adoption criterion(prespecified): adjusted ΔAUC ≥ 0.010 and nested Bootstrap 95% lower bound ≥ 0.005
# output: C:/y03tmp/out/F04_report.txt, F04_*.csv, F04_fig_dAUC.png
# ============================================================================
invisible(try(Sys.setlocale("LC_ALL", "zh_CN.UTF-8"), silent = TRUE))
suppressPackageStartupMessages({ library(data.table); library(pROC); library(splines) })
if (dir.exists("D:/R/library")) .libPaths(c("D:/R/library", .libPaths()))
TMP <- "C:/y03tmp/tmp"; OUT <- "C:/y03tmp/out"
con <- file(file.path(OUT, "F04_report.txt"), open = "wt", encoding = "UTF-8")
say <- function(...) { s <- paste0(..., collapse = ""); cat(s, "\n"); cat(s, "\n", file = con) }
pt  <- function(x) say(paste(capture.output(print(x)), collapse = "\n"))

BETA_OLD <- c("(Intercept)" = -4.19524508778947, age = 0.0492039250378111,
              eduP = -0.290105771000328, eduM = -0.994711800202235, eduH = -1.43588388035287,
              dep = 0.412166484954086, adl = -0.0110119344848014)
LV <- c("Less than primary", "Primary school", "Middle school", "High school or above")

## set F04_SKIP_BOOT=1 can skip Bootstrap(reuse saved F04_bootstrap.csv), regenerate tables and figures only
SKIP_BOOT <- nzchar(Sys.getenv("F04_SKIP_BOOT"))

## ---------------- data ----------------
w <- as.data.table(readRDS(file.path(TMP, "waves_cognition.rds")))
a18 <- w[year == 2018, .(ID = IDn, mem18 = as.numeric(imrc) + as.numeric(dlrc),
                         age18 = as.numeric(age), female = as.numeric(female),
                         rural = as.numeric(rural), edu18 = edu_group,
                         cesd18 = as.numeric(cesd10), adl18 = as.numeric(adl),
                         iadl18 = as.numeric(iadl), iadl2_18 = as.numeric(iadl2),
                         sleep18 = as.numeric(sleep_h), pa18 = as.numeric(pa_any),
                         pain18 = as.numeric(pain_any),
                         hibp = as.numeric(hibp), diab = as.numeric(diab),
                         heart = as.numeric(heart), stroke = as.numeric(stroke),
                         lung = as.numeric(lung), kidney = as.numeric(kidney),
                         dyslip = as.numeric(dyslip),
                         wt18 = as.numeric(weight))]
a20 <- w[year == 2020, .(ID = IDn, mem20 = as.numeric(imrc) + as.numeric(dlrc))]
m <- merge(a18, a20, by = "ID")
m <- m[!is.na(mem18) & !is.na(mem20) & !is.na(age18) & !is.na(edu18) &
         !is.na(cesd18) & !is.na(adl18) & age18 >= 45]

## prespecified variable preparation
m[, dep := as.integer(cesd18 >= 10)]
m[, edu_group := factor(edu18, levels = LV)]
m[, sleep3 := fifelse(sleep18 < 6, "short (<6h)", fifelse(sleep18 > 8, "long (>8h)", "reference (6-8h)"))]
m[!is.na(sleep18) & (sleep18 < 2 | sleep18 > 16), sleep3 := NA_character_]
m[, sleep3 := factor(sleep3, levels = c("reference (6-8h)", "short (<6h)", "long (>8h)"))]
m[, rural01 := fifelse(rural == 1, 1L, fifelse(rural == 0, 0L, NA_integer_))]
ccol <- c("hibp", "diab", "heart", "stroke", "lung", "kidney", "dyslip")
mm <- as.matrix(m[, ..ccol])
m[, n_cc := rowSums(mm, na.rm = TRUE)]
m[rowSums(!is.na(mm)) < 5, n_cc := NA_real_]

say("analysis population n = ", nrow(m))
say("missing candidate variables: sleep ", sum(is.na(m$sleep3)), "｜physical activity ", sum(is.na(m$pa18)),
    "｜IADL ", sum(is.na(m$iadl2_18)), "｜chronic disease burden ", sum(is.na(m$n_cc)),
    "｜urban or rural ", sum(is.na(m$rural01)))

## ---------------- outcome definition(retain F03) ----------------
cut18 <- floor(mean(m$mem18) - sd(m$mem18))
p_anchor <- mean(m$mem18 <= cut18)
cut20 <- as.numeric(quantile(m$mem20, p_anchor, type = 1))
d_all <- copy(m[mem18 > cut18])[, event := as.integer(mem20 <= cut18)]
say("\nD1 definition: 2018 mean=", round(mean(m$mem18), 2), " SD=", round(sd(m$mem18), 2),
    " → cutoff ≤", cut18, "; population at risk n=", nrow(d_all), ", events=", sum(d_all$event),
    "(", round(100*mean(d_all$event), 2), "%)")

## ensure model comparability, all models in"complete common candidate variable sample"fit and evaluate in same participants
need <- c("age18", "edu_group", "dep", "adl18", "cesd18", "sleep3", "pa18",
          "iadl2_18", "iadl18", "pain18", "n_cc", "female", "rural01", "stroke", "diab")
ok <- complete.cases(as.data.frame(d_all[, ..need]))
d <- d_all[ok]
say("common complete candidate variable sample n = ", nrow(d), "(events ", sum(d$event), ", ",
    round(100 * mean(d$event), 2), "%); excluded for missing candidate variables ",
    nrow(d_all) - nrow(d), " participants")
say("(reference: M0 in all ", nrow(d_all), " participant apparent performance AUC see F03 report＝0.706)")

## ---------------- prespecified models ----------------
FORMS <- list(
  M0 = event ~ age18 + edu_group + dep + adl18,
  M1 = event ~ age18 + edu_group + dep + adl18 + sleep3,
  M2 = event ~ age18 + edu_group + dep + adl18 + pa18,
  M3 = event ~ age18 + edu_group + dep + adl18 + iadl2_18,
  M4 = event ~ age18 + edu_group + dep + adl18 + n_cc,
  M5 = event ~ age18 + edu_group + dep + adl18 + female + rural01,
  M6 = event ~ age18 + edu_group + dep + adl18 + stroke + diab,
  M7 = event ~ age18 + edu_group + dep + adl18 + sleep3 + pa18 + iadl2_18 + n_cc,
  M8 = event ~ age18 + edu_group + cesd18 + adl18,
  M9 = event ~ ns(age18, df = 3) + edu_group + dep + adl18,
  M10 = event ~ age18 + edu_group + dep + adl18 + pain18,
  M11 = event ~ age18 + edu_group + dep + adl18 + iadl18
)
LBL <- c(M0 = "existing four items", M1 = "+sleep (three categories)", M2 = "+regular physical activity", M3 = "+IADL",
         M4 = "+chronic disease burden (0-7)", M5 = "+sex+urban or rural", M6 = "+stroke+diabetes",
         M7 = "+sleep+physical activity+IADL+chronic conditions", M8 = "continuous depression variable CES-D", M9 = "spline for age (3df)",
         M10 = "+pain", M11 = "+IADL(complete 6 items)")

fit_one <- function(f, data) tryCatch(glm(f, data = data, family = binomial()),
                                      error = function(e) NULL)
aucv <- function(y, p) {
  if (length(unique(y)) < 2) return(NA_real_)
  as.numeric(auc(roc(y, p, quiet = TRUE, direction = "<")))
}
wt_auc <- function(p, y, wgt) {
  cs <- y == 1; ct <- y == 0
  if (!any(cs) || !any(ct)) return(NA_real_)
  num <- 0
  for (i in which(cs)) num <- num + wgt[i] * sum(wgt[ct] * ((p[i] > p[ct]) + 0.5 * (p[i] == p[ct])))
  num / (sum(wgt[cs]) * sum(wgt[ct]))
}

## ---------------- apparent results ----------------
say("\n========== apparent results(D1 definition)==========")
res <- data.table(); ors <- data.table()
nms <- names(FORMS)
for (nm in nms) {
  fit <- fit_one(FORMS[[nm]], d)
  if (is.null(fit)) { say(nm, " fit failed"); next }
  d[[paste0("p_", nm)]] <- predict(fit, newdata = d, type = "response")
  r <- roc(d$event, d[[paste0("p_", nm)]], quiet = TRUE, direction = "<")
  ci <- as.numeric(ci.auc(r))
  lp <- qlogis(pmin(pmax(d[[paste0("p_", nm)]], 1e-9), 1 - 1e-9))
  cal <- glm(d$event ~ lp, family = binomial())
  res <- rbind(res, data.table(
    model = nm, notes = LBL[[nm]], parameter_count = length(coef(fit)),
    apparent_AUC = round(ci[2], 4), AUC_low = round(ci[1], 4), AUC_high = round(ci[3], 4),
    calibration_slope = round(coef(cal)[2], 3), calibration_intercept = round(coef(cal)[1], 3),
    Brier = round(mean((d[[paste0("p_", nm)]] - d$event)^2), 4),
    weighted_AUC = round(wt_auc(d[!is.na(wt18)][[paste0("p_", nm)]], d[!is.na(wt18)]$event,
                           d[!is.na(wt18)]$wt18), 4)))
  cf <- summary(fit)$coefficients
  base_terms <- c("(Intercept)", "age18", "edu_groupPrimary school", "edu_groupMiddle school",
                  "edu_groupHigh school or above", "dep", "adl18", "cesd18")
  for (t in setdiff(rownames(cf), base_terms)) {
    ors <- rbind(ors, data.table(model = nm, items = t, beta = round(cf[t, 1], 4),
                                 OR = round(exp(cf[t, 1]), 3),
                                 OR_low = round(exp(cf[t, 1] - 1.96*cf[t, 2]), 3),
                                 OR_high = round(exp(cf[t, 1] + 1.96*cf[t, 2]), 3),
                                 p = signif(cf[t, 4], 3)))
  }
}
pt(res)
fwrite(res, file.path(OUT, "F04_apparent.csv"), bom = TRUE)
say("\n additional term coefficients(D1 definition): "); pt(ors)
fwrite(ors, file.path(OUT, "F04_added_terms.csv"), bom = TRUE)

## ---------------- Bootstrap: optimism and ΔAUC ----------------
say("\n========== Bootstrap(B=1000): optimism correction and ΔAUC ==========")
app <- setNames(res$apparent_AUC, res$model)
btfile <- file.path(OUT, "F04_bootstrap.csv")
if (SKIP_BOOT && file.exists(btfile)) {
  bt <- fread(btfile, encoding = "UTF-8")
  say("(skipped Bootstrap, reuse ", basename(btfile), ")")
} else {
  Bset <- 1000; set.seed(20261014)
  bootAUC <- matrix(NA_real_, Bset, length(nms), dimnames = list(NULL, nms))
  optAUC  <- matrix(NA_real_, Bset, length(nms), dimnames = list(NULL, nms))
  for (b in seq_len(Bset)) {
    idx <- sample.int(nrow(d), replace = TRUE); db <- d[idx]
    if (length(unique(db$event)) < 2) next
    for (nm in nms) {
      f <- fit_one(FORMS[[nm]], db); if (is.null(f)) next
      p_b <- tryCatch(predict(f, newdata = db, type = "response"), error = function(e) NULL)
      p_o <- tryCatch(predict(f, newdata = d, type = "response"), error = function(e) NULL)
      if (is.null(p_b) || is.null(p_o)) next
      bootAUC[b, nm] <- aucv(db$event, p_b)
      optAUC[b, nm]  <- aucv(d$event, p_o)
    }
    if (b %% 200 == 0) cat("bootstrap", b, "/", Bset, "\n")
  }
  optimism <- app[nms] - colMeans(optAUC, na.rm = TRUE)
  corr <- app[nms] - optimism
  dA <- sweep(bootAUC, 1, bootAUC[, "M0"], "-")
  bt <- data.table(
    model = nms, notes = unlist(LBL[nms]),
    apparent_AUC = round(app[nms], 4),
    optimism = round(optimism, 4),
    adjusted_AUC = round(corr, 4),
    adjusted_dAUC = round(corr - corr[["M0"]], 4),
    dAUC_boot_median = round(apply(dA[, nms, drop = FALSE], 2, median, na.rm = TRUE), 4),
    dAUC_low = round(apply(dA[, nms, drop = FALSE], 2, quantile, .025, na.rm = TRUE), 4),
    dAUC_high = round(apply(dA[, nms, drop = FALSE], 2, quantile, .975, na.rm = TRUE), 4))
  bt[, criterion_met := fifelse(adjusted_dAUC >= 0.010 & dAUC_low >= 0.005, "yes", "no")]
  fwrite(bt, btfile, bom = TRUE)
}
pt(bt)

## ---------------- comparison with frozen weights ----------------
say("\n========== comparison with frozen weights ==========")
lpo <- BETA_OLD[["(Intercept)"]] + BETA_OLD[["age"]] * d$age18 +
       BETA_OLD[["dep"]] * d$dep + BETA_OLD[["adl"]] * d$adl18 +
       BETA_OLD[["eduP"]] * as.integer(d$edu_group == LV[2]) +
       BETA_OLD[["eduM"]] * as.integer(d$edu_group == LV[3]) +
       BETA_OLD[["eduH"]] * as.integer(d$edu_group == LV[4])
say("frozen(2011)AUC=", round(aucv(d$event, plogis(lpo)), 4),
    "｜M0 reestimate AUC=", round(app[["M0"]], 4),
    "｜best incremental model AUC=", round(max(app[nms]), 4),
    "(", nms[which.max(app[nms])], ")")

## ---------------- figure: ΔAUC ----------------
png(file.path(OUT, "F04_fig_dAUC.png"), width = 1150, height = 620, res = 120)
o <- order(bt$adjusted_dAUC); b <- bt[o]
par(mar = c(4.5, 13, 3.5, 1.5))
plot(b$adjusted_dAUC, seq_len(nrow(b)), xlim = c(-0.01, max(b$dAUC_high, 0.05) + 0.01),
     yaxt = "n", pch = 19, xlab = "Increase in AUC vs current 4-item model (optimism-corrected)",
     ylab = "", main = "Added value of previously screened-out candidate variables", cex = 1.1)
segments(b$dAUC_low, seq_len(nrow(b)), b$dAUC_high, seq_len(nrow(b)))
abline(v = 0, lty = 3); abline(v = 0.01, lty = 2, col = "red")
axis(2, at = seq_len(nrow(b)), labels = paste0(b$model, " ", b$notes), las = 1, cex.axis = 0.75)
legend("bottomright", c("95% nested bootstrap CI", "acceptance threshold 0.010"),
       lty = c(1, 2), col = c("black", "red"), bty = "n", cex = 0.8)
dev.off()

## ---------------- robustness: D2 equipercentile definition ----------------
say("\n========== robustness: D2 equipercentile definition(2020 ≤", cut20, ")==========")
d2 <- copy(d)[, event := as.integer(mem20 <= cut20)]
for (nm in c("M0", "M1", "M2", "M3", "M4", "M7", "M8")) {
  f <- fit_one(FORMS[[nm]], d2)
  if (!is.null(f)) {
    p <- predict(f, newdata = d2, type = "response")
    say(nm, ": AUC=", round(aucv(d2$event, p), 4))
  }
}
say("(D2: n=", nrow(d2), ", events=", sum(d2$event), ")")

close(con)
cat("\n complete ->", file.path(OUT, "F04_report.txt"), "\n")
