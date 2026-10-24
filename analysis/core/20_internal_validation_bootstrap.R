# ============================================================================
# 20_internal_validation_bootstrap.R —— A08: Bootstrap internal validation and optimism correction
# basis: SAP_v3 §8(Bootstrap 1,000 times; optimism correction AUC, calibration slope, Brier)
# model: incident_2015 ~ age + edu_group + dep + adl(options 2 primary model)
# method: Harrell optimism correction(bootstrap within sample performance − bootstrap model performance in original sample)
# generate: 2026-09-26
# ============================================================================
invisible(try(Sys.setlocale("LC_ALL", "zh_CN.UTF-8"), silent = TRUE))
.libPaths(c("D:/R/library", .libPaths()))
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(pROC) })
P <- "D:/BaiduSyncdisk/research/Xu_Zhonghao"; OUT <- file.path(P, "output", "04_longitudinal_cohort")
DIR <- file.path(P, "output", "10_screening_index")
set.seed(20260926)

coh <- readRDS(file.path(OUT, "cohort_change_2026-09-13.rds"))$cohort
w <- coh %>% select(KEY, year, imrc, dlrc) %>% pivot_wider(names_from = year, values_from = c(imrc, dlrc)) %>%
  mutate(mem11 = imrc_2011 + dlrc_2011, mem15 = imrc_2015 + dlrc_2015)
b <- coh %>% filter(year == 2011) %>% transmute(KEY, age, edu_group, dep = as.integer(cesd10 >= 10), adl) %>% left_join(w, by = "KEY")
b$edu_group <- factor(b$edu_group, levels = c("Less than primary","Primary school","Middle school","High school or above"))
d <- b[!is.na(b$mem11) & !is.na(b$mem15) & b$mem11 > 4 & complete.cases(b[, c("age","edu_group","dep","adl")]), ]
d$inc <- as.integer(d$mem15 <= 4)
cat("development sample n =", nrow(d), "｜events =", sum(d$inc), "\n")

f <- inc ~ age + edu_group + dep + adl
perf <- function(data, lp) {                     # lp = linear predictor(calculated by model)
  p <- 1/(1+exp(-lp))
  auc <- as.numeric(auc(roc(data$inc, p, quiet = TRUE)))
  cal <- glm(data$inc ~ lp, family = binomial())$coef          # intercept + slope(slope is number 2 items)
  brier <- mean((p - data$inc)^2)
  c(auc = auc, slope = as.numeric(cal[2]), intercept = as.numeric(cal[1]), brier = brier)
}
m_full <- glm(f, data = d, family = binomial())
lp_full <- predict(m_full, type = "link")
app <- perf(d, lp_full)
cat("\n apparent performance: ", sprintf("AUC=%.4f calibration slope=%.3f Brier=%.4f", app["auc"], app["slope"], app["brier"]), "\n")

B <- 1000; n <- nrow(d); opt <- matrix(NA_real_, B, 4, dimnames = list(NULL, names(app)))
t0 <- Sys.time()
for (i in 1:B) {
  idx <- sample.int(n, n, replace = TRUE)
  db <- d[idx, ]
  mb <- try(glm(f, data = db, family = binomial()), silent = TRUE)
  if (inherits(mb, "try-error")) next
  lp_b_app  <- predict(mb, newdata = db, type = "link")
  lp_b_test <- predict(mb, newdata = d,  type = "link")
  pb_app  <- perf(db, lp_b_app)
  pb_test <- perf(d,  lp_b_test)
  opt[i, ] <- pb_app - pb_test
}
opt_m <- colMeans(opt, na.rm = TRUE)
corr <- app - opt_m
cat("\n optimism(mean, B =", sum(!is.na(opt[,1])), "): ", sprintf("AUC %.4f｜slope %.4f｜Brier %.5f", opt_m["auc"], opt_m["slope"], opt_m["brier"]), "\n")
cat("after optimism correction: ", sprintf("AUC=%.4f  calibration slope=%.3f  Brier=%.4f", corr["auc"], corr["slope"], corr["brier"]), "\n")

## adjusted AUC percentile of CI(using bootstrap test sample AUC distribution approximates)
test_aucs <- app["auc"] - opt[, "auc"]
ci <- quantile(test_aucs, c(0.025, 0.975), na.rm = TRUE)
cat("adjusted AUC 95%CI(percentile method): ", sprintf("%.4f – %.4f", ci[1], ci[2]), "\n")

res <- data.frame(
  metrics = c("AUC", "calibration_slope", "calibration_intercept", "Brier"),
  apparent = round(app[c("auc","slope","intercept","brier")], 4),
  optimism = round(opt_m[c("auc","slope","intercept","brier")], 4),
  adjusted = round(corr[c("auc","slope","intercept","brier")], 4))
print(res, row.names = FALSE)
write.csv(res, file.path(DIR, "A08_bootstrap_optimism_2026-09-26.csv"), row.names = FALSE, fileEncoding = "UTF-8")

md <- c("# A08: Bootstrap internal validation and optimism correction(2026-09-26)", "",
  sprintf("> model: incident_2015 ~ age + edu_group + dep + adl(options 2); development sample n = %d, events = %d(%.1f%%)",
          nrow(d), sum(d$inc), 100*mean(d$inc)),
  sprintf("> Bootstrap = %d times(Harrell optimism correction); random seed 20260926; runtime %.1f minutes", B, as.numeric(difftime(Sys.time(), t0, units = "mins"))), "",
  "## results", "",
  paste(capture.output(print(res, row.names = FALSE)), collapse = "\n"), "",
  sprintf("**adjusted AUC = %.3f(95%%CI %.3f–%.3f)**, calibration slope = %.2f, Brier = %.4f.",
          corr["auc"], ci[1], ci[2], corr["slope"], corr["brier"]), "",
  "## interpretation", "",
  "- optimism(apparent − corrected)depends on model complexity: this model has only 6 parameters, sufficient events(658), optimism should be small; ",
  "- calibration slope < 1 indicates overfitting/prediction distribution too narrow, A09 calibration will use this(isotonic or Platt); ",
  "- adjusted AUC and apparent AUC difference is the optimism to report.")
writeLines(md, file.path(DIR, "A08_bootstrap_optimism_2026-09-26.md"), useBytes = TRUE)
saveRDS(list(apparent = app, optimism = opt_m, corrected = corr, ci_auc = ci, opt_matrix = opt),
        file.path(DIR, "A08_bootstrap_2026-09-26.rds"))
cat("\n saved to", DIR, "\n")
