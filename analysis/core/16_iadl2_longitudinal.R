# ============================================================================
# 16_iadl2_longitudinal.R —— IADL(two item version iadl2)supplementary longitudinal model analysis
# purpose: is A02/A07 decision[IADL included in primary model]provides longitudinal evidence
# specification: and T06(13_weighted_longitudinal.R)consistent——single exposure + age + female + edu_group; 
#       primary sensitivity: ten year change model additionally adjusts baseline ms10(regression to the mean/floor effect)
# generate: 2026-09-26
# output: 
#   output/04_longitudinal_cohort/iadl2_longitudinal_models_2026-09-26.csv
#   output/04_longitudinal_cohort/iadl2_longitudinal_summary_2026-09-26.md
# ============================================================================
invisible(try(Sys.setlocale("LC_ALL", "zh_CN.UTF-8"), silent = TRUE))
.libPaths(c("D:/R/library", .libPaths()))
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(survival); library(survey) })
P <- "D:/BaiduSyncdisk/research/Xu_Zhonghao"; OUT <- file.path(P, "output", "04_longitudinal_cohort")

coh <- readRDS(file.path(OUT, "cohort_change_2026-09-13.rds"))$cohort
inc <- readRDS(file.path(OUT, "interval_incidence_events_2026-09-20.rds"))

b11 <- coh %>% filter(year == 2011) %>%
  transmute(KEY, age, female, edu = edu_group, iadl2, adl, ms10_2011 = ms10, weight, communityID)
delta_tab <- coh %>% select(KEY, year, ms10) %>%
  pivot_wider(names_from = year, values_from = ms10, names_prefix = "y") %>%
  transmute(KEY, delta = y2020 - y2011, ms10_2020 = y2020)

d_inc  <- inc$base %>% inner_join(b11, by = "KEY")
yr <- c("2013" = 2, "2015" = 4, "2018" = 7, "2020" = 9)
im <- cbind(d_inc$im2013, d_inc$im2015, d_inc$im2018, d_inc$im2020)
d_inc$time  <- apply(im, 1, function(z) { i <- which(z == 1); if (length(i)) yr[i[1]] else 9 })
d_inc$event <- as.integer(apply(im, 1, function(z) any(z == 1, na.rm = TRUE)))
d_chg <- b11 %>% inner_join(delta_tab, by = "KEY")

get_est <- function(fit, term, ratio = TRUE) {
  b <- coef(fit)[term]; se <- sqrt(vcov(fit)[term, term])
  if (ratio) list(est = exp(b), lo = exp(b - 1.96*se), hi = exp(b + 1.96*se))
  else list(est = b, lo = b - 1.96*se, hi = b + 1.96*se)
}
fit_all <- function(e, label) {
  rhs  <- paste(e, "+ age + female + edu")
  rhs2 <- paste(rhs, "+ ms10_2011")
  des  <- svydesign(ids = ~communityID, weights = ~weight, data = d_inc)
  des2 <- svydesign(ids = ~communityID, weights = ~weight, data = d_inc)
  des3 <- svydesign(ids = ~communityID, weights = ~weight, data = d_chg)
  out <- list(
    logistic     = c(get_est(svyglm(as.formula(paste("endpoint_2020 ~", rhs)), des, family = quasibinomial()), e),
                     get_est(glm(as.formula(paste("endpoint_2020 ~", rhs)), d_inc, family = binomial()), e)),
    cox          = c(get_est(svycoxph(as.formula(paste("Surv(time, event) ~", rhs)), des2), e),
                     get_est(coxph(as.formula(paste("Surv(time, event) ~", rhs)), d_inc), e)),
    change       = c(get_est(svyglm(as.formula(paste("delta ~", rhs)), des3), e, FALSE),
                     get_est(lm(as.formula(paste("delta ~", rhs)), d_chg), e, FALSE)),
    change_adj   = c(get_est(svyglm(as.formula(paste("delta ~", rhs2)), des3), e, FALSE),
                     get_est(lm(as.formula(paste("delta ~", rhs2)), d_chg), e, FALSE))
  )
  lapply(out, function(z) z)
}
fmt <- function(z) sprintf("%.2f (%.2f-%.2f)", z[[1]], z[[2]], z[[3]])
res <- list()
for (e in c("iadl2", "adl")) {
  r <- fit_all(e, e)
  vars_inc <- c(e, "age", "female", "edu", "endpoint_2020")
  vars_surv <- c(e, "age", "female", "edu", "time", "event")
  vars_chg <- c(e, "age", "female", "edu", "delta")
  n_inc  <- sum(complete.cases(d_inc[, vars_inc]))
  n_surv <- sum(complete.cases(d_inc[, vars_surv]))
  n_chg  <- sum(complete.cases(d_chg[, vars_chg]))
  n_lookup <- c(logistic = n_inc, cox = n_surv, change = n_chg, change_adj = n_chg)
  for (m in names(r)) {
    z <- r[[m]]
    res[[length(res)+1]] <- data.frame(exposure = e, model = m, n = n_lookup[[m]],
      weighted = fmt(list(z[[1]], z[[2]], z[[3]])), unweighted = fmt(list(z[[4]], z[[5]], z[[6]])))
  }
}
res <- do.call(rbind, res)
write.csv(res, file.path(OUT, "iadl2_longitudinal_models_2026-09-26.csv"), row.names = FALSE, fileEncoding = "UTF-8")

desc <- d_chg %>% group_by(iadl2) %>%
  summarise(n = n(), base = round(mean(ms10_2011, na.rm = TRUE), 2), follow = round(mean(ms10_2020, na.rm = TRUE), 2),
            delta = round(mean(delta, na.rm = TRUE), 2), dec1 = sprintf("%.1f%%", 100*mean((ms10_2020 - ms10_2011) <= -1, na.rm = TRUE)), .groups = "drop")
md <- c("# IADL(two item version iadl2)supplementary longitudinal analysis(2026-09-26)", "",
  "> purpose: is[IADL included in primary model]provides longitudinal evidence(previously cross sectional only OR 1.53–2.41).",
  "> specification: and T06 consistent(single exposure + age + sex + education); weighted＝2011 baseline weights + PSU=communityID.", "",
  "## one, descriptive(by 2011 year iadl2)", "", paste(capture.output(print(as.data.frame(desc), row.names = FALSE)), collapse = "\n"), "",
  "## two, model results(weighted / unweighted)", "", paste(capture.output(print(res, row.names = FALSE)), collapse = "\n"), "",
  "## three, interpretation", "",
  "- binary outcome: OR 1.32(0.99–1.75, weighted)/ 1.28(1.00–1.65, unweighted); HR 1.19(1.00–1.41)/ 1.17(1.00–1.38)——**both borderline**, weaker than ADL(OR 1.49, HR 1.20).",
  "- continuous change: without baseline adjustment β=+0.57(apparently[less decline]), **baseline adjustment ms10 after β≈−0.01(p=0.93)** → this[effect]yes**regression to mean artifact from baseline differences**(IADL baseline in limited group ms10 only 6.90 vs 8.19).",
  "- suggested conclusion: **IADL exclude from primary score**; as supplementary/report as secondary analysis(cross sectional evidence + borderline longitudinal association + substantial baseline differences).",
  "", "*scripts: code/16_iadl2_longitudinal.R; data: cohort_change_2026-09-13.rds + interval_incidence_events_2026-09-20.rds*")
writeLines(md, file.path(OUT, "iadl2_longitudinal_summary_2026-09-26.md"), useBytes = TRUE)
cat("\n==== results ====\n"); print(res, row.names = FALSE)
cat("\n==== descriptive ====\n"); print(as.data.frame(desc), row.names = FALSE)
