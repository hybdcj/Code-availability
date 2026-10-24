# ============================================================================
# 05_cohort_risk_factors.R —— fixed cohort: baseline risk factors → ten year cognitive decline
# version: v1.0(2026-09-13)
# input: output/04_longitudinal_cohort/cohort_change_2026-09-13.rds(from code/04_cohort_change.R generate)
# output: output/04_longitudinal_cohort/cohort_risk_2026-09-13.csv(summary of four model results)
# four analyses(all use 2011 baseline exposures as predictors, adjust for age, sex, education): 
#   ① incident impairment(Logistic, OR): unimpaired baseline participants 2011→2020 become impaired
#   ② time to first impairment(Cox, HR): 2013/2015/2018/2020 four follow up waves
#   ③ rate of cognitive decline(linear mixed model, nlme): time×exposure interaction term
#   ④ ten year change score(linear regression, β): 2020 − 2011 score difference
# notes: 2011 unavailable in this wave IADL, therefore baseline exposures exclude IADL; sleep in hours, pain/binary physical activity.
# ============================================================================
invisible(try(Sys.setlocale("LC_ALL", "zh_CN.UTF-8"), silent = TRUE))
if (dir.exists("D:/R/library")) .libPaths(c("D:/R/library", .libPaths()))
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(survival); library(nlme) })

P <- Sys.getenv("CHARLS_PROJ", "D:/BaiduSyncdisk/research/Xu_Zhonghao")
OUT <- file.path(P, "output", "04_longitudinal_cohort")
obj <- readRDS(file.path(OUT, "cohort_change_2026-09-13.rds"))
coh <- obj$cohort

## ---- baseline exposures(2011) ----
base <- coh %>% filter(year == 2011) %>%
  transmute(KEY, age, female, edu = edu_group,
            dep = as.integer(cesd10 >= 10), adl, sleep_h,
            pain = pain_any, pa = pa_any, hibp, diab,
            ms10_2011 = ms10, impaired_2011 = impaired_fixed)

## ---- follow up outcome ----
# ① incident impairment(2011 normal → 2020 impaired)
inc <- base %>% filter(impaired_2011 == 0) %>%
  left_join(coh %>% filter(year == 2020) %>% select(KEY, impaired_2020 = impaired_fixed), by = "KEY")
# ② time to first impairment(year): follow up time point 2013=2, 2015=4, 2018=7, 2020=9
years <- c("2013" = 2, "2015" = 4, "2018" = 7, "2020" = 9)
fu <- coh %>% filter(year != 2011) %>% select(KEY, year, impaired_fixed) %>%
  mutate(t = years[as.character(year)]) %>% arrange(KEY, t)
surv_df <- fu %>% group_by(KEY) %>% summarise(
  time = if (any(impaired_fixed == 1, na.rm = TRUE)) min(t[impaired_fixed == 1], na.rm = TRUE) else 9,
  event = as.integer(any(impaired_fixed == 1, na.rm = TRUE)), .groups = "drop")
surv_df <- base %>% filter(impaired_2011 == 0) %>% left_join(surv_df, by = "KEY")
# ④ ten year change score
delta <- coh %>% select(KEY, year, ms10) %>%
  pivot_wider(names_from = year, values_from = ms10, names_prefix = "y") %>%
  transmute(KEY, delta = y2020 - y2011)
base <- base %>% left_join(delta, by = "KEY")

exposures <- c("dep", "adl", "sleep_h", "pain", "pa", "hibp", "diab")
labels <- c(dep = "depression(CES-D-10 ≥10)", adl = "ADL limited", sleep_h = "sleep duration(per 1 hours)",
            pain = "pain", pa = "regular physical activity", hibp = "hypertension", diab = "diabetes")
adj <- " + age + female + edu"
fmt <- function(est, lo, hi, p) sprintf("%.2f (%.2f-%.2f), p=%s", est, lo, hi,
                                        ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))

res <- do.call(rbind, lapply(exposures, function(e) {
  # ① Logistic
  m1 <- glm(as.formula(paste("impaired_2020 ~", e, adj)), data = inc, family = binomial)
  s1 <- summary(m1)$coefficients[2, ]
  or <- exp(s1[1]); ci1 <- exp(s1[1] - 1.96*s1[2]); ci2 <- exp(s1[1] + 1.96*s1[2])
  # ② Cox
  m2 <- coxph(as.formula(paste("Surv(time, event) ~", e, adj)), data = surv_df)
  s2 <- summary(m2)$coefficients[e, ]; s2c <- summary(m2)$conf.int[e, ]
  # ③ mixed model: time × exposure
  d <- coh %>% select(KEY, year, ms10) %>% mutate(t = c("2011"=0,"2013"=2,"2015"=4,"2018"=7,"2020"=9)[as.character(year)]) %>%
    left_join(base %>% select(KEY, all_of(exposures), age, female, edu), by = "KEY") %>% filter(!is.na(ms10))
  m3 <- try(lme(as.formula(paste("ms10 ~ t *", e, adj)), random = ~ 1 + t | KEY, data = d, na.action = na.omit), silent = TRUE)
  if (inherits(m3, "try-error")) { beta <- NA; blo <- NA; bhi <- NA; bp <- NA } else {
    s3 <- summary(m3)$tTable; rn <- grep(paste0("^t:", e), rownames(s3))
    if (length(rn)) { beta <- s3[rn[1],1]; blo <- beta - 1.96*s3[rn[1],2]; bhi <- beta + 1.96*s3[rn[1],2]; bp <- s3[rn[1],5] }
    else { beta <- NA; blo <- NA; bhi <- NA; bp <- NA } }
  # ④ ten year change score
  m4 <- lm(as.formula(paste("delta ~", e, adj)), data = base)
  s4 <- summary(m4)$coefficients[2, ]
  data.frame(exposure = labels[[e]],
             incident_impairment_OR = fmt(or, ci1, ci2, s1[4]),
             first_impairment_HR = fmt(s2c[1], s2c[3], s2c[4], s2[5]),
             decline_rate_beta = ifelse(is.na(beta), "—", sprintf("%.3f (%.3f-%.3f), p=%s", beta, blo, bhi,
                                                        ifelse(bp < 0.001, "<0.001", sprintf("%.3f", bp)))),
             ten_year_change_beta = sprintf("%.3f (%.3f-%.3f), p=%s", s4[1], s4[1]-1.96*s4[2], s4[1]+1.96*s4[2],
                                     ifelse(s4[4] < 0.001, "<0.001", sprintf("%.3f", s4[4]))),
             row.names = NULL)
}))

cat("=== fixed cohort risk factor analysis(n=", nrow(base), "; incidence analysis n=", nrow(inc), ")===\n", sep = "")
print(res, row.names = FALSE)
write.csv(res, file.path(OUT, "cohort_risk_2026-09-13.csv"), row.names = FALSE, fileEncoding = "UTF-8")
cat("\n event count: cognitive impairment during follow up ", sum(surv_df$event, na.rm = TRUE), " cases(",
    sprintf("%.1f%%", 100*mean(surv_df$event, na.rm = TRUE)), ")\n", sep = "")
cat("saved: output/04_longitudinal_cohort/cohort_risk_2026-09-13.csv\n")
