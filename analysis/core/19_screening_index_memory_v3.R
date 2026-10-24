# ============================================================================
# 19_screening_index_memory_v3.R —— A07(options 2): memory outcome screening index
# outcome: memory score(immediate+delayed recall, 0–20); fixed cutoff = 2011 entire cohort mean−1SD(retain across waves)
# design: longitudinal CHARLS 2011 → 2015(both waves use 1 learning trials, same word list, comparable)
# secondary: 2011 concurrent(cross sectional)same cutoff
# variables: age(per 5 years old), education(4 classes), depression(CES-D-10 ≥10), ADL(any limitation), all 2011 baseline
# unweighted primary analysis; weighted(2011 weights + PSU=communityID)for sensitivity analysis
# generate: 2026-09-26
# ============================================================================
invisible(try(Sys.setlocale("LC_ALL", "zh_CN.UTF-8"), silent = TRUE))
.libPaths(c("D:/R/library", .libPaths()))
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(survival); library(survey); library(pROC) })
P <- "D:/BaiduSyncdisk/research/Xu_Zhonghao"; OUT <- file.path(P, "output", "04_longitudinal_cohort")
DIR <- file.path(P, "output", "10_screening_index"); if (!dir.exists(DIR)) dir.create(DIR, recursive = TRUE)

coh <- readRDS(file.path(OUT, "cohort_change_2026-09-13.rds"))$cohort
w <- coh %>% select(KEY, year, imrc, dlrc) %>%
  pivot_wider(names_from = year, values_from = c(imrc, dlrc)) %>%
  mutate(mem2011 = imrc_2011 + dlrc_2011, mem2015 = imrc_2015 + dlrc_2015)

b <- coh %>% filter(year == 2011) %>%
  transmute(KEY, age, female, edu_group, dep = as.integer(cesd10 >= 10), adl, iadl2, ms10_2011 = ms10,
            weight, communityID) %>%
  left_join(w, by = "KEY")
b$edu_group <- factor(b$edu_group, levels = c("Less than primary","Primary school","Middle school","High school or above"))

mu <- mean(b$mem2011, na.rm = TRUE); sdd <- sd(b$mem2011, na.rm = TRUE); cut <- mu - sdd
cat(sprintf("memory score 2011: mean=%.2f SD=%.2f → fixed cutoff mean−1SD=%.2f(take <= %d)\n", mu, sdd, cut, floor(cut)))
b$impaired_2011 <- as.integer(b$mem2011 <= floor(cut))
# event definition: only among[unimpaired at baseline]assess incident events(baseline impaired participants excluded from risk population)
b$incident_2015 <- ifelse(is.na(b$mem2011) | is.na(b$mem2015) | b$mem2011 <= floor(cut),
                          NA_integer_,
                          as.integer(b$mem2015 <= floor(cut)))
cat(sprintf("concurrent impairment prevalence %.1f%%｜longitudinal incidence %.1f%%(follow up 4 year)\n\n",
    100*mean(b$impaired_2011, na.rm=TRUE), 100*mean(b$incident_2015, na.rm=TRUE)))

## analysis sample(primary: longitudinal)
d <- b %>% filter(!is.na(incident_2015)) %>% filter(complete.cases(age, edu_group, dep, adl))
cat("longitudinal analysis sample n =", nrow(d), "｜events =", sum(d$incident_2015), "\n")
d_conc <- b %>% filter(complete.cases(impaired_2011, age, edu_group, dep, adl))

## ---------- model ----------
f <- incident_2015 ~ age + edu_group + dep + adl
m_u <- glm(f, data = d, family = binomial())
m_w <- svyglm(f, design = svydesign(ids = ~communityID, weights = ~weight, data = d), family = quasibinomial())
cf <- data.frame(term = names(coef(m_u)),
  OR_unweighted = exp(coef(m_u)), ci_u = sprintf("%.2f (%.2f-%.2f)", exp(coef(m_u)),
    exp(coef(m_u) - 1.96*sqrt(diag(vcov(m_u)))), exp(coef(m_u) + 1.96*sqrt(diag(vcov(m_u))))),
  OR_weighted = exp(coef(m_w)[names(coef(m_u))]), ci_w = sprintf("%.2f (%.2f-%.2f)",
    exp(coef(m_w)[names(coef(m_u))]), exp(coef(m_w)[names(coef(m_u))] - 1.96*sqrt(diag(vcov(m_w)))[names(coef(m_u))]),
    exp(coef(m_w)[names(coef(m_u))] + 1.96*sqrt(diag(vcov(m_w)))[names(coef(m_u))])))
cat("\n== primary model(unweighted / weighted)==\n"); print(cf, row.names = FALSE)
write.csv(cf, file.path(DIR, "screening_index_memory_v3_models_2026-09-26.csv"), row.names = FALSE, fileEncoding = "UTF-8")

## ---------- integerization(SAP §7.3)----------
bb <- coef(m_u); b_age5 <- 5*bb[["age"]]
terms_beta <- c(age5 = b_age5, edu_primary = bb[["edu_groupPrimary school"]],
                edu_middle = bb[["edu_groupMiddle school"]], edu_high = bb[["edu_groupHigh school or above"]],
                dep = bb[["dep"]], adl = bb[["adl"]])
# revised rules(SAP v3.1): after scaling use |β| largest as 10 point reference, round remaining values proportionally
# (original[minimum |β| is 1 point reference]can inflate scores when a variable has negligible effect, see output/10_screening_index/ notes)
bmax <- max(abs(terms_beta)); pts <- round(10 * terms_beta / bmax)
score_tab <- data.frame(project = c("age(per 5 years old)","education: primary school","education: middle school","education_high_school_or_above","depression","ADL limited"),
  beta_5_years_equivalent = round(terms_beta, 3), per_score_unit = as.integer(pts),
  coding = c("45–49=0; per increase 5 years old +2 points","vs no schooling","vs no schooling","vs no schooling","CES-D-10 ≥10","any limitation(β≈0 → 0 points, see SAP clause)"))
print(score_tab, row.names = FALSE)
write.csv(score_tab, file.path(DIR, "screening_index_memory_v3_score_table_2026-09-26.csv"), row.names = FALSE, fileEncoding = "UTF-8")

d$score <- as.integer(pts[["age5"]]) * pmax(floor((d$age - 45)/5), 0) +
  as.integer(pts[["edu_primary"]]) * as.integer(d$edu_group == "Primary school") +
  as.integer(pts[["edu_middle"]]) * as.integer(d$edu_group == "Middle school") +
  as.integer(pts[["edu_high"]]) * as.integer(d$edu_group == "High school or above") +
  as.integer(pts[["dep"]]) * d$dep + as.integer(pts[["adl"]]) * d$adl
cat("\n total score range:", range(d$score), "\n")

## ---------- stratification ----------
q <- quantile(d$score, c(1/3, 2/3))
d$stratum <- cut(d$score, c(-Inf, q[1], q[2], Inf), labels = c("low","in","high"), right = TRUE)
strata <- d %>% group_by(stratum) %>% summarise(n = n(), events = sum(incident_2015),
  event_rate = sprintf("%.1f%%", 100*mean(incident_2015)), mean_score = round(mean(score),1), .groups="drop")
print(strata, row.names = FALSE)
write.csv(strata, file.path(DIR, "screening_index_memory_v3_risk_strata_2026-09-26.csv"), row.names = FALSE, fileEncoding = "UTF-8")

## ---------- performance ----------
d$p_hat <- predict(m_u, type = "response")
r1 <- roc(d$incident_2015, d$p_hat, quiet = TRUE); r2 <- roc(d$incident_2015, d$score, quiet = TRUE)
lp_hat <- qlogis(pmin(pmax(d$p_hat, 1e-6), 1 - 1e-6))      # calibration slope must be in logit scale calculation
cal <- glm(incident_2015 ~ lp_hat, data = d, family = binomial())
perf <- data.frame(metrics = c("AUC(continuous predicted probability)","AUC(integer total score)","calibration_intercept","calibration_slope","Brier"),
  value = c(sprintf("%.3f (%.3f-%.3f)", as.numeric(auc(r1)), ci.auc(r1)[1], ci.auc(r1)[3]),
         sprintf("%.3f (%.3f-%.3f)", as.numeric(auc(r2)), ci.auc(r2)[1], ci.auc(r2)[3]),
         sprintf("%.3f (%.3f-%.3f)", coef(cal)[1], coef(cal)[1]-1.96*summary(cal)$coef[1,2], coef(cal)[1]+1.96*summary(cal)$coef[1,2]),
         sprintf("%.2f (%.2f-%.2f)", coef(cal)[2], coef(cal)[2]-1.96*summary(cal)$coef[2,2], coef(cal)[2]+1.96*summary(cal)$coef[2,2]),
         sprintf("%.4f", mean((d$p_hat - d$incident_2015)^2))))
print(perf, row.names = FALSE)
write.csv(perf, file.path(DIR, "screening_index_memory_v3_performance_2026-09-26.csv"), row.names = FALSE, fileEncoding = "UTF-8")

## threshold
thr <- do.call(rbind, lapply(c(0.10, 0.15, 0.20), function(t) {
  pos <- d$p_hat >= t
  data.frame(cutoff_threshold = sprintf("predicted probability ≥%d%%", t*100), positive_count = sum(pos), positive_fraction = sprintf("%.1f%%", 100*mean(pos)),
    event_rate_among_screen_positives = sprintf("%.1f%%", 100*mean(d$incident_2015[pos])),
    sensitivity = sprintf("%.1f%%", 100*sum(d$incident_2015[pos])/sum(d$incident_2015)),
    specificity = sprintf("%.1f%%", 100*mean(d$incident_2015[!pos] == 0)))
}))
print(thr, row.names = FALSE)
write.csv(thr, file.path(DIR, "screening_index_memory_v3_thresholds_2026-09-26.csv"), row.names = FALSE, fileEncoding = "UTF-8")

## secondary: concurrent(2011 cross sectional)
m_c <- glm(impaired_2011 ~ age + edu_group + dep + adl, data = d_conc, family = binomial())
cat("\n== secondary: 2011 concurrent AUC =", sprintf("%.3f", as.numeric(auc(roc(d_conc$impaired_2011,
    predict(m_c, newdata = d_conc, type = "response"), quiet = TRUE)))), "(n =", nrow(d_conc), ")\n")

saveRDS(list(model_u = m_u, model_w = m_w, coef = cf, score_table = score_tab, strata = strata,
             perf = perf, thresholds = thr, cut = floor(cut), data = d, model_concurrent = m_c),
        file.path(DIR, "screening_index_memory_v3_2026-09-26.rds"))
cat("\n saved to", DIR, "\n")
