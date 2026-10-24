# ============================================================================
# 18_screening_index_v3.R —— A07: construct cognitive vulnerability screening index(v3 frozen version)
# basis: SAP_v3_2026-09-26.md §7(primary model 4 items + integer scoring rules + risk stratification)
# variables: age, education(low/in/high), depression(CES-D-10 ≥10), ADL(any limitation), all are 2011 baseline
# outcome: impaired_2020(2011 among unimpaired participants 2020 endpoint assessment); secondary: first impairment Cox
# unweighted primary analysis; weighted(2011 baseline weights + PSU=communityID)for sensitivity analysis
# generate: 2026-09-26
# ============================================================================
invisible(try(Sys.setlocale("LC_ALL", "zh_CN.UTF-8"), silent = TRUE))
.libPaths(c("D:/R/library", .libPaths()))
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(survival); library(survey); library(pROC) })
P <- "D:/BaiduSyncdisk/research/Xu_Zhonghao"; OUT <- file.path(P, "output", "04_longitudinal_cohort")
DIR <- file.path(P, "output", "10_screening_index"); if (!dir.exists(DIR)) dir.create(DIR, recursive = TRUE)

coh <- readRDS(file.path(OUT, "cohort_change_2026-09-13.rds"))$cohort
inc <- readRDS(file.path(OUT, "interval_incidence_events_2026-09-20.rds"))

b11 <- coh %>% filter(year == 2011) %>%
  transmute(KEY, age, female, edu_group, dep = as.integer(cesd10 >= 10), adl,
            ms10_2011 = ms10, weight, communityID)
d <- inc$base %>% inner_join(b11, by = "KEY")
yr <- c("2013" = 2, "2015" = 4, "2018" = 7, "2020" = 9)
im <- cbind(d$im2013, d$im2015, d$im2018, d$im2020)
d$time  <- apply(im, 1, function(z) { i <- which(z == 1); if (length(i)) yr[i[1]] else 9 })
d$event <- as.integer(apply(im, 1, function(z) any(z == 1, na.rm = TRUE)))
d$edu_group <- factor(d$edu_group, levels = c("Less than primary","Primary school","Middle school","High school or above"))

cat("analysis sample:", nrow(d), "｜2020 impairment event:", sum(d$endpoint_2020), "｜first impairment event:", sum(d$event), "\n")
cat("complete cases(4 variables + outcome):", sum(complete.cases(d[, c("age","edu_group","dep","adl","endpoint_2020")])), "\n\n")

## ---------- 1) primary model(unweighted = primary analysis; weighted = sensitivity) ----------
f <- endpoint_2020 ~ age + edu_group + dep + adl
m_u <- glm(f, data = d, family = binomial())
m_w <- svyglm(f, design = svydesign(ids = ~communityID, weights = ~weight, data = d), family = quasibinomial())
cf <- data.frame(
  term = names(coef(m_u)),
  beta_unweighted = as.numeric(coef(m_u)),
  se_unweighted = as.numeric(sqrt(diag(vcov(m_u)))),
  beta_weighted = as.numeric(coef(m_w)[names(coef(m_u))]),
  se_weighted = as.numeric(sqrt(diag(vcov(m_w)))[names(coef(m_u))])
)
cf$OR_unweighted <- exp(cf$beta_unweighted); cf$OR_weighted <- exp(cf$beta_weighted)
cf$ci_u <- sprintf("%.2f (%.2f-%.2f)", cf$OR_unweighted, exp(cf$beta_unweighted - 1.96*cf$se_unweighted), exp(cf$beta_unweighted + 1.96*cf$se_unweighted))
cf$ci_w <- sprintf("%.2f (%.2f-%.2f)", cf$OR_weighted, exp(cf$beta_weighted - 1.96*cf$se_weighted), exp(cf$beta_weighted + 1.96*cf$se_weighted))
print(cf[, c("term","ci_u","ci_w")], row.names = FALSE)
write.csv(cf, file.path(DIR, "screening_index_v3_models_2026-09-26.csv"), row.names = FALSE, fileEncoding = "UTF-8")

## ---------- 2) integerized score(SAP §7.3; based on unweighted primary model) ----------
b <- coef(m_u)
b_age5 <- 5 * b[["age"]]                       # age per 5 years equivalent
terms_beta <- c(age5 = b_age5,
                edu_primary = b[["edu_groupPrimary school"]],
                edu_middle  = b[["edu_groupMiddle school"]],
                edu_high    = b[["edu_groupHigh school or above"]],
                dep = b[["dep"]], adl = b[["adl"]])
names(terms_beta) <- c("age(per 5 years old)","education: primary school(vs no schooling)","education: middle school(vs no schooling)","education: high school or above(vs no schooling)","depression","ADL limited")
bmin <- min(abs(terms_beta))
pts <- round(terms_beta / bmin)
score_tab <- data.frame(project = names(terms_beta), beta_using_5_years_equivalent = round(terms_beta, 4),
                        per_score_unit = as.integer(pts), notes = c(
  "45–49=0 points; per increase 5 years old +1 score units", "relative to Less than primary", "relative to Less than primary", "relative to Less than primary",
  "CES-D-10 ≥10", "any ADL limited"))
print(score_tab, row.names = FALSE)
write.csv(score_tab, file.path(DIR, "screening_index_v3_score_table_2026-09-26.csv"), row.names = FALSE, fileEncoding = "UTF-8")

## calculate individual total score
d$score <- as.integer(pts[["age(per 5 years old)"]]) * pmax(floor((d$age - 45)/5), 0) +
  as.integer(pts[["education: primary school(vs no schooling)"]]) * as.integer(d$edu_group == "Primary school") +
  as.integer(pts[["education: middle school(vs no schooling)"]]) * as.integer(d$edu_group == "Middle school") +
  as.integer(pts[["education: high school or above(vs no schooling)"]]) * as.integer(d$edu_group == "High school or above") +
  as.integer(pts[["depression"]]) * d$dep +
  as.integer(pts[["ADL limited"]]) * d$adl
cat("\n total score range:", range(d$score, na.rm = TRUE), "｜distribution:\n"); print(table(d$score, useNA = "ifany"))

## ---------- 3) risk stratification(tertiles + probability threshold) ----------
q <- quantile(d$score, c(1/3, 2/3), na.rm = TRUE)
d$stratum <- cut(d$score, c(-Inf, q[1], q[2], Inf), labels = c("low", "in", "high"), right = TRUE)
strata <- d %>% group_by(stratum) %>% summarise(n = n(), events = sum(endpoint_2020),
  event_rate = sprintf("%.1f%%", 100*mean(endpoint_2020)), mean_score = round(mean(score), 1), .groups = "drop")
print(strata, row.names = FALSE)

d$p_hat <- predict(m_u, newdata = d, type = "response")
dd <- d[!is.na(d$p_hat), ]
thr <- do.call(rbind, lapply(c(0.10, 0.20, 0.30), function(t) {
  pos <- dd$p_hat >= t
  data.frame(cutoff_threshold = sprintf("predicted probability ≥%d%%", t*100), positive_count = sum(pos), positive_fraction = sprintf("%.1f%%", 100*mean(pos)),
             event_rate_among_screen_positives = sprintf("%.1f%%", 100*mean(dd$endpoint_2020[pos])),
             sensitivity = sprintf("%.1f%%", 100*sum(dd$endpoint_2020[pos])/sum(dd$endpoint_2020)),
             specificity = sprintf("%.1f%%", 100*mean(dd$endpoint_2020[!pos] == 0)))
}))
print(thr, row.names = FALSE)
write.csv(strata, file.path(DIR, "screening_index_v3_risk_strata_2026-09-26.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(thr, file.path(DIR, "screening_index_v3_thresholds_2026-09-26.csv"), row.names = FALSE, fileEncoding = "UTF-8")

## ---------- 4) preliminary performance(formal version in A08/A09) ----------
ok <- complete.cases(d[, c("endpoint_2020","p_hat","score")])
r1 <- roc(d$endpoint_2020[ok], d$p_hat[ok], quiet = TRUE)
r2 <- roc(d$endpoint_2020[ok], d$score[ok], quiet = TRUE)
cal <- glm(endpoint_2020 ~ p_hat, data = d[ok, ], family = binomial())
brier <- mean((d$p_hat[ok] - d$endpoint_2020[ok])^2)
perf <- data.frame(metrics = c("AUC(continuous predicted probability)","AUC(integer total score)","calibration slope(logistic fixed intercept 0? including intercept)","calibration slope β1","Brier"),
                   value = c(sprintf("%.3f (%.3f-%.3f)", as.numeric(auc(r1)), ci.auc(r1)[1], ci.auc(r1)[3]),
                          sprintf("%.3f (%.3f-%.3f)", as.numeric(auc(r2)), ci.auc(r2)[1], ci.auc(r2)[3]),
                          sprintf("%.3f (%.3f-%.3f)", coef(cal)[1], coef(cal)[1]-1.96*summary(cal)$coef[1,2], coef(cal)[1]+1.96*summary(cal)$coef[1,2]),
                          sprintf("%.2f (%.2f-%.2f)", coef(cal)[2], coef(cal)[2]-1.96*summary(cal)$coef[2,2], coef(cal)[2]+1.96*summary(cal)$coef[2,2]),
                          sprintf("%.4f", brier)))
print(perf, row.names = FALSE)
write.csv(perf, file.path(DIR, "screening_index_v3_performance_2026-09-26.csv"), row.names = FALSE, fileEncoding = "UTF-8")

saveRDS(list(model_u = m_u, model_w = m_w, coef_table = cf, score_table = score_tab,
             strata = strata, thresholds = thr, performance = perf, data = d),
        file.path(DIR, "screening_index_v3_2026-09-26.rds"))
cat("\n saved to", DIR, "\n")
