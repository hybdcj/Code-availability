# ============================================================================
# 24_external_validation.R —— A12：外部队列正式外部验证
# 主分析 S0 = baseline_all_mean_sd；正式敏感性 S1 = complete_fu_mean_sd、
# S2 = fixed_charls_le4。
# 队列：HRS（w10→w12）、ELSA（w4→w6）、SHARE（w4→w6）
# 模型：CHARLS 固定模型（年龄 + 教育 + 抑郁 + ADL），不重新拟合。
# 输出：output/11_外部验证/
# 生成：2026-09-26
# ============================================================================
invisible(try(Sys.setlocale("LC_ALL", "zh_CN.UTF-8"), silent = TRUE))
.libPaths(c("D:/R/library", .libPaths()))
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(pROC); library(Hmisc); library(knitr)
})

set.seed(20260926)  # cluster bootstrap 可复现

P    <- "D:/BaiduSyncdisk/科研/徐中豪"
A11  <- file.path(P, "external_data")
OUT  <- file.path(P, "output", "11_外部验证")
FIG  <- file.path(P, "output", "09_图表")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)

a07 <- readRDS(file.path(P, "output/10_筛查指数/screening_index_memory_v3_2026-09-26.rds"))
model_u <- a07$model_u
model_levels <- model_u$xlevels$edu_group

cohorts <- c("HRS","ELSA","SHARE")
cut_defs <- c("baseline_all_mean_sd","complete_fu_mean_sd","fixed_charls_le4")
scenario_lab <- c(
  baseline_all_mean_sd = "S0_baseline_all",
  complete_fu_mean_sd  = "S1_complete_fu",
  fixed_charls_le4     = "S2_fixed_le4"
)
scenario_name <- c(
  baseline_all_mean_sd = "S0 主分析：baseline-all mean−1SD",
  complete_fu_mean_sd  = "S1 敏感性：complete-fu mean−1SD",
  fixed_charls_le4     = "S2 敏感性：固定 CHARLS ≤4"
)
cohort_cols <- c(HRS = "#1f77b4", ELSA = "#2ca02c", SHARE = "#d62728")

# ---------- helpers ----------
wilson_ci <- function(x, n, z = 1.96) {
  if (is.na(n) || n <= 0) return(c(NA_real_, NA_real_))
  p <- x / n
  den <- 1 + z^2 / n
  cen <- (p + z^2 / (2*n)) / den
  half <- z * sqrt(p*(1-p)/n + z^2/(4*n^2)) / den
  c(max(0, cen-half), min(1, cen+half))
}
auc_ci <- function(y, x) {
  ok <- !is.na(y) & !is.na(x)
  y <- y[ok]; x <- x[ok]
  if (length(y) < 10 || length(unique(y)) < 2) return(c(auc=NA, low=NA, high=NA))
  r <- pROC::roc(y, x, direction = "<", quiet = TRUE)
  ci <- as.numeric(pROC::ci.auc(r))
  c(auc = as.numeric(pROC::auc(r)), low = ci[1], high = ci[3])
}
calibration_metrics <- function(y, p) {
  ok <- !is.na(y) & !is.na(p) & p > 0 & p < 1
  y <- y[ok]; p <- p[ok]
  lp <- qlogis(p)
  slope_fit <- glm(y ~ lp, family = binomial())
  slope <- coef(slope_fit)[["lp"]]
  slope_ci <- confint.default(slope_fit)["lp", ]
  int_fit <- glm(y ~ offset(lp), family = binomial())
  int <- coef(int_fit)[["(Intercept)"]]
  int_ci <- confint.default(int_fit)["(Intercept)", ]
  brier <- mean((p - y)^2)
  p0 <- mean(y)
  scaled <- if (p0 > 0 && p0 < 1) 1 - brier/(p0*(1-p0)) else NA_real_
  c(slope = unname(slope), slope_low = unname(slope_ci[1]), slope_high = unname(slope_ci[2]),
    intercept = unname(int), intercept_low = unname(int_ci[1]), intercept_high = unname(int_ci[2]),
    brier = unname(brier), scaled_brier = unname(scaled))
}
add_pred <- function(d, model) {
  d <- as.data.frame(d)
  d$edu_group <- factor(as.character(d$edu_group), levels = model_levels)
  nd <- data.frame(age = d$age, edu_group = d$edu_group, dep = d$dep, adl = d$adl)
  d$p_hat <- as.numeric(predict(model, newdata = nd, type = "response"))
  d$lp <- qlogis(pmin(pmax(d$p_hat, 1e-6), 1 - 1e-6))
  d
}
evaluate_one <- function(d, cohort, cd, model) {
  d <- add_pred(d, model)
  m <- d[!is.na(d$p_hat) & !is.na(d$event) & !is.na(d$score), , drop = FALSE]
  a_model <- auc_ci(m$event, m$p_hat)
  a_score_full <- auc_ci(d$event, d$score)
  a_score_same <- auc_ci(m$event, m$score)
  cal <- calibration_metrics(m$event, m$p_hat)
  strata_full <- d %>% mutate(stratum = cut(score, c(-Inf,-3,4,Inf), labels = c("低","中","高"))) %>%
    group_by(stratum) %>% summarise(n = n(), events = sum(event), rate = mean(event), .groups = "drop")
  strata_model <- m %>% mutate(stratum = cut(score, c(-Inf,-3,4,Inf), labels = c("低","中","高"))) %>%
    group_by(stratum) %>% summarise(n_model = n(), events_model = sum(event),
                                    pred_mean = mean(p_hat), .groups = "drop")
  strata <- full_join(strata_full, strata_model, by = "stratum") %>%
    mutate(oe_model = events_model/(n_model*pred_mean), cohort = cohort, cut_definition = cd) %>%
    relocate(cohort, cut_definition)
  lo <- strata_full$rate[strata_full$stratum == "低"]; hi <- strata_full$rate[strata_full$stratum == "高"]
  high_low <- if (length(lo) == 1 && length(hi) == 1 && !is.na(lo) && lo > 0) hi/lo else NA_real_
  thr <- do.call(rbind, lapply(c(0.10,0.15,0.20), function(t) {
    pos <- m$p_hat >= t
    TP <- sum(pos & m$event == 1); FP <- sum(pos & m$event == 0)
    FN <- sum(!pos & m$event == 1); TN <- sum(!pos & m$event == 0)
    sens <- if ((TP+FN) > 0) TP/(TP+FN) else NA_real_
    spec <- if ((TN+FP) > 0) TN/(TN+FP) else NA_real_
    ppv  <- if ((TP+FP) > 0) TP/(TP+FP) else NA_real_
    npv  <- if ((TN+FN) > 0) TN/(TN+FN) else NA_real_
    prev <- mean(m$event)
    nb_m <- TP/nrow(m) - FP/nrow(m) * t/(1-t)
    nb_a <- prev - (1-prev) * t/(1-t)
    data.frame(cohort = cohort, cut_definition = cd, threshold = t, n_positive = sum(pos),
               pos_rate = mean(pos), sensitivity = sens, specificity = spec,
               ppv = ppv, npv = npv, nb_model = nb_m, nb_all = nb_a)
  }))
  m$decile <- dplyr::ntile(m$p_hat, 10)
  dec <- m %>% group_by(decile) %>% summarise(n = n(), pred_mean = mean(p_hat),
                                              events = sum(event), observed = mean(event), .groups = "drop")
  dec$observed_low  <- mapply(function(x,n) wilson_ci(x,n)[1], dec$events, dec$n)
  dec$observed_high <- mapply(function(x,n) wilson_ci(x,n)[2], dec$events, dec$n)
  dec <- dec %>% mutate(cohort = cohort, cut_definition = cd) %>% relocate(cohort, cut_definition)
  dca <- bind_rows(lapply(seq(0.05, 0.40, 0.01), function(t) {
    pos <- m$p_hat >= t
    TP <- sum(pos & m$event == 1); FP <- sum(pos & m$event == 0)
    nb_m <- TP/nrow(m) - FP/nrow(m) * t/(1-t)
    prev <- mean(m$event)
    nb_a <- prev - (1-prev) * t/(1-t)
    data.frame(cohort = cohort, cut_definition = cd, threshold = t,
               strategy = c("model","all","none"), net_benefit = c(nb_m, nb_a, 0))
  }))
  perf <- data.frame(
    cohort = cohort, cut_definition = cd, scenario = scenario_lab[[cd]],
    cut_int = unique(d$cut_int), n_total = nrow(d), events_total = sum(d$event),
    event_rate_total = mean(d$event), n_model = nrow(m), events_model = sum(m$event),
    event_rate_model = mean(m$event),
    auc_model = a_model["auc"], auc_model_low = a_model["low"], auc_model_high = a_model["high"],
    auc_score_full = a_score_full["auc"], auc_score_full_low = a_score_full["low"], auc_score_full_high = a_score_full["high"],
    auc_score_same = a_score_same["auc"], auc_score_same_low = a_score_same["low"], auc_score_same_high = a_score_same["high"],
    cal_slope = cal["slope"], cal_slope_low = cal["slope_low"], cal_slope_high = cal["slope_high"],
    cal_intercept = cal["intercept"], cal_intercept_low = cal["intercept_low"], cal_intercept_high = cal["intercept_high"],
    brier = cal["brier"], scaled_brier = cal["scaled_brier"], high_low_ratio = high_low
  )
  list(perf = perf, strata = strata, deciles = dec, thresholds = thr, dca = dca, data = m)
}

# ---------- 读取三套风险集 ----------
risk_sets <- list()
for (co in cohorts) {
  f <- file.path(A11, co, "analysis", paste0(co, "_sensitivity_risk_sets_2026-09-26.rds"))
  x <- readRDS(f)
  for (cd in cut_defs) risk_sets[[paste(co, cd, sep = "__")]] <- x[[paste(co, cd, sep = "__")]]
}

# ---------- 主分析：三队列 × S0/S1/S2 ----------
all_perf <- list(); all_strata <- list(); all_dec <- list(); all_thr <- list(); all_dca <- list(); all_data <- list()
for (co in cohorts) {
  for (cd in cut_defs) {
    res <- evaluate_one(risk_sets[[paste(co, cd, sep = "__")]], co, cd, model_u)
    key <- paste(co, cd, sep = "__")
    all_perf[[key]] <- res$perf; all_strata[[key]] <- res$strata; all_dec[[key]] <- res$deciles
    all_thr[[key]] <- res$thresholds; all_dca[[key]] <- res$dca; all_data[[key]] <- res$data
  }
}
perf <- bind_rows(all_perf); strata <- bind_rows(all_strata); deciles <- bind_rows(all_dec)
thresholds <- bind_rows(all_thr); dca <- bind_rows(all_dca)

# ---------- SHARE 国家层面 + cluster bootstrap ----------
share_s0 <- add_pred(risk_sets[["SHARE__baseline_all_mean_sd"]], model_u)
share_ok <- share_s0[!is.na(share_s0$p_hat), ]
share_country <- bind_rows(lapply(sort(unique(share_ok$country_mod)), function(cc) {
  d <- share_ok[share_ok$country_mod == cc, , drop = FALSE]
  a1 <- auc_ci(d$event, d$p_hat); a2 <- auc_ci(d$event, d$score)
  cal <- tryCatch(calibration_metrics(d$event, d$p_hat), error = function(e) rep(NA_real_, 8))
  data.frame(country_mod = cc, country_label = unique(d$country_label), n = nrow(d),
             events = sum(d$event), event_rate = mean(d$event),
             auc_model = a1["auc"], auc_model_low = a1["low"], auc_model_high = a1["high"],
             auc_score = a2["auc"], auc_score_low = a2["low"], auc_score_high = a2["high"],
             cal_slope = cal["slope"], cal_intercept = cal["intercept"], brier = cal["brier"])
}))
cluster_boot_share_metric <- function(d, metric = c("auc_model","auc_score"), B = 1000) {
  metric <- match.arg(metric)
  dd <- if (metric == "auc_model") d[!is.na(d$p_hat), , drop = FALSE] else d
  countries <- unique(dd$country_mod)
  out <- numeric(B)
  for (b in seq_len(B)) {
    samp <- sample(countries, length(countries), replace = TRUE)
    db <- bind_rows(lapply(samp, function(cc) dd[dd$country_mod == cc, , drop = FALSE]))
    out[b] <- if (metric == "auc_model") auc_ci(db$event, db$p_hat)["auc"] else auc_ci(db$event, db$score)["auc"]
  }
  orig <- if (metric == "auc_model") auc_ci(dd$event, dd$p_hat)["auc"] else auc_ci(dd$event, dd$score)["auc"]
  data.frame(metric = metric, orig_auc = as.numeric(orig),
             boot_mean = mean(out, na.rm = TRUE),
             boot_low = as.numeric(quantile(out, 0.025, na.rm = TRUE)),
             boot_high = as.numeric(quantile(out, 0.975, na.rm = TRUE)), B = B)
}
share_boot_s0 <- bind_rows(
  cluster_boot_share_metric(share_s0, "auc_model", B = 1000),
  cluster_boot_share_metric(share_s0, "auc_score", B = 1000)
); share_boot_s0$cut_definition <- "baseline_all_mean_sd"
share_s1 <- add_pred(risk_sets[["SHARE__complete_fu_mean_sd"]], model_u)
share_boot_s1 <- share_boot_s0; share_boot_s1$cut_definition <- "complete_fu_mean_sd"  # SHARE S1 与 S0 风险集/事件相同
share_s2 <- add_pred(risk_sets[["SHARE__fixed_charls_le4"]], model_u)
share_boot_s2 <- bind_rows(
  cluster_boot_share_metric(share_s2, "auc_model", B = 1000),
  cluster_boot_share_metric(share_s2, "auc_score", B = 1000)
); share_boot_s2$cut_definition <- "fixed_charls_le4"
share_boot <- bind_rows(share_boot_s0, share_boot_s1, share_boot_s2) %>%
  relocate(cut_definition, metric, orig_auc, boot_mean, boot_low, boot_high, B)

# ---------- 跨队列随机效应 meta（DL，logit AUC） ----------
meta_dl <- function(auc, low, high, labels) {
  yi <- qlogis(auc); sei <- (qlogis(high) - qlogis(low)) / (2*1.96)
  wi <- 1/sei^2
  fixed <- sum(wi*yi)/sum(wi)
  Q <- sum(wi*(yi-fixed)^2); df <- length(yi)-1
  C <- sum(wi) - sum(wi^2)/sum(wi)
  tau2 <- max(0, (Q-df)/C)
  ws <- 1/(sei^2 + tau2)
  pooled <- sum(ws*yi)/sum(ws); se <- sqrt(1/sum(ws))
  I2 <- if (Q > 0) max(0, (Q-df)/Q)*100 else 0
  data.frame(k = length(yi), pooled_auc = plogis(pooled),
             pooled_low = plogis(pooled - 1.96*se), pooled_high = plogis(pooled + 1.96*se),
             tau2 = tau2, Q = Q, p_Q = pchisq(Q, df, lower.tail = FALSE), I2 = I2,
             cohorts = paste(labels, collapse = ", "))
}
meta_rows <- list()
for (cd in cut_defs) {
  for (metric in c("model","score_full")) {
    sub <- perf[perf$cut_definition == cd, ]
    if (metric == "model") {
      a <- sub$auc_model; lo <- sub$auc_model_low; hi <- sub$auc_model_high
    } else {
      a <- sub$auc_score_full; lo <- sub$auc_score_full_low; hi <- sub$auc_score_full_high
    }
    # SHARE 用与当前切点定义匹配的 cluster bootstrap CI
    sb_metric <- ifelse(metric == "model", "auc_model", "auc_score")
    sb <- share_boot[share_boot$cut_definition == cd & share_boot$metric == sb_metric, ]
    if (nrow(sb) != 1) stop("SHARE cluster bootstrap CI not found for ", cd, "/", sb_metric)
    if (any(sub$cohort == "SHARE")) {
      lo[sub$cohort == "SHARE"] <- sb$boot_low
      hi[sub$cohort == "SHARE"] <- sb$boot_high
    }
    m <- meta_dl(a, lo, hi, sub$cohort)
    m$cut_definition <- cd; m$scenario <- scenario_lab[[cd]]; m$metric <- metric
    meta_rows[[paste(cd, metric)]] <- m
  }
}
meta <- bind_rows(meta_rows)

# ---------- 额外敏感性 ----------
# 1) 抑郁切点：HRS/ELSA CES-D-8 ≥3；SHARE EURO-D ≥5
extra_dep <- bind_rows(lapply(cohorts, function(co) {
  d <- add_pred(risk_sets[[paste(co, "baseline_all_mean_sd", sep = "__")]], model_u)
  d$dep_alt <- as.integer(d$dep_raw >= ifelse(co == "SHARE", 5, 3))
  d$score_dep_alt <- 2*pmax(floor((d$age-45)/5), 0) + d$edu_pts + 3*d$dep_alt
  d$p_alt <- as.numeric(predict(model_u, newdata = data.frame(age = d$age,
    edu_group = factor(as.character(d$edu_group), levels = model_levels), dep = d$dep_alt, adl = d$adl), type = "response"))
  m <- d[!is.na(d$p_alt) & !is.na(d$event), ]
  a1 <- auc_ci(m$event, m$p_alt); a2 <- auc_ci(m$event, m$score_dep_alt)
  cal <- calibration_metrics(m$event, m$p_alt)
  data.frame(cohort = co, sensitivity = "depression_alt_cutoff",
             dep_rule = ifelse(co == "SHARE", "EURO-D ≥5", "CES-D-8 ≥3"),
             n = nrow(m), events = sum(m$event),
             auc_model = a1["auc"], auc_model_low = a1["low"], auc_model_high = a1["high"],
             auc_score = a2["auc"], auc_score_low = a2["low"], auc_score_high = a2["high"],
             cal_slope = cal["slope"], cal_intercept = cal["intercept"], brier = cal["brier"])
}))
# 2) score_alt
extra_score_alt <- bind_rows(lapply(cohorts, function(co) {
  d <- risk_sets[[paste(co, "baseline_all_mean_sd", sep = "__")]]
  a <- auc_ci(d$event, d$score_alt)
  data.frame(cohort = co, sensitivity = "score_alt_positional_education",
             n = nrow(d), events = sum(d$event), auc_score_alt = a["auc"],
             auc_low = a["low"], auc_high = a["high"])
}))
# 3) ELSA age 90 顶码：排除 age=90
extra_elsa_age <- local({
  d <- add_pred(risk_sets[["ELSA__baseline_all_mean_sd"]], model_u)
  d2 <- d[d$age < 90, , drop = FALSE]
  m <- d2[!is.na(d2$p_hat) & !is.na(d2$event), ]
  a1 <- auc_ci(m$event, m$p_hat); a2 <- auc_ci(m$event, m$score)
  cal <- calibration_metrics(m$event, m$p_hat)
  data.frame(cohort = "ELSA", sensitivity = "exclude_age90_topcode",
             n = nrow(m), events = sum(m$event),
             auc_model = a1["auc"], auc_model_low = a1["low"], auc_model_high = a1["high"],
             auc_score = a2["auc"], auc_score_low = a2["low"], auc_score_high = a2["high"],
             cal_slope = cal["slope"], cal_intercept = cal["intercept"], brier = cal["brier"])
})

# ---------- 保存 CSV / RDS ----------
write.csv(perf, file.path(OUT, "A12_主结果_三队列_S0_S1_S2.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(strata, file.path(OUT, "A12_风险分层_S0_S1_S2.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(deciles, file.path(OUT, "A12_校准十分位_S0_S1_S2.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(thresholds, file.path(OUT, "A12_阈值表现_S0_S1_S2.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(dca, file.path(OUT, "A12_DCA_S0_S1_S2.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(share_country, file.path(OUT, "A12_SHARE国家层面.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(share_boot, file.path(OUT, "A12_SHARE国家聚类bootstrap.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(meta, file.path(OUT, "A12_跨队列随机效应.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(extra_dep, file.path(OUT, "A12_敏感性_抑郁切点.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(extra_score_alt, file.path(OUT, "A12_敏感性_score_alt.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(extra_elsa_age, file.path(OUT, "A12_敏感性_ELSA_age90.csv"), row.names = FALSE, fileEncoding = "UTF-8")
saveRDS(list(perf = perf, strata = strata, deciles = deciles, thresholds = thresholds, dca = dca,
             share_country = share_country, share_boot = share_boot, meta = meta,
             extra_dep = extra_dep, extra_score_alt = extra_score_alt, extra_elsa_age = extra_elsa_age),
        file.path(OUT, "A12_results_2026-09-26.rds"))

writeLines(capture.output(sessionInfo()), file.path(OUT, "A12_sessionInfo_2026-09-26.txt"), useBytes = TRUE)

cat("\nA12 主结果（S0）：\n")
print(perf %>% filter(cut_definition == "baseline_all_mean_sd") %>%
        select(cohort, n_model, events_model, auc_model, auc_model_low, auc_model_high,
               auc_score_full, cal_slope, cal_intercept, brier) %>%
        mutate(across(where(is.numeric), ~round(.x, 4))) %>% as.data.frame())
cat("\nA12 敏感性（S1/S2）AUC：\n")
print(perf %>% select(cohort, scenario, auc_model, auc_score_full) %>%
        mutate(across(where(is.numeric), ~round(.x, 4))) %>% as.data.frame())
cat("\nA12 跨队列随机效应：\n")
print(meta %>% select(cut_definition, metric, pooled_auc, pooled_low, pooled_high, I2, p_Q) %>%
        mutate(across(where(is.numeric), ~round(.x, 4))) %>% as.data.frame())
cat("\n已写入：", OUT, "\n")