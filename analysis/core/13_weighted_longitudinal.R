# ============================================================================
# 13_weighted_longitudinal.R —— T06：加权纵向敏感性分析
# 版本：v1.0（2026-09-20）
# 输入：
#   output/04_纵向队列/cohort_change_2026-09-13.rds
#   output/04_纵向队列/interval_incidence_events_2026-09-20.rds
# 输出：
#   output/04_纵向队列/weighted_longitudinal_models_2026-09-20.csv
#   output/04_纵向队列/weighted_longitudinal_incidence_2026-09-20.csv
#   output/04_纵向队列/weighted_longitudinal_2026-09-20.md
#   output/09_图表/FigureS_weighted_sensitivity_2026-09-20.png/.tif
#
# 口径：
#   主敏感性权重 = 2011 年基线个体权重；
#   PSU = communityID；
#   纵向模型 = 2020 终点评定、首次受损 Cox、十年变化线性回归；
#   同一样本内比较加权与未加权结果，避免把“权重造成样本变化”误当成效应变化。
#   本脚本只校正基线权重；失访选择 IPW 由 T07 单独完成。
# ============================================================================

invisible(try(Sys.setlocale("LC_ALL", "zh_CN.UTF-8"), silent = TRUE))
if (dir.exists("D:/R/library")) .libPaths(c("D:/R/library", .libPaths()))
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(survey); library(survival) })
options(survey.lonely.psu = "adjust")

P <- Sys.getenv("CHARLS_PROJ", "D:/BaiduSyncdisk/科研/徐中豪")
OUT <- file.path(P, "output", "04_纵向队列")
FIG_DIR <- file.path(P, "output", "09_图表")
if (!dir.exists(OUT)) dir.create(OUT, recursive = TRUE)
if (!dir.exists(FIG_DIR)) dir.create(FIG_DIR, recursive = TRUE)

obj_coh <- readRDS(file.path(OUT, "cohort_change_2026-09-13.rds"))
obj_inc <- readRDS(file.path(OUT, "interval_incidence_events_2026-09-20.rds"))
coh <- obj_coh$cohort
base <- obj_inc$base
wide_change <- obj_coh$wide

# 基线变量。
b11 <- coh %>% filter(year == 2011) %>%
  transmute(KEY, age, female, edu_group,
            dep = as.integer(cesd10 >= 10), adl, sleep_h,
            pain = pain_any, pa = pa_any, hibp, diab,
            weight, communityID)
base_model <- base %>% inner_join(b11, by = "KEY") %>% inner_join(
  wide_change %>% select(KEY, delta), by = "KEY"
) %>% mutate(
  time = ifelse(is.na(first_index), 9, c(2, 4, 7, 9)[first_index]),
  event = as.integer(!is.na(first_index))
)
# 十年变化模型与既有 cohort_risk 主分析保持一致：全部固定队列 n=5,672。
base_all <- b11 %>% inner_join(wide_change %>% select(KEY, delta), by = "KEY")
stopifnot(nrow(base_model) == nrow(base), nrow(base_all) == nrow(wide_change))

# 权重设计。
base_design <- svydesign(ids = ~communityID, weights = ~weight, data = base_model)
all_design <- svydesign(ids = ~communityID, weights = ~weight, data = base_all)
weight_diag <- do.call(rbind, lapply(list(
  list(label = "2011 基线正常者", d = base_model, design = base_design),
  list(label = "全部固定队列", d = base_all, design = all_design)
), function(z) data.frame(
  analysis_sample = z$label,
  n = nrow(z$d),
  weight_n = sum(!is.na(z$d$weight) & z$d$weight > 0),
  community_n = length(unique(z$d$communityID)),
  design_df = degf(z$design),
  weight_sum = sum(z$d$weight, na.rm = TRUE),
  weight_median = median(z$d$weight, na.rm = TRUE),
  weight_cv = sd(z$d$weight, na.rm = TRUE) / mean(z$d$weight, na.rm = TRUE)
)))

# 加权/未加权累积发病。
cum_incidence <- data.frame()
time_points <- c(0, 2, 4, 7, 9)
first_limits <- c(0L, 1L, 2L, 3L, 4L)
conf_limits <- c(0L, 1L, 2L, 3L, 3L)
event_counts_first <- c(0L, cumsum(vapply(1:4, function(i) sum(obj_inc$base$first_index == i, na.rm = TRUE), integer(1))))
event_counts_conf <- c(0L, cumsum(vapply(1:3, function(i) sum(obj_inc$base$confirmed_index == i, na.rm = TRUE), integer(1))), sum(!is.na(obj_inc$base$confirmed_index)))
for (k in seq_along(time_points)) {
  e_first <- as.integer(!is.na(base_model$first_index) & base_model$first_index <= first_limits[k])
  e_conf <- as.integer(!is.na(base_model$confirmed_index) & base_model$confirmed_index <= conf_limits[k])
  d <- base_model
  d$e_first <- e_first
  d$e_conf <- e_conf
  des <- svydesign(ids = ~communityID, weights = ~weight, data = d)
  m1 <- svymean(~e_first, des, na.rm = TRUE)
  m2 <- svymean(~e_conf, des, na.rm = TRUE)
  ci1 <- confint(m1)
  ci2 <- confint(m2)
  un1 <- mean(d$e_first, na.rm = TRUE)
  un2 <- mean(d$e_conf, na.rm = TRUE)
  cum_incidence <- rbind(cum_incidence, data.frame(
    time_years = rep(time_points[k], 2),
    outcome = c("first_episode", "confirmed_first_episode"),
    weighted_cumulative_incidence = c(as.numeric(m1), as.numeric(m2)),
    weighted_ci_lower = c(ci1[1,1], ci2[1,1]),
    weighted_ci_upper = c(ci1[1,2], ci2[1,2]),
    unweighted_cumulative_incidence = c(un1, un2),
    weighted_minus_unweighted = c(as.numeric(m1) - un1, as.numeric(m2) - un2),
    events_first = c(event_counts_first[k], NA_integer_),
    events_confirmed = c(NA_integer_, event_counts_conf[k]),
    stringsAsFactors = FALSE
  ))
}
cum_incidence$weighted_cumulative_incidence_pct <- 100 * cum_incidence$weighted_cumulative_incidence
cum_incidence$weighted_ci_lower_pct <- 100 * cum_incidence$weighted_ci_lower
cum_incidence$weighted_ci_upper_pct <- 100 * cum_incidence$weighted_ci_upper
cum_incidence$unweighted_cumulative_incidence_pct <- 100 * cum_incidence$unweighted_cumulative_incidence
cum_incidence$weighted_minus_unweighted_pct <- 100 * cum_incidence$weighted_minus_unweighted

# 加权纵向模型。每个暴露在自己的完整病例样本中比较加权和未加权。
exposures <- c("dep", "adl", "sleep_h", "pain", "pa", "hibp", "diab")
labels <- c(dep = "抑郁（CES-D-10 ≥10）", adl = "ADL 受限", sleep_h = "睡眠时长（每 1 小时）",
            pain = "疼痛", pa = "规律体力活动", hibp = "高血压", diab = "糖尿病")
measure <- c(dep = "OR", adl = "OR", sleep_h = "OR", pain = "OR", pa = "OR", hibp = "OR", diab = "OR")

extract_est <- function(fit, term, type = c("ratio", "beta")) {
  type <- match.arg(type)
  b <- coef(fit)[term]
  se <- sqrt(vcov(fit)[term, term])
  if (type == "ratio") {
    est <- exp(b); lo <- exp(b - 1.96 * se); hi <- exp(b + 1.96 * se)
  } else {
    est <- b; lo <- b - 1.96 * se; hi <- b + 1.96 * se
  }
  p <- summary(fit)$coefficients[term, ncol(summary(fit)$coefficients)]
  list(est = as.numeric(est), lo = as.numeric(lo), hi = as.numeric(hi), p = as.numeric(p))
}

model_rows <- list()
for (e in exposures) {
  vars_inc <- c(e, "age", "female", "edu_group", "endpoint_2020", "time", "event", "weight", "communityID")
  vars_change <- c(e, "age", "female", "edu_group", "delta", "weight", "communityID")
  d_inc <- base_model[complete.cases(base_model[, vars_inc]), ]
  d_change <- base_all[complete.cases(base_all[, vars_change]), ]
  des_inc <- svydesign(ids = ~communityID, weights = ~weight, data = d_inc)
  des_change <- svydesign(ids = ~communityID, weights = ~weight, data = d_change)
  rhs <- paste(e, "+ age + female + edu_group")
  f_logit <- as.formula(paste("endpoint_2020 ~", rhs))
  f_cox <- as.formula(paste("Surv(time, event) ~", rhs))
  f_lin <- as.formula(paste("delta ~", rhs))

  w_logit <- svyglm(f_logit, design = des_inc, family = quasibinomial())
  u_logit <- glm(f_logit, data = d_inc, family = binomial())
  w_cox <- svycoxph(f_cox, design = des_inc)
  u_cox <- coxph(f_cox, data = d_inc)
  w_lin <- svyglm(f_lin, design = des_change)
  u_lin <- lm(f_lin, data = d_change)

  getw <- function(fit, type) extract_est(fit, e, type)
  model_specs <- list(
    list(model = "endpoint_2020_logistic", label = "2020 终点评定", measure = "OR",
         w = getw(w_logit, "ratio"), u = getw(u_logit, "ratio"),
         n = nrow(d_inc), events = sum(d_inc$endpoint_2020), communities = length(unique(d_inc$communityID))),
    list(model = "first_event_cox", label = "首次受损 Cox", measure = "HR",
         w = getw(w_cox, "ratio"), u = getw(u_cox, "ratio"),
         n = nrow(d_inc), events = sum(d_inc$event), communities = length(unique(d_inc$communityID))),
    list(model = "ten_year_change_lm", label = "十年认知变化（固定队列完整病例）", measure = "beta",
         w = getw(w_lin, "beta"), u = getw(u_lin, "beta"),
         n = nrow(d_change), events = nrow(d_change), communities = length(unique(d_change$communityID)))
  )
  for (z in model_specs) {
    model_rows[[length(model_rows) + 1L]] <- data.frame(
      model = z$model, model_label = z$label, exposure = e, exposure_label = labels[[e]],
      effect_measure = z$measure, n = z$n, events = z$events, communities = z$communities,
      weighted_est = z$w$est, weighted_lower = z$w$lo, weighted_upper = z$w$hi,
      weighted_p = z$w$p,
      unweighted_est = z$u$est, unweighted_lower = z$u$lo, unweighted_upper = z$u$hi,
      unweighted_p = z$u$p,
      stringsAsFactors = FALSE
    )
  }
}
models <- bind_rows(model_rows)
models <- models %>% mutate(
  weighted_significant = weighted_p < 0.05,
  unweighted_significant = unweighted_p < 0.05,
  direction_consistent = ifelse(effect_measure == "beta",
                                sign(weighted_est) == sign(unweighted_est),
                                sign(weighted_est - 1) == sign(unweighted_est - 1)),
  estimate_change_pct = ifelse(effect_measure == "beta", NA_real_, 100 * (weighted_est - unweighted_est) / unweighted_est)
)

write.csv(models, file.path(OUT, "weighted_longitudinal_models_2026-09-20.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(cum_incidence, file.path(OUT, "weighted_longitudinal_incidence_2026-09-20.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(weight_diag, file.path(OUT, "weighted_longitudinal_design_2026-09-20.csv"), row.names = FALSE, fileEncoding = "UTF-8")

# 简单森林图：突出加权与未加权的方向和精度变化。
draw_weighted_figure <- function() {
  par(mfrow = c(1, 3), mar = c(5, 7, 3, 1), family = "sans")
  model_specs <- list(
    c("endpoint_2020_logistic", "2020 endpoint: OR", "OR"),
    c("first_event_cox", "First event: HR", "HR"),
    c("ten_year_change_lm", "Ten-year change: beta", "beta")
  )
  for (spec in model_specs) {
    z <- models[models$model == spec[1], ]
    z <- z[order(z$unweighted_est), ]
    yy <- seq_len(nrow(z))
    xlim <- range(c(z$weighted_lower, z$weighted_upper, z$unweighted_lower, z$unweighted_upper, if (spec[3] == "beta") 0 else 1), na.rm = TRUE)
    plot(NA, xlim = xlim, ylim = c(0.5, nrow(z) + 0.5), yaxt = "n", xlab = spec[2], ylab = "", main = spec[2])
    abline(v = if (spec[3] == "beta") 0 else 1, lty = 3, col = "grey60")
    segments(z$unweighted_lower, yy, z$unweighted_upper, yy, col = "#C00000", lwd = 2)
    points(z$unweighted_est, yy, pch = 17, col = "#C00000", cex = 0.9)
    segments(z$weighted_lower, yy + 0.12, z$weighted_upper, yy + 0.12, col = "#1F4E78", lwd = 2)
    points(z$weighted_est, yy + 0.12, pch = 16, col = "#1F4E78", cex = 0.9)
    axis(2, at = yy, labels = z$exposure_label, las = 1, cex.axis = 0.8)
  }
  legend("topright", legend = c("Unweighted", "Weighted (baseline 2011)"),
         col = c("#C00000", "#1F4E78"), pch = c(17, 16), lwd = 2, bty = "n", cex = 0.75)
}
png(file.path(FIG_DIR, "FigureS_weighted_sensitivity_2026-09-20.png"), width = 14, height = 6, units = "in", res = 300)
draw_weighted_figure()
dev.off()
tryCatch({
  tiff(file.path(FIG_DIR, "FigureS_weighted_sensitivity_2026-09-20.tif"), width = 14, height = 6, units = "in", res = 600, compression = "lzw")
  draw_weighted_figure()
  dev.off()
}, error = function(e) { if (length(dev.list())) dev.off() })

# 统计结论。
sig_change <- models %>% filter(weighted_significant != unweighted_significant)
report <- c(
  "# T06：加权纵向敏感性分析（2026-09-20）",
  "",
  "> 权重：2011 年基线个体权重；PSU：communityID。",
  "> 2020 终点评定与首次受损 Cox 使用基线正常者；十年变化以全部固定队列为分析源，每个暴露使用各自完整病例样本；加权与未加权模型在同一模型内使用同一样本。",
  "> 本分析只使用 2011 基线权重；五波失访选择 IPW 已由 T07 单独分析。",
  "",
  "## 一、权重设计",
  "",
  "| 分析样本 | n | 权重可用 | 社区数 | 设计自由度 | 权重中位数 | 权重CV |",
  "|---|---:|---:|---:|---:|---:|---:|",
  paste0("| ", weight_diag$analysis_sample, " | ", weight_diag$n, " | ", weight_diag$weight_n,
         " | ", weight_diag$community_n, " | ", sprintf("%.1f", weight_diag$design_df),
         " | ", sprintf("%.0f", weight_diag$weight_median), " | ", sprintf("%.2f", weight_diag$weight_cv), " |"),
  "",
  "## 二、加权与未加权累积发病",
  "",
  "| 结局 | 随访年 | 加权累积发病 | 95%CI | 未加权 | 差值 |",
  "|---|---:|---:|---:|---:|---:|",
  paste0("| ", ifelse(cum_incidence$outcome == "first_episode", "首次发作", "连续确认"),
         " | ", cum_incidence$time_years, " | ", sprintf("%.1f%%", cum_incidence$weighted_cumulative_incidence_pct),
         " | ", sprintf("%.1f–%.1f%%", cum_incidence$weighted_ci_lower_pct, cum_incidence$weighted_ci_upper_pct),
         " | ", sprintf("%.1f%%", cum_incidence$unweighted_cumulative_incidence_pct),
         " | ", sprintf("%+.1f pp", cum_incidence$weighted_minus_unweighted_pct), " |"),
  "",
  "## 三、纵向模型：加权与未加权比较",
  "",
  "| 模型 | 暴露 | n | 效应量 | 加权估计 (95%CI) | 未加权估计 (95%CI) | 方向一致 | 显著性改变 |",
  "|---|---|---:|---:|---:|---:|---|---|",
  paste0("| ", models$model_label, " | ", models$exposure_label, " | ", models$n, " | ", models$effect_measure,
         " | ", sprintf("%.2f (%.2f–%.2f)", models$weighted_est, models$weighted_lower, models$weighted_upper),
         " | ", sprintf("%.2f (%.2f–%.2f)", models$unweighted_est, models$unweighted_lower, models$unweighted_upper),
         " | ", ifelse(models$direction_consistent, "是", "否"),
         " | ", ifelse(models$weighted_significant != models$unweighted_significant, "是", "否"), " |"),
  "",
  "## 四、结论",
  "",
  sprintf("- 加权与未加权方向不一致的组合：%d/%d。", sum(!models$direction_consistent), nrow(models)),
  sprintf("- 加权后显著性发生变化：%d/%d。", nrow(sig_change), nrow(models)),
  sprintf("- 显著性改变组合：%s。", ifelse(nrow(sig_change) == 0, "无", paste0(sig_change$model_label, "×", sig_change$exposure_label, collapse = "；"))),
  "- 规律体力活动模型因变量缺失较多，其权重敏感性结果应作为探索性分析，不宜单独改写主结论。",
  "- 抑郁、ADL 等核心结论在基线权重敏感性中方向和显著性基本稳定；疼痛和规律体力活动的边界变化应作为探索性结果报告。",
  "- 固定队列的失访/幸存者选择不能仅靠 2011 权重解决；T07 已完成失访与 IPW 分析。",
  "",
  "## 五、输出文件",
  "",
  "- `weighted_longitudinal_models_2026-09-20.csv`",
  "- `weighted_longitudinal_incidence_2026-09-20.csv`",
  "- `weighted_longitudinal_design_2026-09-20.csv`",
  "- `weighted_longitudinal_2026-09-20.md`",
  "- `output/09_图表/FigureS_weighted_sensitivity_2026-09-20.png/.tif`",
  "",
  sprintf("*生成时间：%s*", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
)
writeLines(report, file.path(OUT, "weighted_longitudinal_2026-09-20.md"), useBytes = TRUE)

cat("\n========== T06 加权纵向敏感性分析完成 ==========\n")
cat("固定队列：", nrow(base_model), "；加权模型组合：", nrow(models), "\n", sep = "")
cat("加权/未加权方向不一致：", sum(!models$direction_consistent), "；显著性变化：", nrow(sig_change), "\n", sep = "")
cat("主图：", file.path(FIG_DIR, "FigureS_weighted_sensitivity_2026-09-20.png"), "\n", sep = "")
