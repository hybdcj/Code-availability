# ============================================================================
# 14_attrition_selection.R —— T07：失访、固定队列选择偏倚与 IPW 敏感性分析
# 版本：v1.0（2026-09-20）
# 输入：output/02_五波数据与质控/waves_cognition_2026-09-13.rds
# 输出：
#   output/04_纵向队列/attrition_pattern_2026-09-20.csv
#   output/04_纵向队列/attrition_balance_2026-09-20.csv
#   output/04_纵向队列/selection_model_2026-09-20.csv
#   output/04_纵向队列/ipw_diagnostics_2026-09-20.csv
#   output/04_纵向队列/ipw_longitudinal_models_2026-09-20.csv
#   output/04_纵向队列/ipw_cumulative_incidence_2026-09-20.csv
#   output/04_纵向队列/attrition_selection_2026-09-20.md
#   output/09_图表/FigureS_attrition_selection_2026-09-20.png/.tif
#
# 重要口径：
#   目标人群 = 2011 年年龄 ≥45 岁且基线 ms10（定向+连续减7）可计分者。
#   满足该条件且五波 ms10 均可计分者为年龄合格固定队列 n=5,672。
#   旧 longitudinal 分析曾使用 5,884 人，其中含 212 名年龄<45或年龄缺失者；
#   2026-09-20 已将 T04-T06 统一为年龄合格队列口径。
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

w <- readRDS(file.path(P, "output", "02_五波数据与质控", "waves_cognition_2026-09-13.rds"))
w$KEY <- ifelse(is.na(w$IDn), w$ID, w$IDn)
w$ms10 <- as.numeric(w$orient) + as.numeric(w$ser7)
stopifnot(!anyDuplicated(w[, c("KEY", "year")]))

# 五波认知可计分次数。
nwaves <- w %>% filter(!is.na(ms10)) %>% count(KEY, name = "n_waves")
b2011 <- w %>% filter(year == 2011)
eligible_all <- b2011 %>% filter(!is.na(age), age >= 45)
eligible_ms10 <- eligible_all %>% filter(!is.na(ms10)) %>% left_join(nwaves, by = "KEY") %>%
  mutate(fixed = as.integer(n_waves == 5))
eligible_all <- eligible_all %>% left_join(nwaves, by = "KEY") %>% mutate(fixed = as.integer(!is.na(n_waves) & n_waves == 5))

stopifnot(nrow(eligible_all) == 16831)
stopifnot(nrow(eligible_ms10) == 13710, sum(eligible_ms10$fixed) == 5672)

# ------------------------- 1. 失访模式与各波保留率 -------------------------
pattern_tab <- eligible_ms10 %>% count(n_waves, name = "n") %>%
  mutate(percent = 100 * n / sum(n)) %>%
  arrange(n_waves)
colnames(pattern_tab) <- c("waves_with_ms10", "n", "percent")

ids <- eligible_ms10$KEY
retention <- w %>% filter(KEY %in% ids) %>%
  group_by(year) %>%
  summarise(interviewed_n = n(), ms10_n = sum(!is.na(ms10)), .groups = "drop") %>%
  mutate(retained_percent = 100 * ms10_n / nrow(eligible_ms10))
write.csv(pattern_tab, file.path(OUT, "attrition_pattern_2026-09-20.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(retention, file.path(OUT, "attrition_retention_2026-09-20.csv"), row.names = FALSE, fileEncoding = "UTF-8")

# ------------------------- 2. 固定队列与未纳入者基线比较 -------------------------
balance_vars <- c("age", "female", "rural", "ms10", "cesd10", "adl", "sleep_h", "pain_any", "hibp", "diab")
balance_labels <- c(age = "年龄（岁）", female = "女性", rural = "农村", ms10 = "基线 ms10",
                    cesd10 = "CES-D-10 得分", adl = "ADL 受限", sleep_h = "睡眠时长（小时）",
                    pain_any = "疼痛", hibp = "高血压", diab = "糖尿病")

smd_cont <- function(x1, x2) {
  x1 <- x1[!is.na(x1)]; x2 <- x2[!is.na(x2)]
  if (length(x1) < 2 || length(x2) < 2) return(NA_real_)
  s <- sqrt((var(x1) + var(x2)) / 2)
  if (s == 0) return(0)
  (mean(x1) - mean(x2)) / s
}
wmean <- function(x, w) sum(x * w, na.rm = TRUE) / sum(w[!is.na(x)], na.rm = TRUE)
wvar <- function(x, w) {
  keep <- !is.na(x) & !is.na(w)
  x <- x[keep]; w <- w[keep]
  mu <- sum(w * x) / sum(w)
  sum(w * (x - mu)^2) / sum(w)
}
smd_weighted <- function(x1, w1, x2) {
  keep1 <- !is.na(x1) & !is.na(w1); keep2 <- !is.na(x2)
  x1 <- x1[keep1]; w1 <- w1[keep1]; x2 <- x2[keep2]
  if (length(x1) < 2 || length(x2) < 2) return(NA_real_)
  s <- sqrt((wvar(x1, w1) + var(x2)) / 2)
  if (s == 0) return(0)
  (wmean(x1, w1) - mean(x2)) / s
}

fixed <- eligible_ms10 %>% filter(fixed == 1)
nonfixed <- eligible_ms10 %>% filter(fixed == 0)
balance <- do.call(rbind, lapply(balance_vars, function(v) {
  x1 <- fixed[[v]]; x2 <- nonfixed[[v]]
  data.frame(
    variable = v, label = balance_labels[[v]],
    fixed_n = sum(!is.na(x1)), fixed_mean = mean(x1, na.rm = TRUE), fixed_sd = sd(x1, na.rm = TRUE),
    nonfixed_n = sum(!is.na(x2)), nonfixed_mean = mean(x2, na.rm = TRUE), nonfixed_sd = sd(x2, na.rm = TRUE),
    smd = smd_cont(x1, x2),
    stringsAsFactors = FALSE
  )
}))
write.csv(balance, file.path(OUT, "attrition_balance_2026-09-20.csv"), row.names = FALSE, fileEncoding = "UTF-8")

# ------------------------- 3. 选择概率模型 -------------------------
pred_vars <- c("age", "female", "edu_group", "rural", "ms10", "cesd10", "adl", "sleep_h", "pain_any", "hibp", "diab")
ps_data <- eligible_ms10 %>% filter(complete.cases(.[, pred_vars])) %>% mutate(fixed_int = as.integer(fixed))
ps_formula <- as.formula(paste("fixed_int ~", paste(c("age", "female", "edu_group", "rural", "ms10", "cesd10", "adl", "sleep_h", "pain_any", "hibp", "diab"), collapse = " + ")))
ps_fit <- glm(ps_formula, data = ps_data, family = binomial())
ps_data$p_hat <- pmin(pmax(predict(ps_fit, type = "response"), 0.02), 0.98)
ps_data$ipw_stab <- mean(ps_data$fixed_int) / ps_data$p_hat
ps_fixed <- ps_data %>% filter(fixed_int == 1)

selection_coef <- summary(ps_fit)$coefficients
selection_tab <- data.frame(
  term = rownames(selection_coef), estimate = selection_coef[, 1], se = selection_coef[, 2],
  p = selection_coef[, ncol(selection_coef)],
  OR = exp(selection_coef[, 1]),
  OR_lower = exp(selection_coef[, 1] - 1.96 * selection_coef[, 2]),
  OR_upper = exp(selection_coef[, 1] + 1.96 * selection_coef[, 2]),
  stringsAsFactors = FALSE
)
write.csv(selection_tab, file.path(OUT, "selection_model_2026-09-20.csv"), row.names = FALSE, fileEncoding = "UTF-8")

# IPW 诊断；有效样本量常用 sum(w)^2 / sum(w^2)。
ipw_weight <- ps_fixed$ipw_stab
ipw_diag <- data.frame(
  selection_model_n = nrow(ps_data), fixed_n = nrow(ps_fixed), selection_fraction = mean(ps_data$fixed_int),
  p_min = min(ps_data$p_hat), p_max = max(ps_data$p_hat), p_median = median(ps_data$p_hat),
  ipw_min = min(ipw_weight), ipw_median = median(ipw_weight), ipw_max = max(ipw_weight),
  ipw_mean = mean(ipw_weight), effective_n = sum(ipw_weight)^2 / sum(ipw_weight^2),
  stringsAsFactors = FALSE
)
write.csv(ipw_diag, file.path(OUT, "ipw_diagnostics_2026-09-20.csv"), row.names = FALSE, fileEncoding = "UTF-8")

# 用 IPW 重新加权固定队列，与基线完整病例池比较。
target <- ps_data %>% filter(fixed_int == 0 | fixed_int == 1)
balance_ipw <- do.call(rbind, lapply(balance_vars, function(v) {
  # target mean is the full eligible complete-case pool; IPW reweights fixed cases to it.
  x_pool <- target[[v]]
  x_fixed <- ps_fixed[[v]]
  w_ipw <- ps_fixed$ipw_stab
  data.frame(
    variable = v, label = balance_labels[[v]],
    fixed_unweighted = mean(x_fixed, na.rm = TRUE),
    pool_unweighted = mean(x_pool, na.rm = TRUE),
    fixed_ipw_weighted = wmean(x_fixed, w_ipw),
    smd_before = smd_cont(x_fixed, x_pool),
    smd_after = smd_weighted(x_fixed, w_ipw, x_pool),
    stringsAsFactors = FALSE
  )
}))
balance_ipw$smd_abs_before <- abs(balance_ipw$smd_before)
balance_ipw$smd_abs_after <- abs(balance_ipw$smd_after)
write.csv(balance_ipw, file.path(OUT, "attrition_balance_ipw_2026-09-20.csv"), row.names = FALSE, fileEncoding = "UTF-8")

# ------------------------- 4. 构建年龄合格固定队列结局数据 -------------------------
fixed_ids <- eligible_ms10$KEY[eligible_ms10$fixed == 1]
# 连续认知分数矩阵。
score_wide <- w %>% filter(KEY %in% fixed_ids) %>%
  select(KEY, year, ms10) %>%
  pivot_wider(names_from = year, values_from = ms10, names_prefix = "ms") %>%
  arrange(KEY)
# 2011 年教育分层固定切点。
base11 <- w %>% filter(KEY %in% fixed_ids, year == 2011) %>% select(KEY, edu_group)
thr <- w %>% filter(KEY %in% fixed_ids, year == 2011) %>%
  group_by(edu_group) %>% summarise(anchor = mean(ms10, na.rm = TRUE) - sd(ms10, na.rm = TRUE), .groups = "drop")
# 受损状态矩阵。
status_wide <- w %>% filter(KEY %in% fixed_ids) %>%
  select(KEY, year, ms10, edu_group) %>%
  left_join(thr, by = "edu_group") %>%
  mutate(impaired = as.integer(ms10 < anchor)) %>%
  select(KEY, year, impaired) %>%
  pivot_wider(names_from = year, values_from = impaired, names_prefix = "im") %>%
  arrange(KEY)
wide <- score_wide %>% left_join(status_wide, by = "KEY")

# 基线暴露和十年认知变化。
base_cov <- ps_data %>% select(KEY, ipw_stab, p_hat, fixed_int)
fixed_base <- w %>% filter(KEY %in% fixed_ids, year == 2011) %>%
  transmute(KEY, age, female, edu_group, rural,
            dep = as.integer(cesd10 >= 10), adl, sleep_h, pain = pain_any, pa = pa_any,
            hibp, diab, weight, communityID)
fixed_base <- fixed_base %>% left_join(base_cov, by = "KEY") %>% left_join(
  wide %>% transmute(
    KEY,
    delta = ms2020 - ms2011,
    first_index = case_when(
      im2013 == 1 ~ 1L, im2015 == 1 ~ 2L, im2018 == 1 ~ 3L, im2020 == 1 ~ 4L,
      TRUE ~ NA_integer_
    ),
    confirmed_index = case_when(
      im2013 == 1 & im2015 == 1 ~ 1L,
      im2015 == 1 & im2018 == 1 ~ 2L,
      im2018 == 1 & im2020 == 1 ~ 3L,
      TRUE ~ NA_integer_
    ),
    endpoint_2020 = im2020,
    ms10_2011 = ms2011
  ), by = "KEY"
) %>% mutate(
  time = ifelse(is.na(first_index), 9, c(2, 4, 7, 9)[first_index]),
  event = as.integer(!is.na(first_index)),
  ipw_combined = weight * ipw_stab
)
# 在固定队列内标准化 IPW 组合权重；回归中的常数缩放不改变估计。
fixed_base$ipw_combined <- fixed_base$ipw_combined / mean(fixed_base$ipw_combined, na.rm = TRUE)

# ------------------------- 5. IPW 纵向模型 -------------------------
exposures <- c("dep", "adl", "sleep_h", "pain", "pa", "hibp", "diab")
labels <- c(dep = "抑郁（CES-D-10 ≥10）", adl = "ADL 受限", sleep_h = "睡眠时长（每 1 小时）",
            pain = "疼痛", pa = "规律体力活动", hibp = "高血压", diab = "糖尿病")
extract_est <- function(fit, term, type = c("ratio", "beta")) {
  type <- match.arg(type); b <- coef(fit)[term]; se <- sqrt(vcov(fit)[term, term])
  if (type == "ratio") c(exp(b), exp(b - 1.96 * se), exp(b + 1.96 * se))
  else c(b, b - 1.96 * se, b + 1.96 * se)
}
ipw_rows <- list()
for (e in exposures) {
  vars_inc <- c(e, "age", "female", "edu_group", "endpoint_2020", "time", "event", "ipw_combined", "communityID")
  vars_change <- c(e, "age", "female", "edu_group", "delta", "ipw_combined", "communityID")
  d_inc <- fixed_base[complete.cases(fixed_base[, vars_inc]), ]
  d_change <- fixed_base[complete.cases(fixed_base[, vars_change]), ]
  des_inc <- svydesign(ids = ~communityID, weights = ~ipw_combined, data = d_inc)
  des_change <- svydesign(ids = ~communityID, weights = ~ipw_combined, data = d_change)
  rhs <- paste(e, "+ age + female + edu_group")
  f1 <- as.formula(paste("endpoint_2020 ~", rhs))
  f2 <- as.formula(paste("Surv(time, event) ~", rhs))
  f3 <- as.formula(paste("delta ~", rhs))
  wl <- svyglm(f1, design = des_inc, family = quasibinomial())
  wc <- svycoxph(f2, design = des_inc)
  wch <- svyglm(f3, design = des_change)
  u1 <- glm(f1, data = d_inc, family = binomial())
  u2 <- coxph(f2, data = d_inc)
  u3 <- lm(f3, data = d_change)
  specs <- list(
    list(model = "endpoint_2020_logistic", label = "2020 终点评定", measure = "OR",
         w = extract_est(wl, e, "ratio"), u = extract_est(u1, e, "ratio"),
         n = nrow(d_inc), events = sum(d_inc$endpoint_2020), communities = length(unique(d_inc$communityID))),
    list(model = "first_event_cox", label = "首次受损 Cox", measure = "HR",
         w = extract_est(wc, e, "ratio"), u = extract_est(u2, e, "ratio"),
         n = nrow(d_inc), events = sum(d_inc$event), communities = length(unique(d_inc$communityID))),
    list(model = "ten_year_change_lm", label = "十年变化（固定队列完整病例）", measure = "beta",
         w = extract_est(wch, e, "beta"), u = extract_est(u3, e, "beta"),
         n = nrow(d_change), events = nrow(d_change), communities = length(unique(d_change$communityID)))
  )
  for (z in specs) {
    ipw_rows[[length(ipw_rows) + 1L]] <- data.frame(
      model = z$model, model_label = z$label, exposure = e, exposure_label = labels[[e]],
      effect_measure = z$measure, n = z$n, events = z$events, communities = z$communities,
      ipw_est = z$w[1], ipw_lower = z$w[2], ipw_upper = z$w[3],
      unweighted_est = z$u[1], unweighted_lower = z$u[2], unweighted_upper = z$u[3],
      stringsAsFactors = FALSE
    )
  }
}
ipw_models <- bind_rows(ipw_rows)
write.csv(ipw_models, file.path(OUT, "ipw_longitudinal_models_2026-09-20.csv"), row.names = FALSE, fileEncoding = "UTF-8")

# ------------------------- 6. IPW 累积发病 -------------------------
time_points <- c(0,2,4,7,9); first_limits <- c(0L,1L,2L,3L,4L); conf_limits <- c(0L,1L,2L,3L,3L)
ipw_inc <- list()
for (k in seq_along(time_points)) {
  d <- fixed_base
  d$e_first <- as.integer(!is.na(d$first_index) & d$first_index <= first_limits[k])
  d$e_conf <- as.integer(!is.na(d$confirmed_index) & d$confirmed_index <= conf_limits[k])
  des <- svydesign(ids = ~communityID, weights = ~ipw_combined, data = d[!is.na(d$ipw_combined), ])
  m1 <- svymean(~e_first, des, na.rm = TRUE); m2 <- svymean(~e_conf, des, na.rm = TRUE)
  ipw_inc[[k]] <- data.frame(
    time_years = rep(time_points[k], 2), outcome = c("first_episode", "confirmed_first_episode"),
    ipw_cumulative_incidence = c(as.numeric(m1), as.numeric(m2)),
    unweighted_cumulative_incidence = c(mean(d$e_first, na.rm = TRUE), mean(d$e_conf, na.rm = TRUE)),
    stringsAsFactors = FALSE
  )
}
ipw_inc <- bind_rows(ipw_inc)
write.csv(ipw_inc, file.path(OUT, "ipw_cumulative_incidence_2026-09-20.csv"), row.names = FALSE, fileEncoding = "UTF-8")

# ------------------------- 7. 图表：保留率 + SMD + IPW 森林图 -------------------------
draw_attrition_figure <- function() {
  par(mfrow = c(1, 3), mar = c(5, 7, 3, 1), family = "sans")
  # Retention.
  plot(retention$year, retention$retained_percent, type = "b", pch = 16, lwd = 2, col = "#1F4E78",
       ylim = c(0, 100), xaxt = "n", ylab = "ms10 available (%)", xlab = "Year", main = "Follow-up availability")
  axis(1, at = retention$year); grid(col = "grey88")
  # SMD before/after.
  z <- balance_ipw[order(abs(balance_ipw$smd_before)), ]
  yy <- seq_len(nrow(z)); xlim <- range(c(z$smd_before, z$smd_after, -0.1, 0.1), na.rm = TRUE)
  plot(NA, xlim = xlim, ylim = c(0.5, nrow(z) + 0.5), yaxt = "n", xlab = "Standardized mean difference", ylab = "", main = "Covariate balance")
  abline(v = c(-0.1, 0, 0.1), lty = c(2,1,2), col = c("grey60","grey40","grey60"))
  segments(z$smd_before, yy - 0.10, z$smd_after, yy + 0.10, col = "grey80")
  points(z$smd_before, yy - 0.10, pch = 17, col = "#C00000", cex = 0.9)
  points(z$smd_after, yy + 0.10, pch = 16, col = "#1F4E78", cex = 0.9)
  axis(2, at = yy, labels = z$label, las = 1, cex.axis = 0.75)
  legend("topright", legend = c("Before IPW", "After IPW"), col = c("#C00000","#1F4E78"), pch = c(17,16), bty = "n", cex = 0.75)
  # IPW vs unweighted Cox.
  z2 <- ipw_models[ipw_models$model == "first_event_cox", ]
  z2 <- z2[order(z2$ipw_est), ]; yy <- seq_len(nrow(z2))
  xlim2 <- range(c(z2$ipw_lower, z2$ipw_upper, z2$unweighted_lower, z2$unweighted_upper, 1), na.rm = TRUE)
  plot(NA, xlim = xlim2, ylim = c(0.5, nrow(z2)+0.5), yaxt = "n", xlab = "Hazard ratio", ylab = "", main = "IPW sensitivity: first event")
  abline(v = 1, lty = 3, col = "grey60")
  segments(z2$unweighted_lower, yy - 0.1, z2$unweighted_upper, yy - 0.1, col = "#C00000", lwd = 2)
  points(z2$unweighted_est, yy - 0.1, pch = 17, col = "#C00000")
  segments(z2$ipw_lower, yy + 0.1, z2$ipw_upper, yy + 0.1, col = "#1F4E78", lwd = 2)
  points(z2$ipw_est, yy + 0.1, pch = 16, col = "#1F4E78")
  axis(2, at = yy, labels = z2$exposure_label, las = 1, cex.axis = 0.75)
  legend("topright", legend = c("Unweighted", "IPW (selection-adjusted)"), col = c("#C00000","#1F4E78"), pch = c(17,16), bty = "n", cex = 0.72)
}
png(file.path(FIG_DIR, "FigureS_attrition_selection_2026-09-20.png"), width = 15, height = 6, units = "in", res = 300)
draw_attrition_figure(); dev.off()
tryCatch({ tiff(file.path(FIG_DIR, "FigureS_attrition_selection_2026-09-20.tif"), width = 15, height = 6, units = "in", res = 600, compression = "lzw"); draw_attrition_figure(); dev.off() }, error = function(e) { if (length(dev.list())) dev.off() })

# ------------------------- 8. 报告 -------------------------
fmt_pct <- function(x) sprintf("%.1f%%", 100*x)
report <- c(
  "# T07：失访、固定队列选择偏倚与 IPW 敏感性分析（2026-09-20）",
  "",
  "> 目标人群：2011 年年龄 ≥45 岁且基线 ms10 可计分者。",
  sprintf("> 基线可分析者 n=%d；年龄合格固定队列 n=%d；旧 5,884 人队列含 212 名年龄<45或年龄缺失者。", nrow(eligible_ms10), nrow(fixed_base)),
  "",
  "## 一、失访模式",
  "",
  "| 五波中可计分波数 | 人数 | 比例 |",
  "|---:|---:|---:|",
  paste0("| ", pattern_tab[[1]], " | ", pattern_tab[[2]], " | ", sprintf("%.1f%%", pattern_tab[[3]]), " |"),
  "",
  "| 年份 | 实际访问人数 | ms10 可计分 | 相对 2011 保留 |",
  "|---:|---:|---:|---:|",
  paste0("| ", retention$year, " | ", retention$interviewed_n, " | ", retention$ms10_n, " | ", fmt_pct(retention$retained_percent/100), " |"),
  "",
  "## 二、固定队列与未纳入者基线差异",
  "",
  "| 指标 | 固定队列 | 未纳入 | SMD |",
  "|---|---:|---:|---:|",
  paste0("| ", balance$label, " | ", sprintf("%.2f", balance$fixed_mean), " | ", sprintf("%.2f", balance$nonfixed_mean), " | ", sprintf("%.2f", balance$smd), " |"),
  "",
  "> SMD 绝对值越接近 0 越均衡；|SMD|≥0.10 提示有意义的基线差异。",
  "",
  "## 三、选择概率与 IPW",
  "",
  sprintf("- 选择模型完整病例：%d 人；其中固定队列 %d 人。", ipw_diag$selection_model_n, ipw_diag$fixed_n),
  sprintf("- 截断后的预测选择概率：中位数 %.2f，范围 %.2f–%.2f。", ipw_diag$p_median, ipw_diag$p_min, ipw_diag$p_max),
  sprintf("- 稳定化 IPW：中位数 %.2f，最大 %.2f，有效样本量 %.0f。", ipw_diag$ipw_median, ipw_diag$ipw_max, ipw_diag$effective_n),
  "",
  "| 指标 | 固定队列未加权 | 全池未加权 | IPW 重加权 | SMD 前 | SMD 后 |",
  "|---|---:|---:|---:|---:|---:|",
  paste0("| ", balance_ipw$label, " | ", sprintf("%.2f", balance_ipw$fixed_unweighted), " | ", sprintf("%.2f", balance_ipw$pool_unweighted), " | ", sprintf("%.2f", balance_ipw$fixed_ipw_weighted), " | ", sprintf("%.2f", balance_ipw$smd_before), " | ", sprintf("%.2f", balance_ipw$smd_after), " |"),
  "",
  "## 四、IPW 纵向模型",
  "",
  "| 模型 | 暴露 | n | IPW 估计 | 未加权 |",
  "|---|---|---:|---:|---:|",
  paste0("| ", ipw_models$model_label, " | ", ipw_models$exposure_label, " | ", ipw_models$n, " | ",
         sprintf("%.2f (%.2f–%.2f)", ipw_models$ipw_est, ipw_models$ipw_lower, ipw_models$ipw_upper), " | ",
         sprintf("%.2f (%.2f–%.2f)", ipw_models$unweighted_est, ipw_models$unweighted_lower, ipw_models$unweighted_upper), " |"),
  "",
  "## 五、结论",
  "",
  sprintf("- 年龄合格固定队列占基线可分析者 %.1f%%。", 100 * nrow(fixed_base) / nrow(eligible_ms10)),
  sprintf("- IPW 后 SMD 平均绝对值由 %.2f 变为 %.2f。", mean(abs(balance_ipw$smd_before)), mean(abs(balance_ipw$smd_after))),
  "- 固定队列明显更年轻、更健康、基线认知更高，选择偏倚真实存在。",
  "- IPW 只能校正可观测基线因素，不能校正未测量选择、死亡和认知测量缺失。",
  "- T04–T06 已按同一年龄合格口径重跑；本报告中的 5,672 人与当前纵向主分析一致。",
  "",
  "## 六、输出文件",
  "",
  "- `attrition_pattern_2026-09-20.csv`",
  "- `attrition_retention_2026-09-20.csv`",
  "- `attrition_balance_2026-09-20.csv`",
  "- `attrition_balance_ipw_2026-09-20.csv`",
  "- `selection_model_2026-09-20.csv`",
  "- `ipw_diagnostics_2026-09-20.csv`",
  "- `ipw_longitudinal_models_2026-09-20.csv`",
  "- `ipw_cumulative_incidence_2026-09-20.csv`",
  "- `output/09_图表/FigureS_attrition_selection_2026-09-20.png/.tif`",
  "",
  sprintf("*生成时间：%s*", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
)
writeLines(report, file.path(OUT, "attrition_selection_2026-09-20.md"), useBytes = TRUE)

cat("\n========== T07 失访与选择偏倚分析完成 ==========\n")
cat("基线≥45岁且ms10可分析：", nrow(eligible_ms10), "；年龄合格固定队列：", nrow(fixed_base), "\n", sep = "")
cat("选择模型完整病例：", nrow(ps_data), "；固定队列：", nrow(ps_fixed), "\n", sep = "")
cat("IPW平均|SMD|：", sprintf("%.2f -> %.2f", mean(abs(balance_ipw$smd_before)), mean(abs(balance_ipw$smd_after))), "\n", sep = "")
cat("主图：", file.path(FIG_DIR, "FigureS_attrition_selection_2026-09-20.png"), "\n", sep = "")
