# ============================================================
# 血液统计第三步：PLT/WBC/CRP 与十年认知受损/变化关联
# 主分析：
#   1) 2011 PLT/WBC/logCRP 预测 2020 新发认知受损
#   2) 2011 PLT/WBC/logCRP 与 2011→2020 认知变化
#   3) 2011→2015 血液变化预测 2015→2020 认知变化（探索）
# 说明：
#   - 连续暴露按 1 SD 标准化；
#   - CRP 使用 log(CRP+1)，四分位仍按原始值排序；
#   - 未加权为主，血液权重用于后续敏感性；
#   - 极端值采用 1%/99% winsorization 敏感性分析。
# ============================================================

invisible(try(Sys.setlocale("LC_ALL", "zh_CN.UTF-8"), silent = TRUE))
.libPaths(c("D:/R/library", .libPaths()))
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(openxlsx)
})

in_analysis <- "output/07_血液指标/analysis_blood_two_waves_fixedcohort_2026-09-19.rds"
in_cohort   <- "output/04_纵向队列/cohort_change_2026-09-13.rds"
out_dir     <- "output/07_血液指标"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

norm_id <- function(x) {
  x <- as.character(x)
  ifelse(nchar(x) == 11, paste0(substr(x, 1, 9), sprintf("%03d", as.integer(substr(x, 10, 11)))), x)
}

a <- readRDS(in_analysis)
cohort_obj <- readRDS(in_cohort)
cohort_long <- cohort_obj$cohort %>% mutate(KEY = norm_id(ifelse(!is.na(IDn), IDn, ID)))
base_cov <- cohort_long %>%
  filter(year == 2011) %>%
  transmute(
    KEY,
    rural = as.numeric(rural),
    adl = as.numeric(adl),
    cesd10 = as.numeric(cesd10),
    hibp = as.numeric(hibp),
    diab = as.numeric(diab),
    heart = as.numeric(heart),
    stroke = as.numeric(stroke),
    lung = as.numeric(lung),
    kidney = as.numeric(kidney),
    dyslip = as.numeric(dyslip),
    sleep_h = as.numeric(sleep_h),
    pain_any = as.numeric(pain_any),
    communityID = as.character(communityID)
  )

d <- a %>% left_join(base_cov, by = "KEY") %>%
  mutate(
    edu_group = factor(edu_group),
    edu_group = relevel(edu_group, ref = "Less than primary"),
    y_change_10 = y2020 - y2011,
    y_change_5 = y2020 - y2015,
    logcrp_2011 = log1p(pmax(crp_2011, 0)),
    logcrp_2015 = log1p(pmax(crp_2015, 0)),
    z_plt_2011 = as.numeric(scale(plt_2011)),
    z_wbc_2011 = as.numeric(scale(wbc_2011)),
    z_logcrp_2011 = as.numeric(scale(logcrp_2011)),
    z_dplt = as.numeric(scale(d_plt)),
    z_dwbc = as.numeric(scale(d_wbc)),
    z_dlogcrp = as.numeric(scale(d_logcrp)),
    q_plt_2011 = ntile(plt_2011, 4),
    q_wbc_2011 = ntile(wbc_2011, 4),
    q_crp_2011 = ntile(crp_2011, 4),
    q_dplt = ntile(d_plt, 4),
    q_dwbc = ntile(d_wbc, 4),
    q_dcrp = ntile(d_logcrp, 4),
    q_plt_2011 = factor(q_plt_2011, levels = 1:4, labels = paste0("Q", 1:4)),
    q_wbc_2011 = factor(q_wbc_2011, levels = 1:4, labels = paste0("Q", 1:4)),
    q_crp_2011 = factor(q_crp_2011, levels = 1:4, labels = paste0("Q", 1:4)),
    q_dplt = factor(q_dplt, levels = 1:4, labels = paste0("Q", 1:4)),
    q_dwbc = factor(q_dwbc, levels = 1:4, labels = paste0("Q", 1:4)),
    q_dcrp = factor(q_dcrp, levels = 1:4, labels = paste0("Q", 1:4))
  )

# 1%/99% winsorization 敏感性变量
wins <- function(x) {
  lo <- quantile(x, 0.01, na.rm = TRUE, names = FALSE)
  hi <- quantile(x, 0.99, na.rm = TRUE, names = FALSE)
  pmin(pmax(x, lo), hi)
}
d <- d %>%
  mutate(
    plt_2011_w = wins(plt_2011),
    wbc_2011_w = wins(wbc_2011),
    crp_2011_w = wins(crp_2011),
    logcrp_2011_w = log1p(pmax(crp_2011_w, 0)),
    z_plt_2011_w = as.numeric(scale(plt_2011_w)),
    z_wbc_2011_w = as.numeric(scale(wbc_2011_w)),
    z_logcrp_2011_w = as.numeric(scale(logcrp_2011_w))
  )

# 基础模型与Fully-adjusted模型
base_covars <- c("age2011", "female", "edu_group", "rural", "y2011")
full_covars <- c(base_covars, "cesd10", "adl", "hibp", "diab", "heart", "stroke", "lung", "kidney", "dyslip", "sleep_h", "pain_any")

extract_glm <- function(fit, exposure_label, model_label, outcome_label) {
  co <- coef(summary(fit))
  ci <- confint.default(fit)
  keep <- grepl(paste0("^", attr(terms(fit), "term.labels")[1]), rownames(co))
  data.frame(
    outcome = outcome_label,
    model = model_label,
    exposure = exposure_label,
    term = rownames(co)[keep],
    n = length(fit$fitted.values),
    events = if (!is.null(fit$y)) sum(fit$y == 1, na.rm = TRUE) else NA_integer_,
    estimate = co[keep, 1],
    se = co[keep, 2],
    p = co[keep, 4],
    OR = exp(co[keep, 1]),
    lower = exp(ci[keep, 1]),
    upper = exp(ci[keep, 2]),
    stringsAsFactors = FALSE
  )
}

extract_lm <- function(fit, exposure_label, model_label, outcome_label) {
  co <- coef(summary(fit))
  ci <- confint.default(fit)
  keep <- grepl(paste0("^", attr(terms(fit), "term.labels")[1]), rownames(co))
  data.frame(
    outcome = outcome_label,
    model = model_label,
    exposure = exposure_label,
    term = rownames(co)[keep],
    n = length(fit$fitted.values),
    events = NA_integer_,
    estimate = co[keep, 1],
    se = co[keep, 2],
    p = co[keep, 4],
    beta = co[keep, 1],
    lower = ci[keep, 1],
    upper = ci[keep, 2],
    stringsAsFactors = FALSE
  )
}

fit_logit <- function(data, exposure, covars) {
  vars <- c("incident_2020", exposure, covars)
  dd <- data %>% select(all_of(vars)) %>% drop_na()
  # 新发分析限定基线未受损
  dd <- dd[data$impaired_2011[match(rownames(dd), rownames(data))] == 0, , drop = FALSE]
  f <- as.formula(paste("incident_2020 ~", paste(c(exposure, covars), collapse = " + ")))
  list(fit = glm(f, data = dd, family = binomial()), data = dd)
}

# 新发模型：更稳妥地单独处理基线受损限制
make_incident_data <- function(data, exposure, covars) {
  vars <- unique(c("incident_2020", "impaired_2011", exposure, covars))
  data %>% select(all_of(vars)) %>% filter(impaired_2011 == 0) %>% drop_na()
}
fit_logit2 <- function(data, exposure, covars) {
  dd <- make_incident_data(data, exposure, covars)
  f <- as.formula(paste("incident_2020 ~", paste(c(exposure, covars), collapse = " + ")))
  list(fit = glm(f, data = dd, family = binomial()), data = dd)
}
fit_lm2 <- function(data, outcome, exposure, covars) {
  vars <- unique(c(outcome, exposure, covars))
  dd <- data %>% select(all_of(vars)) %>% drop_na()
  f <- as.formula(paste(outcome, "~", paste(c(exposure, covars), collapse = " + ")))
  list(fit = lm(f, data = dd), data = dd)
}

cont_exposures <- c("z_plt_2011", "z_wbc_2011", "z_logcrp_2011")
cont_labels <- c("PLT per 1 SD", "WBC per 1 SD", "log(CRP+1) per 1 SD")
quart_exposures <- c("q_plt_2011", "q_wbc_2011", "q_crp_2011")
quart_labels <- c("PLT quartile", "WBC quartile", "CRP quartile")

main_cont <- bind_rows(lapply(seq_along(cont_exposures), function(i) {
  ex <- cont_exposures[i]
  lab <- cont_labels[i]
  m1 <- fit_logit2(d, ex, base_covars)$fit
  m2 <- fit_logit2(d, ex, full_covars)$fit
  bind_rows(
    extract_glm(m1, lab, "Base-adjusted", "Incident cognitive impairment 2020"),
    extract_glm(m2, lab, "Fully-adjusted", "Incident cognitive impairment 2020")
  )
}))
main_cont <- bind_rows(
  main_cont,
  bind_rows(lapply(seq_along(cont_exposures), function(i) {
    ex <- cont_exposures[i]; lab <- cont_labels[i]
    m1 <- fit_lm2(d, "y_change_10", ex, base_covars)$fit
    m2 <- fit_lm2(d, "y_change_10", ex, full_covars)$fit
    bind_rows(
      extract_lm(m1, lab, "Base-adjusted", "Cognitive change 2011-2020"),
      extract_lm(m2, lab, "Fully-adjusted", "Cognitive change 2011-2020")
    )
  }))
)

main_quart <- bind_rows(lapply(seq_along(quart_exposures), function(i) {
  ex <- quart_exposures[i]; lab <- quart_labels[i]
  m1 <- fit_logit2(d, ex, base_covars)$fit
  m2 <- fit_logit2(d, ex, full_covars)$fit
  bind_rows(
    extract_glm(m1, lab, "Base-adjusted", "Incident cognitive impairment 2020"),
    extract_glm(m2, lab, "Fully-adjusted", "Incident cognitive impairment 2020")
  )
}))
main_quart <- bind_rows(
  main_quart,
  bind_rows(lapply(seq_along(quart_exposures), function(i) {
    ex <- quart_exposures[i]; lab <- quart_labels[i]
    m1 <- fit_lm2(d, "y_change_10", ex, base_covars)$fit
    m2 <- fit_lm2(d, "y_change_10", ex, full_covars)$fit
    bind_rows(
      extract_lm(m1, lab, "Base-adjusted", "Cognitive change 2011-2020"),
      extract_lm(m2, lab, "Fully-adjusted", "Cognitive change 2011-2020")
    )
  }))
)

# 变化暴露：以 2015 认知为基线，检验 2011→2015 血液变化
change_exposures <- c("z_dplt", "z_dwbc", "z_dlogcrp")
change_labels <- c("delta PLT per 1 SD", "delta WBC per 1 SD", "delta log(CRP+1) per 1 SD")
change_covars <- c("age2011", "female", "edu_group", "rural", "y2015")
change_models <- bind_rows(lapply(seq_along(change_exposures), function(i) {
  ex <- change_exposures[i]; lab <- change_labels[i]
  m1 <- fit_lm2(d, "y_change_5", ex, change_covars)$fit
  m2 <- fit_lm2(d, "y_change_5", ex, c(change_covars, "cesd10", "adl", "hibp", "diab", "heart", "stroke", "lung", "kidney", "dyslip", "sleep_h", "pain_any"))$fit
  bind_rows(
    extract_lm(m1, lab, "Base-adjusted", "Cognitive change 2015-2020"),
    extract_lm(m2, lab, "Fully-adjusted", "Cognitive change 2015-2020")
  )
}))

# 极端值敏感性：1%/99% winsorization
sens_cont <- bind_rows(lapply(seq_along(cont_exposures), function(i) {
  ex <- c("z_plt_2011_w", "z_wbc_2011_w", "z_logcrp_2011_w")[i]
  lab <- cont_labels[i]
  m1 <- fit_logit2(d, ex, base_covars)$fit
  m2 <- fit_lm2(d, "y_change_10", ex, base_covars)$fit
  bind_rows(
    extract_glm(m1, lab, "Winsorized base-adjusted", "Incident cognitive impairment 2020"),
    extract_lm(m2, lab, "Winsorized base-adjusted", "Cognitive change 2011-2020")
  )
}))

# FDR：主连续暴露，按结局/模型家族分别校正
add_fdr <- function(df) {
  out <- df
  out$q <- NA_real_
  keys <- unique(out[, c("outcome", "model")])
  for (i in seq_len(nrow(keys))) {
    idx <- which(out$outcome == keys$outcome[i] & out$model == keys$model[i] & grepl("per 1 SD", out$exposure))
    if (length(idx)) out$q[idx] <- p.adjust(out$p[idx], method = "BH")
  }
  out
}
main_cont <- add_fdr(main_cont)
change_models <- add_fdr(change_models)
sens_cont <- add_fdr(sens_cont)

# 输出 CSV（数据标签全部使用 ASCII，避免中文编码问题）
write.csv(main_cont, file.path(out_dir, "blood_main_continuous_2026-09-19.csv"), row.names = FALSE)
write.csv(main_quart, file.path(out_dir, "blood_main_quartile_2026-09-19.csv"), row.names = FALSE)
write.csv(change_models, file.path(out_dir, "blood_change_models_2026-09-19.csv"), row.names = FALSE)
write.csv(sens_cont, file.path(out_dir, "blood_winsor_sensitivity_2026-09-19.csv"), row.names = FALSE)

fmt <- function(x, d = 3) formatC(x, format = "f", digits = d)
ci_or <- function(r) paste0(fmt(r$OR), " (", fmt(r$lower), "-", fmt(r$upper), ")")
ci_b <- function(r) paste0(fmt(r$beta), " (", fmt(r$lower, 3), "-", fmt(r$upper, 3), ")")
exposure_cn <- c(
  "PLT per 1 SD" = "PLT（每1 SD）",
  "WBC per 1 SD" = "WBC（每1 SD）",
  "log(CRP+1) per 1 SD" = "log(CRP+1)（每1 SD）",
  "PLT quartile" = "PLT四分位",
  "WBC quartile" = "WBC四分位",
  "CRP quartile" = "CRP四分位",
  "delta PLT per 1 SD" = "ΔPLT（每1 SD）",
  "delta WBC per 1 SD" = "ΔWBC（每1 SD）",
  "delta log(CRP+1) per 1 SD" = "Δlog(CRP+1)（每1 SD）"
)
tr <- function(x) ifelse(x %in% names(exposure_cn), unname(exposure_cn[x]), x)

lines <- c(
  "# 血液统计第三步 · PLT/WBC/CRP 与认知结局",
  "",
  "> 生成日期：2026-09-20（年龄口径修正后重算）",
  "",
  "## 一、主模型结果摘要",
  "",
  "### 2020 新发认知受损（Base-adjusted）",
  "",
  "| 暴露 | OR (95%CI) | p | q | n |",
  "|---|---:|---:|---:|---:|"
)
tmp <- main_cont %>% filter(outcome == "Incident cognitive impairment 2020", model == "Base-adjusted")
for (i in seq_len(nrow(tmp))) lines <- c(lines, paste0("| ", tr(tmp$exposure[i]), " | ", ci_or(tmp[i, ]), " | ", fmt(tmp$p[i], 4), " | ", fmt(tmp$q[i], 4), " | ", tmp$n[i], " |"))
lines <- c(lines, "", "### 2011→2020 认知变化（Base-adjusted）", "", "| 暴露 | β (95%CI) | p | q | n |", "|---|---:|---:|---:|---:|")
tmp <- main_cont %>% filter(outcome == "Cognitive change 2011-2020", model == "Base-adjusted")
for (i in seq_len(nrow(tmp))) lines <- c(lines, paste0("| ", tr(tmp$exposure[i]), " | ", ci_b(tmp[i, ]), " | ", fmt(tmp$p[i], 4), " | ", fmt(tmp$q[i], 4), " | ", tmp$n[i], " |"))
lines <- c(lines, "", "### 2011→2015 血液变化 → 2015→2020 认知变化（Base-adjusted）", "", "| 暴露 | β (95%CI) | p | q | n |", "|---|---:|---:|---:|---:|")
tmp <- change_models %>% filter(model == "Base-adjusted")
for (i in seq_len(nrow(tmp))) lines <- c(lines, paste0("| ", tr(tmp$exposure[i]), " | ", ci_b(tmp[i, ]), " | ", fmt(tmp$p[i], 4), " | ", fmt(tmp$q[i], 4), " | ", tmp$n[i], " |"))
lines <- c(lines, "", "## 二、四分位结果", "", "四分位模型见 Excel 的 `Main_quartile`；Q1 为参照。", "", "## 三、敏感性分析", "", "1%/99% winsorization 结果见 `Winsorized`。", "", "## 四、解读边界", "", "- 血液指标为观察性关联，不是因果效应。", "- 血样子样本存在健康参与者/存活者选择，需在论文中披露。", "- Neu/NLR/PLR 不可从 CHARLS 已释放数据计算。", "- 下一步应补：血液权重、纵向缺失、极端值排除方案及与主论文 Table 5/Figure 6 的对接。")
writeLines(lines, file.path(out_dir, "blood_analysis_report_2026-09-19.md"), useBytes = TRUE)

cat("完成：blood model CSV outputs\\n")
cat("主连续模型：\n"); print(main_cont %>% filter(model == "Base-adjusted") %>% select(outcome, exposure, OR, lower, upper, p, q, n))
cat("变化模型：\n"); print(change_models %>% filter(model == "Base-adjusted") %>% select(exposure, beta, lower, upper, p, q, n))
