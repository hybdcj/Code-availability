# F39 TRIPOD补充分析（2026-10-17）
# 仓库适配：由 tools/prepare_workspace.py 替换项目路径及 R 库路径。
# 来自本地F39脚本；统计表达式保持不变，不再依赖临时目录链接。
options(warn = 1)
.libPaths(c("D:/R/library", .libPaths()))
suppressPackageStartupMessages({ library(dplyr); library(pROC) })
options(stringsAsFactors = FALSE)
set.seed(20261017)

OUT <- "D:/BaiduSyncdisk/科研/徐中豪/16.论文修订第三轮/F39_TRIPOD补充分析_2026-10-17"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
logf <- file.path(OUT, "F39_run_log.txt")
log_con <- file(logf, open = "wt", encoding = "UTF-8")
logmsg <- function(...) {
  msg <- paste0(...)
  cat(msg, file = stdout())
  cat(msg, file = log_con)
}
logmsg("=== F39 TRIPOD supplemental ===\n")
logmsg("start: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

fixed_score <- function(d) {
  2 * pmax(floor((d$age - 45) / 5), 0) +
    ifelse(d$edu_group == "Primary school", -2,
      ifelse(d$edu_group == "Middle school", -7,
        ifelse(d$edu_group == "High school or above", -10, 0))) +
    3 * d$dep
}
b <- c(
  `(Intercept)` = -4.19524508778947,
  age = 0.0492039250378111,
  `edu_groupPrimary school` = -0.290105771000328,
  `edu_groupMiddle school` = -0.994711800202235,
  `edu_groupHigh school or above` = -1.43588388035287,
  dep = 0.412166484954086,
  adl = -0.0110119344848014
)
fixed_lp_manual <- function(d, b) {
  e <- as.character(d$edu_group)
  b[["(Intercept)"]] + b[["age"]] * d$age +
    b[["edu_groupPrimary school"]] * (e == "Primary school") +
    b[["edu_groupMiddle school"]] * (e == "Middle school") +
    b[["edu_groupHigh school or above"]] * (e == "High school or above") +
    b[["dep"]] * d$dep + b[["adl"]] * d$adl
}
auc_ci <- function(y, x) {
  ok <- !is.na(y) & !is.na(x)
  y <- y[ok]; x <- x[ok]
  if (length(y) < 10 || length(unique(y)) < 2)
    return(c(auc = NA_real_, low = NA_real_, high = NA_real_))
  r <- pROC::roc(response = y, predictor = x, direction = "<", quiet = TRUE)
  ci <- as.numeric(pROC::ci.auc(r))
  c(auc = as.numeric(pROC::auc(r)), low = ci[1], high = ci[3])
}
cal_metrics <- function(y, p) {
  ok <- !is.na(y) & !is.na(p) & p > 0 & p < 1
  y <- y[ok]; p <- p[ok]
  if (length(y) < 10 || length(unique(y)) < 2)
    return(c(intercept = NA, intercept_low = NA, intercept_high = NA,
             slope = NA, slope_low = NA, slope_high = NA, brier = NA))
  lp <- qlogis(p)
  fi <- glm(y ~ offset(lp), family = binomial())
  cii <- confint.default(fi)["(Intercept)", ]
  fs <- glm(y ~ lp, family = binomial())
  cis <- confint.default(fs)["lp", ]
  c(intercept = unname(coef(fi)[["(Intercept)"]]),
    intercept_low = unname(cii[1]), intercept_high = unname(cii[2]),
    slope = unname(coef(fs)[["lp"]]),
    slope_low = unname(cis[1]), slope_high = unname(cis[2]),
    brier = mean((p - y)^2))
}

logmsg("frozen coef from A14 CSV / A07 model_u:\n"); print(b)

a07 <- readRDS("D:/BaiduSyncdisk/科研/徐中豪/output/10_筛查指数/screening_index_memory_v3_2026-09-26.rds")
b07 <- coef(a07$model_u)
logmsg("max |b - b07| = ", max(abs(b - b07[names(b)])), "\n")
stopifnot(max(abs(b - b07[names(b)])) < 1e-8)

w <- readRDS("D:/BaiduSyncdisk/科研/徐中豪/output/02_五波数据与质控/waves_cognition_2026-09-13.rds")
w$KEY <- ifelse(is.na(w$IDn), w$ID, w$IDn)
w$ms10 <- as.numeric(w$orient) + as.numeric(w$ser7)

base <- w %>% filter(year == 2011, !is.na(age), age >= 45) %>%
  transmute(
    KEY, age, female,
    edu_group = factor(edu_group,
      levels = c("Less than primary", "Primary school", "Middle school", "High school or above")),
    rural, communityID, weight, cesd10, dep = as.integer(cesd10 >= 10),
    adl, mem2011 = imrc + dlrc
  )
m15 <- w %>% filter(year == 2015) %>% transmute(KEY, mem2015 = imrc + dlrc)
n_ms10 <- w %>% filter(!is.na(ms10)) %>% count(KEY, name = "n_ms10")
base <- base %>% left_join(m15, by = "KEY") %>% left_join(n_ms10, by = "KEY") %>%
  mutate(fixed5 = as.integer(!is.na(n_ms10) & n_ms10 == 5))
rm(w); invisible(gc(verbose = FALSE))

pool <- base %>% filter(!is.na(mem2011), mem2011 > 4)
main <- pool %>%
  filter(fixed5 == 1, !is.na(mem2015),
         complete.cases(.[, c("age", "edu_group", "dep", "adl")])) %>%
  transmute(KEY, age, edu_group, dep, adl, mem2011, mem2015,
            event = as.integer(mem2015 <= 4), female, rural, cesd10)
broad <- pool %>%
  filter(!is.na(mem2015), complete.cases(.[, c("age", "edu_group", "dep", "adl")])) %>%
  transmute(KEY, age, edu_group, dep, adl, mem2011, mem2015,
            event = as.integer(mem2015 <= 4), female, rural, cesd10, fixed5)

logmsg("flow: base=", nrow(base), " mem2011=", sum(!is.na(base$mem2011)),
       " pool=", nrow(pool), " mem2015=", sum(!is.na(pool$mem2015)),
       " broad=", nrow(broad), " main=", nrow(main), " events=", sum(main$event), "\n")
stopifnot(nrow(main) == 4441, sum(main$event) == 658, nrow(broad) == 8296)
stopifnot(setequal(main$KEY, a07$data$KEY))
main$score <- fixed_score(main)
stopifnot(sum(main$score != a07$data$score[match(main$KEY, a07$data$KEY)], na.rm = TRUE) == 0)
main$lp <- fixed_lp_manual(main, b)
main$p_hat <- plogis(main$lp)
rm(a07); invisible(gc(verbose = FALSE))

ov_m <- auc_ci(main$event, main$p_hat)
ov_s <- auc_ci(main$event, main$score)
logmsg(sprintf("overall model AUC=%.6f (%.6f-%.6f)\n", ov_m["auc"], ov_m["low"], ov_m["high"]))
logmsg(sprintf("overall score AUC=%.6f (%.6f-%.6f)\n", ov_s["auc"], ov_s["low"], ov_s["high"]))
if (is.na(ov_m["auc"]) || abs(ov_m["auc"] - 0.702) > 0.005) {
  logmsg("ERROR: overall AUC mismatch; stop.\n"); close(log_con); quit(save = "no", status = 1)
}

one_row <- function(d, cohort, subgroup, level) {
  y <- d$event; p <- d$p_hat; sc <- d$score
  a_m <- auc_ci(y, p); a_s <- auc_ci(y, sc); cal <- cal_metrics(y, p)
  data.frame(
    cohort = cohort, subgroup = subgroup, level = as.character(level),
    n = nrow(d), events = sum(y, na.rm = TRUE), event_rate = mean(y, na.rm = TRUE),
    auc_model = unname(a_m["auc"]), auc_model_low = unname(a_m["low"]), auc_model_high = unname(a_m["high"]),
    auc_score = unname(a_s["auc"]), auc_score_low = unname(a_s["low"]), auc_score_high = unname(a_s["high"]),
    cal_intercept = unname(cal["intercept"]), cal_intercept_low = unname(cal["intercept_low"]),
    cal_intercept_high = unname(cal["intercept_high"]),
    cal_slope = unname(cal["slope"]), cal_slope_low = unname(cal["slope_low"]), cal_slope_high = unname(cal["slope_high"]),
    brier = unname(cal["brier"]), stringsAsFactors = FALSE, row.names = NULL
  )
}
run_subgroups <- function(d, cohort) {
  rows <- list(one_row(d, cohort, "overall", "all"))
  if ("female" %in% names(d) && any(!is.na(d$female))) {
    for (v in sort(unique(as.numeric(d$female[!is.na(d$female)])))) {
      lab <- if (v == 1) "female" else if (v == 0) "male" else as.character(v)
      rows[[length(rows) + 1]] <- one_row(d[!is.na(d$female) & as.numeric(d$female) == v, ], cohort, "sex", lab)
    }
  } else logmsg(cohort, ": sex unavailable\n")
  if ("edu_group" %in% names(d)) {
    for (lv in levels(factor(d$edu_group))) {
      dd <- d[!is.na(d$edu_group) & as.character(d$edu_group) == lv, ]
      if (nrow(dd) > 0) rows[[length(rows) + 1]] <- one_row(dd, cohort, "education", lv)
    }
  }
  if ("age" %in% names(d)) {
    ag1 <- cut(d$age, breaks = c(44, 59, 69, Inf), labels = c("45-59", "60-69", ">=70"), right = TRUE)
    for (lv in levels(ag1)) {
      dd <- d[!is.na(ag1) & ag1 == lv, ]
      if (nrow(dd) > 0) rows[[length(rows) + 1]] <- one_row(dd, cohort, "age_band", as.character(lv))
    }
    ag2 <- ifelse(d$age < 65, "<65", ">=65")
    for (lv in c("<65", ">=65")) {
      dd <- d[!is.na(d$age) & ag2 == lv, ]
      if (nrow(dd) > 0) rows[[length(rows) + 1]] <- one_row(dd, cohort, "age_65", lv)
    }
  }
  if ("rural" %in% names(d) && any(!is.na(d$rural))) {
    for (v in sort(unique(as.numeric(d$rural[!is.na(d$rural)])))) {
      lab <- if (v == 1) "rural" else if (v == 0) "urban" else as.character(v)
      rows[[length(rows) + 1]] <- one_row(d[!is.na(d$rural) & as.numeric(d$rural) == v, ], cohort, "rural_urban", lab)
    }
  } else logmsg(cohort, ": rural unavailable\n")
  bind_rows(rows)
}

logmsg("\n=== (1) CHARLS development ===\n")
print(table(main$female, useNA = "ifany"))
print(table(main$rural, useNA = "ifany"))
dev_tab <- run_subgroups(main, "CHARLS_dev")
m_row <- dev_tab[dev_tab$subgroup == "sex" & dev_tab$level == "male", ]
f_row <- dev_tab[dev_tab$subgroup == "sex" & dev_tab$level == "female", ]
logmsg(sprintf("R07 check male AUC=%.6f (expect ~0.710)\n", m_row$auc_model))
logmsg(sprintf("R07 check female AUC=%.6f (expect ~0.702)\n", f_row$auc_model))
write.csv(dev_tab, file.path(OUT, "F39_23a_CHARLS_dev_subgroups.csv"), row.names = FALSE, fileEncoding = "UTF-8")

logmsg("\n=== (2) External S0 ===\n")
ext_all <- list()
expected <- c(HRS = 0.696, ELSA = 0.763, SHARE = 0.751)
ext_paths <- c(
  HRS = "D:/BaiduSyncdisk/科研/徐中豪/external_data/HRS/analysis/HRS_analysis_2026-09-26.rds",
  ELSA = "D:/BaiduSyncdisk/科研/徐中豪/external_data/ELSA/analysis/ELSA_analysis_2026-09-26.rds",
  SHARE = "D:/BaiduSyncdisk/科研/徐中豪/external_data/SHARE/analysis/SHARE_analysis_2026-09-26.rds"
)
for (co in c("HRS", "ELSA", "SHARE")) {
  d <- readRDS(ext_paths[[co]])
  d$edu_group <- factor(as.character(d$edu_group),
    levels = c("Less than primary", "Primary school", "Middle school", "High school or above"))
  d$lp <- fixed_lp_manual(d, b)
  d$p_hat <- plogis(d$lp)
  if (!"score" %in% names(d)) d$score <- fixed_score(d)
  logmsg(co, " nrow=", nrow(d), " events=", sum(d$event),
         " ADL_NA=", sum(is.na(d$adl)),
         " female=", "female" %in% names(d), "\n")
  d_m <- d[!is.na(d$adl) & !is.na(d$p_hat), ]
  ov <- auc_ci(d_m$event, d_m$p_hat)
  logmsg(sprintf("  modelComplete n=%d events=%d AUC=%.6f (%.6f-%.6f) expect~%.3f\n",
                 nrow(d_m), sum(d_m$event), ov["auc"], ov["low"], ov["high"], expected[[co]]))
  if (abs(ov["auc"] - expected[[co]]) > 0.005)
    logmsg("  WARNING: external overall AUC differs from manuscript by >0.005\n")
  ext_all[[length(ext_all) + 1]] <- run_subgroups(d_m, paste0(co, "_S0_modelComplete"))
  ext_all[[length(ext_all) + 1]] <- run_subgroups(d, paste0(co, "_S0_scoreFull"))
  rm(d, d_m); invisible(gc(verbose = FALSE))
}
ext_tab <- bind_rows(ext_all)
write.csv(ext_tab, file.path(OUT, "F39_23a_external_S0_subgroups.csv"), row.names = FALSE, fileEncoding = "UTF-8")

logmsg("\n=== (3) Missingness ===\n")
miss_one <- function(df, stage, vars) {
  n <- nrow(df)
  bind_rows(lapply(vars, function(v) {
    if (!v %in% names(df)) {
      data.frame(stage = stage, variable = v, n_total = n,
                 n_missing = NA_integer_, pct_missing = NA_real_, note = "variable absent",
                 stringsAsFactors = FALSE)
    } else {
      nm <- sum(is.na(df[[v]]))
      data.frame(stage = stage, variable = v, n_total = n,
                 n_missing = nm, pct_missing = 100 * nm / max(n, 1), note = "",
                 stringsAsFactors = FALSE)
    }
  }))
}
vars_base <- c("age", "female", "edu_group", "cesd10", "dep", "adl", "mem2011", "mem2015", "rural")
s1 <- base
s2 <- base %>% filter(!is.na(mem2011))
s3 <- base %>% filter(!is.na(mem2011), mem2011 > 4)
s4 <- pool %>% filter(!is.na(mem2015))
s4$miss_any4 <- is.na(s4$age) | is.na(s4$edu_group) | is.na(s4$dep) | is.na(s4$adl)
logmsg("s4 n=", nrow(s4), " any4miss=", sum(s4$miss_any4), " complete4=", sum(!s4$miss_any4), "\n")
logmsg("s4 miss age/edu/dep/adl/cesd10/female=",
       sum(is.na(s4$age)), "/", sum(is.na(s4$edu_group)), "/", sum(is.na(s4$dep)), "/",
       sum(is.na(s4$adl)), "/", sum(is.na(s4$cesd10)), "/", sum(is.na(s4$female)), "\n")

miss_tab <- bind_rows(
  miss_one(s1, "1_2011_age>=45", vars_base),
  miss_one(s2, "2_mem2011_measurable", vars_base),
  miss_one(s3, "3_mem2011>4", vars_base),
  miss_one(s4, "4_mem2015_measurable", vars_base)
)
write.csv(miss_tab, file.path(OUT, "F39_20b_CHARLS_missing_by_var.csv"), row.names = FALSE, fileEncoding = "UTF-8")

ext_miss <- list()
for (co in c("HRS", "ELSA", "SHARE")) {
  d <- readRDS(ext_paths[[co]])
  vars_e <- c("age", "edu_group", "dep", "adl", "mem_base", "mem_fu", "event", "score", "female", "sex", "rural")
  ext_miss[[co]] <- miss_one(d, paste0(co, "_S0_riskset"), vars_e)
  rm(d)
}
ext_miss_tab <- bind_rows(ext_miss)
write.csv(ext_miss_tab, file.path(OUT, "F39_20b_external_missing_by_var.csv"), row.names = FALSE, fileEncoding = "UTF-8")

all_perf <- bind_rows(dev_tab, ext_tab)
write.csv(all_perf, file.path(OUT, "F39_23a_all_subgroups.csv"), row.names = FALSE, fileEncoding = "UTF-8")

# Also write xlsx if openxlsx available
if (requireNamespace("openxlsx", quietly = TRUE)) {
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, "CHARLS_dev"); openxlsx::writeData(wb, "CHARLS_dev", dev_tab)
  openxlsx::addWorksheet(wb, "external_S0"); openxlsx::writeData(wb, "external_S0", ext_tab)
  openxlsx::addWorksheet(wb, "CHARLS_missing"); openxlsx::writeData(wb, "CHARLS_missing", miss_tab)
  openxlsx::addWorksheet(wb, "external_missing"); openxlsx::writeData(wb, "external_missing", ext_miss_tab)
  openxlsx::saveWorkbook(wb, file.path(OUT, "F39_TRIPOD_tables.xlsx"), overwrite = TRUE)
  logmsg("wrote xlsx\n")
} else logmsg("openxlsx not installed; skipped xlsx\n")

fmt3 <- function(x) ifelse(is.na(x), "NA", sprintf("%.3f", x))
fmt4 <- function(x) ifelse(is.na(x), "NA", sprintf("%.4f", x))
fmt_ci <- function(e, lo, hi) ifelse(is.na(e), "NA", sprintf("%.3f (%.3f\u2013%.3f)", e, lo, hi))

md <- c(
  "# F39 TRIPOD \u8865\u5145\u5206\u6790\u62a5\u544a\uff082026-10-17\uff09",
  "",
  "> \u53ea\u8bfb\u65e2\u6709\u6570\u636e\u4e0e\u51bb\u7ed3\u6a21\u578b\uff1b\u672a\u6539\u52a8 code/output\u3002\u8f93\u51fa\u4ec5\u5728\u672c\u76ee\u5f55\u3002",
  "> \u51bb\u7ed3\u7cfb\u6570\uff1aA14_\u6700\u7ec8\u6a21\u578b\u7cfb\u6570.csv\uff08\u4e0e A07 model_u \u4e00\u81f4\uff09\u3002",
  "> \u5f00\u53d1\u6837\u672c\uff1awaves \u6309 A15 \u7b2c1\u8282\u903b\u8f91\u91cd\u5efa\uff0cKEY/score \u4e0e A07 \u5bf9\u9f50\u3002",
  "> \u79cd\u5b50\uff1a20261017\u3002",
  "",
  "## 1. \u590d\u73b0\u6838\u5bf9",
  "",
  sprintf("- \u5f00\u53d1\u6837\u672c n=%d\uff0c\u4e8b\u4ef6=%d\uff08\u8981\u6c42 4441 / 658\uff09", nrow(main), sum(main$event)),
  sprintf("- \u603b\u4f53\u8fde\u7eed\u6a21\u578b AUC=%s\uff08\u7a3f\u7ea6\u8868\u89c2 0.702\uff09", fmt_ci(ov_m["auc"], ov_m["low"], ov_m["high"])),
  sprintf("- \u603b\u4f53\u6574\u6570\u8bc4\u5206 AUC=%s\uff08\u7a3f\u7ea6\u8868\u89c2 0.700\uff09", fmt_ci(ov_s["auc"], ov_s["low"], ov_s["high"])),
  sprintf("- R07 \u7537\u6027 AUC=%s\uff08\u671f\u671b \u22480.710\uff09", fmt_ci(m_row$auc_model, m_row$auc_model_low, m_row$auc_model_high)),
  sprintf("- R07 \u5973\u6027 AUC=%s\uff08\u671f\u671b \u22480.702\uff09", fmt_ci(f_row$auc_model, f_row$auc_model_low, f_row$auc_model_high)),
  "",
  "## 2. TRIPOD 23a\uff1aCHARLS \u5f00\u53d1\u6837\u672c\u4e9a\u7ec4",
  "",
  "| \u4e9a\u7ec4 | \u6c34\u5e73 | n | \u4e8b\u4ef6 | \u8fde\u7eedAUC (95%CI) | \u6574\u6570AUC (95%CI) | \u6821\u51c6\u622a\u8ddd | \u6821\u51c6\u659c\u7387 | Brier |",
  "|---|---|---:|---:|---|---|---|---|---:|"
)
for (i in seq_len(nrow(dev_tab))) {
  r <- dev_tab[i, ]
  md <- c(md, sprintf("| %s | %s | %d | %d | %s | %s | %s | %s | %s |",
    r$subgroup, r$level, r$n, r$events,
    fmt_ci(r$auc_model, r$auc_model_low, r$auc_model_high),
    fmt_ci(r$auc_score, r$auc_score_low, r$auc_score_high),
    fmt_ci(r$cal_intercept, r$cal_intercept_low, r$cal_intercept_high),
    fmt_ci(r$cal_slope, r$cal_slope_low, r$cal_slope_high),
    fmt4(r$brier)))
}
md <- c(md, "",
  "> \u8bf4\u660e\uff1a\u6559\u80b2\u5c42\u5185\u533a\u5206\u5ea6\u4e3b\u8981\u6765\u81ea\u5c42\u5185\u5e74\u9f84\u4e0e\u6291\u90c1\u7684\u6b8b\u5dee\u53d8\u5f02\uff1b\u5e74\u9f84\u662f\u8bc4\u5206\u6784\u6210\u9879\uff0c\u5e74\u9f84\u7ec4\u5185 AUC \u901a\u5e38\u4f4e\u4e8e\u603b\u4f53\u3002",
  "",
  "## 3. TRIPOD 23a\uff1a\u5916\u90e8 S0",
  "",
  "- \u5916\u90e8 analysis RDS **\u4e0d\u542b sex/female \u6216 rural**\uff0c\u5916\u90e8\u6027\u522b/\u57ce\u4e61\u4e9a\u7ec4\u4e0d\u53ef\u7528\uff08\u672a\u91cd\u5efa\u534f\u8c03\u5316\uff09\u3002",
  "- \u4e0b\u8868\uff1aADL \u5b8c\u6574\u75c5\u4f8b\uff08\u8fde\u7eed\u6a21\u578b\u53e3\u5f84\uff09\u3002\u6574\u6570\u5206\u5168\u98ce\u9669\u96c6\u89c1 CSV `*_scoreFull`\u3002",
  "",
  "| \u961f\u5217 | \u4e9a\u7ec4 | \u6c34\u5e73 | n | \u4e8b\u4ef6 | \u8fde\u7eedAUC | \u6574\u6570AUC | \u622a\u8ddd | \u659c\u7387 | Brier |",
  "|---|---|---|---:|---:|---|---|---|---|---:|"
)
ext_show <- ext_tab[grepl("_modelComplete$", ext_tab$cohort), ]
for (i in seq_len(nrow(ext_show))) {
  r <- ext_show[i, ]
  md <- c(md, sprintf("| %s | %s | %s | %d | %d | %s | %s | %s | %s | %s |",
    sub("_S0_modelComplete", "", r$cohort), r$subgroup, r$level, r$n, r$events,
    fmt_ci(r$auc_model, r$auc_model_low, r$auc_model_high),
    fmt_ci(r$auc_score, r$auc_score_low, r$auc_score_high),
    fmt3(r$cal_intercept), fmt3(r$cal_slope), fmt4(r$brier)))
}
md <- c(md, "",
  "## 4. TRIPOD 20b\uff1aCHARLS \u9010\u53d8\u91cf\u7f3a\u5931",
  "",
  "| \u9636\u6bb5 | \u53d8\u91cf | n | \u7f3a\u5931n | \u7f3a\u5931% |",
  "|---|---|---:|---:|---:|"
)
for (i in seq_len(nrow(miss_tab))) {
  r <- miss_tab[i, ]
  md <- c(md, sprintf("| %s | %s | %d | %s | %s |",
    r$stage, r$variable, r$n_total,
    ifelse(is.na(r$n_missing), "\u2014", as.character(r$n_missing)),
    ifelse(is.na(r$pct_missing), "\u2014", sprintf("%.2f", r$pct_missing))))
}
md <- c(md, "",
  sprintf("- \u6838\u5bf9\uff1a2015 \u8bb0\u5fc6\u53ef\u6d4b n=%d\uff1b\u5e74\u9f84/\u6559\u80b2/\u6291\u90c1/ADL \u4efb\u4e00\u7f3a\u5931 %d \u4eba \u2192 \u56db\u53d8\u91cf\u5b8c\u6574 %d\uff08\u8981\u6c42 8296\uff09\u3002",
          nrow(s4), sum(s4$miss_any4), sum(!s4$miss_any4)),
  "",
  "## 5. \u5916\u90e8\u7f3a\u5931\uff08S0 \u98ce\u9669\u96c6\uff09",
  "",
  "| \u961f\u5217\u9636\u6bb5 | \u53d8\u91cf | n | \u7f3a\u5931n | \u7f3a\u5931% |",
  "|---|---|---:|---:|---:|"
)
for (i in seq_len(nrow(ext_miss_tab))) {
  r <- ext_miss_tab[i, ]
  md <- c(md, sprintf("| %s | %s | %d | %s | %s |",
    r$stage, r$variable, r$n_total,
    ifelse(is.na(r$n_missing), "\u53d8\u91cf\u4e0d\u5b58\u5728", as.character(r$n_missing)),
    ifelse(is.na(r$pct_missing), "\u2014", sprintf("%.2f", r$pct_missing))))
}
md <- c(md, "",
  "## 6. \u5c40\u9650",
  "",
  "1. \u4e9a\u7ec4\u6821\u51c6\u4e0e AUC \u5747\u4e3a\u8868\u89c2\u4f30\u8ba1\uff0c\u672a\u505a\u7ec4\u5185\u4e50\u89c2\u5ea6\u6821\u6b63\u3002",
  "2. \u6559\u80b2/\u5e74\u9f84\u4e9a\u7ec4\u7684\u533a\u5206\u5ea6\u53d7\u8bc4\u5206\u6784\u6210\u9879\u538b\u7f29\uff0c\u5c42\u5185 AUC \u504f\u4f4e\u5c5e\u9884\u671f\u73b0\u8c61\u3002",
  "3. \u5916\u90e8\u534f\u8c03\u5316\u6587\u4ef6\u65e0\u6027\u522b\u3001\u57ce\u4e61\uff0c\u672a\u91cd\u5efa\u6570\u636e\u3002",
  "4. \u6574\u6570\u8bc4\u5206\u65e0\u6982\u7387\u6620\u5c04\uff0c\u6821\u51c6\u6307\u6807\u6765\u81ea\u540c\u4e00\u6279\u4eba\u7684\u8fde\u7eed\u6a21\u578b\u9884\u6d4b\u6982\u7387\uff08\u82e5 ADL \u5b8c\u6574\uff09\u3002",
  ""
)
writeLines(md, file.path(OUT, "F39_TRIPOD_report_zh.md"), useBytes = TRUE)

# Copy script into OUT for provenance
# script copied separately by runner

logmsg("\n=== sessionInfo ===\n")
si <- capture.output(print(sessionInfo()))
writeLines(si, log_con); writeLines(si, stdout())
logmsg("\ndone: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
close(log_con)
