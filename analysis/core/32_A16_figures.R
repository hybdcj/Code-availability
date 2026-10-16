# ============================================================================
# 32_A16_figures.R —— A16：Figure 1–5 中英文投稿版（PNG 300 dpi + TIF 600 dpi）
# 输入：A07/A08/A09/A12/A13/A14/A15 冻结结果
# 输出：output/14_A16_图表/{中文,English}/Figure*_ZH/EN.png/.tif
# 生成：2026-09-27
# ============================================================================
invisible(try(Sys.setlocale("LC_ALL", "zh_CN.UTF-8"), silent = TRUE))
.libPaths(c("D:/R/library", .libPaths()))
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(pROC); library(readr) })
P <- "D:/BaiduSyncdisk/科研/徐中豪"
OUT <- file.path(P, "output", "14_A16_图表")
ZH <- file.path(OUT, "中文"); EN <- file.path(OUT, "English"); DATA <- file.path(OUT, "data")
for (d in c(OUT, ZH, EN, DATA)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

save_both <- function(draw_fun, basename, lang) {
  dir <- if (lang == "ZH") ZH else EN
  for (ext in c("png", "tif")) {
    path <- file.path(dir, paste0(basename, ".", ext))
    ok <- FALSE
    for (i in seq_len(10)) {
      ok <- tryCatch({
        if (ext == "png") {
          png(path, width = 12, height = 7, units = "in", res = 300)
        } else {
          tiff(path, width = 12, height = 7, units = "in", res = 600, compression = "lzw")
        }
        draw_fun(lang)
        dev.off()
        TRUE
      }, error = function(e) {
        if (length(dev.list())) dev.off()
        FALSE
      })
      if (ok) break
      Sys.sleep(1)
    }
    if (!ok) stop("无法写入图件（重试10次）: ", path)
  }
}
safe_write_csv <- function(x, file, ...) {
  for (i in seq_len(10)) {
    ok <- tryCatch({ utils::write.csv(x, file, ...); TRUE }, error = function(e) FALSE)
    if (ok) return(invisible(TRUE))
    Sys.sleep(1)
  }
  stop("无法写入文件（重试10次）: ", file)
}
parse_leading_num <- function(x) suppressWarnings(as.numeric(sub(" .*$", "", trimws(x))))
parse_est <- function(x) {
  x <- sub(", p=.*$", "", trimws(x))
  m <- regmatches(x, regexec("^([-0-9.]+) \\(([-0-9.]+)-([-0-9.]+)\\)$", x))[[1]]
  if (length(m) == 4) return(as.numeric(m[2:4]))
  c(NA_real_, NA_real_, NA_real_)
}

# ---------------------------------------------------------------------------
# Figure 1：研究流程图
# ---------------------------------------------------------------------------
draw_flow <- function(lang) {
  zh <- lang == "ZH"
  par(mar = c(1, 1, 1, 1), xpd = NA, family = "sans")
  plot(NA, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE, xlab = "", ylab = "")
  box <- function(x, y, w, h, label, col = "white", border = "#2C3E50", lwd = 1.5, cex = 0.85) {
    rect(x - w/2, y - h/2, x + w/2, y + h/2, col = col, border = border, lwd = lwd)
    text(x, y, label, cex = cex, font = 1)
  }
  ar <- function(x1, y1, x2, y2, label = NULL, col = "#2C3E50") {
    arrows(x1, y1, x2, y2, length = 0.08, lwd = 1.5, col = col)
    if (!is.null(label)) text((x1+x2)/2 + 0.01, (y1+y2)/2, label, cex = 0.72, col = col, srt = 0)
  }
  title_zh <- "CHARLS 开发队列与外部验证流程"
  title_en <- "Development and external validation flow"
  text(0.5, 0.98, if (zh) title_zh else title_en, cex = 1.3, font = 2)
  # development column
  text(0.25, 0.925, if (zh) "CHARLS 开发队列" else "CHARLS development cohort", cex = 1.0, font = 2)
  boxes <- list(
    c(0.25, 0.845, if (zh) "2011 年年龄≥45 岁\nn=16,831" else "2011 age ≥45 years\nn=16,831"),
    c(0.25, 0.725, if (zh) "2011 年记忆分可测\nn=13,444" else "2011 memory measured\nn=13,444"),
    c(0.25, 0.605, if (zh) "基线未受损（记忆分>4）\nn=10,422" else "Baseline unimpaired (memory >4)\nn=10,422"),
    c(0.25, 0.485, if (zh) "2015 年记忆分可测\nn=8,371" else "2015 memory measured\nn=8,371"),
    c(0.25, 0.365, if (zh) "四变量完整\nn=8,296" else "Four variables complete\nn=8,296"),
    c(0.25, 0.225, if (zh) "最终开发样本（五波 ms10 固定队列）\nn=4,441；事件=658（14.8%）" else "Final development sample (five-wave ms10 fixed cohort)\nn=4,441; events=658 (14.8%)")
  )
  for (i in seq_along(boxes)) {
    b <- boxes[[i]]
    box(as.numeric(b[1]), as.numeric(b[2]), 0.42, 0.095,
        b[3], col = if (i == length(boxes)) "#D6EAF8" else "white", cex = 0.82)
    if (i < length(boxes)) ar(as.numeric(b[1]), as.numeric(b[2]) - 0.048,
                              as.numeric(boxes[[i+1]][1]), as.numeric(boxes[[i+1]][2]) + 0.048)
  }
  # external validation
  text(0.76, 0.925, if (zh) "外部验证队列（S0）" else "External validation cohorts (S0)", cex = 1.0, font = 2)
  ext <- list(
    c(0.76, 0.82, if (zh) "HRS w10→w12\nn=14,051；事件=1,842（13.1%）" else "HRS w10→w12\nn=14,051; events=1,842 (13.1%)"),
    c(0.76, 0.68, if (zh) "ELSA w4→w6\nn=6,725；事件=554（8.2%）" else "ELSA w4→w6\nn=6,725; events=554 (8.2%)"),
    c(0.76, 0.54, if (zh) "SHARE w4→w6\nn=26,879；事件=2,930（10.9%）" else "SHARE w4→w6\nn=26,879; events=2,930 (10.9%)")
  )
  for (b in ext) box(as.numeric(b[1]), as.numeric(b[2]), 0.42, 0.105, b[3], col = "#FDEBD0", cex = 0.82)
  ar(0.47, 0.225, 0.55, 0.70, if (zh) "冻结 A14 模型\n不重新拟合" else "Frozen A14 model\nno refitting")
  box(0.76, 0.36, 0.42, 0.13,
      if (zh) "S0：各队列基线全样本 mean−1SD\nS1/S2：敏感性分析\n外部 AUC 0.692–0.763" else
        "S0: baseline-all mean−1SD\nS1/S2: sensitivity analyses\nExternal AUC 0.692–0.763",
      col = "#D5F5E3", cex = 0.80)
  text(0.25, 0.135, if (zh) "固定五波 ms10 完整限制排除 3,855 人" else "3,855 excluded by five-wave ms10 completeness", cex = 0.72, col = "#7F8C8D")
  text(0.5, 0.055,
       if (zh) "注：CHARLS 结局为 2011→2015 记忆分跨阈值（≤4）；外部队列为 4 年随访跨队列特异性阈值。ms10（定向 + 连续减 7）仅用于定义五波固定队列。" else
         "Note: CHARLS outcome is 2011→2015 memory cross-threshold (≤4); external cohorts use 4-year follow-up and cohort-specific thresholds. ms10 (orientation + serial 7s) was used only to define the five-wave fixed cohort.",
       cex = 0.66, adj = 0.5)
}
save_both(draw_flow, "Figure1_study_flow_ZH", "ZH")
save_both(draw_flow, "Figure1_study_flow_EN", "EN")

# ---------------------------------------------------------------------------
# Figure 2：候选因素稳定性
# ---------------------------------------------------------------------------
orm <- read.csv(file.path(P, "output/03_五波横断面/五波分波Logistic_OR矩阵_2026-09-19.csv"),
                check.names = FALSE, fileEncoding = "UTF-8")
risk <- read.csv(file.path(P, "output/04_纵向队列/cohort_risk_2026-09-13.csv"),
                 check.names = FALSE, fileEncoding = "UTF-8")
iadl2 <- read.csv(file.path(P, "output/04_纵向队列/iadl2_longitudinal_models_2026-09-26.csv"),
                  check.names = FALSE, fileEncoding = "UTF-8")
or_stats <- function(pattern) {
  idx <- grep(pattern, orm[[1]])
  if (length(idx) == 0) return(c(NA, NA, NA))
  vals <- unlist(orm[idx[1], 3:7], use.names = FALSE)
  nums <- parse_leading_num(vals); nums <- nums[!is.na(nums)]
  if (length(nums) == 0) return(c(NA, NA, NA))
  c(min(nums), median(nums), max(nums))
}
risk_est <- function(rowname) {
  idx <- which(risk[[1]] == rowname)
  if (length(idx) == 0) return(c(NA, NA, NA))
  parse_est(risk[idx[1], 2])
}
risk_hr <- function(rowname) {
  idx <- which(risk[[1]] == rowname)
  if (length(idx) == 0) return(c(NA, NA, NA))
  parse_est(risk[idx[1], 3])
}
iadl_est <- function(model) {
  x <- iadl2$unweighted[iadl2$exposure == "iadl2" & iadl2$model == model]
  if (length(x) == 0) return(c(NA, NA, NA))
  parse_est(x[1])
}
plot_candidate_stability <- function(lang) {
  zh <- lang == "ZH"
  par(mfrow = c(1, 2), mar = c(4.5, 6, 3, 2), family = "sans")
  rows <- list(
    list(name_zh = "抑郁", name_en = "Depression", role = "Core", or = or_stats("抑郁"), lr = risk_est("抑郁（CES-D-10 ≥10）"), hr = risk_hr("抑郁（CES-D-10 ≥10）")),
    list(name_zh = "ADL 受限", name_en = "ADL limitation", role = "Core", or = or_stats("ADL 受限"), lr = risk_est("ADL 受限"), hr = risk_hr("ADL 受限")),
    list(name_zh = "IADL 受限", name_en = "IADL limitation", role = "Supplementary", or = or_stats("IADL 受限"), lr = iadl_est("logistic"), hr = iadl_est("cox")),
    list(name_zh = "睡眠 <6 h", name_en = "Sleep <6 h", role = "Incremental", or = or_stats("睡眠 <6"), lr = risk_est("睡眠时长（每 1 小时）"), hr = risk_hr("睡眠时长（每 1 小时）")),
    list(name_zh = "睡眠 >8 h", name_en = "Sleep >8 h", role = "Incremental", or = or_stats("睡眠 >8"), lr = c(NA, NA, NA), hr = c(NA, NA, NA)),
    list(name_zh = "规律体力活动", name_en = "Physical activity", role = "Incremental", or = or_stats("规律体力活动"), lr = risk_est("规律体力活动"), hr = risk_hr("规律体力活动")),
    list(name_zh = "高血压", name_en = "Hypertension", role = "Background", or = or_stats("高血压"), lr = risk_est("高血压"), hr = risk_hr("高血压")),
    list(name_zh = "糖尿病", name_en = "Diabetes", role = "Incremental", or = or_stats("糖尿病"), lr = risk_est("糖尿病"), hr = risk_hr("糖尿病")),
    list(name_zh = "卒中", name_en = "Stroke", role = "Incremental", or = or_stats("卒中"), lr = c(NA, NA, NA), hr = c(NA, NA, NA)),
    list(name_zh = "疼痛", name_en = "Pain", role = "Excluded", or = or_stats("躯体疼痛"), lr = risk_est("疼痛"), hr = risk_hr("疼痛"))
  )
  cols <- c(Core = "#1F4E78", Supplementary = "#2E7D32", Incremental = "#E67E22", Background = "#7D3C98", Excluded = "#C0392B")
  # Panel A：cross-sectional OR range
  yy <- seq_along(rows)
  plot(NA, xlim = c(0.45, 2.6), ylim = c(0.5, length(rows) + 0.5), log = "x",
       yaxt = "n", xlab = if (zh) "OR（对数刻度）" else "OR (log scale)", ylab = "",
       main = if (zh) "(a) 五波横断面 OR 范围" else "(a) Cross-sectional OR range across five waves", bty = "n")
  abline(v = 1, lty = 2, col = "grey60")
  for (i in seq_along(rows)) {
    z <- rows[[i]]; if (all(is.na(z$or))) next
    segments(z$or[1], i, z$or[3], i, col = cols[z$role], lwd = 3)
    points(z$or[2], i, pch = 18, cex = 1.3, col = cols[z$role])
  }
  axis(2, at = yy, labels = if (zh) vapply(rows, `[[`, "", "name_zh") else vapply(rows, `[[`, "", "name_en"),
       las = 1, cex.axis = 0.75)
  legend("bottomright", if (zh) c("核心", "补充", "增量", "背景", "排除") else names(cols), col = cols, lwd = 3, pch = 18, bty = "n", cex = 0.7)
  # Panel B：longitudinal OR/HR forest (only factors with longitudinal estimates)
  rows_b <- rows[!vapply(rows, function(z) all(is.na(z$lr)) && all(is.na(z$hr)), logical(1))]
  lab <- if (zh) vapply(rows_b, `[[`, "", "name_zh") else vapply(rows_b, `[[`, "", "name_en")
  y2 <- rev(seq_along(rows_b))
  plot(NA, xlim = c(0.45, 2.6), ylim = c(0.5, length(rows_b) + 0.5), log = "x",
       yaxt = "n", xlab = if (zh) "OR / HR（对数刻度）" else "OR / HR (log scale)", ylab = "",
       main = if (zh) "(b) 纵向效应量" else "(b) Longitudinal effect estimates", bty = "n")
  abline(v = 1, lty = 2, col = "grey60")
  for (i in seq_along(rows_b)) {
    z <- rows_b[[i]]; y <- y2[i]
    if (!all(is.na(z$lr))) {
      segments(z$lr[2], y + 0.12, z$lr[3], y + 0.12, col = "#1F4E78", lwd = 2)
      points(z$lr[1], y + 0.12, pch = 19, col = "#1F4E78", cex = 1.1)
    }
    if (!all(is.na(z$hr))) {
      segments(z$hr[2], y - 0.12, z$hr[3], y - 0.12, col = "#C0392B", lwd = 2)
      points(z$hr[1], y - 0.12, pch = 17, col = "#C0392B", cex = 1.1)
    }
  }
  axis(2, at = y2, labels = lab, las = 1, cex.axis = 0.75)
  legend("bottomright", c(if (zh) "终点评定 OR" else "Endpoint OR", if (zh) "首次受损 HR" else "First impairment HR"),
         col = c("#1F4E78", "#C0392B"), pch = c(19, 17), lwd = 2, bty = "n", cex = 0.75)
}
save_both(plot_candidate_stability, "Figure2_candidate_stability_ZH", "ZH")
save_both(plot_candidate_stability, "Figure2_candidate_stability_EN", "EN")

safe_write_csv(data.frame(pattern = c("抑郁", "ADL 受限", "IADL 受限", "睡眠 <6", "睡眠 >8", "规律体力活动", "高血压", "糖尿病", "卒中", "躯体疼痛")),
               file.path(DATA, "Figure2_candidate_patterns.csv"), row.names = FALSE, fileEncoding = "UTF-8")

# ---------------------------------------------------------------------------
# Figure 3：风险分层
# ---------------------------------------------------------------------------
a14_strata <- read.csv(file.path(P, "output/12_模型定稿/A14_风险分层.csv"),
                       check.names = FALSE, fileEncoding = "UTF-8")
a09_strata <- read.csv(file.path(P, "output/10_筛查指数/A09_strata_observed_vs_predicted_2026-09-26.csv"),
                       check.names = FALSE, fileEncoding = "UTF-8")
draw_risk_strata <- function(lang) {
  zh <- lang == "ZH"
  par(mar = c(5, 5, 3, 2), family = "sans", cex.lab = 1.2)
  bp <- barplot(a09_strata$观测率, names.arg = if (zh) c("低危", "中危", "高危") else c("Low", "Intermediate", "High"),
                col = c("#C6DBEF", "#6BAED6", "#2171B5"), ylim = c(0, max(a09_strata$观测率) * 1.35),
                ylab = if (zh) "4 年新发认知受损率" else "4-year incident cognitive impairment",
                xlab = if (zh) "筛查指数风险层" else "Screening index risk stratum",
                main = if (zh) "筛查指数风险分层（CHARLS 开发队列）" else "Screening index risk stratification (CHARLS development cohort)")
  text(bp, a09_strata$观测率, sprintf("%.1f%%\n(n=%s)", 100*a09_strata$观测率, format(a09_strata$n, big.mark=",")), pos = 3, cex = 0.95)
  points(bp, a09_strata$预测率, pch = 18, cex = 1.6, col = "#C0392B")
  legend("topleft", if (zh) c("观测事件率", "模型预测率") else c("Observed event rate", "Model-predicted risk"),
         pch = c(15, 18), col = c("#2171B5", "#C0392B"), bty = "n", cex = 0.9)
}
save_both(draw_risk_strata, "Figure3_risk_stratification_ZH", "ZH")
save_both(draw_risk_strata, "Figure3_risk_stratification_EN", "EN")

# ---------------------------------------------------------------------------
# Figure 4：内部验证 ROC / 校准 / DCA
# ---------------------------------------------------------------------------
a07 <- readRDS(file.path(P, "output/10_筛查指数/screening_index_memory_v3_2026-09-26.rds"))
a09_cal <- read.csv(file.path(P, "output/10_筛查指数/A09_calibration_deciles_2026-09-26.csv"),
                    check.names = FALSE, fileEncoding = "UTF-8")
a09_dca <- read.csv(file.path(P, "output/10_筛查指数/A09_dca_2026-09-26.csv"),
                    check.names = FALSE, fileEncoding = "UTF-8")
d <- a07$data; m <- a07$model_u
d$p_hat <- predict(m, type = "response")
r <- roc(d$incident_2015, d$p_hat, direction = "<", quiet = TRUE)
auc_v <- as.numeric(auc(r)); ci_v <- as.numeric(ci.auc(r))
cal_slope <- glm(incident_2015 ~ qlogis(pmin(pmax(p_hat, 1e-6), 1 - 1e-6)), data = d, family = binomial())$coefficients[[2]]
cal_int <- glm(incident_2015 ~ offset(qlogis(pmin(pmax(p_hat, 1e-6), 1 - 1e-6))), data = d, family = binomial())$coefficients[[1]]
draw_internal <- function(lang) {
  zh <- lang == "ZH"
  par(mfrow = c(1, 3), mar = c(5, 5, 3, 2), family = "sans", cex.lab = 1.0, cex.axis = 0.9)
  plot(1 - r$specificities, r$sensitivities, type = "l", lwd = 2.5, col = "#2C7FB8",
       xlab = if (zh) "1 − 特异度" else "1 − Specificity", ylab = if (zh) "灵敏度" else "Sensitivity",
       main = if (zh) "(a) ROC 曲线" else "(a) ROC curve", xlim = c(0, 1), ylim = c(0, 1))
  abline(0, 1, lty = 3, col = "grey60")
  text(0.60, 0.18, sprintf("AUC = %.3f\n(95%% CI %.3f–%.3f)", auc_v, ci_v[1], ci_v[3]), cex = 0.9)
  plot(a09_cal$预测均值, a09_cal$观测率, pch = 19, cex = 1.2, col = "#31A354",
       xlim = c(0, max(a09_cal$观测CI高) * 1.05), ylim = c(0, max(a09_cal$观测CI高) * 1.05),
       xlab = if (zh) "预测概率" else "Predicted probability", ylab = if (zh) "观测事件率" else "Observed event rate",
       main = if (zh) "(b) 校准图（十分位）" else "(b) Calibration plot (deciles)")
  segments(a09_cal$预测均值, a09_cal$观测CI低, a09_cal$预测均值, a09_cal$观测CI高, col = "#31A354")
  abline(0, 1, lty = 3, col = "grey60")
  lines(lowess(a09_cal$预测均值, a09_cal$观测率), col = "#756BB1", lwd = 1.5)
  text(max(a09_cal$预测均值) * 0.42, max(a09_cal$观测CI高) * 0.95,
       sprintf(if (zh) "校准斜率 = %.2f\n截距 = %.2f" else "Calibration slope = %.2f\nIntercept = %.2f", cal_slope, cal_int), cex = 0.9, adj = 0)
  plot(a09_dca$阈值, a09_dca$净获益_模型, type = "l", lwd = 2.5, col = "#2C7FB8",
       ylim = c(min(a09_dca$净获益_模型, -0.02), max(a09_dca$净获益_模型, mean(d$incident_2015)) * 1.05),
       xlab = if (zh) "阈值概率" else "Threshold probability", ylab = if (zh) "净获益" else "Net benefit",
       main = if (zh) "(c) 决策曲线分析" else "(c) Decision curve analysis")
  lines(a09_dca$阈值, a09_dca$净获益_全部干预, lty = 2, col = "#E6550D", lwd = 2)
  abline(h = 0, lty = 3, col = "grey60")
  legend("topright", if (zh) c("筛查指数", "全部干预", "不干预") else c("Screening index", "Treat all", "Treat none"),
         lty = c(1, 2, 3), lwd = c(2.5, 2, 1), col = c("#2C7FB8", "#E6550D", "grey60"), bty = "n", cex = 0.8)
}
save_both(draw_internal, "Figure4_internal_validation_ZH", "ZH")
save_both(draw_internal, "Figure4_internal_validation_EN", "EN")

# ---------------------------------------------------------------------------
# Figure 5：外部验证 ROC / 校准 / AUC forest
# ---------------------------------------------------------------------------
a12 <- readRDS(file.path(P, "output/11_外部验证/A12_results_2026-09-26.rds"))
model_u <- a07$model_u; model_levels <- model_u$xlevels$edu_group
cohorts <- c("HRS", "ELSA", "SHARE")
cohort_cols <- c(HRS = "#1F77B4", ELSA = "#2CA02C", SHARE = "#D62728")
risk_sets <- list()
for (co in cohorts) {
  x <- readRDS(file.path(P, "external_data", co, "analysis",
                         paste0(co, "_sensitivity_risk_sets_2026-09-26.rds")))[[paste0(co, "__baseline_all_mean_sd")]]
  x$edu_group <- factor(as.character(x$edu_group), levels = model_levels)
  x$p_hat <- as.numeric(predict(model_u, newdata = data.frame(age = x$age, edu_group = x$edu_group,
                                                              dep = x$dep, adl = x$adl), type = "response"))
  risk_sets[[co]] <- x
}
draw_external <- function(lang) {
  zh <- lang == "ZH"
  par(mfrow = c(1, 3), mar = c(4.8, 4.8, 3, 1.5), family = "sans", cex.lab = 1.0, cex.axis = 0.85)
  # (a) ROC
  plot(NA, xlim = c(0, 1), ylim = c(0, 1), xlab = if (zh) "1 − 特异度" else "1 − Specificity",
       ylab = if (zh) "灵敏度" else "Sensitivity", main = if (zh) "(a) 外部队列 ROC（S0）" else "(a) External ROC (S0 model)", bty = "n")
  abline(0, 1, lty = 2, col = "grey70")
  for (co in cohorts) {
    d <- risk_sets[[co]]; d <- d[!is.na(d$p_hat), ]
    r <- roc(d$event, d$p_hat, direction = "<", quiet = TRUE)
    lines(1 - r$specificities, r$sensitivities, col = cohort_cols[co], lwd = 2)
  }
  a <- a12$perf[a12$perf$cut_definition == "baseline_all_mean_sd", ]
  a <- a[match(cohorts, a$cohort), ]
  legend("bottomright", sprintf("%s AUC %.3f (%.3f–%.3f)", a$cohort, a$auc_model, a$auc_model_low, a$auc_model_high),
         col = cohort_cols[a$cohort], lwd = 2, bty = "n", cex = 0.70)
  # (b) Calibration
  dec <- a12$deciles[a12$deciles$cut_definition == "baseline_all_mean_sd", ]
  ymax <- max(c(dec$observed_high, dec$pred_mean), na.rm = TRUE) * 1.05
  plot(NA, xlim = c(0, ymax), ylim = c(0, ymax), xlab = if (zh) "预测风险" else "Predicted risk",
       ylab = if (zh) "观测风险" else "Observed risk", main = if (zh) "(b) 校准图（S0）" else "(b) Calibration by decile (S0 model)", bty = "n")
  abline(0, 1, lty = 2, col = "grey70")
  for (co in cohorts) {
    d <- dec[dec$cohort == co, ]
    arrows(d$pred_mean, d$observed_low, d$pred_mean, d$observed_high, angle = 90, code = 3, length = 0.03, col = cohort_cols[co])
    lines(d$pred_mean, d$observed, type = "b", pch = 19, col = cohort_cols[co], lwd = 2)
  }
  legend("bottomright", cohorts, col = cohort_cols[cohorts], lwd = 2, pch = 19, bty = "n", cex = 0.75)
  # (c) AUC forest
  m <- a12$meta[a12$meta$cut_definition == "baseline_all_mean_sd" & a12$meta$metric == "model", ]
  labels <- c(cohorts, if (zh) "随机效应合并" else "Pooled RE")
  est <- c(a$auc_model, m$pooled_auc); lo <- c(a$auc_model_low, m$pooled_low); hi <- c(a$auc_model_high, m$pooled_high)
  y <- rev(seq_along(labels))
  plot(NA, xlim = c(0.62, 0.86), ylim = c(0.5, length(labels) + 0.5), yaxt = "n",
       xlab = if (zh) "AUC（95% CI）" else "AUC (95% CI)", ylab = "", main = if (zh) "(c) AUC 森林图（S0）" else "(c) AUC forest (S0 model)", bty = "n")
  abline(v = 0.5, lty = 3, col = "grey70")
  segments(lo, y, hi, y, col = c(cohort_cols[a$cohort], "black"), lwd = 2)
  points(est, y, pch = c(19, 19, 19, 18), col = c(cohort_cols[a$cohort], "black"), cex = c(1, 1, 1, 1.4))
  axis(2, at = y, labels = labels, las = 1, cex.axis = 0.8)
}
save_both(draw_external, "Figure5_external_validation_ZH", "ZH")
save_both(draw_external, "Figure5_external_validation_EN", "EN")

a_s0_csv <- a12$perf[a12$perf$cut_definition == "baseline_all_mean_sd", ]
a_s0_csv <- a_s0_csv[match(cohorts, a_s0_csv$cohort), ]
safe_write_csv(data.frame(cohort = cohorts, n = a_s0_csv$n_total, events = a_s0_csv$events_total,
                          auc = a_s0_csv$auc_model),
               file.path(DATA, "Figure5_external_summary.csv"), row.names = FALSE, fileEncoding = "UTF-8")
cat("A16 figures written.\n")
