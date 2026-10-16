# ============================================================================
# 12_cumulative_incidence.R —— T05：固定队列累积发病曲线与风险人数表
# 版本：v1.0（2026-09-20）
# 输入：output/04_纵向队列/interval_incidence_events_2026-09-20.rds
# 输出：
#   output/04_纵向队列/cumulative_incidence_2026-09-20.csv
#   output/04_纵向队列/cumulative_incidence_risk_table_2026-09-20.csv
#   output/04_纵向队列/cumulative_incidence_km_sensitivity_2026-09-20.csv
#   output/04_纵向队列/cumulative_incidence_2026-09-20.md
#   output/09_图表/Figure4_cumulative_incidence_2026-09-20.png
#   output/09_图表/Figure4_cumulative_incidence_2026-09-20.tif
#
# 说明：
#   主曲线：按计划随访时点 2、4、7、9 年计算离散累积发病；
#   敏感性曲线：连续两波确认受损；
#   KM：事件时间未知，使用区间中点近似 1、3、5.5、8 年，仅作敏感性分析。
# ============================================================================

invisible(try(Sys.setlocale("LC_ALL", "zh_CN.UTF-8"), silent = TRUE))
if (dir.exists("D:/R/library")) .libPaths(c("D:/R/library", .libPaths()))
suppressPackageStartupMessages(library(survival))

P <- Sys.getenv("CHARLS_PROJ", "D:/BaiduSyncdisk/科研/徐中豪")
OUT <- file.path(P, "output", "04_纵向队列")
FIG_DIR <- file.path(P, "output", "09_图表")
if (!dir.exists(OUT)) dir.create(OUT, recursive = TRUE)
if (!dir.exists(FIG_DIR)) dir.create(FIG_DIR, recursive = TRUE)

obj <- readRDS(file.path(OUT, "interval_incidence_events_2026-09-20.rds"))
base <- obj$base
n_base <- nrow(base)
stopifnot(n_base > 0, all(c("first_index", "confirmed_index", "endpoint_2020") %in% names(base)))

times <- c(0, 2, 4, 7, 9)
interval_len <- c(2, 2, 3, 2)
midpoints <- c(1, 3, 5.5, 8)

binomial_ci <- function(events, n) {
  if (n <= 0) return(c(NA_real_, NA_real_))
  as.numeric(binom.test(events, n)$conf.int)
}
fmt_pct <- function(x) ifelse(is.na(x), "—", sprintf("%.1f%%", 100 * x))

# 主曲线：首次发作。
events_first <- vapply(1:4, function(i) sum(base$first_index == i, na.rm = TRUE), integer(1))
cum_first <- c(0L, cumsum(events_first))
ci_first <- t(vapply(cum_first, function(x) binomial_ci(x, n_base), numeric(2)))
risk_first <- n_base - cum_first
# t=9 的风险人数是末次随访前尚未发生首次事件者；之后统一为研究结束。
risk_first[5] <- n_base - cum_first[5]

# 敏感性曲线：连续两波确认。2018-2020 无下一波，无法确认新事件。
events_conf <- vapply(1:3, function(i) sum(base$confirmed_index == i, na.rm = TRUE), integer(1))
cum_conf <- c(0L, cumsum(events_conf), sum(events_conf))
ci_conf <- t(vapply(cum_conf, function(x) binomial_ci(x, n_base), numeric(2)))
risk_conf <- c(n_base, n_base - cumsum(events_conf), NA_integer_)
risk_conf_status <- c("可估计", "可估计", "可估计", "可估计", "不可确认")

curve_tab <- data.frame(
  time_years = rep(times, 2),
  outcome = rep(c("first_episode", "confirmed_first_episode"), each = length(times)),
  outcome_label = rep(c("九年内首次受损", "连续两波确认受损"), each = length(times)),
  cumulative_events = c(cum_first, cum_conf),
  cumulative_incidence = c(cum_first / n_base, cum_conf / n_base),
  ci_lower = c(ci_first[, 1], ci_conf[, 1]),
  ci_upper = c(ci_first[, 2], ci_conf[, 2]),
  event_free_at_time = c(risk_first, risk_conf),
  risk_status = c(rep("可估计", 5), risk_conf_status),
  stringsAsFactors = FALSE
)
curve_tab$cumulative_incidence_pct <- 100 * curve_tab$cumulative_incidence
curve_tab$ci_lower_pct <- 100 * curve_tab$ci_lower
curve_tab$ci_upper_pct <- 100 * curve_tab$ci_upper

# 风险人数表：只保留论文图需要的两行。
risk_tab <- data.frame(
  time_years = times,
  at_risk_first_episode = risk_first,
  cumulative_events_first_episode = cum_first,
  at_risk_confirmed = risk_conf,
  cumulative_events_confirmed = cum_conf,
  confirmed_risk_status = risk_conf_status,
  stringsAsFactors = FALSE
)

# KM 敏感性：事件时间放在区间中点。
km_time <- ifelse(is.na(base$first_index), 9, midpoints[base$first_index])
km_event <- as.integer(!is.na(base$first_index))
km_fit <- survfit(Surv(km_time, km_event) ~ 1)
km_s <- summary(km_fit, times = c(2, 4, 7, 9), extend = TRUE)
km_tab <- data.frame(
  time_years = km_s$time,
  n_risk = km_s$n.risk,
  n_event = km_s$n.event,
  cumulative_incidence_km = 1 - km_s$surv,
  ci_lower = 1 - km_s$upper,
  ci_upper = 1 - km_s$lower,
  stringsAsFactors = FALSE
)
km_tab$cumulative_incidence_pct <- 100 * km_tab$cumulative_incidence_km
km_tab$ci_lower_pct <- 100 * km_tab$ci_lower
km_tab$ci_upper_pct <- 100 * km_tab$ci_upper

# 输出数据。
write.csv(curve_tab, file.path(OUT, "cumulative_incidence_2026-09-20.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(risk_tab, file.path(OUT, "cumulative_incidence_risk_table_2026-09-20.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(km_tab, file.path(OUT, "cumulative_incidence_km_sensitivity_2026-09-20.csv"), row.names = FALSE, fileEncoding = "UTF-8")
saveRDS(list(curve = curve_tab, risk = risk_tab, km = km_tab, fit = km_fit),
        file.path(OUT, "cumulative_incidence_2026-09-20.rds"))

# 主图：英文标签，便于后续投稿。risk table 放在 x 轴下方。
draw_figure <- function() {
  layout(matrix(c(1, 2), ncol = 1), heights = c(4.1, 1.7))
  # Upper panel: cumulative incidence curves.
  par(mar = c(2.8, 8.0, 3.2, 1.2), family = "sans")
  ylim <- c(0, 50)
  plot(times, 100 * cum_first / n_base, type = "n", xaxt = "n", yaxt = "n",
       xlim = c(0, 9.4), ylim = ylim,
       xlab = "", ylab = "Cumulative incidence (%)",
       main = "Cumulative Incidence of Cognitive Impairment (Baseline Normal)")
  axis(1, at = times, labels = FALSE, tcl = -0.25)
  axis(2, at = seq(0, 50, 10), las = 1)
  grid(col = "grey88", lty = 3)
  # First episode: 95% CI and line.
  segments(times, 100 * ci_first[, 1], times, 100 * ci_first[, 2], col = "#1F4E78", lwd = 1.5)
  lines(times, 100 * cum_first / n_base, type = "s", lwd = 2, col = "#1F4E78")
  points(times, 100 * cum_first / n_base, pch = 16, cex = 0.9, col = "#1F4E78")
  # Confirmed episode: last time is carry-forward, not a new confirmable event.
  ctimes <- c(0, 2, 4, 7)
  segments(ctimes, 100 * ci_conf[1:4, 1], ctimes, 100 * ci_conf[1:4, 2], col = "#C00000", lwd = 1.5)
  lines(ctimes, 100 * cum_conf[1:4] / n_base, type = "s", lwd = 2, col = "#C00000")
  points(ctimes, 100 * cum_conf[1:4] / n_base, pch = 17, cex = 0.9, col = "#C00000")
  points(9, 100 * cum_conf[5] / n_base, pch = 1, cex = 1.35, col = "#C00000", lwd = 2)
  legend("topleft",
         legend = c("First impairment episode (any wave)",
                    "Confirmed in two consecutive waves",
                    "Not confirmable at 9 years"),
         col = c("#1F4E78", "#C00000", "#C00000"),
         pch = c(16, 17, 1), lwd = c(2, 2, NA), bty = "n", cex = 0.82)
  # Lower panel: risk set and cumulative events table.
  par(mar = c(1.0, 8.0, 0.2, 1.2), family = "sans")
  plot(NA, xlim = c(0, 9.4), ylim = c(0, 5.6), axes = FALSE, xlab = "", ylab = "")
  mtext("No. at risk", side = 2, line = 5.2, at = 5.35, las = 1, adj = 0, font = 2, cex = 0.82)
  mtext("First episode", side = 2, line = 5.2, at = 4.15, las = 1, adj = 0, col = "#1F4E78", font = 2, cex = 0.80)
  mtext("Confirmed", side = 2, line = 5.2, at = 3.00, las = 1, adj = 0, col = "#C00000", font = 2, cex = 0.80)
  mtext("First-event count", side = 2, line = 5.2, at = 1.85, las = 1, adj = 0, font = 2, cex = 0.80)
  text(times, rep(4.15, length(times)), risk_first, cex = 0.78, col = "#1F4E78")
  conf_risk_display <- c(as.character(risk_conf[1:4]), "NE")
  text(times, rep(3.00, length(times)), conf_risk_display, cex = 0.78, col = "#C00000")
  text(times, rep(1.85, length(times)), cum_first, cex = 0.78, col = "black")
  text(times, rep(0.25, length(times)), times, cex = 0.82, font = 2)
  mtext("Follow-up (years)", side = 1, line = 1.85, cex = 0.82)
}

png(file.path(FIG_DIR, "Figure4_cumulative_incidence_2026-09-20.png"), width = 10, height = 7.2, units = "in", res = 300)
draw_figure()
dev.off()

tiff_ok <- TRUE
tryCatch({
  tiff(file.path(FIG_DIR, "Figure4_cumulative_incidence_2026-09-20.tif"), width = 10, height = 7.2, units = "in", res = 600, compression = "lzw")
  draw_figure()
  dev.off()
}, error = function(e) {
  tiff_ok <<- FALSE
  if (length(dev.list())) dev.off()
})

# 中文报告。
fmt_int <- function(x) ifelse(is.na(x), "—", format(as.integer(x), big.mark = ",", scientific = FALSE, trim = TRUE))
report <- c(
  "# T05：固定队列累积发病曲线（2026-09-20）",
  "",
  "> 输入：T04 的 `interval_incidence_events_2026-09-20.rds`",
  sprintf("> 主结局：2011 基线正常者（n=%s）首次发生 `impaired_fixed=1`", format(n_base, big.mark = ",")),
  "> **注意：** 本曲线是 2011 年≥45岁且五波 ms10 均可计分固定队列的条件累积发病，不能外推为全 CHARLS 人群发病率。",
  "> 主曲线按计划随访时点 2、4、7、9 年展示离散累积发病；风险人数为该时点仍无对应事件的参与者。KM 使用事件时间区间中点进行一致性检查。",
  "",
  "## 一、主曲线数据",
  "",
  "| 年份 | 随访年 | 累积事件 | 累积发病 | 95%CI | 该时点仍无事件人数 |",
  "|---|---:|---:|---:|---:|---:|",
  paste0("| ", c("2011","2013","2015","2018","2020"), " | ", times, " | ", cum_first, " | ",
         fmt_pct(cum_first / n_base), " | ",
         ifelse(is.na(ci_first[,1]), "—", sprintf("%.1f–%.1f%%", 100*ci_first[,1], 100*ci_first[,2])),
         " | ", fmt_int(risk_first), " |"),
  "",
  "## 二、连续两波确认累积发病（敏感性）",
  "",
  "| 年份 | 随访年 | 累积确认事件 | 累积发病 | 95%CI | 该时点仍未确认人数 |",
  "|---|---:|---:|---:|---:|---:|",
  paste0("| ", c("2011","2013","2015","2018","2020"), " | ", times, " | ", cum_conf, " | ",
         fmt_pct(cum_conf / n_base), " | ",
         ifelse(is.na(ci_conf[,1]), "—", sprintf("%.1f–%.1f%%", 100*ci_conf[,1], 100*ci_conf[,2])),
         " | ", ifelse(is.na(risk_conf), "不可确认", fmt_int(risk_conf)), " |"),
  "",
  "> 2018–2020 区间无法确认新的连续两波事件，因为需要 2022 年下一波数据；该时点标为“不可确认”，不是 0 事件。",
  "> 连续确认事件的时间记为连续两波中第一波受损的时间，确认则依赖后一波结果。",
  "",
  "## 三、KM 一致性检查（固定队列无删失）",
  "",
  "KM 使用区间中点 1、3、5.5、8 年作为事件时间。由于本固定队列在 9 年内没有中途删失，KM 估计与离散累积发病数一致；该表用于确认计算一致性，不构成独立敏感性结论。",
  "",
  "| 随访年 | 风险人数 | KM 累积发病 | 95%CI |",
  "|---:|---:|---:|---:|",
  paste0("| ", km_tab$time_years, " | ", km_tab$n_risk, " | ", fmt_pct(km_tab$cumulative_incidence_km), " | ",
         sprintf("%.1f–%.1f%%", km_tab$ci_lower_pct, km_tab$ci_upper_pct), " |"),
  "",
  "## 四、结论",
  "",
  sprintf("- 首次发作累积发病：2 年 %s，4 年 %s，7 年 %s，9 年 %s。",
          fmt_pct(cum_first[2]/n_base), fmt_pct(cum_first[3]/n_base),
          fmt_pct(cum_first[4]/n_base), fmt_pct(cum_first[5]/n_base)),
  sprintf("- 连续两波确认累积发病：2 年 %s，4 年 %s，7 年 %s；9 年不能新增确认。",
          fmt_pct(cum_conf[2]/n_base), fmt_pct(cum_conf[3]/n_base), fmt_pct(cum_conf[4]/n_base)),
  "- 主图适合用于 Figure 4 的累积发病面板，但必须同时展示风险人数。",
  sprintf("- 单次阈值跨越比例较高，连续确认曲线更适合作为稳健性展示；不要在论文中把首次发作 %s 写成持续认知障碍率。", fmt_pct(cum_first[5] / n_base)),
  "",
  "## 五、输出文件",
  "",
  "- `cumulative_incidence_2026-09-20.csv`",
  "- `cumulative_incidence_risk_table_2026-09-20.csv`",
  "- `cumulative_incidence_km_sensitivity_2026-09-20.csv`",
  "- `cumulative_incidence_2026-09-20.rds`",
  "- `output/09_图表/Figure4_cumulative_incidence_2026-09-20.png`",
  "- `output/09_图表/Figure4_cumulative_incidence_2026-09-20.tif`",
  "",
  sprintf("*生成时间：%s*", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
)
writeLines(report, file.path(OUT, "cumulative_incidence_2026-09-20.md"), useBytes = TRUE)

cat("\n========== T05 累积发病曲线完成 ==========\n")
cat("首次发作累积发病：", paste(fmt_pct(cum_first / n_base), collapse = " | "), "\n", sep = "")
cat("连续确认累积发病：", paste(fmt_pct(cum_conf / n_base), collapse = " | "), "\n", sep = "")
cat("主图：", file.path(FIG_DIR, "Figure4_cumulative_incidence_2026-09-20.png"), "\n", sep = "")
cat("TIFF：", ifelse(tiff_ok, file.path(FIG_DIR, "Figure4_cumulative_incidence_2026-09-20.tif"), "生成失败"), "\n", sep = "")
