# ============================================================================
# 25_a12_figures_report.R —— A12 图表与结果报告
# 输入：output/11_外部验证/A12_results_2026-09-26.rds
# 输出：Figure5、FigureS1、FigureS2（PNG 300 dpi + TIF 600 dpi）
#       A12_外部验证结果_2026-09-26.md
# ============================================================================
invisible(try(Sys.setlocale("LC_ALL", "zh_CN.UTF-8"), silent = TRUE))
.libPaths(c("D:/R/library", .libPaths()))
suppressPackageStartupMessages({ library(dplyr); library(pROC); library(knitr) })

P   <- "D:/BaiduSyncdisk/科研/徐中豪"
OUT <- file.path(P, "output", "11_外部验证")
FIG <- file.path(P, "output", "09_图表")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)

res <- readRDS(file.path(OUT, "A12_results_2026-09-26.rds"))
perf <- res$perf; strata <- res$strata; deciles <- res$deciles; thresholds <- res$thresholds
dca <- res$dca; share_country <- res$share_country; share_boot <- res$share_boot; meta <- res$meta
extra_dep <- res$extra_dep; extra_score_alt <- res$extra_score_alt; extra_elsa_age <- res$extra_elsa_age

a07 <- readRDS(file.path(P, "output/10_筛查指数/screening_index_memory_v3_2026-09-26.rds"))
model_u <- a07$model_u; model_levels <- model_u$xlevels$edu_group
cohorts <- c("HRS","ELSA","SHARE")
cohort_cols <- c(HRS = "#1f77b4", ELSA = "#2ca02c", SHARE = "#d62728")
scenario_lab <- c(baseline_all_mean_sd = "S0_baseline_all", complete_fu_mean_sd = "S1_complete_fu", fixed_charls_le4 = "S2_fixed_le4")

add_pred <- function(d, model) {
  d <- as.data.frame(d)
  d$edu_group <- factor(as.character(d$edu_group), levels = model_levels)
  d$p_hat <- as.numeric(predict(model, newdata = data.frame(age = d$age, edu_group = d$edu_group,
                                                            dep = d$dep, adl = d$adl), type = "response"))
  d$lp <- qlogis(pmin(pmax(d$p_hat, 1e-6), 1 - 1e-6))
  d
}
risk_sets <- list()
for (co in cohorts) {
  x <- readRDS(file.path(P, "external_data", co, "analysis", paste0(co, "_sensitivity_risk_sets_2026-09-26.rds")))
  risk_sets[[co]] <- add_pred(x[[paste(co, "baseline_all_mean_sd", sep = "__")]], model_u)
}

# ---------- 绘图函数 ----------
draw_main <- function() {
  par(mfrow = c(1,3), mar = c(4.2,4.2,3,1), oma = c(0,0,0,0))
  # (a) ROC
  plot(NA, xlim = c(0,1), ylim = c(0,1), xlab = "1 - Specificity", ylab = "Sensitivity",
       main = "(a) External ROC (S0 model)", bty = "n")
  abline(0,1,lty = 2,col = "grey70")
  for (co in cohorts) {
    d <- risk_sets[[co]]; d <- d[!is.na(d$p_hat), ]
    r <- roc(d$event, d$p_hat, direction = "<", quiet = TRUE)
    lines(1 - r$specificities, r$sensitivities, col = cohort_cols[co], lwd = 2)
  }
  a <- perf %>% filter(cut_definition == "baseline_all_mean_sd")
  leg <- sprintf("%s AUC %.3f (%.3f-%.3f)", a$cohort, a$auc_model, a$auc_model_low, a$auc_model_high)
  legend("bottomright", legend = leg, col = cohort_cols[a$cohort], lwd = 2, bty = "n", cex = 0.72)
  # (b) Calibration deciles
  dec <- deciles %>% filter(cut_definition == "baseline_all_mean_sd")
  ymax <- max(c(dec$observed_high, dec$pred_mean), na.rm = TRUE) * 1.05
  plot(NA, xlim = c(0, ymax), ylim = c(0, ymax), xlab = "Predicted risk", ylab = "Observed risk",
       main = "(b) Calibration by decile (S0 model)", bty = "n")
  abline(0,1,lty = 2,col = "grey70")
  for (co in cohorts) {
    d <- dec[dec$cohort == co, ]
    arrows(d$pred_mean, d$observed_low, d$pred_mean, d$observed_high, angle = 90, code = 3,
           length = 0.03, col = cohort_cols[co])
    lines(d$pred_mean, d$observed, type = "b", pch = 19, col = cohort_cols[co], lwd = 2)
  }
  legend("bottomright", legend = cohorts, col = cohort_cols[cohorts], lwd = 2, pch = 19, bty = "n", cex = 0.8)
  # (c) Forest plot
  a <- perf %>% filter(cut_definition == "baseline_all_mean_sd") %>% arrange(match(cohort, cohorts))
  m <- meta %>% filter(cut_definition == "baseline_all_mean_sd", metric == "model")
  labels <- c(a$cohort, "Pooled RE")
  est <- c(a$auc_model, m$pooled_auc); lo <- c(a$auc_model_low, m$pooled_low); hi <- c(a$auc_model_high, m$pooled_high)
  y <- rev(seq_along(labels))
  plot(NA, xlim = c(0.62, 0.86), ylim = c(0.5, length(labels)+0.5), yaxt = "n",
       xlab = "AUC (95% CI)", ylab = "", main = "(c) AUC forest (S0 model)", bty = "n")
  abline(v = 0.5, lty = 3, col = "grey70")
  segments(lo, y, hi, y, col = c(cohort_cols[a$cohort], "black"), lwd = 2)
  points(est, y, pch = c(19,19,19,18), col = c(cohort_cols[a$cohort], "black"), cex = c(1,1,1,1.4))
  axis(2, at = y, labels = labels, las = 1, cex.axis = 0.85)
}
draw_sens <- function() {
  par(mfrow = c(1,3), mar = c(4,3.2,3,1))
  for (cd in c("baseline_all_mean_sd","complete_fu_mean_sd","fixed_charls_le4")) {
    a <- perf %>% filter(cut_definition == cd) %>% arrange(match(cohort, cohorts))
    m <- meta %>% filter(cut_definition == cd, metric == "model")
    est <- c(a$auc_model, m$pooled_auc); lo <- c(a$auc_model_low, m$pooled_low); hi <- c(a$auc_model_high, m$pooled_high)
    labels <- c(a$cohort, "Pooled")
    y <- rev(seq_along(labels))
    plot(NA, xlim = c(0.62, 0.88), ylim = c(0.5, length(labels)+0.5), yaxt = "n",
         xlab = "AUC", ylab = "", main = scenario_lab[[cd]], bty = "n")
    abline(v = 0.5, lty = 3, col = "grey70")
    segments(lo, y, hi, y, col = c(cohort_cols[a$cohort], "black"), lwd = 2)
    points(est, y, pch = 19, col = c(cohort_cols[a$cohort], "black"), cex = c(1,1,1,1.3))
    axis(2, at = y, labels = labels, las = 1, cex.axis = 0.8)
  }
}
draw_share_country <- function() {
  d <- share_country %>% arrange(auc_model)
  y <- seq_len(nrow(d))
  par(mar = c(4,7,3,1))
  plot(NA, xlim = c(0.6, 0.9), ylim = c(0.5, nrow(d)+0.5), yaxt = "n",
       xlab = "AUC (95% CI)", ylab = "", main = "SHARE country-level AUC (S0 model)", bty = "n")
  abline(v = 0.5, lty = 3, col = "grey70")
  segments(d$auc_model_low, y, d$auc_model_high, y, col = "#d62728", lwd = 2)
  points(d$auc_model, y, pch = 19, col = "#d62728", cex = 0.9)
  axis(2, at = y, labels = d$country_label, las = 1, cex.axis = 0.78)
}
draw_cal_sens <- function() {
  par(mfrow = c(1,3), mar = c(4,3.5,3,1))
  for (cd in c("baseline_all_mean_sd","complete_fu_mean_sd","fixed_charls_le4")) {
    dec <- deciles %>% filter(cut_definition == cd)
    ymax <- max(c(dec$observed_high, dec$pred_mean), na.rm = TRUE) * 1.05
    plot(NA, xlim = c(0, ymax), ylim = c(0, ymax), xlab = "Predicted risk", ylab = "Observed risk",
         main = scenario_lab[[cd]], bty = "n")
    abline(0,1,lty = 2,col = "grey70")
    for (co in cohorts) {
      d <- dec[dec$cohort == co, ]
      arrows(d$pred_mean, d$observed_low, d$pred_mean, d$observed_high, angle = 90, code = 3,
             length = 0.03, col = cohort_cols[co])
      lines(d$pred_mean, d$observed, type = "b", pch = 19, col = cohort_cols[co], lwd = 1.7)
    }
    if (cd == "baseline_all_mean_sd") legend("bottomright", legend = cohorts, col = cohort_cols[cohorts], lwd = 2, pch = 19, bty = "n", cex = 0.75)
  }
}
draw_dca <- function() {
  par(mfrow = c(1,3), mar = c(4,4,3,1))
  for (co in cohorts) {
    d <- dca %>% filter(cohort == co, cut_definition == "baseline_all_mean_sd")
    ylim <- range(d$net_benefit, na.rm = TRUE)
    plot(NA, xlim = c(0.05,0.40), ylim = ylim, xlab = "Threshold probability", ylab = "Net benefit",
         main = paste0(co, " (S0)"), bty = "n")
    abline(h = 0, lty = 3, col = "grey70")
    for (st in c("model","all")) {
      dd <- d[d$strategy == st, ]
      lines(dd$threshold, dd$net_benefit, col = ifelse(st == "model", cohort_cols[co], "grey40"),
            lwd = 2, lty = ifelse(st == "model", 1, 2))
    }
    legend("topright", legend = c("Model","Treat all"), col = c(cohort_cols[co],"grey40"),
           lwd = 2, lty = c(1,2), bty = "n", cex = 0.75)
  }
}
save_both <- function(base, width, height, draw) {
  png(paste0(base, ".png"), width = width, height = height, units = "in", res = 300)
  draw(); dev.off()
  tiff(paste0(base, ".tif"), width = width, height = height, units = "in", res = 600, compression = "lzw")
  draw(); dev.off()
}
save_both(file.path(FIG, "Figure5_external_validation_2026-09-26"), 13, 4.5, draw_main)
save_both(file.path(FIG, "FigureS1_external_sensitivity_forest_2026-09-26"), 12, 4.2, draw_sens)
save_both(file.path(FIG, "FigureS2_SHARE_country_forest_2026-09-26"), 8, 5.5, draw_share_country)
save_both(file.path(FIG, "FigureS3_external_calibration_S0_S1_S2_2026-09-26"), 12, 4.2, draw_cal_sens)
save_both(file.path(FIG, "FigureS4_external_DCA_S0_2026-09-26"), 12, 4.2, draw_dca)

# ---------- 表格格式化 ----------
fmt3 <- function(x) ifelse(is.na(x), "NA", sprintf("%.3f", x))
fmt2 <- function(x) ifelse(is.na(x), "NA", sprintf("%.2f", x))
fmt_pct <- function(x) ifelse(is.na(x), "NA", sprintf("%.1f%%", 100*x))
rng <- function(a,b) paste0(fmt3(a), " (", fmt3(b), ")", ifelse(is.na(a) | is.na(b), "", ""))
rng2 <- function(a,b) paste0(fmt2(a), " (", fmt2(b), ")")

perf_s0 <- perf %>% filter(cut_definition == "baseline_all_mean_sd") %>% arrange(match(cohort, cohorts))
tab_main <- data.frame(
  队列 = perf_s0$cohort,
  模型n = perf_s0$n_model,
  事件 = perf_s0$events_model,
  连续模型AUC = paste0(fmt3(perf_s0$auc_model), " (", fmt3(perf_s0$auc_model_low), "–", fmt3(perf_s0$auc_model_high), ")"),
  整数分AUC = paste0(fmt3(perf_s0$auc_score_full), " (", fmt3(perf_s0$auc_score_full_low), "–", fmt3(perf_s0$auc_score_full_high), ")"),
  校准斜率 = paste0(fmt2(perf_s0$cal_slope), " (", fmt2(perf_s0$cal_slope_low), " to ", fmt2(perf_s0$cal_slope_high), ")"),
  校准截距 = paste0(fmt2(perf_s0$cal_intercept), " (", fmt2(perf_s0$cal_intercept_low), " to ", fmt2(perf_s0$cal_intercept_high), ")"),
  Brier = fmt3(perf_s0$brier), 标化Brier = fmt3(perf_s0$scaled_brier)
)
tab_strata <- strata %>% filter(cut_definition == "baseline_all_mean_sd") %>%
  mutate(分层 = as.character(stratum),
         观察事件率 = fmt_pct(rate), 预测风险 = fmt_pct(pred_mean), OE = fmt2(oe_model)) %>%
  select(队列 = cohort, 分层, n, 事件 = events, 观察事件率, 预测风险, OE)
tab_thr <- thresholds %>% filter(cut_definition == "baseline_all_mean_sd") %>%
  mutate(阈值 = paste0(threshold*100, "%"), 灵敏度 = fmt_pct(sensitivity), 特异度 = fmt_pct(specificity),
         PPV = fmt_pct(ppv), NPV = fmt_pct(npv), 净获益 = fmt3(nb_model), 全部干预净获益 = fmt3(nb_all)) %>%
  mutate(阳性率 = fmt_pct(pos_rate)) %>% select(队列 = cohort, 阈值, 阳性数 = n_positive, 阳性率, 灵敏度, 特异度, PPV, NPV, 净获益, 全部干预净获益)
tab_meta <- meta %>% mutate(cut = scenario_lab[cut_definition], metric = ifelse(metric == "model", "连续模型", "整数分"),
                            pooled = paste0(fmt3(pooled_auc), " (", fmt3(pooled_low), "–", fmt3(pooled_high), ")"),
                            I2 = sprintf("%.1f%%", I2)) %>%
  select(切点 = cut, 指标 = metric, 合并AUC = pooled, I2, k)
tab_sens <- perf %>% filter(cut_definition != "baseline_all_mean_sd") %>%
  mutate(切点 = scenario_lab[cut_definition],
         连续模型AUC = paste0(fmt3(auc_model), " (", fmt3(auc_model_low), "–", fmt3(auc_model_high), ")"),
         整数分AUC = paste0(fmt3(auc_score_full), " (", fmt3(auc_score_full_low), "–", fmt3(auc_score_full_high), ")"),
         校准斜率 = paste0(fmt2(cal_slope), " (", fmt2(cal_slope_low), " to ", fmt2(cal_slope_high), ")"), 校准截距 = paste0(fmt2(cal_intercept), " (", fmt2(cal_intercept_low), " to ", fmt2(cal_intercept_high), ")"), Brier = fmt3(brier)) %>%
  select(队列 = cohort, 切点, n = n_model, 事件 = events_model, 连续模型AUC, 整数分AUC, 校准斜率, 校准截距, Brier)
tab_share <- share_country %>% arrange(desc(n)) %>%
  mutate(AUC = paste0(fmt3(auc_model), " (", fmt3(auc_model_low), "–", fmt3(auc_model_high), ")"),
         事件率 = fmt_pct(event_rate), 校准斜率 = fmt2(cal_slope), 校准截距 = fmt2(cal_intercept)) %>%
  select(国家 = country_label, n, 事件 = events, 事件率, AUC, 校准斜率, 校准截距)
tab_extra <- bind_rows(
  extra_dep %>% transmute(队列 = cohort, 敏感性 = "抑郁切点", 规则 = dep_rule, n, 事件 = events,
                          连续模型AUC = paste0(fmt3(auc_model), " (", fmt3(auc_model_low), "–", fmt3(auc_model_high), ")"),
                          整数分AUC = paste0(fmt3(auc_score), " (", fmt3(auc_score_low), "–", fmt3(auc_score_high), ")")),
  extra_score_alt %>% transmute(队列 = cohort, 敏感性 = "score_alt", 规则 = "位置式教育映射", n, 事件 = events,
                                连续模型AUC = "—", 整数分AUC = paste0(fmt3(auc_score_alt), " (", fmt3(auc_low), "–", fmt3(auc_high), ")")),
  extra_elsa_age %>% transmute(队列 = cohort, 敏感性 = "ELSA 排除age=90", 规则 = "年龄顶码敏感性", n, 事件 = events,
                               连续模型AUC = paste0(fmt3(auc_model), " (", fmt3(auc_model_low), "–", fmt3(auc_model_high), ")"),
                               整数分AUC = paste0(fmt3(auc_score), " (", fmt3(auc_score_low), "–", fmt3(auc_score_high), ")"))
)
tab_dca <- dca %>% filter(cut_definition == "baseline_all_mean_sd", strategy %in% c("model","all"),
                          round(threshold, 2) %in% c(0.10,0.15,0.20,0.25,0.30)) %>%
  mutate(阈值 = paste0(threshold*100, "%")) %>%
  select(队列 = cohort, 阈值, 策略 = strategy, 净获益 = net_benefit) %>%
  mutate(净获益 = fmt3(净获益))

# ---------- 报告 ----------
m_s0_model <- meta %>% filter(cut_definition == "baseline_all_mean_sd", metric == "model")
m_s1_model <- meta %>% filter(cut_definition == "complete_fu_mean_sd", metric == "model")
m_s2_model <- meta %>% filter(cut_definition == "fixed_charls_le4", metric == "model")
meta_text <- sprintf(
  "> 随机效应采用 DerSimonian–Laird 方法，在 logit(AUC) 尺度合并。S0 连续模型合并 AUC = %.3f（%.3f–%.3f），但 **I² = %.1f%%**；S1 合并 AUC = %.3f（%.3f–%.3f），I² = %.1f%%；S2 合并 AUC = %.3f（%.3f–%.3f），I² = %.1f%%。异质性很高，合并值仅作描述，论文应重点报告各队列结果与异质性来源。",
  m_s0_model$pooled_auc, m_s0_model$pooled_low, m_s0_model$pooled_high, m_s0_model$I2,
  m_s1_model$pooled_auc, m_s1_model$pooled_low, m_s1_model$pooled_high, m_s1_model$I2,
  m_s2_model$pooled_auc, m_s2_model$pooled_low, m_s2_model$pooled_high, m_s2_model$I2)
sb_model <- share_boot[share_boot$cut_definition == "baseline_all_mean_sd" & share_boot$metric == "auc_model", ]
sb_score <- share_boot[share_boot$cut_definition == "baseline_all_mean_sd" & share_boot$metric == "auc_score", ]
share_boot_text <- sprintf(
  "> 按国家 cluster bootstrap（%d 次）：连续模型 AUC **%.3f（%.3f–%.3f）**，整数分 AUC **%.3f（%.3f–%.3f）**。点估计与经典 AUC 接近；cluster bootstrap CI 更宽，反映国家聚类不确定性。国家层面差异应在论文中重点说明。",
  sb_model$B, sb_model$orig_auc, sb_model$boot_low, sb_model$boot_high,
  sb_score$orig_auc, sb_score$boot_low, sb_score$boot_high)
share_s0_fig <- risk_sets[["SHARE"]]
n_missing_adl_share <- sum(is.na(share_s0_fig$p_hat))
events_missing_adl_share <- sum(share_s0_fig$event[is.na(share_s0_fig$p_hat)])
share_auc_text <- paste0(
  sprintf("> S0 SHARE 国家异质性明显：国家 AUC 范围约 %.3f–%.3f，事件率范围 %.1f%%–%.1f%%。", 
          min(share_country$auc_model), max(share_country$auc_model),
          100*min(share_country$event_rate), 100*max(share_country$event_rate)),
  "国家层面校准斜率/截距样本量较小，属探索性结果，不作为主要结论。
",
  sprintf("> 国家层面表为连续模型完整样本 n=%s；全 SHARE 风险集 n=%s（%d 人 ADL 缺失，%d 个事件）。",
          format(sum(!is.na(share_s0_fig$p_hat)), big.mark = ","), format(nrow(share_s0_fig), big.mark = ","),
          n_missing_adl_share, events_missing_adl_share))
md <- c(
"# A12：外部队列正式外部验证结果（2026-09-26）",
"",
"> 主分析：**S0 = 各队列基线全部符合条件者 mean−1SD**；正式敏感性：**S1 = complete-fu mean−1SD**、**S2 = 固定 CHARLS ≤4**。",
"> 模型：CHARLS 固定模型（年龄 + 教育 + 抑郁 + ADL），**不重新拟合**。",
"> 队列：HRS（w10→w12）、ELSA（w4→w6）、SHARE（w4→w6）。",
"",
"## 方法概要",
"",
"- 使用 CHARLS 固定模型（年龄 + 教育 + 抑郁 + ADL），在外部队列中**不重新拟合、不重新选变量**；",
"- 主分析 S0 = 各队列基线全部符合条件者 mean−1SD；S1 = complete-fu mean−1SD；S2 = 固定 CHARLS ≤4；",
"- AUC 用 DeLong 法计算 95%CI；校准斜率 = glm(event ~ logit(p)) 中 logit(p) 的系数；校准截距 = glm(event ~ offset(logit(p))) 的截距（calibration-in-the-large）；Brier 与标化 Brier 按模型完整样本计算；",
"- 跨队列合并使用 DerSimonian–Laird 随机效应模型（logit AUC 尺度），SE 由各队列 95%CI 反推；SHARE 在合并时使用国家 cluster bootstrap CI（固定种子 20260926），HRS/ELSA 使用 DeLong CI；",
"- 风险分层使用 CHARLS 固定切点 低 ≤−3｜中 −2～4｜高 ≥5。",
"",
"## 一、主分析（S0）：连续模型外部表现",
"",
kable(tab_main, format = "pipe", row.names = FALSE),
"",
"> 连续模型分析使用 ADL 完整案例：HRS 13,972（缺失 79 人，其中事件 9）、ELSA 6,725（缺失 0）、SHARE 26,874（缺失 5，事件 0）。整数分 AUC 使用全部风险人群。
> 在模型完整样本中，整数分 AUC 分别为 HRS 0.693、ELSA 0.759、SHARE 0.752，与连续模型几乎一致；上表整数分 AUC 使用全部风险人群。主表 SHARE 的 CI 为 DeLong CI，未考虑国家聚类；按国家 cluster bootstrap 的 CI 见第六节。",
"> **主要发现**：连续模型 AUC 为 HRS 0.696、ELSA 0.763、SHARE 0.751；整数分 AUC 为 0.692、0.759、0.752。点估计较接近（约 0.70–0.76），但跨队列统计异质性高（I² 85%–95%，见第四节）。",
"> **校准异质性明显**：HRS 校准斜率 1.23、截距 +0.33（系统低估风险）；ELSA 斜率 1.82、截距 −0.19；SHARE 斜率 1.09、截距 −0.31（系统高估风险）。",
"",
"## 二、S0 风险分层（固定评分切点）",
"",
kable(tab_strata, format = "pipe", row.names = FALSE),
"",
"> 分层切点为 CHARLS 固定切点 低 ≤−3｜中 −2～4｜高 ≥5。观察事件率按全部风险人群计算；预测风险与 O/E 按模型完整样本计算。",
"> HRS 高危组事件率 36.3%、ELSA 31.7%、SHARE 27.1%；高危/低危比分别约 **4.9×、10.5×、7.3×**。",
"",
"## 三、S0 阈值表现（预测概率）",
"",
kable(tab_thr, format = "pipe", row.names = FALSE),
"",
"> 阈值表现基于模型完整样本；对应的 DCA 见 FigureS4。",
"",
"## 四、跨队列随机效应汇总",
"",
kable(tab_meta, format = "pipe", row.names = FALSE),
"",
meta_text,
"",
"## 五、正式敏感性分析（S1/S2）",
"",
kable(tab_sens, format = "pipe", row.names = FALSE),
"",
"> S1/S2 表的 n 与事件为**连续模型完整样本**；对应全风险集事件数为：S2 HRS 945、ELSA 324、SHARE 2,256，S1 ELSA 706，HRS/SHARE 与 S0 相同。
> ELSA 是唯一被切点定义实质影响的队列：S0 ≤6 vs S1 ≤7；S2 固定 ≤4 在三个队列均提高 AUC，但事件率明显降低（全风险集 HRS 6.1%、ELSA 4.5%、SHARE 7.9%），且依赖跨队列记忆分绝对可比性。三套结果均应报告。",
"",
"## 六、SHARE 国家层面与聚类",
"",
kable(tab_share, format = "pipe", row.names = FALSE),
"",
share_auc_text,
share_boot_text,
"",
"## 七、额外敏感性分析",
"",
kable(tab_extra, format = "pipe", row.names = FALSE),
"",
"> 抑郁切点与 ELSA age90 行为模型完整样本；score_alt 行为全风险集。抑郁切点改用 CES-D-8 ≥3 / EURO-D ≥5 后，AUC 基本不变；位置式教育映射 score_alt 的 AUC 为 HRS 0.707、ELSA 0.755、SHARE 0.756；ELSA 排除年龄顶码 90 后 AUC 轻微下降（模型 0.757、整数分 0.753）。",
"",
"## 八、决策曲线（S0，净获益节选）",
"",
kable(tab_dca, format = "pipe", row.names = FALSE),
"",
"> DCA 基于模型完整样本；在 10%–30% 阈值范围内，S0 连续模型的净获益均高于“全部干预”策略；约 35% 后逐渐接近 0。",
"",
"## 九、解释与限制",
"",
"1. **定位**：该工具是术前相关认知脆弱性**筛查/分流**工具，不预测术后谵妄、PND、ICU 或住院结局。",
"2. **结局口径**：外部队列结局是 4 年随访“跨阈值”，含基线位置、测量误差与回归均值成分，不等同于临床痴呆/MCI 发病率。",
"3. **区分度**：连续模型 AUC 约 0.70–0.76，属中等；不能声称高精度预测。",
"4. **校准**：HRS 系统性低估风险、SHARE 系统性高估风险、ELSA 斜率明显 >1；这些校准差异必须如实报告，不能只报告 AUC。",
"5. **异质性**：跨队列合并 I² 85%–95%，提示队列间差异大；SHARE 国家异质性也很大。",
"6. **测量**：ELSA 年龄在 90 岁顶码；跨队列记忆词表、语言与测试实施不同；S2 固定 ≤4 假设绝对可比，需谨慎。",
"7. **缺失**：连续模型分析排除了 ADL 缺失的 HRS 79 人、SHARE 5 人；整数分分析不受影响。",
"8. **临床验证**：尚未在真实手术队列验证，不能声称改善围术期决策或临床结局。",
"9. **随访选择**：分析限定于 4 年后记忆分可测者；死亡、失访或记忆分缺失者被排除，可能存在选择偏倚，且未将死亡作为竞争风险处理。",
"10. **阈值含义**：S0/S1 使用各队列内 mean−1SD，与 CHARLS 固定 ≤4 不是同一绝对结局；S2 固定 ≤4 回答“同一绝对阈值”问题，但假设跨队列记忆分绝对可比。三套结果互补，不能互相替代。",
"",
"## 十、输出文件",
"",
"- `A12_主结果_三队列_S0_S1_S2.csv`",
"- `A12_风险分层_S0_S1_S2.csv`",
"- `A12_校准十分位_S0_S1_S2.csv`",
"- `A12_阈值表现_S0_S1_S2.csv`",
"- `A12_DCA_S0_S1_S2.csv`",
"- `A12_SHARE国家层面.csv`、`A12_SHARE国家聚类bootstrap.csv`",
"- `A12_跨队列随机效应.csv`",
"- `A12_敏感性_抑郁切点.csv`、`A12_敏感性_score_alt.csv`、`A12_敏感性_ELSA_age90.csv`",
"- `A12_sessionInfo_2026-09-26.txt`",
"- `A12_独立复核_2026-09-26.txt`",
"- `Figure5_external_validation_2026-09-26.png/.tif`",
"- `FigureS1_external_sensitivity_forest_2026-09-26.png/.tif`",
"- `FigureS2_SHARE_country_forest_2026-09-26.png/.tif`",
"- `FigureS3_external_calibration_S0_S1_S2_2026-09-26.png/.tif`",
"- `FigureS4_external_DCA_S0_2026-09-26.png/.tif`"
)
writeLines(md, file.path(OUT, "A12_外部验证结果_2026-09-26.md"), useBytes = TRUE)
cat("Figures and report written.\n")