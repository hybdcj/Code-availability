# ============================================================================
# F03 完整版：时段验证（A）＋ 同批人两套权重正面对比（B）
# 人群：CHARLS 2018 基线 → 2020 结局；2018 未受损者进风险人群
# 结局两种口径：
#   D1 绝对切点：2018 记忆分 mean−1SD 取整（≤3），2020 沿用同一绝对分
#   D2 等百分位锚定：取 D1 切点在 2018 分布的百分位 p，2020 用自身分布同一百分位
# 两套权重：冻结（2011→2015，A14）vs 重估（2018→2020 现拟合）
# 乐观偏倚：Harrell bootstrap 校正 ＋ 嵌套 bootstrap 的 ΔAUC 区间
# 输出：C:/y03tmp/out/F03_report.txt、F03_*.csv、F03_*.png
# ============================================================================
invisible(try(Sys.setlocale("LC_ALL", "zh_CN.UTF-8"), silent = TRUE))
suppressPackageStartupMessages({ library(data.table); library(pROC) })
if (dir.exists("D:/R/library")) .libPaths(c("D:/R/library", .libPaths()))
TMP <- "C:/y03tmp/tmp"; OUT <- "C:/y03tmp/out"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
con <- file(file.path(OUT, "F03_report.txt"), open = "wt", encoding = "UTF-8")
say <- function(...) { s <- paste0(..., collapse = ""); cat(s, "\n"); cat(s, "\n", file = con) }
pt  <- function(x) say(paste(capture.output(print(x)), collapse = "\n"))

BETA_OLD <- c("(Intercept)" = -4.19524508778947, age = 0.0492039250378111,
              eduP = -0.290105771000328, eduM = -0.994711800202235, eduH = -1.43588388035287,
              dep = 0.412166484954086, adl = -0.0110119344848014)
LV <- c("Less than primary", "Primary school", "Middle school", "High school or above")

## ---------------- 数据集 ----------------
w <- as.data.table(readRDS(file.path(TMP, "waves_cognition.rds")))
a18 <- w[year == 2018, .(ID = IDn, mem18 = as.numeric(imrc) + as.numeric(dlrc),
                         age18 = as.numeric(age), edu18 = edu_group,
                         cesd18 = as.numeric(cesd10), adl18 = as.numeric(adl),
                         wt18 = as.numeric(weight), comm18 = as.character(communityID))]
a20 <- w[year == 2020, .(ID = IDn, mem20 = as.numeric(imrc) + as.numeric(dlrc),
                         wt20 = as.numeric(weight))]
m <- merge(a18, a20, by = "ID")
m <- m[!is.na(mem18) & !is.na(mem20) & !is.na(age18) & !is.na(edu18) &
         !is.na(cesd18) & !is.na(adl18)]
m[, dep := as.integer(cesd18 >= 10)]
m[, edu_group := factor(edu18, levels = LV)]
say("分析人群 n = ", nrow(m), "（其中 2018 年 <45 岁：", sum(m$age18 < 45), "）")
m <- m[age18 >= 45]
say("限定 2018 年 ≥45 岁后 n = ", nrow(m))

## ---------------- 结局口径 ----------------
cut18 <- floor(mean(m$mem18) - sd(m$mem18))
p_anchor <- mean(m$mem18 <= cut18)
cut20 <- as.numeric(quantile(m$mem20, p_anchor, type = 1))
say("\n2018 记忆分 mean=", round(mean(m$mem18), 2), " SD=", round(sd(m$mem18), 2),
    " → D1 切点 = ≤", cut18, "（该切点在 2018 分布的百分位 = ", round(100 * p_anchor, 1), "%）")
say("2020 记忆分 mean=", round(mean(m$mem20), 2), " SD=", round(sd(m$mem20), 2),
    " → D2 等百分位切点 = ≤", cut20)

defs <- list(D1 = list(k18 = cut18, k20 = cut18, lab = paste0("D1 绝对切点 ≤", cut18)),
             D2 = list(k18 = cut18, k20 = cut20, lab = paste0("D2 等百分位锚定（2020 ≤", cut20, "）")))

## ---------------- 工具函数 ----------------
pred_old <- function(d) {
  plogis(BETA_OLD[["(Intercept)"]] + BETA_OLD[["age"]] * d$age18 + BETA_OLD[["dep"]] * d$dep +
           BETA_OLD[["adl"]] * d$adl18 +
           BETA_OLD[["eduP"]] * as.integer(d$edu_group == LV[2]) +
           BETA_OLD[["eduM"]] * as.integer(d$edu_group == LV[3]) +
           BETA_OLD[["eduH"]] * as.integer(d$edu_group == LV[4]))
}
fit_new <- function(d) glm(event ~ age18 + edu_group + dep + adl18, data = d, family = binomial())
calc <- function(d, p, tag = "") {
  r <- roc(d$event, p, quiet = TRUE, direction = "<")
  ci <- as.numeric(ci.auc(r))
  lp <- qlogis(pmin(pmax(p, 1e-9), 1 - 1e-9))
  cal <- glm(d$event ~ lp, family = binomial())
  data.table(标签 = tag, n = nrow(d), 事件 = sum(d$event),
             AUC = round(ci[2], 3), AUC低 = round(ci[1], 3), AUC高 = round(ci[3], 3),
             校准斜率 = round(coef(cal)[2], 3), 校准截距 = round(coef(cal)[1], 3),
             Brier = round(mean((p - d$event)^2), 4),
             OE比 = round(sum(d$event) / sum(p), 2))
}
wt_auc <- function(p, y, wgt) {
  cs <- y == 1; ct <- y == 0
  num <- 0
  for (i in which(cs)) num <- num + wgt[i] * sum(wgt[ct] * ((p[i] > p[ct]) + 0.5 * (p[i] == p[ct])))
  num / (sum(wgt[cs]) * sum(wgt[ct]))
}
nb_curve <- function(p, y, th) sapply(th, function(t) sum(p >= t & y == 1)/length(y) -
                                        sum(p >= t & y == 0)/length(y) * (t/(1-t)))

## ---------------- A：冻结模型时段验证 ----------------
say("\n========== A. 冻结模型在 2018→2020 的时段验证 ==========")
A <- data.table(); stA <- NULL; d1 <- NULL
for (nm in names(defs)) {
  dfn <- defs[[nm]]
  d <- copy(m[mem18 > dfn$k18])[, event := as.integer(mem20 <= dfn$k20)]
  d[, phat := pred_old(d)]
  A <- rbind(A, calc(d, d$phat, dfn$lab))
  if (nm == "D1") {
    d1 <- copy(d)
    d1[, score := 2L * pmax(floor((age18 - 45)/5), 0) +
          (-2L) * as.integer(edu_group == LV[2]) + (-7L) * as.integer(edu_group == LV[3]) +
          (-10L) * as.integer(edu_group == LV[4]) + 3L * dep]
    stA <- d1[, .(n = .N, 事件 = sum(event), 事件率 = paste0(round(100*mean(event), 1), "%"),
                  平均分 = round(mean(score), 1)),
              by = .(风险层 = cut(score, c(-Inf, -3.5, 4.5, Inf), labels = c("低(≤-3)", "中(-2~+4)", "高(≥+5)")))][order(风险层)]
  }
}
pt(A)
say("\n冻结评分风险分层（", defs$D1$lab, "）："); pt(stA)
say("（对照：开发样本表观 AUC 0.702、Bootstrap 校正 0.698；本时段 n=", nrow(d1),
    "，事件=", sum(d1$event), "，事件率=", round(100*mean(d1$event), 2), "%）")
fwrite(A, file.path(OUT, "F03_A_temporal_validation.csv"), bom = TRUE)
fwrite(stA, file.path(OUT, "F03_A_frozen_strata.csv"), bom = TRUE)

## ---------------- B：同批人两套权重正面对比 ----------------
say("\n========== B. 同批人两套权重正面对比 ==========")
Bset <- 2000; set.seed(20261014)
B <- data.table(); DCAtab <- data.table(); stAll <- NULL; dtab <- NULL
for (nm in names(defs)) {
  dfn <- defs[[nm]]
  d <- copy(m[mem18 > dfn$k18])[, event := as.integer(mem20 <= dfn$k20)]
  d[, phat_old := pred_old(d)]
  fit <- fit_new(d)
  d[, phat_new := predict(fit, type = "response")]
  say("\n---- ", dfn$lab, "（n=", nrow(d), "，事件=", sum(d$event), "）----")
  say("重估系数（2018→2020）："); pt(round(coef(fit), 4))

  rowA <- calc(d, d$phat_old, "旧权重（2011 拟合）")
  rowB <- calc(d, d$phat_new, "新权重（2018→2020 拟合，表观）")
  r_old <- roc(d$event, d$phat_old, quiet = TRUE, direction = "<")
  r_new <- roc(d$event, d$phat_new, quiet = TRUE, direction = "<")
  pdel <- tryCatch(roc.test(r_old, r_new, paired = TRUE)$p.value, error = function(e) NA_real_)

  opt <- rep(NA_real_, Bset); dA <- rep(NA_real_, Bset)
  for (b in seq_len(Bset)) {
    idx <- sample.int(nrow(d), replace = TRUE); db <- d[idx]
    f <- tryCatch(glm(event ~ age18 + edu_group + dep + adl18, data = db, family = binomial()),
                  error = function(e) NULL)
    if (is.null(f) || length(unique(db$event)) < 2) next
    p_orig <- predict(f, newdata = d, type = "response")
    p_boot <- predict(f, type = "response")
    opt[b] <- as.numeric(auc(roc(db$event, p_boot, quiet = TRUE, direction = "<"))) -
              as.numeric(auc(roc(d$event, p_orig, quiet = TRUE, direction = "<")))
    dA[b]  <- as.numeric(auc(roc(db$event, p_boot, quiet = TRUE, direction = "<"))) -
              as.numeric(auc(roc(db$event, pred_old(db), quiet = TRUE, direction = "<")))
  }
  auc_new_corrected <- as.numeric(auc(r_new)) - mean(opt, na.rm = TRUE)
  say("\n旧权重 AUC=", round(rowA$AUC, 3), " (", round(rowA$AUC低, 3), "-", round(rowA$AUC高, 3),
      ")｜校准斜率=", rowA$校准斜率, " 截距=", rowA$校准截距, "｜Brier=", rowA$Brier)
  say("新权重 AUC=", round(rowB$AUC, 3), " (", round(rowB$AUC低, 3), "-", round(rowB$AUC高, 3),
      ")｜校准斜率=", rowB$校准斜率, " 截距=", rowB$校准截距, "｜Brier=", rowB$Brier)
  say("新权重 Harrell 校正 AUC=", round(auc_new_corrected, 3), "（乐观度 ", round(mean(opt, na.rm = TRUE), 4), "）")
  say("表观 ΔAUC=", round(rowB$AUC - rowA$AUC, 4), "｜配对 DeLong p=", format(pdel, digits = 3))
  say("嵌套 bootstrap ΔAUC 中位数=", round(median(dA, na.rm = TRUE), 4),
      "，95% 区间=[", round(quantile(dA, .025, na.rm = TRUE), 4), ", ",
      round(quantile(dA, .975, na.rm = TRUE), 4), "]")
  yo <- function(r) { co <- coords(r, "best", best.method = "youden",
                                   ret = c("sensitivity", "specificity"))
                      data.table(敏感度 = round(co$sensitivity, 3), 特异度 = round(co$specificity, 3),
                                 约登 = round(co$sensitivity + co$specificity - 1, 3)) }
  say("\n约登切点：旧权重"); pt(yo(r_old)); say("约登切点：新权重"); pt(yo(r_new))

  ## 新权重整数化（沿用 A14 规则）
  bb <- coef(fit)
  tb <- c(age5 = 5 * bb[["age18"]], eduP = bb[["edu_groupPrimary school"]],
          eduM = bb[["edu_groupMiddle school"]], eduH = bb[["edu_groupHigh school or above"]],
          dep = bb[["dep"]], adl = bb[["adl18"]])
  pts <- round(10 * tb / max(abs(tb)))
  say("\n新权重折算 β："); pt(round(tb, 3)); say("新权重整数化分值："); pt(pts)
  d[, score_new := as.integer(pts[["age5"]]) * pmax(floor((d$age18 - 45)/5), 0) +
        as.integer(pts[["eduP"]]) * as.integer(d$edu_group == LV[2]) +
        as.integer(pts[["eduM"]]) * as.integer(d$edu_group == LV[3]) +
        as.integer(pts[["eduH"]]) * as.integer(d$edu_group == LV[4]) +
        as.integer(pts[["dep"]]) * d$dep]
  d[, score_old := 2L * pmax(floor((d$age18 - 45)/5), 0) +
        (-2L) * as.integer(d$edu_group == LV[2]) + (-7L) * as.integer(d$edu_group == LV[3]) +
        (-10L) * as.integer(d$edu_group == LV[4]) + 3L * d$dep]
  say("\n整数分 AUC：旧 ", round(as.numeric(auc(roc(d$event, d$score_old, quiet = TRUE, direction = "<"))), 3),
      "｜新 ", round(as.numeric(auc(roc(d$event, d$score_new, quiet = TRUE, direction = "<"))), 3))
  LVL <- c("低(≤-3)", "中(-2~+4)", "高(≥+5)")
  d[, 层_旧 := cut(score_old, c(-Inf, -3.5, 4.5, Inf), labels = LVL)]
  d[, 层_新 := cut(score_new, c(-Inf, -3.5, 4.5, Inf), labels = LVL)]
  s_old <- d[, .(n = .N, 事件 = sum(event), 事件率 = round(100*mean(event), 1)), by = .(风险层 = 层_旧)]
  s_old[, 评分 := "旧(2011)"]
  s_new <- d[, .(n = .N, 事件 = sum(event), 事件率 = round(100*mean(event), 1)), by = .(风险层 = 层_新)]
  s_new[, 评分 := "新(2018-2020)"]
  stk <- rbind(s_old[, .(评分, 风险层, n, 事件, 事件率)][order(风险层)],
               s_new[, .(评分, 风险层, n, 事件, 事件率)][order(风险层)])
  pt(stk)
  say("\n按“高危(≥+5)”判阳性：")
  for (sc in c("score_old", "score_new")) {
    pos <- d[[sc]] >= 5
    say("  ", sc, "：阳性 n=", sum(pos), "（", round(100*mean(pos), 1), "%）",
        "｜敏感度=", round(100*sum(d$event[pos])/sum(d$event), 1), "%",
        "｜特异度=", round(100*sum(!pos & d$event == 0)/sum(d$event == 0), 1), "%",
        "｜阳性者事件率=", round(100*mean(d$event[pos]), 1), "%")
  }
  ## 模型更新阶梯：① 原样 ② 仅重校准截距（1 参数）③ 斜率＋截距重校准（2 参数）④ 全重估
  lp_old <- qlogis(pmin(pmax(d$phat_old, 1e-9), 1 - 1e-9))
  cal_int <- glm(event ~ 1, offset = lp_old, data = d, family = binomial())
  d[, phat_int := plogis(lp_old + coef(cal_int)[1])]
  rowC <- calc(d, d$phat_int, "冻结＋仅重校准截距")
  cal_sl <- glm(event ~ lp_old, data = d, family = binomial())
  d[, phat_slope := predict(cal_sl, type = "response")]
  rowD <- calc(d, d$phat_slope, "冻结＋斜率与截距重校准")
  say("\n仅重校准截距：平移量=", round(coef(cal_int)[1], 3), "｜Brier=", rowC$Brier,
      "｜AUC 不变（单调变换）=", rowC$AUC)
  say("斜率＋截距重校准：斜率=", round(coef(cal_sl)[2], 3), " 截距=", round(coef(cal_sl)[1], 3),
      "｜Brier=", rowD$Brier, "｜AUC（单调变换）=", rowD$AUC)

  dw <- d[!is.na(wt18)]
  w_old <- wt_auc(dw$phat_old, dw$event, dw$wt18)
  w_int <- wt_auc(dw$phat_int, dw$event, dw$wt18)
  w_slope <- wt_auc(dw$phat_slope, dw$event, dw$wt18)
  w_new <- wt_auc(dw$phat_new, dw$event, dw$wt18)
  say("\n加权敏感性（2018 基线权重，可用 n=", nrow(dw), "）：原样=", round(w_old, 3),
      "｜仅截距重校准=", round(w_int, 3), "｜斜率截距重校准=", round(w_slope, 3),
      "｜全重估=", round(w_new, 3))

  th <- seq(0.02, 0.30, 0.02)
  DCAtab <- rbind(DCAtab, data.table(口径 = dfn$lab, 阈值 = th,
                                     冻结原样 = nb_curve(d$phat_old, d$event, th),
                                     冻结仅截距 = nb_curve(d$phat_int, d$event, th),
                                     冻结斜率截距 = nb_curve(d$phat_slope, d$event, th),
                                     全重估 = nb_curve(d$phat_new, d$event, th)))
  stAll <- rbind(stAll, stk)
  nr <- data.table(
    口径 = dfn$lab,
    权重 = c("① 冻结原样（2011 拟合）", "② 冻结＋仅重校准截距", "③ 冻结＋斜率与截距重校准",
             "④ 全重估（2018-2020 拟合，表观）", "④ 全重估（Harrell 校正）"),
    AUC = c(rowA$AUC, rowC$AUC, rowD$AUC, rowB$AUC, round(auc_new_corrected, 3)),
    校准斜率 = c(rowA$校准斜率, rowC$校准斜率, rowD$校准斜率, rowB$校准斜率, NA_real_),
    校准截距 = c(rowA$校准截距, rowC$校准截距, rowD$校准截距, rowB$校准截距, NA_real_),
    Brier = c(rowA$Brier, rowC$Brier, rowD$Brier, rowB$Brier, NA_real_),
    加权AUC = c(round(w_old, 3), round(w_int, 3), round(w_slope, 3), round(w_new, 3), NA_real_),
    dAUC = c(NA_real_, 0, round(rowD$AUC - rowA$AUC, 4), round(rowB$AUC - rowA$AUC, 4), NA_real_),
    DeLong_p = c(NA_real_, NA_real_, NA_real_, round(pdel, 4), NA_real_),
    dAUC_boot低 = c(NA_real_, NA_real_, NA_real_, round(quantile(dA, .025, na.rm = TRUE), 4), NA_real_),
    dAUC_boot高 = c(NA_real_, NA_real_, NA_real_, round(quantile(dA, .975, na.rm = TRUE), 4), NA_real_))
  dtab <- rbind(dtab, nr)
}
fwrite(dtab, file.path(OUT, "F03_B_headtohead.csv"), bom = TRUE)
fwrite(stAll, file.path(OUT, "F03_B_integer_strata.csv"), bom = TRUE)
fwrite(DCAtab, file.path(OUT, "F03_B_dca.csv"), bom = TRUE)

## ---------------- 图 ----------------
d <- copy(m[mem18 > defs$D1$k18])[, event := as.integer(mem20 <= defs$D1$k20)]
d[, phat_old := pred_old(d)]
d[, phat_new := predict(fit_new(d), type = "response")]
d[, lp_old := qlogis(pmin(pmax(phat_old, 1e-9), 1 - 1e-9))]
d[, phat_int := plogis(lp_old + coef(glm(event ~ 1, offset = lp_old, family = binomial()))[1])]

png(file.path(OUT, "F03_fig_DCA.png"), width = 1100, height = 850, res = 120)
th <- seq(0, 0.4, 0.01)
plot(th, nb_curve(d$phat_old, d$event, th), type = "l", lwd = 2, col = "grey40", lty = 2,
     xlab = "Threshold probability", ylab = "Net benefit", ylim = c(-0.02, 0.10),
     main = "DCA: frozen (2011) vs refit (2018-2020) weights")
lines(th, nb_curve(d$phat_int, d$event, th), lwd = 2, col = "black")
lines(th, nb_curve(d$phat_new, d$event, th), lwd = 2, col = "red", lty = 1)
abline(h = 0, lty = 3)
legend("topright", c("Frozen 2011 weights (as published)", "Frozen + intercept recalibration",
                     "Refit 2018-2020 weights", "Treat none"),
       col = c("grey40", "black", "red", "black"), lty = c(2, 1, 1, 3), lwd = c(2, 2, 2, 1), bty = "n")
dev.off()

png(file.path(OUT, "F03_fig_calibration.png"), width = 1500, height = 520, res = 120)
par(mfrow = c(1, 3))
PVS  <- c(phat_old = "1) Frozen 2011 (as published)", phat_int = "2) Frozen + intercept recalibration",
          phat_new = "4) Refit 2018-2020")
for (pv in names(PVS)) {
  br <- unique(quantile(d[[pv]], seq(0, 1, 0.1)))
  d[, grp := cut(d[[pv]], br, include.lowest = TRUE)]
  agg <- d[, .(pred = mean(d[[pv]]), obs = mean(event), n = .N), by = grp][order(pred)]
  plot(agg$pred, agg$obs, pch = 19, cex = 1.2, xlim = c(0, 0.20), ylim = c(0, 0.20),
       xlab = "Predicted risk", ylab = "Observed risk", main = PVS[[pv]])
  abline(0, 1, lty = 3)
  text(agg$pred, agg$obs, labels = agg$n, pos = 3, cex = 0.7)
}
dev.off()

close(con)
cat("\n完成 ->", file.path(OUT, "F03_report.txt"), "\n")
