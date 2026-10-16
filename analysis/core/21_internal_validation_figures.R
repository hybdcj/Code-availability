# ============================================================================
# 21_internal_validation_figures.R —— A09：内部验证图表（ROC/校准/DCA/分层）
# 依据：SAP_v3 §8（校准图、DCA 阈值 5–40%）与 §11（Figure 3/4）
# 模型：incident_2015 ~ age + edu_group + dep + adl（选项 2）
# 生成：2026-09-26 ｜ 输出：output/09_图表/Figure3_*、Figure4_*（png 300dpi + tif 600dpi）
#        附表：output/10_筛查指数/A09_*（校准分位表、DCA 表）
# ============================================================================
invisible(try(Sys.setlocale("LC_ALL", "zh_CN.UTF-8"), silent = TRUE))
.libPaths(c("D:/R/library", .libPaths()))
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(pROC) })
P <- "D:/BaiduSyncdisk/科研/徐中豪"; OUT <- file.path(P, "output", "04_纵向队列")
DIR <- file.path(P, "output", "10_筛查指数"); FIG <- file.path(P, "output", "09_图表")
if (!dir.exists(FIG)) dir.create(FIG, recursive = TRUE)

coh <- readRDS(file.path(OUT, "cohort_change_2026-09-13.rds"))$cohort
w <- coh %>% select(KEY, year, imrc, dlrc) %>% pivot_wider(names_from = year, values_from = c(imrc, dlrc)) %>%
  mutate(mem11 = imrc_2011 + dlrc_2011, mem15 = imrc_2015 + dlrc_2015)
b <- coh %>% filter(year == 2011) %>% transmute(KEY, age, edu_group, dep = as.integer(cesd10 >= 10), adl) %>% left_join(w, by = "KEY")
b$edu_group <- factor(b$edu_group, levels = c("Less than primary","Primary school","Middle school","High school or above"))
d <- b[!is.na(b$mem11) & !is.na(b$mem15) & b$mem11 > 4 & complete.cases(b[, c("age","edu_group","dep","adl")]), ]
d$inc <- as.integer(d$mem15 <= 4)
m <- glm(inc ~ age + edu_group + dep + adl, data = d, family = binomial())
d$p <- predict(m, type = "response")
lp <- qlogis(pmin(pmax(d$p, 1e-6), 1 - 1e-6))
cal <- glm(inc ~ lp, data = d, family = binomial())
slope <- coef(cal)[2]; intercept <- coef(cal)[1]
auc_v <- as.numeric(auc(roc(d$inc, d$p, quiet = TRUE)))

## ---------- 校准分位表（10 分位）----------
d$dec <- cut(d$p, quantile(d$p, seq(0, 1, 0.1)), include.lowest = TRUE, labels = FALSE)
cal_tab <- d %>% group_by(dec) %>% summarise(n = n(), 预测均值 = mean(p), 事件 = sum(inc),
  观测率 = mean(inc), .groups = "drop") %>%
  mutate(观测CI低 = mapply(function(k, n) binom.test(k, n)$conf.int[1], 事件, n),
         观测CI高 = mapply(function(k, n) binom.test(k, n)$conf.int[2], 事件, n))
write.csv(cal_tab, file.path(DIR, "A09_calibration_deciles_2026-09-26.csv"), row.names = FALSE, fileEncoding = "UTF-8")

## ---------- DCA（净获益，阈值 5–40%）----------
prev <- mean(d$inc); n <- nrow(d)
dca_tab <- do.call(rbind, lapply(seq(0.05, 0.40, 0.01), function(pt) {
  pos <- d$p >= pt
  TP <- sum(pos & d$inc == 1); FP <- sum(pos & d$inc == 0)
  nb_model <- TP/n - FP/n * (pt/(1-pt))
  nb_all   <- prev - (1-prev) * (pt/(1-pt))
  data.frame(阈值 = pt, 阳性率 = mean(pos), 净获益_模型 = nb_model, 净获益_全部干预 = nb_all, 净获益_不干预 = 0)
}))
write.csv(dca_tab, file.path(DIR, "A09_dca_2026-09-26.csv"), row.names = FALSE, fileEncoding = "UTF-8")

## ---------- 风险分层（**整数总分**，与工具一致：SAP v3.1 规则）----------
bb <- coef(m); b_age5 <- 5*bb[["age"]]
tb <- c(年龄每5岁 = b_age5, 教育小学 = bb[["edu_groupPrimary school"]], 教育初中 = bb[["edu_groupMiddle school"]],
        教育高中 = bb[["edu_groupHigh school or above"]], 抑郁 = bb[["dep"]], ADL = bb[["adl"]])
pts <- round(10*tb/max(abs(tb)))
d$score <- as.integer(pts[["年龄每5岁"]])*pmax(floor((d$age-45)/5),0) +
  as.integer(pts[["教育小学"]])*as.integer(d$edu_group=="Primary school") +
  as.integer(pts[["教育初中"]])*as.integer(d$edu_group=="Middle school") +
  as.integer(pts[["教育高中"]])*as.integer(d$edu_group=="High school or above") +
  as.integer(pts[["抑郁"]])*d$dep + as.integer(pts[["ADL"]])*d$adl
q <- quantile(d$score, c(1/3, 2/3))   # 整数总分三分位（与 A07 一致）
d$stratum <- cut(d$score, c(-Inf, q[1], q[2], Inf), labels = c("低","中","高"), right = TRUE)
strata <- d %>% group_by(stratum) %>% summarise(n = n(), 事件 = sum(inc), 观测率 = mean(inc),
  预测率 = mean(p), 平均分 = round(mean(score),1), .groups = "drop")
cat("整数总分三分位切点：", q, "\n")
write.csv(strata, file.path(DIR, "A09_strata_observed_vs_predicted_2026-09-26.csv"), row.names = FALSE, fileEncoding = "UTF-8")

## ---------- Figure 4：ROC / 校准 / DCA ----------
draw4 <- function() {
  par(mfrow = c(1, 3), mar = c(5, 5, 3, 2), cex.lab = 1.2, cex.axis = 1.0)
  # (a) ROC
  ro <- roc(d$inc, d$p, quiet = TRUE)
  plot(1 - ro$specificities, ro$sensitivities, type = "l", lwd = 2.5, col = "#2C7FB8",
       xlab = "1 − 特异度", ylab = "灵敏度", main = "(a) ROC 曲线", xlim = c(0,1), ylim = c(0,1))
  abline(0, 1, lty = 3, col = "grey50")
  text(0.62, 0.18, sprintf("AUC = %.3f\n(95%%CI %.3f–%.3f)", auc_v,
       as.numeric(ci.auc(ro)[1]), as.numeric(ci.auc(ro)[3])), cex = 1.0)
  # (b) 校准图
  plot(cal_tab$预测均值, cal_tab$观测率, pch = 19, cex = 1.2, col = "#31A354",
       xlim = c(0, max(cal_tab$预测均值, cal_tab$观测CI高)*1.05),
       ylim = c(0, max(cal_tab$观测CI高)*1.05),
       xlab = "预测概率", ylab = "观测事件率", main = "(b) 校准图（十分位）")
  segments(cal_tab$预测均值, cal_tab$观测CI低, cal_tab$预测均值, cal_tab$观测CI高, col = "#31A354")
  abline(0, 1, lty = 3, col = "grey50")
  lines(lowess(cal_tab$预测均值, cal_tab$观测率), col = "#756BB1", lwd = 1.5)
  text(max(cal_tab$预测均值)*0.45, max(cal_tab$观测CI高)*0.95,
       sprintf("校准斜率 = %.2f\n截距 = %.2f", slope, intercept), cex = 1.0, adj = 0)
  # (c) DCA
  plot(dca_tab$阈值, dca_tab$净获益_模型, type = "l", lwd = 2.5, col = "#2C7FB8",
       ylim = c(min(dca_tab$净获益_模型, -0.02), max(dca_tab$净获益_模型, prev)*1.05),
       xlab = "阈值概率", ylab = "净获益", main = "(c) 决策曲线分析")
  lines(dca_tab$阈值, dca_tab$净获益_全部干预, lty = 2, col = "#E6550D", lwd = 2)
  abline(h = 0, lty = 3, col = "grey50")
  legend("topright", c("筛查指数", "全部干预", "不干预"), lty = c(1,2,3), lwd = c(2.5,2,1),
         col = c("#2C7FB8","#E6550D","grey50"), bty = "n", cex = 0.95)
}
png(file.path(FIG, "Figure4_internal_validation_2026-09-26.png"), width = 18, height = 6, units = "in", res = 300)
draw4(); dev.off()
tiff(file.path(FIG, "Figure4_internal_validation_2026-09-26.tif"), width = 18, height = 6, units = "in", res = 600, compression = "lzw")
draw4(); dev.off()

## ---------- Figure 3：风险分层 ----------
draw3 <- function() {
  par(mar = c(5, 5, 3, 2), cex.lab = 1.2)
  bp <- barplot(strata$观测率, names.arg = c("低危","中危","高危"), col = c("#C6DBEF","#6BAED6","#2171B5"),
                ylim = c(0, max(strata$观测率)*1.3), ylab = "4 年内新发认知受损率", xlab = "筛查指数风险层",
                main = "筛查指数风险分层（开发队列）")
  text(bp, strata$观测率, sprintf("%.1f%%\n(n=%d)", 100*strata$观测率, strata$n), pos = 3, cex = 1.0)
  points(bp, strata$预测率, pch = 18, cex = 1.6, col = "red")
  legend("topleft", c("观测事件率", "模型预测率"), pch = c(15, 18), col = c("#2171B5","red"), bty = "n")
}
png(file.path(FIG, "Figure3_risk_stratification_2026-09-26.png"), width = 8, height = 6, units = "in", res = 300)
draw3(); dev.off()
tiff(file.path(FIG, "Figure3_risk_stratification_2026-09-26.tif"), width = 8, height = 6, units = "in", res = 600, compression = "lzw")
draw3(); dev.off()

cat("== 校准十分位表 ==\n"); print(as.data.frame(cal_tab), row.names = FALSE)
cat("\n== DCA 关键阈值 ==\n"); print(dca_tab[seq(1, 36, 5), ], row.names = FALSE)
cat("\n== 分层（观测 vs 预测）==\n"); print(as.data.frame(strata), row.names = FALSE)
cat("\n图已保存：", FIG, "\n")
