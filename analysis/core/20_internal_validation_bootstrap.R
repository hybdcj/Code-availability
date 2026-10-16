# ============================================================================
# 20_internal_validation_bootstrap.R —— A08：Bootstrap 内部验证与乐观度校正
# 依据：SAP_v3 §8（Bootstrap 1,000 次；乐观度校正 AUC、校准斜率、Brier）
# 模型：incident_2015 ~ age + edu_group + dep + adl（选项 2 主模型）
# 方法：Harrell 乐观度校正（bootstrap 内表现 − bootstrap 模型在原始样本表现）
# 生成：2026-09-26
# ============================================================================
invisible(try(Sys.setlocale("LC_ALL", "zh_CN.UTF-8"), silent = TRUE))
.libPaths(c("D:/R/library", .libPaths()))
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(pROC) })
P <- "D:/BaiduSyncdisk/科研/徐中豪"; OUT <- file.path(P, "output", "04_纵向队列")
DIR <- file.path(P, "output", "10_筛查指数")
set.seed(20260926)

coh <- readRDS(file.path(OUT, "cohort_change_2026-09-13.rds"))$cohort
w <- coh %>% select(KEY, year, imrc, dlrc) %>% pivot_wider(names_from = year, values_from = c(imrc, dlrc)) %>%
  mutate(mem11 = imrc_2011 + dlrc_2011, mem15 = imrc_2015 + dlrc_2015)
b <- coh %>% filter(year == 2011) %>% transmute(KEY, age, edu_group, dep = as.integer(cesd10 >= 10), adl) %>% left_join(w, by = "KEY")
b$edu_group <- factor(b$edu_group, levels = c("Less than primary","Primary school","Middle school","High school or above"))
d <- b[!is.na(b$mem11) & !is.na(b$mem15) & b$mem11 > 4 & complete.cases(b[, c("age","edu_group","dep","adl")]), ]
d$inc <- as.integer(d$mem15 <= 4)
cat("开发样本 n =", nrow(d), "｜事件 =", sum(d$inc), "\n")

f <- inc ~ age + edu_group + dep + adl
perf <- function(data, lp) {                     # lp = 线性预测子（已按模型计算）
  p <- 1/(1+exp(-lp))
  auc <- as.numeric(auc(roc(data$inc, p, quiet = TRUE)))
  cal <- glm(data$inc ~ lp, family = binomial())$coef          # 截距 + 斜率（斜率为第 2 个）
  brier <- mean((p - data$inc)^2)
  c(auc = auc, slope = as.numeric(cal[2]), intercept = as.numeric(cal[1]), brier = brier)
}
m_full <- glm(f, data = d, family = binomial())
lp_full <- predict(m_full, type = "link")
app <- perf(d, lp_full)
cat("\n表观表现：", sprintf("AUC=%.4f 校准斜率=%.3f Brier=%.4f", app["auc"], app["slope"], app["brier"]), "\n")

B <- 1000; n <- nrow(d); opt <- matrix(NA_real_, B, 4, dimnames = list(NULL, names(app)))
t0 <- Sys.time()
for (i in 1:B) {
  idx <- sample.int(n, n, replace = TRUE)
  db <- d[idx, ]
  mb <- try(glm(f, data = db, family = binomial()), silent = TRUE)
  if (inherits(mb, "try-error")) next
  lp_b_app  <- predict(mb, newdata = db, type = "link")
  lp_b_test <- predict(mb, newdata = d,  type = "link")
  pb_app  <- perf(db, lp_b_app)
  pb_test <- perf(d,  lp_b_test)
  opt[i, ] <- pb_app - pb_test
}
opt_m <- colMeans(opt, na.rm = TRUE)
corr <- app - opt_m
cat("\n乐观度（平均，B =", sum(!is.na(opt[,1])), "）：", sprintf("AUC %.4f｜斜率 %.4f｜Brier %.5f", opt_m["auc"], opt_m["slope"], opt_m["brier"]), "\n")
cat("乐观度校正后：", sprintf("AUC=%.4f  校准斜率=%.3f  Brier=%.4f", corr["auc"], corr["slope"], corr["brier"]), "\n")

## 校正后 AUC 的百分位 CI（用 bootstrap 检验样本 AUC 的分布近似）
test_aucs <- app["auc"] - opt[, "auc"]
ci <- quantile(test_aucs, c(0.025, 0.975), na.rm = TRUE)
cat("校正后 AUC 95%CI（百分位法）：", sprintf("%.4f – %.4f", ci[1], ci[2]), "\n")

res <- data.frame(
  指标 = c("AUC", "校准斜率", "校准截距", "Brier"),
  表观 = round(app[c("auc","slope","intercept","brier")], 4),
  乐观度 = round(opt_m[c("auc","slope","intercept","brier")], 4),
  校正后 = round(corr[c("auc","slope","intercept","brier")], 4))
print(res, row.names = FALSE)
write.csv(res, file.path(DIR, "A08_bootstrap_optimism_2026-09-26.csv"), row.names = FALSE, fileEncoding = "UTF-8")

md <- c("# A08：Bootstrap 内部验证与乐观度校正（2026-09-26）", "",
  sprintf("> 模型：incident_2015 ~ age + edu_group + dep + adl（选项 2）；开发样本 n = %d，事件 = %d（%.1f%%）",
          nrow(d), sum(d$inc), 100*mean(d$inc)),
  sprintf("> Bootstrap = %d 次（Harrell 乐观度校正）；随机种子 20260926；耗时 %.1f 分钟", B, as.numeric(difftime(Sys.time(), t0, units = "mins"))), "",
  "## 结果", "",
  paste(capture.output(print(res, row.names = FALSE)), collapse = "\n"), "",
  sprintf("**校正后 AUC = %.3f（95%%CI %.3f–%.3f）**，校准斜率 = %.2f，Brier = %.4f。",
          corr["auc"], ci[1], ci[2], corr["slope"], corr["brier"]), "",
  "## 解读", "",
  "- 乐观度（apparent − corrected）与模型复杂度相关：本模型仅 6 个参数、事件数充足（658），乐观度应较小；",
  "- 校准斜率 < 1 提示存在过拟合/预测分布过窄，A09 将据此做校准（等渗或 Platt）；",
  "- 校正后 AUC 与表观 AUC 的差值即为需要报告的乐观度。")
writeLines(md, file.path(DIR, "A08_bootstrap_optimism_2026-09-26.md"), useBytes = TRUE)
saveRDS(list(apparent = app, optimism = opt_m, corrected = corr, ci_auc = ci, opt_matrix = opt),
        file.path(DIR, "A08_bootstrap_2026-09-26.rds"))
cat("\n已保存至", DIR, "\n")
