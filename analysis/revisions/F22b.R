# F26 LEGACY: BDSI实现/输入为历史版本；修正分析请运行F26_bdsi_audit.R。
# 原结果数值保留追溯，不可直接作为原版BDSI的有效比较；Tu限制另见F25。
# F22b：两套阈值下三臂对比的 AUC / 配对检验 / 分层梯度 / 检出效率
invisible(try(Sys.setlocale("LC_ALL", "zh_CN.UTF-8"), silent = TRUE))
suppressPackageStartupMessages({ library(data.table); library(pROC) })
if (dir.exists("D:/R/library")) .libPaths(c("D:/R/library", .libPaths()))
TMP <- "C:/y03tmp/tmp"; OUT <- "C:/y03tmp/out"
con <- file(file.path(OUT, "F22_report.txt"), open = "wt", encoding = "UTF-8")
say <- function(...) { s <- paste0(..., collapse = ""); cat(s, "\n"); cat(s, "\n", file = con) }
pt <- function(x) say(paste(capture.output(print(x)), collapse = "\n"))

KFIX <- 4L
ARMS <- c(score_us = "ben", score_bdsi = "BDSI", lp_tu = "Tu2026")
d <- as.data.table(readRDS(file.path(TMP, "F22_d.rds")))
say("三臂数据 n =", nrow(d))

d1 <- copy(d)
d1[, event := as.integer(mem18 < cutoff)]
d1 <- d1[mem15 >= cutoff]
d2 <- copy(d)
d2[, event := as.integer(mem18 <= KFIX)]
d2 <- d2[mem15 > KFIX]

runset <- function(x, lab) {
  rocs <- list(); res <- data.table(); cm <- data.table()
  for (v in names(ARMS)) {
    r <- roc(x$event, x[[v]], quiet = TRUE, direction = "<")
    rocs[[v]] <- r
    ci <- as.numeric(ci.auc(r))
    res <- rbind(res, data.table(阈值 = lab, 模型 = ARMS[[v]], n = nrow(x),
                                 事件 = sum(x$event), 事件率 = round(100 * mean(x$event), 2),
                                 AUC = round(ci[2], 4), 低 = round(ci[1], 4), 高 = round(ci[3], 4)))
  }
  for (v in c("score_bdsi", "lp_tu")) {
    tt <- roc.test(rocs[["score_us"]], rocs[[v]], paired = TRUE)
    cm <- rbind(cm, data.table(阈值 = lab, 对照 = ARMS[[v]],
                本研究AUC = round(as.numeric(auc(rocs[["score_us"]])), 4),
                对照AUC = round(as.numeric(auc(rocs[[v]])), 4),
                dAUC = round(as.numeric(auc(rocs[["score_us"]])) - as.numeric(auc(rocs[[v]])), 4),
                DeLong_P = signif(tt$p.value, 3)))
  }
  list(res = res, cmp = cm)
}
r1 <- runset(d1, "T1分层")
r2 <- runset(d2, "T2固定4")

say("\n===== 1. 各臂 AUC =====")
pt(rbind(r1$res, r2$res))
fwrite(rbind(r1$res, r2$res), file.path(OUT, "F22_auc.csv"), bom = TRUE)

say("\n===== 2. 配对对比 =====")
pt(rbind(r1$cmp, r2$cmp))
fwrite(rbind(r1$cmp, r2$cmp), file.path(OUT, "F22_pairwise.csv"), bom = TRUE)

say("\n===== 3. 四分位分层梯度 =====")
for (nm in c("T1分层", "T2固定4")) {
  x <- if (nm == "T1分层") d1 else d2
  for (v in names(ARMS)) {
    q <- cut(x[[v]], quantile(x[[v]], seq(0, 1, .25)), include.lowest = TRUE, labels = 1:4)
    g1 <- x[, .(r = paste0(round(100 * mean(event), 1), "%")), by = q][order(q)]$r
    say("  ", nm, " | ", ARMS[[v]], " : ", paste(g1, collapse = " -> "))
  }
}

say("\n===== 4. 目标筛查30%的阈值效率（含同分者） =====")
op <- data.table()
for (nm in c("T1分层", "T2固定4")) {
  x <- if (nm == "T1分层") d1 else d2
  for (v in names(ARMS)) {
    s <- x[[v]]; ev <- x$event; n <- length(ev)
    k <- floor(0.30 * n); th <- sort(s, decreasing = TRUE)[k]; pos <- s >= th
    op <- rbind(op, data.table(阈值 = nm, 模型 = ARMS[[v]], 切点 = th, 目标筛查比例 = "30%", n = length(ev), 事件 = sum(ev), 检出人数 = sum(pos), 实际筛查比例 = round(100 * mean(pos), 2), TP = sum(pos & ev == 1), FP = sum(pos & ev == 0), TN = sum(!pos & ev == 0), FN = sum(!pos & ev == 1),
                 灵敏度 = round(100 * sum(ev[pos]) / sum(ev), 1),
                 特异度 = round(100 * sum(!pos & ev == 0) / sum(ev == 0), 1),
                 阴性预测值 = round(100 * mean(ev[!pos] == 0), 1), 阳性预测值 = round(100 * mean(ev[pos]), 1)))
  }
}
pt(op); fwrite(op, file.path(OUT, "F22_operating.csv"), bom = TRUE)
close(con)
cat("\n完成\n")
