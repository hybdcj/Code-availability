invisible(try(Sys.setlocale("LC_ALL", "zh_CN.UTF-8"), silent = TRUE))
OUT <- "C:/y03tmp/out"
a1 <- read.csv(file.path(OUT, "F15_auc.csv"), fileEncoding = "UTF-8")
a2 <- read.csv(file.path(OUT, "F16_auc.csv"), fileEncoding = "UTF-8")
sel1U <- a1[["口径"]] == "U_全人口" & a1[["预测子"]] == "三项评分卡"
sel1N <- a1[["口径"]] == "N_基线未受损" & a1[["预测子"]] == "三项评分卡"
sel2U <- a2[["口径"]] == "U_全人口" & a2[["预测子"]] == "三项评分卡"
sel2N <- a2[["口径"]] == "N_基线未受损" & a2[["预测子"]] == "三项评分卡"
p <- rbind(
  data.frame(cohort = a1[["窗口"]][sel1U],
             U = a1[["AUC"]][sel1U], N = a1[["AUC"]][sel1N]),
  data.frame(cohort = a2[["队列"]][sel2U],
             U = a2[["AUC"]][sel2U], N = a2[["AUC"]][sel2N]))
lab <- c("CHARLS 2011-2015", "CHARLS 2015-2018", "HRS", "ELSA", "SHARE")
p$cohort <- lab[seq_len(nrow(p))]
rownames(p) <- NULL
write.csv(p, file.path(OUT, "F16_fig_data.csv"), row.names = FALSE, fileEncoding = "UTF-8")
png(file.path(OUT, "F16_fig_wide_vs_restricted.png"), width = 1500, height = 780, res = 150)
par(mar = c(5.5, 5, 3.5, 1.5))
m <- rbind(p$N, p$U)
bp <- barplot(m, beside = TRUE, names.arg = p$cohort, ylim = c(0, 1), xpd = FALSE,
              col = c("grey65", "steelblue"), border = NA,
              ylab = "AUC（4 年记忆低分）", xlab = "",
              main = "全人群口径 vs 现行（基线未受损）口径：四个队列方向一致")
text(bp, m, labels = sprintf("%.3f", m), pos = 3, cex = 0.85, xpd = NA)
abline(h = c(0.5, 0.7), lty = c(3, 2), col = c("black", "firebrick"))
legend("bottomright", legend = c("现行口径：基线未受损", "全人群口径：不设基线记忆限制"),
       fill = c("grey65", "steelblue"), border = NA, bty = "n", cex = 0.95)
dev.off()
print(p)
