invisible(try(Sys.setlocale("LC_ALL", "zh_CN.UTF-8"), silent = TRUE))
OUT <- "C:/y03tmp/out"
a1 <- read.csv(file.path(OUT, "F15_auc.csv"), fileEncoding = "UTF-8")
a2 <- read.csv(file.path(OUT, "F16_auc.csv"), fileEncoding = "UTF-8")
sel1U <- a1[["definition"]] == "U_full_population" & a1[["predictors"]] == "three item scorecard"
sel1N <- a1[["definition"]] == "N_unimpaired_at_baseline" & a1[["predictors"]] == "three item scorecard"
sel2U <- a2[["definition"]] == "U_full_population" & a2[["predictors"]] == "three item scorecard"
sel2N <- a2[["definition"]] == "N_unimpaired_at_baseline" & a2[["predictors"]] == "three item scorecard"
p <- rbind(
  data.frame(cohort = a1[["window"]][sel1U],
             U = a1[["AUC"]][sel1U], N = a1[["AUC"]][sel1N]),
  data.frame(cohort = a2[["cohort"]][sel2U],
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
              ylab = "AUC(4 year low memory score)", xlab = "",
              main = "broader population definition vs current(unimpaired at baseline)definition: consistent direction in four cohorts")
text(bp, m, labels = sprintf("%.3f", m), pos = 3, cex = 0.85, xpd = NA)
abline(h = c(0.5, 0.7), lty = c(3, 2), col = c("black", "firebrick"))
legend("bottomright", legend = c("current definition: unimpaired at baseline", "broader population definition: no baseline memory restriction"),
       fill = c("grey65", "steelblue"), border = NA, bty = "n", cex = 0.95)
dev.off()
print(p)
