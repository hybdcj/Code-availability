# F26 LEGACY: BDSI implementation/input is a historical version; for corrected analyses run F26_bdsi_audit.R.
# original numerical results retained for traceability, cannot directly serve as original model BDSI valid comparison; Tu see also limitations in F25.
# F22b: three arm comparisons under two thresholds AUC / paired test / stratification gradient / detection efficiency
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
say("three arm data n =", nrow(d))

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
    res <- rbind(res, data.table(threshold = lab, model = ARMS[[v]], n = nrow(x),
                                 events = sum(x$event), event_rate = round(100 * mean(x$event), 2),
                                 AUC = round(ci[2], 4), low = round(ci[1], 4), high = round(ci[3], 4)))
  }
  for (v in c("score_bdsi", "lp_tu")) {
    tt <- roc.test(rocs[["score_us"]], rocs[[v]], paired = TRUE)
    cm <- rbind(cm, data.table(threshold = lab, comparison = ARMS[[v]],
                this_study_AUC = round(as.numeric(auc(rocs[["score_us"]])), 4),
                comparison_AUC = round(as.numeric(auc(rocs[[v]])), 4),
                dAUC = round(as.numeric(auc(rocs[["score_us"]])) - as.numeric(auc(rocs[[v]])), 4),
                DeLong_P = signif(tt$p.value, 3)))
  }
  list(res = res, cmp = cm)
}
r1 <- runset(d1, "T1 stratification")
r2 <- runset(d2, "T2 fixed 4")

say("\n===== 1. each arm AUC =====")
pt(rbind(r1$res, r2$res))
fwrite(rbind(r1$res, r2$res), file.path(OUT, "F22_auc.csv"), bom = TRUE)

say("\n===== 2. paired comparison =====")
pt(rbind(r1$cmp, r2$cmp))
fwrite(rbind(r1$cmp, r2$cmp), file.path(OUT, "F22_pairwise.csv"), bom = TRUE)

say("\n===== 3. quartile stratification gradient =====")
for (nm in c("T1 stratification", "T2 fixed 4")) {
  x <- if (nm == "T1 stratification") d1 else d2
  for (v in names(ARMS)) {
    q <- cut(x[[v]], quantile(x[[v]], seq(0, 1, .25)), include.lowest = TRUE, labels = 1:4)
    g1 <- x[, .(r = paste0(round(100 * mean(event), 1), "%")), by = q][order(q)]$r
    say("  ", nm, " | ", ARMS[[v]], " : ", paste(g1, collapse = " -> "))
  }
}

say("\n===== 4. target screening 30%threshold efficiency of(including ties) =====")
op <- data.table()
for (nm in c("T1 stratification", "T2 fixed 4")) {
  x <- if (nm == "T1 stratification") d1 else d2
  for (v in names(ARMS)) {
    s <- x[[v]]; ev <- x$event; n <- length(ev)
    k <- floor(0.30 * n); th <- sort(s, decreasing = TRUE)[k]; pos <- s >= th
    op <- rbind(op, data.table(threshold = nm, model = ARMS[[v]], cutoff = th, target_screening_fraction = "30%", n = length(ev), events = sum(ev), cases_detected = sum(pos), actual_screening_fraction = round(100 * mean(pos), 2), TP = sum(pos & ev == 1), FP = sum(pos & ev == 0), TN = sum(!pos & ev == 0), FN = sum(!pos & ev == 1),
                 sensitivity = round(100 * sum(ev[pos]) / sum(ev), 1),
                 specificity = round(100 * sum(!pos & ev == 0) / sum(ev == 0), 1),
                 negative_predictive_value = round(100 * mean(ev[!pos] == 0), 1), positive_predictive_value = round(100 * mean(ev[pos]), 1)))
  }
}
pt(op); fwrite(op, file.path(OUT, "F22_operating.csv"), bom = TRUE)
close(con)
cat("\n complete\n")
