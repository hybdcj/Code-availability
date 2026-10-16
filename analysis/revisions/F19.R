# ============================================================================
# F19：外部队列全人群口径 × 外部S2敏感性阈值（固定 CHARLS <=4）
# 目的：补充外部固定阈值敏感性及全人群事后分析；不替换S0主分析
# 校验：N 口径（基线 >4）应复现 A11 S2 —— HRS 15521/945/0.7475、
#       ELSA 7203/324/0.8025、SHARE 28391/2256/0.7769
# ============================================================================
invisible(try(Sys.setlocale("LC_ALL", "zh_CN.UTF-8"), silent = TRUE))
suppressPackageStartupMessages({ library(data.table); library(pROC) })
if (dir.exists("D:/R/library")) .libPaths(c("D:/R/library", .libPaths()))
TMP <- "C:/y03tmp/tmp"; OUT <- "C:/y03tmp/out"
con <- file(file.path(OUT, "F19_report.txt"), open = "wt", encoding = "UTF-8")
say <- function(...) { s <- paste0(..., collapse = ""); cat(s, "\n"); cat(s, "\n", file = con) }
pt <- function(x) say(paste(capture.output(print(x)), collapse = "\n"))

B <- c(intercept = -4.19524508778947, age = 0.0492039250378111,
       eduPri = -0.290105771000328, eduMid = -0.994711800202235,
       eduHigh = -1.43588388035287, dep = 0.412166484954086, adl = -0.0110119344848014)
CO <- c("HRS", "ELSA", "SHARE")
K <- 4L

prep <- function(co) {
  d <- as.data.table(readRDS(file.path(TMP, paste0("F16_base_", co, ".rds"))))
  d <- d[!is.na(age) & age >= 45]
  d[, score := 2L * pmax(floor((age - 45) / 5), 0) + edu_pts + 3L * dep]
  d[, lp := B[["intercept"]] + B[["age"]] * age + B[["dep"]] * dep + B[["adl"]] * adl +
      B[["eduPri"]] * as.integer(edu_pts == -2) +
      B[["eduMid"]] * as.integer(edu_pts == -7) +
      B[["eduHigh"]] * as.integer(edu_pts == -10)]
  d[, phat := plogis(lp)]
  d[, event := as.integer(mem_fu <= K)]
  d[, low0 := as.integer(mem_base <= K)]
  d[]
}
D <- lapply(CO, prep); names(D) <- CO

say("阈值：固定 CHARLS <= ", K, "（外部S2敏感性口径；S0仍为主分析）")
say("\n========== 1. 管线校验（复现 A11 S2 的 N 口径）==========")
say("期望：HRS 15521 / 945 / AUC 0.7475；ELSA 7203 / 324 / 0.8025；SHARE 28391 / 2256 / 0.7769")
ver <- rbindlist(lapply(CO, function(co) {
  d <- D[[co]][!is.na(mem_fu) & !is.na(mem_base) & mem_base > K]
  data.table(队列 = co, n = nrow(d), 事件 = sum(d$event),
             AUC_评分卡 = round(as.numeric(auc(roc(d$event, d$score, quiet = TRUE, direction = "<"))), 4))
}))
pt(ver); fwrite(ver, file.path(OUT, "F19_verify_S2.csv"), bom = TRUE)

say("\n========== 2. 各口径 AUC（固定 <=4）==========")
aucrow <- function(d, v, co, arm, vlab) {
  ok <- !is.na(d[[v]]) & !is.na(d$event)
  if (sum(ok) < 30 || length(unique(d$event[ok])) < 2) return(NULL)
  r <- roc(d$event[ok], d[[v]][ok], quiet = TRUE, direction = "<")
  ci <- as.numeric(ci.auc(r))
  data.table(队列 = co, 口径 = arm, 预测子 = vlab, n = sum(ok), 事件 = sum(d$event[ok]),
             事件率 = round(100 * mean(d$event[ok]), 2),
             AUC = round(ci[2], 4), AUC低 = round(ci[1], 4), AUC高 = round(ci[3], 4))
}
mk <- function(d) list(U = d, N = d[!is.na(mem_base) & mem_base > K],
                       I = d[!is.na(mem_base) & mem_base <= K], M = d[is.na(mem_base)])
LB <- c(U = "U_全人口", N = "N_基线未受损", I = "I_基线已受损", M = "M_基线记忆缺失")
aucT <- rbindlist(lapply(CO, function(co) {
  st <- mk(D[[co]][!is.na(mem_fu)])
  rbindlist(lapply(names(st), function(s) {
    if (nrow(st[[s]]) == 0) return(NULL)
    rbind(aucrow(st[[s]], "score", co, LB[[s]], "三项评分卡"),
          aucrow(st[[s]], "phat", co, LB[[s]], "连续模型"))
  }))
}))
pt(aucT); fwrite(aucT, file.path(OUT, "F19_auc.csv"), bom = TRUE)

say("\n========== 3. 分解 ==========")
dec <- rbindlist(lapply(CO, function(co) {
  d <- D[[co]][!is.na(mem_fu) & !is.na(mem_base)]
  r_low <- roc(d$event, d$low0, quiet = TRUE, direction = "<")
  r_scr <- roc(d$event, d$score, quiet = TRUE, direction = "<")
  fit <- glm(event ~ score + low0, data = d, family = binomial())
  r_cmb <- roc(d$event, predict(fit, type = "response"), quiet = TRUE, direction = "<")
  data.table(队列 = co, n = nrow(d), 事件 = sum(d$event),
             事件中基线已低分占比 = round(100 * mean(d$low0[d$event == 1]), 1),
             基线已低分者事件率 = round(100 * mean(d[low0 == 1]$event), 1),
             基线未受损者事件率 = round(100 * mean(d[low0 == 0]$event), 1),
             AUC_仅基线已低分 = round(as.numeric(auc(r_low)), 4),
             AUC_仅评分卡 = round(as.numeric(auc(r_scr)), 4),
             AUC_合并 = round(as.numeric(auc(r_cmb)), 4),
             P_合并优于评分卡 = signif(roc.test(r_cmb, r_scr, paired = TRUE)$p.value, 3))
}))
pt(dec); fwrite(dec, file.path(OUT, "F19_decomposition.csv"), bom = TRUE)

say("\n========== 4. 全人群筛查模拟 ==========")
scr <- rbindlist(lapply(CO, function(co) {
  d <- D[[co]][!is.na(mem_fu)]
  s <- d$score; ev <- d$event; n <- length(ev)
  rbindlist(lapply(c(0.10, 0.20, 0.30), function(fr) {
    k <- floor(fr * n); th <- sort(s, decreasing = TRUE)[k]; pos <- s >= th
    sub <- d[pos]; yld <- sum(sub$event == 1 & sub$low0 == 0, na.rm = TRUE)
    data.table(队列 = co, 目标筛查比例 = paste0(fr * 100, "%"), 阈值 = th, n = length(ev), 事件 = sum(ev), 检出人数 = sum(pos), 实际筛查比例 = round(100 * mean(pos), 2), TP = sum(pos & ev == 1), FP = sum(pos & ev == 0), TN = sum(!pos & ev == 0), FN = sum(!pos & ev == 1),
               灵敏度 = round(100 * sum(ev[pos]) / sum(ev), 1),
               特异度 = round(100 * sum(!pos & ev == 0) / sum(ev == 0), 1),
               阴性预测值 = round(100 * mean(ev[!pos] == 0), 1), 阳性预测值 = round(100 * mean(ev[pos]), 1),
               检出者中基线已低分占比 = round(100 * mean(sub$low0 == 1, na.rm = TRUE), 1),
               新发真事件数 = yld, 每检出1例新发需筛人数 = if (yld > 0) round(n / yld, 1) else NA_real_)
  }))
}))
pt(scr); fwrite(scr, file.path(OUT, "F19_screening.csv"), bom = TRUE)

say("\n========== 5. 与主稿口径的衔接总表 ==========")
pt(rbindlist(lapply(CO, function(co) {
  d0 <- D[[co]][!is.na(mem_fu)]; d1 <- d0[!is.na(mem_base) & mem_base > K]
  data.table(队列 = co,
             N_基线未受损_n = nrow(d1), N_事件 = sum(d1$event),
             N_AUC = round(as.numeric(auc(roc(d1$event, d1$score, quiet = TRUE, direction = "<"))), 4),
             U_全人口_n = nrow(d0), U_事件 = sum(d0$event),
             U_AUC = round(as.numeric(auc(roc(d0$event, d0$score, quiet = TRUE, direction = "<"))), 4))
})))

close(con)
cat("\n完成 ->", file.path(OUT, "F19_report.txt"), "\n")
