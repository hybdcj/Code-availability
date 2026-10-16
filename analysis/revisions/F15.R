# F32（2026-10-16）：以下保留旧共同完整样本的计算以便追溯；CHARLS整数评分自身样本已改用同目录 F32_全人群入组与分母修正.R。不要仅删除lp限制后直接运行旧脚本，因为下游连续模型、分母和表头也需分别处理。
# ============================================================================
# F15  broad-population application of the frozen 3-item memory risk score
# 冻结模型 A14（2011 -> 2015 开发）：三项评分卡 + 连续模型
# 结局：随访波记忆分 < 教育分层 mean-1SD 切点（切点由 2011 波定义）
# 口径：U 全人口（不限基线记忆）| N 基线未受损（现行）| I 基线已受损
# 窗口：W1 = 2011 -> 2015（开发队列）| W2 = 2015 -> 2018（独立样本）
# ============================================================================
invisible(try(Sys.setlocale("LC_ALL", "zh_CN.UTF-8"), silent = TRUE))
suppressPackageStartupMessages({ library(data.table); library(pROC) })
if (dir.exists("D:/R/library")) .libPaths(c("D:/R/library", .libPaths()))
TMP <- "C:/y03tmp/tmp"; OUT <- "C:/y03tmp/out"
con <- file(file.path(OUT, "F15_report.txt"), open = "wt", encoding = "UTF-8")
say <- function(...) { s <- paste0(..., collapse = ""); cat(s, "\n"); cat(s, "\n", file = con) }
pt  <- function(x) say(paste(capture.output(print(x)), collapse = "\n"))

B <- c(intercept = -4.19524508778947, age = 0.0492039250378111,
       eduPri = -0.290105771000328, eduMid = -0.994711800202235,
       eduHigh = -1.43588388035287, dep = 0.412166484954086, adl = -0.0110119344848014)
LV <- c("Less than primary", "Primary school", "Middle school", "High school or above")
G3 <- c("A_none", "B_primid", "C_high")

## ---------------- 1. 数据 ---------------- 
w <- as.data.table(readRDS(file.path(TMP, "waves_cognition.rds")))
g <- w[year %in% c(2011, 2015, 2018), .(ID = IDn, year = as.numeric(year),
        mem = as.numeric(imrc) + as.numeric(dlrc),
        age = as.numeric(age), edu = edu_group, cesd = as.numeric(cesd10),
        adl = as.numeric(adl), female = as.numeric(female), wt = as.numeric(weight))]
say("五波数据可用年份：", paste(sort(unique(w$year)), collapse = ", "))
say("本次使用年份：", paste(sort(unique(g$year)), collapse = ", "),
    "| 行数 ", nrow(g))

## ---------------- 2. 教育分层切点（2011 定义，同 F14）---------------- 
d11 <- g[year == 2011 & !is.na(mem) & !is.na(edu) & age >= 50]
d11[, g3 := fcase(edu == LV[1], G3[1], edu %in% LV[2:3], G3[2], edu == LV[4], G3[3])]
cut3 <- d11[!is.na(g3), .(n = .N, mean = mean(mem), sd = sd(mem)), by = g3][order(g3)]
cut3[, cutoff := mean - sd]
say("\n========== 1. 教育分层切点（2011 波，mean-1SD）==========")
pt(cut3[, .(层 = g3, n, 均值 = round(mean, 2), SD = round(sd, 2), 切点 = round(cutoff, 2))])
fwrite(cut3, file.path(OUT, "F15_cutoffs.csv"), bom = TRUE)

## ---------------- 3. 构造窗口数据 ---------------- 
build <- function(y0, y1) {
  a <- g[year == y0, .(ID, mem0 = mem, age0 = age, edu0 = edu, cesd0 = cesd,
                       adl0 = adl, wt0 = wt, female0 = female)]
  b <- g[year == y1, .(ID, mem1 = mem)]
  d <- merge(a, b, by = "ID", all = FALSE)
  d[, g3 := fcase(edu0 == LV[1], G3[1], edu0 %in% LV[2:3], G3[2], edu0 == LV[4], G3[3])]
  d <- merge(d, cut3[, .(g3, cutoff)], by = "g3", all.x = TRUE)
  d[, edu_f := factor(edu0, levels = LV)]
  d[, score := 2L * pmax(floor((age0 - 45) / 5), 0) +
        (-2L) * as.integer(edu_f == LV[2]) + (-7L) * as.integer(edu_f == LV[3]) +
        (-10L) * as.integer(edu_f == LV[4]) + 3L * as.integer(cesd0 >= 10)]
  d[, lp := B[["intercept"]] + B[["age"]] * age0 + B[["dep"]] * as.integer(cesd0 >= 10) +
        B[["adl"]] * adl0 + B[["eduPri"]] * as.integer(edu_f == LV[2]) +
        B[["eduMid"]] * as.integer(edu_f == LV[3]) + B[["eduHigh"]] * as.integer(edu_f == LV[4])]
  d[, phat := plogis(lp)]
  d[, event := as.integer(mem1 < cutoff)]
  d[, low0 := as.integer(mem0 < cutoff)]
  d[]
}
w1 <- build(2011, 2015)
w2 <- build(2015, 2018)

# 纳入条件同 y04 全国加权分布：年龄 >= 45（原始文件含配偶/子女行，年龄可低至 10 岁）
core <- function(d) !is.na(d$score) & !is.na(d$lp) & !is.na(d$event) & !is.na(d$age0) & d$age0 >= 45
getc <- function(nm) { d <- get(nm, envir = .GlobalEnv); d[core(d)] }

## ---------------- 4. 样本流程 ---------------- 
flow <- data.table()
for (nm in c("w1", "w2")) {
  d <- get(nm)
  lab <- ifelse(nm == "w1", "2011 -> 2015（开发队列）", "2015 -> 2018（独立样本）")
  ok <- core(d)
  flow <- rbind(flow, data.table(
    窗口 = lab,
    两波ID可匹配 = nrow(d),
    随访记忆可计分 = sum(!is.na(d$mem1)),
    随访可计分且有切点 = sum(!is.na(d$mem1) & !is.na(d$cutoff)),
    三项评分可计算 = sum(ok),
    因年龄小于45岁排除 = sum(!is.na(d$score) & !is.na(d$lp) & !is.na(d$event) &
                             !is.na(d$age0) & d$age0 < 45),
    年龄缺失排除 = sum(!is.na(d$score) & !is.na(d$lp) & !is.na(d$event) & is.na(d$age0)),
    基线记忆可计分 = sum(!is.na(d$mem0)),
    U_全人口 = sum(ok),
    U_其中基线记忆缺失 = sum(ok & is.na(d$mem0)),
    N_基线未受损 = sum(ok & d$low0 == 0, na.rm = TRUE),
    I_基线已受损 = sum(ok & d$low0 == 1, na.rm = TRUE)))
}
say("\n========== 2. 样本流程 ==========")
pt(flow)
fwrite(flow, file.path(OUT, "F15_flow.csv"), bom = TRUE)

## ---------------- 5. AUC 主表 ---------------- 
aucrow <- function(d, v, win, arm, vlab) {
  ok <- !is.na(d[[v]]) & !is.na(d$event)
  r <- roc(d$event[ok], d[[v]][ok], quiet = TRUE, direction = "<")
  ci <- as.numeric(ci.auc(r))
  ev <- sum(d$event[ok] == 1); n <- sum(ok)
  data.table(窗口 = win, 口径 = arm, 预测子 = vlab, n = n, 事件 = ev,
             事件率 = round(100 * ev / n, 2),
             AUC = round(ci[2], 4), AUC低 = round(ci[1], 4), AUC高 = round(ci[3], 4))
}
armset <- function(d, win, out) {
  ok <- core(d)
  U <- d[ok]
  N <- d[ok & !is.na(low0) & low0 == 0]
  I <- d[ok & !is.na(low0) & low0 == 1]
  Bm <- d[ok & !is.na(low0)]
  Mu <- d[ok & is.na(low0)]
  sets <- list(U_全人口 = U, N_基线未受损 = N, I_基线已受损 = I,
               M_基线记忆缺失 = Mu, B_基线记忆可得 = Bm)
  for (s in names(sets)) {
    if (nrow(sets[[s]]) == 0) next
    out <- rbind(out, aucrow(sets[[s]], "score", win, s, "三项评分卡"))
    out <- rbind(out, aucrow(sets[[s]], "phat", win, s, "连续模型"))
  }
  out
}
W1L <- "2011 -> 2015（开发队列）"; W2L <- "2015 -> 2018（独立样本）"
aucT <- data.table()
aucT <- armset(w1, W1L, aucT)
aucT <- armset(w2, W2L, aucT)
say("\n========== 3. 各口径 AUC（冻结模型，四项预测子）==========")
pt(aucT)
fwrite(aucT, file.path(OUT, "F15_auc.csv"), bom = TRUE)

## ---------------- 6. 分解：基线已低分本身有多强的预测力 ---------------- 
say("\n========== 4. 分解：基线已低分指示 vs 评分卡 ==========")
decT <- data.table()
for (nm in c("w1", "w2")) {
  d <- getc(nm)[!is.na(low0)]
  lab <- ifelse(nm == "w1", W1L, W2L)
  r_low <- roc(d$event, d$low0, quiet = TRUE, direction = "<")
  r_scr <- roc(d$event, d$score, quiet = TRUE, direction = "<")
  fit <- glm(event ~ score + low0, data = d, family = binomial())
  pp <- predict(fit, type = "response")
  r_cmb <- roc(d$event, pp, quiet = TRUE, direction = "<")
  dt1 <- roc.test(r_scr, r_low, paired = TRUE)
  dt2 <- roc.test(r_cmb, r_scr, paired = TRUE)
  decT <- rbind(decT, data.table(
    窗口 = lab, n = nrow(d), 事件 = sum(d$event),
    基线已低分构成的事件占比 = round(100 * mean(d$low0[d$event == 1]), 1),
    AUC_仅基线已低分 = round(as.numeric(auc(r_low)), 4),
    AUC_仅评分卡 = round(as.numeric(auc(r_scr)), 4),
    AUC_评分卡加基线已低分 = round(as.numeric(auc(r_cmb)), 4),
    P_基线指示优于评分卡 = signif(dt1$p.value, 3),
    P_合并优于评分卡 = signif(dt2$p.value, 3)))
}
pt(decT)
fwrite(decT, file.path(OUT, "F15_decomposition.csv"), bom = TRUE)

## ---------------- 7. 全人口风险分层 ---------------- 
say("\n========== 5. 全人口风险分层（评分为整分值）==========")
strT <- data.table()
for (nm in c("w1", "w2")) {
  d <- copy(getc(nm))
  lab <- ifelse(nm == "w1", W1L, W2L)
  d[, strt := cut(score, breaks = quantile(score, seq(0, 1, .25)), include.lowest = TRUE,
                  labels = c("Q1最低", "Q2", "Q3", "Q4最高"))]
  tab <- d[, .(n = .N, 事件率 = round(100 * mean(event), 1),
               其中基线已低分 = round(100 * mean(low0 == 1, na.rm = TRUE), 1)),
           by = strt][order(strt)]
  say("-- ", lab, " --"); pt(tab)
  strT <- rbind(strT, cbind(窗口 = lab, tab))
}
fwrite(strT, file.path(OUT, "F15_strata.csv"), bom = TRUE)

## ---------------- 8. 全人口筛查模拟 ---------------- 
say("\n========== 6. 全人口筛查模拟（按评分切分数比例）==========")
scrT <- data.table()
for (nm in c("w1", "w2")) {
  d <- copy(getc(nm))
  lab <- ifelse(nm == "w1", W1L, W2L)
  s <- d$score; ev <- d$event; n <- length(ev)
  for (fr in c(0.05, 0.10, 0.20, 0.30)) {
    k <- floor(fr * n); th <- sort(s, decreasing = TRUE)[k]; pos <- s >= th
    sens <- sum(ev[pos]) / sum(ev); spec <- sum(!pos & ev == 0) / sum(ev == 0)
    ppv <- ifelse(sum(pos) > 0, mean(ev[pos]), NA_real_)
    npv <- mean(ev[!pos] == 0)
    sub <- d[pos]
    yld_new <- sum(sub$event == 1 & sub$low0 == 0, na.rm = TRUE)
    scrT <- rbind(scrT, data.table(
      窗口 = lab, 目标筛查比例 = paste0(fr * 100, "%"), 阈值 = th, n = length(ev), 事件 = sum(ev), 检出人数 = sum(pos), 实际筛查比例 = round(100 * mean(pos), 2), TP = sum(pos & ev == 1), FP = sum(pos & ev == 0), TN = sum(!pos & ev == 0), FN = sum(!pos & ev == 1),
      灵敏度 = round(100 * sens, 1), 特异度 = round(100 * spec, 1),
      阳性预测值 = round(100 * ppv, 1), 阴性预测值 = round(100 * npv, 1),
      检出者中基线已低分占比 = round(100 * mean(sub$low0 == 1, na.rm = TRUE), 1),
      新发真事件数 = yld_new, 每检出1例新发需筛人数 = if (yld_new > 0) round(n / yld_new, 1) else NA_real_))
  }
  co <- coords(roc(ev, s, quiet = TRUE, direction = "<"), "best", best.method = "youden",
               ret = c("threshold", "sensitivity", "specificity"))
  co <- co[which.max(co$specificity), , drop = FALSE]
  pos <- s >= co$threshold
  sub <- d[pos]
  yld_new <- sum(sub$event == 1 & sub$low0 == 0, na.rm = TRUE)
  scrT <- rbind(scrT, data.table(
    窗口 = lab, 目标筛查比例 = "约登最优", 阈值 = co$threshold, n = length(ev), 事件 = sum(ev), 检出人数 = sum(pos), 实际筛查比例 = round(100 * mean(pos), 2), TP = sum(pos & ev == 1), FP = sum(pos & ev == 0), TN = sum(!pos & ev == 0), FN = sum(!pos & ev == 1),
    灵敏度 = round(100 * co$sensitivity, 1), 特异度 = round(100 * co$specificity, 1),
    阳性预测值 = round(100 * mean(ev[pos]), 1), 阴性预测值 = round(100 * mean(ev[!pos] == 0), 1),
    检出者中基线已低分占比 = round(100 * mean(sub$low0 == 1, na.rm = TRUE), 1),
    新发真事件数 = yld_new, 每检出1例新发需筛人数 = if (yld_new > 0) round(n / yld_new, 1) else NA_real_))
}
pt(scrT)
fwrite(scrT, file.path(OUT, "F15_screening.csv"), bom = TRUE)

## ---------------- 9. 年龄组稳定性（全人口）---------------- 
say("\n========== 7. 年龄组稳定性（全人口口径）==========")
ageT <- data.table()
for (nm in c("w1", "w2")) {
  d <- copy(getc(nm))
  lab <- ifelse(nm == "w1", W1L, W2L)
  d[, ag := fcase(age0 < 60, "45-59", age0 < 70, "60-69", age0 < 80, "70-79",
                  default = "80+")]
  for (a in c("45-59", "60-69", "70-79", "80+")) {
    s <- d[ag == a]
    if (nrow(s) < 30 || length(unique(s$event)) < 2) next
    r <- roc(s$event, s$score, quiet = TRUE, direction = "<")
    ageT <- rbind(ageT, data.table(窗口 = lab, 年龄组 = a, n = nrow(s),
                                   事件 = sum(s$event),
                                   事件率 = round(100 * mean(s$event), 1),
                                   AUC_三项评分卡 = round(as.numeric(auc(r)), 4)))
  }
}
pt(ageT)
fwrite(ageT, file.path(OUT, "F15_agegroup.csv"), bom = TRUE)

## ---------------- 10. 校准（全人口 vs 未受损）---------------- 
say("\n========== 8. 校准 ==========")
calT <- data.table()
for (nm in c("w1", "w2")) {
  d0 <- copy(getc(nm))
  lab <- ifelse(nm == "w1", W1L, W2L)
  for (s in c("U_全人口", "N_基线未受损")) {
    d <- if (s == "U_全人口") d0 else d0[!is.na(low0) & low0 == 0]
    fit <- glm(event ~ lp, data = d, family = binomial())
    calT <- rbind(calT, data.table(窗口 = lab, 口径 = s, n = nrow(d),
                                   期望事件率 = round(100 * mean(d$phat), 1),
                                   实际事件率 = round(100 * mean(d$event), 1),
                                   校准斜率 = round(coef(fit)[2], 3),
                                   校准截距 = round(coef(fit)[1], 3),
                                   Brier = round(mean((d$phat - d$event)^2), 4)))
  }
}
pt(calT)
fwrite(calT, file.path(OUT, "F15_calibration.csv"), bom = TRUE)

## ---------------- 11. 加权敏感性（全人口）---------------- 
say("\n========== 9. 加权敏感性（基线权重，全人口口径）==========")
wt_auc <- function(p, y, wgt) {
  cs <- y == 1; ct <- y == 0
  num <- 0
  for (i in which(cs)) num <- num + wgt[i] * sum(wgt[ct] * ((p[i] > p[ct]) + 0.5 * (p[i] == p[ct])))
  num / (sum(wgt[cs]) * sum(wgt[ct]))
}
wtT <- data.table()
for (nm in c("w1", "w2")) {
  d <- getc(nm)[!is.na(wt0)]
  lab <- ifelse(nm == "w1", W1L, W2L)
  wtT <- rbind(wtT, data.table(窗口 = lab, n = nrow(d),
                               加权AUC_评分卡 = round(wt_auc(d$score, d$event, d$wt0), 4),
                               未加权AUC_评分卡 =
                                 round(as.numeric(auc(roc(d$event, d$score, quiet = TRUE, direction = "<"))), 4)))
}
pt(wtT)
fwrite(wtT, file.path(OUT, "F15_weighted.csv"), bom = TRUE)

## ---------------- 12. 图 ---------------- 
## ---------------- 12b. 分层稳定性（全人口口径内部）---------------- 
say("\n========== 10. 层内稳定性（全人口口径）==========")
say("（说明：年龄与教育本身是评分的构成项，层内 AUC 必然偏低，此处用于核对方向一致性）")
strataT <- data.table()
for (nm in c("w1", "w2")) {
  d <- copy(getc(nm))
  lab <- ifelse(nm == "w1", W1L, W2L)
  d[, sex := fcase(female0 == 1, "女性", female0 == 0, "男性", default = NA_character_)]
  d[, edulab := fcase(g3 == G3[1], "未上学", g3 == G3[2], "小学/初中",
                      g3 == G3[3], "高中及以上", default = NA_character_)]
  for (sv in c("sex", "edulab")) {
    for (lv in na.omit(unique(d[[sv]]))) {
      s <- d[get(sv) == lv]
      if (nrow(s) < 30 || length(unique(s$event)) < 2) next
      r <- roc(s$event, s$score, quiet = TRUE, direction = "<")
      strataT <- rbind(strataT, data.table(
        窗口 = lab, 分层变量 = ifelse(sv == "sex", "性别", "教育"),
        层 = lv, n = nrow(s), 事件率 = round(100 * mean(s$event), 1),
        AUC_三项评分卡 = round(as.numeric(auc(r)), 4)))
    }
  }
  s <- d[age0 < 50]
  if (nrow(s) >= 30 && length(unique(s$event)) >= 2) {
    r <- roc(s$event, s$score, quiet = TRUE, direction = "<")
    strataT <- rbind(strataT, data.table(窗口 = lab, 分层变量 = "年龄下限",
                                         层 = "45-49", n = nrow(s),
                                          事件率 = round(100 * mean(s$event), 1),
                                          AUC_三项评分卡 = round(as.numeric(auc(r)), 4)))
  }
  say("-- ", lab, "：年龄范围 ", min(d$age0), "-", max(d$age0),
      "，其中 45-49 岁 ", sum(d$age0 < 50), " 人 --")
}
pt(strataT)
fwrite(strataT, file.path(OUT, "F15_strata_within.csv"), bom = TRUE)

## ---------------- 13. 图 ---------------- 
png(file.path(OUT, "F15_fig_auc.png"), width = 1800, height = 720, res = 140)
par(mfrow = c(1, 2), mar = c(5, 13, 3.5, 2))
for (nm in c("w1", "w2")) {
  lab <- ifelse(nm == "w1", "W1 2011->2015", "W2 2015->2018")
  sub <- aucT[窗口 == ifelse(nm == "w1", W1L, W2L) & 预测子 == "三项评分卡"]
  sub <- sub[match(c("U_全人口", "N_基线未受损", "I_基线已受损", "M_基线记忆缺失"), sub$口径)]
  bp <- barplot(rev(sub$AUC), names.arg = rev(sub$口径), horiz = TRUE, las = 1,
                xlim = c(0.4, 1), col = "steelblue", main = lab, xlab = "AUC")
  text(rev(sub$AUC), bp, labels = sprintf("%.3f", rev(sub$AUC)), pos = 4, cex = 0.9)
  abline(v = 0.5, lty = 3)
}
dev.off()

close(con)
cat("\n完成 ->", file.path(OUT, "F15_report.txt"), "\n")
