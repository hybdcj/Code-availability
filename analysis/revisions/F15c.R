# F32（2026-10-16）：以下保留旧共同完整样本的计算以便追溯；CHARLS整数评分自身样本已改用同目录 F32_全人群入组与分母修正.R。不要仅删除lp限制后直接运行旧脚本，因为下游连续模型、分母和表头也需分别处理。
# ============================================================================
# F15c：把 F15 的"全人口口径"换成**主稿原有阈值（全局 ≤4）**，与第一轮 R02
#        的宽样本（n=8,296、连续模型 AUC 0.726）对齐，解决口径不可比问题
# 窗口：CHARLS 2011 -> 2015；结局：2015 记忆分 <= 4（基线 > 4）
# ============================================================================
invisible(try(Sys.setlocale("LC_ALL", "zh_CN.UTF-8"), silent = TRUE))
suppressPackageStartupMessages({ library(data.table); library(pROC) })
if (dir.exists("D:/R/library")) .libPaths(c("D:/R/library", .libPaths()))
TMP <- "C:/y03tmp/tmp"; OUT <- "C:/y03tmp/out"
con <- file(file.path(OUT, "F15c_report.txt"), open = "wt", encoding = "UTF-8")
say <- function(...) { s <- paste0(..., collapse = ""); cat(s, "\n"); cat(s, "\n", file = con) }
pt <- function(x) say(paste(capture.output(print(x)), collapse = "\n"))

B <- c(intercept = -4.19524508778947, age = 0.0492039250378111,
       eduPri = -0.290105771000328, eduMid = -0.994711800202235,
       eduHigh = -1.43588388035287, dep = 0.412166484954086, adl = -0.0110119344848014)
LV <- c("Less than primary", "Primary school", "Middle school", "High school or above")
K <- 4L   # 主稿阈值：2011 五波固定队列 mean-1SD = 4.7647 -> <= 4

w <- as.data.table(readRDS(file.path(TMP, "waves_cognition.rds")))
g <- w[year %in% c(2011, 2015), .(ID = IDn, year = as.numeric(year),
        mem = as.numeric(imrc) + as.numeric(dlrc),
        age = as.numeric(age), edu = edu_group, cesd = as.numeric(cesd10),
        adl = as.numeric(adl), female = as.numeric(female), wt = as.numeric(weight))]
a <- g[year == 2011, .(ID, mem0 = mem, age0 = age, edu0 = edu, cesd0 = cesd,
                       adl0 = adl, wt0 = wt, female0 = female)]
b <- g[year == 2015, .(ID, mem1 = mem)]
d <- merge(a, b, by = "ID")
d[, edu_f := factor(edu0, levels = LV)]
d[, score := 2L * pmax(floor((age0 - 45) / 5), 0) +
      (-2L) * as.integer(edu_f == LV[2]) + (-7L) * as.integer(edu_f == LV[3]) +
      (-10L) * as.integer(edu_f == LV[4]) + 3L * as.integer(cesd0 >= 10)]
d[, lp := B[["intercept"]] + B[["age"]] * age0 + B[["dep"]] * as.integer(cesd0 >= 10) +
      B[["adl"]] * adl0 + B[["eduPri"]] * as.integer(edu_f == LV[2]) +
      B[["eduMid"]] * as.integer(edu_f == LV[3]) + B[["eduHigh"]] * as.integer(edu_f == LV[4])]
d[, phat := plogis(lp)]
d[, event := as.integer(mem1 <= K)]
d[, low0 := as.integer(mem0 <= K)]
ok <- !is.na(d$score) & !is.na(d$lp) & !is.na(d$event) & !is.na(d$age0) & d$age0 >= 45
d <- d[ok]

say("阈值：全局 <= ", K, "（主稿口径；源自 2011 五波固定队列 mean-1SD = 4.7647）")
say("样本流程：两波 ID 可匹配 ", nrow(merge(a, b, by = "ID")),
    " -> 年龄/评分/lp/结局齐全且 >=45 岁：", nrow(d))
say("A11/R02 记录：2011 年 >=45 岁且基线 >4 者 10,422 人；2015 记忆可测 8,371；")
say("             四变量完整 8,296 人（= 宽样本）；其中事件 1,889（22.8%）")

arm <- function(x, lab) {
  if (nrow(x) < 30 || length(unique(x$event)) < 2) return(NULL)
  r1 <- roc(x$event, x$score, quiet = TRUE, direction = "<")
  r2 <- roc(x$event, x$phat, quiet = TRUE, direction = "<")
  c1 <- as.numeric(ci.auc(r1)); c2 <- as.numeric(ci.auc(r2))
  data.table(口径 = lab, n = nrow(x), 事件 = sum(x$event),
             事件率 = round(100 * mean(x$event), 2),
             AUC_评分卡 = round(c1[2], 4), 卡低 = round(c1[1], 4), 卡高 = round(c1[3], 4),
             AUC_连续 = round(c2[2], 4), 连低 = round(c2[1], 4), 连高 = round(c2[3], 4))
}
res <- rbindlist(list(
  arm(d, "U 全人口（不限基线记忆）"),
  arm(d[!is.na(mem0) & mem0 > K], "N 基线未受损（主稿宽样本口径）"),
  arm(d[!is.na(mem0) & mem0 <= K], "I 基线已受损"),
  arm(d[is.na(mem0)], "M 基线记忆缺失")), fill = TRUE)
say("\n========== 全局 <=4 阈值下各口径 ==========")
pt(res)
fwrite(res, file.path(OUT, "F15c_auc_globalcut.csv"), bom = TRUE)

say("\n========== 与第一轮 R02 宽样本的对照 ==========")
pt(data.table(
  指标 = c("n", "事件", "事件率", "连续模型 AUC", "评分卡 AUC"),
  R02宽样本 = c("8,296", "1,889", "22.8%", "0.726", "未报告"),
  本次N口径 = c(format(nrow(d[!is.na(mem0) & mem0 > K]), big.mark = ","),
                 format(sum(d[!is.na(mem0) & mem0 > K]$event), big.mark = ","),
                 paste0(round(100 * mean(d[!is.na(mem0) & mem0 > K]$event), 1), "%"),
                 as.character(res[口径 == "N 基线未受损（主稿宽样本口径）"]$AUC_连续),
                 as.character(res[口径 == "N 基线未受损（主稿宽样本口径）"]$AUC_评分卡))))

say("\n========== 放宽阶梯（2011 -> 2015 窗口，全局 <=4 阈值）==========")
uq <- d[!is.na(mem0)]
say("① 主稿原开发样本（五波共同题集完整）：4,441 人 / 658 事件 / 表观 AUC 0.702")
say("② 宽样本（取消五波要求，基线 >4）：", nrow(uq[mem0 > K]), " 人 / ",
    sum(uq[mem0 > K]$event), " 事件 / AUC 0.726（R02 记录 8,296 / 1,889）")
say("③ 全人口（再取消基线 >4 要求，此为 F15 的口径）：", nrow(d), " 人 / ",
    sum(d$event), " 事件")
say("④ 全人口 + 不限基线记忆：" , nrow(d), " 人 / ", sum(d$event), " 事件")

say("\n========== 分解：全人口事件里基线已低分的占比 ==========")
dd <- d[!is.na(mem0)]
say("事件中基线 <=4 者占比 = ", round(100 * mean(dd$low0[dd$event == 1]), 1), "%")
say("其中：基线已低分者事件率 = ", round(100 * mean(dd[low0 == 1]$event), 1),
    "%｜基线未受损者事件率 = ", round(100 * mean(dd[low0 == 0]$event), 1), "%")

close(con)
cat("\n完成 ->", file.path(OUT, "F15c_report.txt"), "\n")
