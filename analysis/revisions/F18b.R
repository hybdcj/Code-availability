# F32（2026-10-16）：以下保留旧共同完整样本的计算以便追溯；CHARLS整数评分自身样本已改用同目录 F32_全人群入组与分母修正.R。不要仅删除lp限制后直接运行旧脚本，因为下游连续模型、分母和表头也需分别处理。
# ============================================================================
# F18b：补描述"阶梯第 ③ 级比第 ② 级新增的 3,457 人"
# 口径：主稿阈值 <=4，CHARLS 2011 -> 2015
# ============================================================================
invisible(try(Sys.setlocale("LC_ALL", "zh_CN.UTF-8"), silent = TRUE))
suppressPackageStartupMessages({ library(data.table) })
if (dir.exists("D:/R/library")) .libPaths(c("D:/R/library", .libPaths()))
TMP <- "C:/y03tmp/tmp"; OUT <- "C:/y03tmp/out"
con <- file(file.path(OUT, "F18b_report.txt"), open = "wt", encoding = "UTF-8")
say <- function(...) { s <- paste0(..., collapse = ""); cat(s, "\n"); cat(s, "\n", file = con) }
pt <- function(x) say(paste(capture.output(print(x)), collapse = "\n"))
LV <- c("Less than primary", "Primary school", "Middle school", "High school or above")
B <- c(intercept = -4.19524508778947, age = 0.0492039250378111,
       eduPri = -0.290105771000328, eduMid = -0.994711800202235,
       eduHigh = -1.43588388035287, dep = 0.412166484954086, adl = -0.0110119344848014)
K <- 4L

w <- as.data.table(readRDS(file.path(TMP, "waves_cognition.rds")))
g <- w[year %in% c(2011, 2015), .(ID = IDn, year = as.numeric(year),
        mem = as.numeric(imrc) + as.numeric(dlrc), age = as.numeric(age),
        edu = edu_group, cesd = as.numeric(cesd10), adl = as.numeric(adl))]
a <- g[year == 2011, .(ID, mem0 = mem, age0 = age, edu0 = edu, cesd0 = cesd, adl0 = adl)]
b <- g[year == 2015, .(ID, mem1 = mem)]
d <- merge(a, b, by = "ID")
d[, edu_f := factor(edu0, levels = LV)]
d[, score := 2L * pmax(floor((age0 - 45) / 5), 0) +
      (-2L) * as.integer(edu_f == LV[2]) + (-7L) * as.integer(edu_f == LV[3]) +
      (-10L) * as.integer(edu_f == LV[4]) + 3L * as.integer(cesd0 >= 10)]
d[, lp := B[["intercept"]] + B[["age"]] * age0 + B[["dep"]] * as.integer(cesd0 >= 10) +
      B[["adl"]] * adl0 + B[["eduPri"]] * as.integer(edu_f == LV[2]) +
      B[["eduMid"]] * as.integer(edu_f == LV[3]) + B[["eduHigh"]] * as.integer(edu_f == LV[4])]
d <- d[!is.na(score) & !is.na(lp) & !is.na(age0) & age0 >= 45 & !is.na(mem1)]
d[, event := as.integer(mem1 <= K)]

grp <- function(x, lab) data.table(
  组 = lab, n = nrow(x),
  平均年龄 = round(mean(x$age0), 1),
  未上学比例 = round(100 * mean(x$edu0 == LV[1]), 1),
  基线记忆分可测 = sum(!is.na(x$mem0)),
  基线记忆分均值 = round(mean(x$mem0, na.rm = TRUE), 2),
  事件 = sum(x$event), 事件率 = round(100 * mean(x$event), 1),
  评分均值 = round(mean(x$score), 2))

T <- rbindlist(list(
  grp(d[!is.na(mem0) & mem0 > K], "① 宽样本（基线>4，R02 第②级）"),
  grp(d[!is.na(mem0) & mem0 <= K], "② 基线已低分（新增）"),
  grp(d[is.na(mem0)], "③ 基线记忆不可测（新增）"),
  grp(d, "合计（全人口，第③级）")))
say("按主稿阈值 <=4、CHARLS 2011->2015。")
say("第③级 11,753 人 − 第②级 8,296 人 = 新增 ", nrow(d) - nrow(d[!is.na(mem0) & mem0 > K]), " 人")
say("校验：② + ③ = ", nrow(d[!is.na(mem0) & mem0 <= K]) + nrow(d[is.na(mem0)]))
pt(T)
fwrite(T, file.path(OUT, "F18b_newcomers.csv"), bom = TRUE)
close(con)
cat("\n完成\n")
