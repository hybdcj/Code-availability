# F32（2026-10-16）：以下保留旧共同完整样本的计算以便追溯；CHARLS整数评分自身样本已改用同目录 F32_全人群入组与分母修正.R。不要仅删除lp限制后直接运行旧脚本，因为下游连续模型、分母和表头也需分别处理。
# ============================================================================
# F20：高龄（>=90 岁）人数在**全人群口径**下的重算（更新第一轮 R09）
# 四队列：CHARLS 2011->2015、2015->2018；HRS、ELSA、SHARE
# 同时给出"主稿阈值 <=4"与"各队列自身 mean-1SD"两套口径下的事件数
# ============================================================================
invisible(try(Sys.setlocale("LC_ALL", "zh_CN.UTF-8"), silent = TRUE))
suppressPackageStartupMessages({ library(data.table) })
if (dir.exists("D:/R/library")) .libPaths(c("D:/R/library", .libPaths()))
TMP <- "C:/y03tmp/tmp"; OUT <- "C:/y03tmp/out"
con <- file(file.path(OUT, "F20_report.txt"), open = "wt", encoding = "UTF-8")
say <- function(...) { s <- paste0(..., collapse = ""); cat(s, "\n"); cat(s, "\n", file = con) }
pt <- function(x) say(paste(capture.output(print(x)), collapse = "\n"))
LV <- c("Less than primary", "Primary school", "Middle school", "High school or above")
G3 <- c("A_none", "B_primid", "C_high")
K <- 4L
B <- c(intercept = -4.19524508778947, age = 0.0492039250378111,
       eduPri = -0.290105771000328, eduMid = -0.994711800202235,
       eduHigh = -1.43588388035287, dep = 0.412166484954086, adl = -0.0110119344848014)

## ---------- CHARLS：构造两个窗口的宽表 ----------
w <- as.data.table(readRDS(file.path(TMP, "waves_cognition.rds")))
g <- w[year %in% c(2011, 2015, 2018), .(ID = IDn, year = as.numeric(year),
        mem = as.numeric(imrc) + as.numeric(dlrc), age = as.numeric(age),
        edu = edu_group, cesd = as.numeric(cesd10), adl = as.numeric(adl))]
d11 <- g[year == 2011 & !is.na(mem) & !is.na(edu) & age >= 50]
d11[, g3 := fcase(edu == LV[1], G3[1], edu %in% LV[2:3], G3[2], edu == LV[4], G3[3])]
cut3 <- d11[!is.na(g3), .(cutoff = mean(mem) - sd(mem)), by = g3]

make_win <- function(y0, y1) {
  a <- g[year == y0, .(ID, mem0 = mem, age0 = age, edu0 = edu, cesd0 = cesd, adl0 = adl)]
  b <- g[year == y1, .(ID, mem1 = mem)]
  d <- merge(a, b, by = "ID")
  d[, g3 := fcase(edu0 == LV[1], G3[1], edu0 %in% LV[2:3], G3[2], edu0 == LV[4], G3[3])]
  d <- merge(d, cut3[, .(g3, cutoff)], by = "g3", all.x = TRUE)
  d[, edu_f := factor(edu0, levels = LV)]
  d[, score := 2L * pmax(floor((age0 - 45) / 5), 0) +
        (-2L) * as.integer(edu_f == LV[2]) + (-7L) * as.integer(edu_f == LV[3]) +
        (-10L) * as.integer(edu_f == LV[4]) + 3L * as.integer(cesd0 >= 10)]
  d[, lp := B[["intercept"]] + B[["age"]] * age0 + B[["dep"]] * as.integer(cesd0 >= 10) +
      B[["adl"]] * adl0 + B[["eduPri"]] * as.integer(edu_f == LV[2]) +
      B[["eduMid"]] * as.integer(edu_f == LV[3]) + B[["eduHigh"]] * as.integer(edu_f == LV[4])]
  d <- d[!is.na(score) & !is.na(lp) & !is.na(age0) & age0 >= 45 & !is.na(mem1) & !is.na(cutoff)]
  d[]
}
W1 <- make_win(2011, 2015); W2 <- make_win(2015, 2018)

rowU <- function(d, co, win) {
  data.table(队列 = co, 窗口 = win, 口径 = "U 全人口", n = nrow(d),
             n90 = sum(d$age0 >= 90), pct90 = round(100 * mean(d$age0 >= 90), 2),
             事件_主稿阈值4 = sum(d$mem1 <= K), 事件_教育分层 = sum(d$mem1 < d$cutoff),
             年龄范围 = paste0(min(d$age0), "-", max(d$age0)))
}
rowN <- function(d, co, win) {
  s <- d[!is.na(mem0) & mem0 >= cutoff]
  data.table(队列 = co, 窗口 = win, 口径 = "N 基线未受损", n = nrow(s),
             n90 = sum(s$age0 >= 90), pct90 = round(100 * mean(s$age0 >= 90), 2),
             事件_主稿阈值4 = sum(s$mem1 <= K), 事件_教育分层 = sum(s$mem1 < s$cutoff),
             年龄范围 = paste0(min(s$age0), "-", max(s$age0)))
}
ch <- rbind(rowU(W1, "CHARLS", "2011->2015"), rowN(W1, "CHARLS", "2011->2015"),
            rowU(W2, "CHARLS", "2015->2018"), rowN(W2, "CHARLS", "2015->2018"))
say("========== CHARLS ==========")
pt(ch); fwrite(ch, file.path(OUT, "F20_charls.csv"), bom = TRUE)

## ---------- 外部队列 ----------
mk_ext <- function(co) {
  d <- as.data.table(readRDS(file.path(TMP, paste0("F16_base_", co, ".rds"))))
  d <- d[!is.na(age) & age >= 45]
  cs <- d[!is.na(mem_base)]
  d[, ci := floor(mean(cs$mem_base) - sd(cs$mem_base))]
  d[, score := 2L * pmax(floor((age - 45) / 5), 0) + edu_pts + 3L * dep]
  d[, lp := B[["intercept"]] + B[["age"]] * age + B[["dep"]] * dep + B[["adl"]] * adl +
      B[["eduPri"]] * as.integer(edu_pts == -2) + B[["eduMid"]] * as.integer(edu_pts == -7) +
      B[["eduHigh"]] * as.integer(edu_pts == -10)]
  # 与 F16 的 U 口径一致：只要求评分可算与随访记忆可得（不要求 ADL）
  d[!is.na(score) & !is.na(mem_fu)]
}
ext_row <- function(d, co, lab) {
  data.table(队列 = co, 窗口 = "基线->随访", 口径 = lab, n = nrow(d),
             n90 = sum(d$age >= 90), pct90 = round(100 * mean(d$age >= 90), 2),
             事件_主稿阈值4 = sum(d$mem_fu <= K), 事件_自身切点 = sum(d$mem_fu <= d$ci),
             年龄范围 = paste0(min(d$age), "-", max(d$age)))
}
ext <- rbindlist(lapply(c("HRS", "ELSA", "SHARE"), function(co) {
  d <- mk_ext(co)
  rbind(ext_row(d, co, "U 全人口"),
        ext_row(d[!is.na(mem_base) & mem_base > ci], co, "N 基线未受损"))
}))
say("\n========== 外部队列 ==========")
pt(ext); fwrite(ext, file.path(OUT, "F20_external.csv"), bom = TRUE)

## ---------- 高龄亚组（>=90）自身的事件 ----------
sub90 <- function(d, memvar, agevar, co) {
  s <- d[!is.na(d[[agevar]]) & d[[agevar]] >= 90 & !is.na(d[[memvar]])]
  data.table(队列 = co, n = nrow(s), 事件_主稿阈值4 = sum(s[[memvar]] <= K),
             事件率 = round(100 * mean(s[[memvar]] <= K), 1),
             年龄范围 = paste0(min(s[[agevar]]), "-", max(s[[agevar]])))
}
col1 <- sub90(W1, "mem1", "age0", "CHARLS 2011->2015")
col2 <- sub90(W2, "mem1", "age0", "CHARLS 2015->2018")
col3 <- rbindlist(lapply(c("HRS", "ELSA", "SHARE"), function(co) {
  sub90(mk_ext(co), "mem_fu", "age", co)
}))
hh <- rbind(col1, col2, col3)
say("\n========== 高龄亚组（>=90 岁）自身 ==========")
pt(hh); fwrite(hh, file.path(OUT, "F20_90plus_subgroup.csv"), bom = TRUE)

close(con)
cat("\n完成 ->", file.path(OUT, "F20_report.txt"), "\n")
