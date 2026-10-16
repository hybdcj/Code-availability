# F26 LEGACY: BDSI实现/输入为历史版本；修正分析请运行F26_bdsi_audit.R。
# 原结果数值保留追溯，不可直接作为原版BDSI的有效比较；Tu限制另见F25。
# ============================================================================
# F14：与主流模型的正面对比（含敏感性、显著性、灵敏度/特异度）
# 三臂：① 本研究三项评分卡  ② BDSI（Barnes 2014，7 项）
#       ③ Tu 2026 六因子模型（含记忆＋定向＋握力）
# 窗口：CHARLS 2015 → 2018；结局：记忆分跌至「教育分层 mean−1SD 切点」以下
# 输出：C:/y03tmp/out/F14_report.txt、F14_*.csv
# ============================================================================
invisible(try(Sys.setlocale("LC_ALL", "zh_CN.UTF-8"), silent = TRUE))
suppressPackageStartupMessages({ library(data.table); library(pROC); library(haven) })
if (dir.exists("D:/R/library")) .libPaths(c("D:/R/library", .libPaths()))
TMP <- "C:/y03tmp/tmp"; OUT <- "C:/y03tmp/out"
con <- file(file.path(OUT, "F14_report.txt"), open = "wt", encoding = "UTF-8")
say <- function(...) { s <- paste0(..., collapse = ""); cat(s, "\n"); cat(s, "\n", file = con) }
pt  <- function(x) say(paste(capture.output(print(x)), collapse = "\n"))

B_TU <- c(age_lt65 = -0.53, age_ge75 = 0.26, male = 0.38,
          edu_primid = 0.45, edu_none = 0.03, grip = -0.16, memory = -0.33, orient = -0.15)
LV <- c("Less than primary", "Primary school", "Middle school", "High school or above")

## ---------------- 1. 认知与协变量 ----------------
w <- as.data.table(readRDS(file.path(TMP, "waves_cognition.rds")))
g <- w[year %in% c(2011, 2015, 2018), .(ID = IDn, year,
                                        mem = as.numeric(imrc) + as.numeric(dlrc),
                                        orient = as.numeric(orient), age = as.numeric(age),
                                        female = as.numeric(female), edu_group = edu_group,
                                        cesd10 = as.numeric(cesd10), adl = as.numeric(adl),
                                        iadl2 = as.numeric(iadl2),
                                        stroke = as.numeric(stroke), diab = as.numeric(diab),
                                        weight = as.numeric(weight))]
d15 <- g[year == 2015, .(ID, mem15 = mem, orient15 = orient, age15 = age, female15 = female,
                         edu15 = edu_group, cesd15 = cesd10, adl15 = adl, iadl2_15 = iadl2,
                         stroke15 = stroke, diab15 = diab, wt15 = weight)]
d18 <- g[year == 2018, .(ID, mem18 = mem)]

## ---------------- 2. 体格与握力（含哨兵值清洗）----------------
b <- as.data.table(read_dta(file.path(TMP, "w3_biomarker.dta"),
                            col_select = any_of(c("ID", "qi002", "ql002", "qc003", "qc004", "qc005", "qc006"))))
b[, `:=`(height = fifelse(qi002 >= 120 & qi002 <= 220, qi002, NA_real_),
         weight = fifelse(ql002 >= 25 & ql002 <= 200, ql002, NA_real_))]
for (v in c("qc003", "qc004", "qc005", "qc006"))
  b[[v]] <- fifelse(b[[v]] > 0 & b[[v]] <= 100, b[[v]], NA_real_)
b[, `:=`(bmi = weight / (height / 100)^2,
         grip_max = pmax(qc003, qc004, qc005, qc006, na.rm = TRUE))]
b[, grip_alt := rowMeans(cbind(pmax(qc003, qc004, na.rm = TRUE), pmax(qc005, qc006, na.rm = TRUE)),
                         na.rm = TRUE)]
b <- b[, .(ID = as.character(ID), bmi, grip_max, grip_alt)]
say("体格数据：身高可用 ", sum(!is.na(b$bmi)), " 人（BMI 中位 ",
    round(median(b$bmi, na.rm = TRUE), 2), "）")

## ---------------- 3. 教育分层切点（2011 重建，同 F07）----------------
d11 <- g[year == 2011 & !is.na(mem) & !is.na(edu_group) & age >= 50]
d11[, edu3 := fcase(edu_group == LV[1], "No formal schooling",
                    edu_group %in% LV[2:3], "Primary/middle",
                    edu_group == LV[4], "High school and above")]
cut3 <- d11[, .(n = .N, cutoff = mean(mem) - sd(mem)), by = edu3]

## ---------------- 4. 分析样本 ----------------
mm <- merge(merge(d15, d18, by = "ID"), b, by = "ID", all.x = TRUE)
mm[, edu3 := fcase(edu15 == LV[1], "No formal schooling",
                   edu15 %in% LV[2:3], "Primary/middle",
                   edu15 == LV[4], "High school and above")]
mm <- merge(mm, cut3[, .(edu3, cutoff)], by = "edu3", all.x = TRUE)
say("\n样本流程：2015+2018 记忆均可计分 ", nrow(mm[!is.na(mem15) & !is.na(mem18)]),
    " → ≥50 岁 ", nrow(mm[!is.na(mem15) & !is.na(mem18) & age15 >= 50]),
    " → 三臂变量齐全 ", nrow(mm[!is.na(mem15) & !is.na(mem18) & age15 >= 50 & !is.na(cutoff) &
                                   !is.na(bmi) & !is.na(grip_max) & !is.na(orient15) &
                                   !is.na(cesd15) & !is.na(iadl2_15) & !is.na(stroke15) & !is.na(diab15)]))
d <- mm[!is.na(mem15) & !is.na(mem18) & age15 >= 50 & !is.na(cutoff) & !is.na(bmi) &
          !is.na(grip_max) & !is.na(orient15) & !is.na(cesd15) & !is.na(iadl2_15) &
          !is.na(stroke15) & !is.na(diab15)]
d <- d[mem15 >= cutoff]                        # 2015 记忆正常
d[, event := as.integer(mem18 < cutoff)]
say("分析样本 n = ", nrow(d), "｜事件 = ", sum(d$event), "（", round(100 * mean(d$event), 2), "%）")

## ---------------- 5. 三臂评分 ----------------
d[, edu_f := factor(edu15, levels = LV)]
d[, male := as.integer(female15 == 0)]
## ① 本研究三项评分卡
d[, score_us := 2L * pmax(floor((age15 - 45) / 5), 0) +
      (-2L) * as.integer(edu_f == LV[2]) + (-7L) * as.integer(edu_f == LV[3]) +
      (-10L) * as.integer(edu_f == LV[4]) + 3L * as.integer(cesd15 >= 10)]
## ② BDSI（Barnes 2014）：年龄 1 分/岁(65–79)、<12 年教育 9、卒中 6、糖尿病 3、BMI<18.5 8、需协助理财/服药 10、抑郁 6
d[, `:=`(b_age = pmax(0, pmin(age15, 79) - 64),
         b_edu = 9L * as.integer(edu_f != LV[4]),      # <12 年 ≈ 高中以下
         b_stroke = 6L * as.integer(stroke15 == 1),
         b_diab = 3L * as.integer(diab15 == 1),
         b_bmi = 8L * as.integer(bmi < 18.5),
         b_func = 10L * as.integer(iadl2_15 == 1),
         b_dep = 6L * as.integer(cesd15 >= 10))]
d[, score_bdsi := b_age + b_edu + b_stroke + b_diab + b_bmi + b_func + b_dep]
## ③ Tu 2026 六因子线性预测子
d[, lp_tu := B_TU[["age_lt65"]] * as.integer(age15 < 65) + B_TU[["age_ge75"]] * as.integer(age15 >= 75) +
      B_TU[["male"]] * male + B_TU[["edu_primid"]] * as.integer(edu3 == "Primary/middle") +
      B_TU[["edu_none"]] * as.integer(edu3 == "No formal schooling") +
      B_TU[["grip"]] * grip_max + B_TU[["memory"]] * mem15 + B_TU[["orient"]] * orient15]

say("\nBDSI 得分分布：均值 ", round(mean(d$score_bdsi), 1), " 中位 ", median(d$score_bdsi),
    " 范围 ", min(d$score_bdsi), "-", max(d$score_bdsi))
say("BDSI 各构成项均值：年龄 ", round(mean(d$b_age), 2), "｜教育 ", round(mean(d$b_edu), 2),
    "｜卒中 ", round(mean(d$b_stroke), 2), "｜糖尿病 ", round(mean(d$b_diab), 2),
    "｜BMI ", round(mean(d$b_bmi), 2), "｜功能 ", round(mean(d$b_func), 2),
    "｜抑郁 ", round(mean(d$b_dep), 2))
say("本研究三项卡：均值 ", round(mean(d$score_us), 1), " 范围 ", min(d$score_us), "-", max(d$score_us))

## ---------------- 6. 主结果：AUC 与两两显著性 ----------------
arms <- c(score_us = "本研究三项评分卡", score_bdsi = "BDSI（7 项）", lp_tu = "Tu 2026（6 因子，含认知测验）")
res <- data.table()
rocs <- list()
for (v in names(arms)) {
  r <- roc(d$event, d[[v]], quiet = TRUE, direction = "<")
  rocs[[v]] <- r
  ci <- as.numeric(ci.auc(r))
  res <- rbind(res, data.table(模型 = arms[[v]], AUC = round(ci[2], 4),
                               低 = round(ci[1], 4), 高 = round(ci[3], 4)))
}
say("\n========== AUC ==========")
pt(res)
fwrite(res, file.path(OUT, "F14_auc.csv"), bom = TRUE)

say("\n========== 与本研究的三项卡作配对比较 ==========")
cmp <- data.table()
for (v in c("score_bdsi", "lp_tu")) {
  dt <- roc.test(rocs[["score_us"]], rocs[[v]], paired = TRUE)
  set.seed(20261015); B <- 2000; dv <- rep(NA_real_, B)
  for (i in seq_len(B)) {
    idx <- sample.int(nrow(d), replace = TRUE); db <- d[idx]
    if (length(unique(db$event)) < 2) next
    dv[i] <- as.numeric(auc(roc(db$event, db$score_us, quiet = TRUE, direction = "<"))) -
             as.numeric(auc(roc(db$event, db[[v]],  quiet = TRUE, direction = "<")))
  }
  cmp <- rbind(cmp, data.table(对照 = arms[[v]],
                               本研究AUC = round(as.numeric(auc(rocs[["score_us"]])), 4),
                               对照AUC = round(as.numeric(auc(rocs[[v]])), 4),
                               ΔAUC = round(as.numeric(auc(rocs[["score_us"]])) - as.numeric(auc(rocs[[v]])), 4),
                               DeLong_P = signif(dt$p.value, 3),
                               boot_低 = round(quantile(dv, .025, na.rm = TRUE), 4),
                               boot_高 = round(quantile(dv, .975, na.rm = TRUE), 4)))
}
pt(cmp)
fwrite(cmp, file.path(OUT, "F14_pairwise.csv"), bom = TRUE)

## ---------------- 7. 灵敏度／特异度（目标筛查比例（含同分者，报告实际比例）+ 约登）----------------
oper <- function(v, lab) {
  s <- d[[v]]; ev <- d$event; n <- length(ev)
  out <- data.table()
  for (fr in c(0.10, 0.20, 0.30)) {
    k <- floor(fr * n); th <- sort(s, decreasing = TRUE)[k]
    pos <- s >= th
    out <- rbind(out, data.table(模型 = lab, 口径 = paste0("目标筛查 ", fr * 100, "%"),
                                 阈值 = th, n = length(ev), 事件 = sum(ev), 检出人数 = sum(pos), 实际筛查比例 = round(100 * mean(pos), 2), TP = sum(pos & ev == 1), FP = sum(pos & ev == 0), TN = sum(!pos & ev == 0), FN = sum(!pos & ev == 1), 灵敏度 = round(100 * sum(ev[pos]) / sum(ev), 1),
                                 特异度 = round(100 * sum(!pos & ev == 0) / sum(ev == 0), 1),
                                 阴性预测值 = round(100 * mean(ev[!pos] == 0), 1), 阳性者事件率 = round(100 * mean(ev[pos]), 1)))
  }
  co <- coords(roc(ev, s, quiet = TRUE, direction = "<"), "best", best.method = "youden",
               ret = c("threshold", "sensitivity", "specificity"))
  co <- co[which.max(co$specificity), , drop = FALSE]
  pos <- s >= co$threshold
  out <- rbind(out, data.table(模型 = lab, 口径 = "约登最优",
                               阈值 = round(co$threshold, 2), n = length(ev), 事件 = sum(ev), 检出人数 = sum(pos), 实际筛查比例 = round(100 * mean(pos), 2), TP = sum(pos & ev == 1), FP = sum(pos & ev == 0), TN = sum(!pos & ev == 0), FN = sum(!pos & ev == 1),
                               灵敏度 = round(100 * co$sensitivity, 1),
                               特异度 = round(100 * co$specificity, 1),
                               阴性预测值 = round(100 * mean(ev[!pos] == 0), 1), 阳性者事件率 = round(100 * mean(ev[pos]), 1)))
  out
}
op <- rbind(oper("score_us", "本研究三项评分卡"), oper("score_bdsi", "BDSI（7 项）"),
            oper("lp_tu", "Tu 2026（6 因子）"))
say("\n========== 灵敏度／特异度 ==========")
pt(op)
fwrite(op, file.path(OUT, "F14_operating.csv"), bom = TRUE)

## ---------------- 8. 风险分层（各自四分位）----------------
say("\n========== 风险分层事件率（各自四分位）==========")
for (v in names(arms)) {
  q <- cut(d[[v]], quantile(d[[v]], seq(0, 1, .25)), include.lowest = TRUE, labels = 1:4)
  cat("", arms[[v]], "：", paste(d[, .(r = paste0(round(100 * mean(event), 1), "%")), by = q][order(q)]$r,
                                collapse = "  →  "), "\n")
  say("  ", arms[[v]], "：", paste(d[, .(r = paste0(round(100 * mean(event), 1), "%")), by = q][order(q)]$r,
                                   collapse = "  →  "))
}

## ---------------- 9. 敏感性分析 ----------------
say("\n========== 敏感性分析 ==========")
## (a) 结局改用本研究全局切点（≤3）
mu15 <- mean(mm$mem15, na.rm = TRUE); sd15 <- sd(mm$mem15, na.rm = TRUE); k15 <- floor(mu15 - sd15)
d2 <- copy(d)[, event2 := as.integer(mem18 <= k15)]
say("(a) 结局改用本研究全局切点 ≤", k15, "（事件 ", sum(d2$event2), "，",
    round(100 * mean(d2$event2), 2), "%）")
for (v in names(arms))
  say("    ", arms[[v]], " AUC = ",
      round(as.numeric(auc(roc(d2$event2, d2[[v]], quiet = TRUE, direction = "<"))), 4))
## (b) BDSI 抑郁项改严格阈值 CES-D≥16
d3 <- copy(d)[, score_bdsi_b := b_age + b_edu + b_stroke + b_diab + b_bmi + b_func +
                 6L * as.integer(cesd15 >= 16)]
say("(b) BDSI 抑郁项改 CES-D≥16：AUC = ",
    round(as.numeric(auc(roc(d3$event, d3$score_bdsi_b, quiet = TRUE, direction = "<"))), 4),
    "（原 CES-D≥10 为 ", round(as.numeric(auc(rocs[["score_bdsi"]])), 4), "）")
## (c) 握力改用两手最好值均值（影响 Tu 臂）
d4 <- copy(d)[, lp_tu_b := B_TU[["age_lt65"]] * as.integer(age15 < 65) + B_TU[["age_ge75"]] * as.integer(age15 >= 75) +
                 B_TU[["male"]] * male + B_TU[["edu_primid"]] * as.integer(edu3 == "Primary/middle") +
                 B_TU[["edu_none"]] * as.integer(edu3 == "No formal schooling") +
                 B_TU[["grip"]] * grip_alt + B_TU[["memory"]] * mem15 + B_TU[["orient"]] * orient15]
say("(c) Tu 臂握力改两手最好值均值：AUC = ",
    round(as.numeric(auc(roc(d4$event, d4$lp_tu_b, quiet = TRUE, direction = "<"))), 4))
## (d) 加权（2015 基线权重）
dw <- d[!is.na(wt15)]
wt_auc <- function(p, y, wgt) {
  cs <- y == 1; ct <- y == 0
  num <- 0
  for (i in which(cs)) num <- num + wgt[i] * sum(wgt[ct] * ((p[i] > p[ct]) + 0.5 * (p[i] == p[ct])))
  num / (sum(wgt[cs]) * sum(wgt[ct]))
}
say("(d) 加权敏感性（2015 基线权重，可用 n = ", nrow(dw), "）：")
for (v in names(arms))
  say("    ", arms[[v]], " 加权 AUC = ", round(wt_auc(dw[[v]], dw$event, dw$wt15), 4))

close(con)
cat("\n完成 ->", file.path(OUT, "F14_report.txt"), "\n")
