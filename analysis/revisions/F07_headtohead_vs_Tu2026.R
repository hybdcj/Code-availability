# ============================================================================
# F07：与 Tu 等（GeroScience 2026）模型的正面对比
# 窗口：CHARLS 2015 基线 → 2018 结局（正是他们的时段验证窗口）
# 结局：采用**他们的口径**——记忆分跌到"教育分层 mean−1SD 切点"以下
# 对比：他们的 6 因子 CoxPH 线性预测子（用其论文 Table 2 已发表系数）
#       vs 我们的 4 因子冻结评分（A14，2011→2015 拟合，直接搬运）
# 输出：C:/y03tmp/out/F07_report.txt、F07_*.csv
# 事前声明：本次为事后正面对比；我们的评分是 4 年模型搬到 3 年窗口，
#           绝对值会失校准，只比较排序类指标（AUC）与相对结构。
# ============================================================================
invisible(try(Sys.setlocale("LC_ALL", "zh_CN.UTF-8"), silent = TRUE))
suppressPackageStartupMessages({ library(data.table); library(pROC); library(haven) })
if (dir.exists("D:/R/library")) .libPaths(c("D:/R/library", .libPaths()))
TMP <- "C:/y03tmp/tmp"; OUT <- "C:/y03tmp/out"
con <- file(file.path(OUT, "F07_report.txt"), open = "wt", encoding = "UTF-8")
say <- function(...) { s <- paste0(..., collapse = ""); cat(s, "\n"); cat(s, "\n", file = con) }
pt  <- function(x) say(paste(capture.output(print(x)), collapse = "\n"))

## ---------- 他们的系数（论文 Table 2，CoxPH）----------
B_TU <- c(age_lt65 = -0.53, age_ge75 = 0.26, male = 0.38,
          edu_primid = 0.45, edu_none = 0.03,
          grip = -0.16, memory = -0.33, orient = -0.15)
## ---------- 我们的冻结模型（A14）----------
B_US <- c("(Intercept)" = -4.19524508778947, age = 0.0492039250378111,
          eduP = -0.290105771000328, eduM = -0.994711800202235, eduH = -1.43588388035287,
          dep = 0.412166484954086, adl = -0.0110119344848014)
LV <- c("Less than primary", "Primary school", "Middle school", "High school or above")

## ---------------- 1. 认知数据 ----------------
w <- as.data.table(readRDS(file.path(TMP, "waves_cognition.rds")))
g <- w[year %in% c(2011, 2015, 2018), .(ID = IDn, year,
                                        mem = as.numeric(imrc) + as.numeric(dlrc),
                                        orient = as.numeric(orient),
                                        age = as.numeric(age), female = as.numeric(female),
                                        edu_group = edu_group,
                                        cesd10 = as.numeric(cesd10), adl = as.numeric(adl))]

## ---------------- 2. 2015 握力 ----------------
b <- as.data.table(read_dta(file.path(TMP, "w3_biomarker.dta"),
                            col_select = any_of(c("ID", "qc003", "qc004", "qc005", "qc006"))))
say("生物标记文件 ID 列：", paste(setdiff(names(b), c("qc003","qc004","qc005","qc006")), collapse = ", "))
grip <- b[, .(ID = as.character(ID),
              grip_max = do.call(pmax, c(.SD, na.rm = TRUE)),
              grip_n   = rowSums(!is.na(.SD))),
          .SDcols = c("qc003", "qc004", "qc005", "qc006")]
say("2015 握力：非缺失 ", sum(!is.na(grip$grip_max)), " 人；中位数 ",
    round(median(grip$grip_max, na.rm = TRUE), 1), " kg")

## ---------------- 3. 教育分层切点（在 2011 波上重建）----------------
d11 <- g[year == 2011 & !is.na(mem) & !is.na(edu_group) & age >= 50]
d11[, edu3 := fcase(edu_group == LV[1], "No formal schooling",
                    edu_group %in% LV[2:3], "Primary/middle",
                    edu_group == LV[4], "High school and above")]
cut3 <- d11[, .(n = .N, mean = mean(mem), sd = sd(mem)), by = edu3]
cut3[, cutoff := mean - sd]
say("\n教育分层切点（在 CHARLS 2011、≥50 岁、记忆可计分者上重建）：")
pt(cut3[, .(edu3, n, 均值 = round(mean, 2), SD = round(sd, 2), 切点 = round(cutoff, 2))])

## ---------------- 4. 构建 2015→2018 分析样本 ----------------
d15 <- g[year == 2015, .(ID, mem15 = mem, orient15 = orient, age15 = age,
                         female15 = female, edu_group15 = edu_group,
                         cesd10_15 = cesd10, adl15 = adl)]
d18 <- g[year == 2018, .(ID, mem18 = mem, orient18 = orient)]
mm <- merge(merge(d15, d18, by = "ID"), grip, by = "ID", all.x = TRUE)
mm[, edu3 := fcase(edu_group15 == LV[1], "No formal schooling",
                   edu_group15 %in% LV[2:3], "Primary/middle",
                   edu_group15 == LV[4], "High school and above")]
mm <- merge(mm, cut3[, .(edu3, cutoff)], by = "edu3", all.x = TRUE)

say("\n样本流程：")
say("  2015 与 2018 记忆均可计分：", nrow(mm[!is.na(mem15) & !is.na(mem18)]))
say("  ＋2015 ≥50 岁：", nrow(mm[!is.na(mem15) & !is.na(mem18) & age15 >= 50]))
say("  ＋有切点：", nrow(mm[!is.na(mem15) & !is.na(mem18) & age15 >= 50 & !is.na(cutoff)]))
say("  ＋握力可得：", nrow(mm[!is.na(mem15) & !is.na(mem18) & age15 >= 50 & !is.na(cutoff) & !is.na(grip_max)]))
say("  ＋定向分可得：", nrow(mm[!is.na(mem15) & !is.na(mem18) & age15 >= 50 & !is.na(cutoff) &
                                 !is.na(grip_max) & !is.na(orient15)]))
d <- mm[!is.na(mem15) & !is.na(mem18) & age15 >= 50 & !is.na(cutoff) &
           !is.na(grip_max) & !is.na(orient15)]
d <- d[!is.na(cesd10_15) & !is.na(adl15)]
say("  ＋我们的协变量齐全：", nrow(d))

## 基线记忆正常者入组（他们的定义：≥切点为 intact）
d <- d[mem15 >= cutoff]
say("  ＋2015 记忆正常（≥教育分层切点）：", nrow(d))
d[, event := as.integer(mem18 < cutoff)]
d[, dep := as.integer(cesd10_15 >= 10)]
d[, edu_group := factor(edu_group15, levels = LV)]
d[, male := as.integer(female15 == 0)]
d[, agecat := fcase(age15 < 65, "lt65", age15 < 75, "6574", default = "ge75")]
say("\n分析样本 n = ", nrow(d), "｜结局（他们的口径）事件 = ", sum(d$event),
    "（", round(100 * mean(d$event), 2), "%）")
say("（对照：他们报告的时段验证队列 n=1,333，事件 210，15.8%）")

## ---------------- 5. 两套模型 ----------------
d[, lp_tu := B_TU[["age_lt65"]] * as.integer(agecat == "lt65") +
             B_TU[["age_ge75"]] * as.integer(agecat == "ge75") +
             B_TU[["male"]] * male +
             B_TU[["edu_primid"]] * as.integer(edu3 == "Primary/middle") +
             B_TU[["edu_none"]] * as.integer(edu3 == "No formal schooling") +
             B_TU[["grip"]] * grip_max + B_TU[["memory"]] * mem15 +
             B_TU[["orient"]] * orient15]
d[, lp_us := B_US[["(Intercept)"]] + B_US[["age"]] * age15 + B_US[["dep"]] * dep +
             B_US[["adl"]] * adl15 +
             B_US[["eduP"]] * as.integer(edu_group == LV[2]) +
             B_US[["eduM"]] * as.integer(edu_group == LV[3]) +
             B_US[["eduH"]] * as.integer(edu_group == LV[4])]
d[, p_tu := plogis(scale(lp_tu)[, 1])][, p_us := plogis(lp_us)]
## 简化评分卡：仅 3 项计分（年龄 +2/5 岁、教育 0/−2/−7/−10、抑郁 +3）；ADL 保留在连续模型但计 0 分
d[, score3 := 2L * pmax(floor((age15 - 45)/5), 0) +
      (-2L) * as.integer(edu_group == LV[2]) +
      (-7L) * as.integer(edu_group == LV[3]) +
      (-10L) * as.integer(edu_group == LV[4]) +
      3L * dep]

## 我们的**简化评分卡**（三项：年龄＋教育＋抑郁；ADL 计 0 分，等同不含）
d[, score_us := 2L * pmax(floor((age15 - 45)/5), 0) +
        (-2L) * as.integer(edu_group == LV[2]) + (-7L) * as.integer(edu_group == LV[3]) +
        (-10L) * as.integer(edu_group == LV[4]) + 3L * dep]

## ---------------- 6. 性能 ----------------
say("\n========== 主要对比（结局＝他们的口径）==========")
r_tu <- roc(d$event, d$lp_tu, quiet = TRUE, direction = "<")
r_us <- roc(d$event, d$lp_us, quiet = TRUE, direction = "<")
ci_tu <- as.numeric(ci.auc(r_tu)); ci_us <- as.numeric(ci.auc(r_us))
dt <- roc.test(r_us, r_tu, paired = TRUE)
say("Tu 6 因子：AUC = ", round(ci_tu[2], 3), "（", round(ci_tu[1], 3), "–", round(ci_tu[3], 3), "）")
say("我们 4 因子：AUC = ", round(ci_us[2], 3), "（", round(ci_us[1], 3), "–", round(ci_us[3], 3), "）")
say("ΔAUC（我们 − 他们）= ", round(ci_us[2] - ci_tu[2], 4),
    "｜配对 DeLong p = ", format(dt$p.value, digits = 3))
r_s3 <- roc(d$event, d$score3, quiet = TRUE, direction = "<")
ci_s3 <- as.numeric(ci.auc(r_s3))
say("我们简化评分卡（**3 项计分**，ADL 计 0）：AUC = ", round(ci_s3[2], 3),
    "（", round(ci_s3[1], 3), "–", round(ci_s3[3], 3), "）",
    "｜ΔAUC（3 项卡 − 他们）= ", round(ci_s3[2] - ci_tu[2], 4),
    "｜配对 DeLong p = ", format(roc.test(r_s3, r_tu, paired = TRUE)$p.value, digits = 3))

## 简化评分卡（三项，无 ADL）作为第三个臂
r_sc <- roc(d$event, d$score_us, quiet = TRUE, direction = "<")
ci_sc <- as.numeric(ci.auc(r_sc))
dt_sc <- roc.test(r_sc, r_tu, paired = TRUE)
say("我们 3 项评分卡（无 ADL）：AUC = ", round(ci_sc[2], 3), "（", round(ci_sc[1], 3), "–",
    round(ci_sc[3], 3), "）｜相对他们 ΔAUC = ", round(ci_sc[2] - ci_tu[2], 4),
    "｜配对 DeLong p = ", format(dt_sc$p.value, digits = 3))
say("（连续模型 vs 评分卡：ΔAUC = ", round(ci_us[2] - ci_sc[2], 4), "，即 ADL 的实际贡献）")

## 配对 Bootstrap（重抽样个体）
set.seed(20261014); B <- 2000; dv <- rep(NA_real_, B); au <- av <- rep(NA_real_, B)
for (i in seq_len(B)) {
  idx <- sample.int(nrow(d), replace = TRUE); db <- d[idx]
  if (length(unique(db$event)) < 2) next
  au[i] <- as.numeric(auc(roc(db$event, db$lp_tu, quiet = TRUE, direction = "<")))
  av[i] <- as.numeric(auc(roc(db$event, db$lp_us, quiet = TRUE, direction = "<")))
  dv[i] <- av[i] - au[i]
}
say("配对 Bootstrap 95% 区间：他们 [", round(quantile(au, .025, na.rm = TRUE), 3), ", ",
    round(quantile(au, .975, na.rm = TRUE), 3), "]；我们 [",
    round(quantile(av, .025, na.rm = TRUE), 3), ", ", round(quantile(av, .975, na.rm = TRUE), 3), "]")
say("ΔAUC 的 95% 区间：[", round(quantile(dv, .025, na.rm = TRUE), 4), ", ",
    round(quantile(dv, .975, na.rm = TRUE), 4), "]；P(ΔAUC>0) = ",
    round(mean(dv > 0, na.rm = TRUE), 3))

## 校准与 Brier 不作比较：他们的 CoxPH 线性预测子不含基线风险函数，
## 与我们的 logistic 绝对概率不在同一标尺上，比较校准会失真。
say("\n（校准与 Brier 不做对比：他们的 CoxPH 线性预测子没有基线风险函数，")
say("  两者绝对风险不在同一标尺上，只比较排序类指标。）")
cal_ok <- function(lp, y) coef(glm(y ~ as.numeric(scale(lp)), family = binomial()))
say("  参考：把各自 LP 标准化后的校准斜率——他们 ", round(cal_ok(d$lp_tu, d$event)[2], 3),
    "；我们 ", round(cal_ok(d$lp_us, d$event)[2], 3), "（越大者风险分散度越大）")

## 风险分层（各自四分位）
say("\n风险分层事件率（各自四分位）：")
for (v in c("lp_tu", "lp_us")) {
  d[, q := cut(get(v), quantile(get(v), seq(0, 1, .25)), include.lowest = TRUE, labels = 1:4)]
  say(paste0("  ", v, "：", paste(d[, .(事件率 = paste0(round(100*mean(event),1), "%"), n = .N),
                                    by = q][order(q)]$事件率, collapse = "  →  ")))
}

## ---------------- 7. 次要：换成我们的结局口径 ----------------
say("\n========== 次要：换成我们的口径（全局 mean−1SD 取整切点）==========")
mu15 <- mean(mm$mem15, na.rm = TRUE); sd15 <- sd(mm$mem15, na.rm = TRUE)
k15 <- floor(mu15 - sd15)
dd <- copy(d)[, event2 := as.integer(mem18 <= k15)]
say("2015 全局切点 = ≤", k15, "；事件 = ", sum(dd$event2), "（", round(100*mean(dd$event2),2), "%）")
r2t <- roc(dd$event2, dd$lp_tu, quiet = TRUE, direction = "<")
r2u <- roc(dd$event2, dd$lp_us, quiet = TRUE, direction = "<")
say("Tu 6 因子 AUC = ", round(as.numeric(auc(r2t)), 3),
    "；我们 4 因子 AUC = ", round(as.numeric(auc(r2u)), 3),
    "；ΔAUC = ", round(as.numeric(auc(r2u)) - as.numeric(auc(r2t)), 4),
    "｜配对 DeLong p = ", format(roc.test(r2u, r2t, paired = TRUE)$p.value, digits = 3))
r2s <- roc(dd$event2, dd$score3, quiet = TRUE, direction = "<")
say("我们 3 项评分卡 AUC = ", round(as.numeric(auc(r2s)), 3),
    "；ΔAUC（3 项卡 − 他们）= ", round(as.numeric(auc(r2s)) - as.numeric(auc(r2t)), 4),
    "｜配对 DeLong p = ", format(roc.test(r2s, r2t, paired = TRUE)$p.value, digits = 3))
r2s <- roc(dd$event2, dd$score_us, quiet = TRUE, direction = "<")
say("Tu 6 因子 AUC = ", round(as.numeric(auc(r2t)), 3),
    "；我们 3 项评分卡 AUC = ", round(as.numeric(auc(r2s)), 3),
    "；ΔAUC = ", round(as.numeric(auc(r2s)) - as.numeric(auc(r2t)), 4),
    "｜配对 DeLong p = ", format(roc.test(r2s, r2t, paired = TRUE)$p.value, digits = 3))

## ---------------- 8. 敏感性：握力定义 ----------------

## ---------------- 9. 每个模型"优于随机"的检验 ----------------
say("\n========== 附：每个模型自身 vs 随机（AUC=0.5）的检验 ==========")
auc_vs_chance <- function(y, p, tag) {
  r <- roc(y, p, quiet = TRUE, direction = "<")
  ci <- as.numeric(ci.auc(r))
  se <- (ci[3] - ci[1]) / (2 * 1.96)
  z <- (ci[2] - 0.5) / se
  data.table(模型 = tag, AUC = round(ci[2], 3), SE = round(se, 4),
             z = round(z, 1), p_vs_0.5 = format(2 * pnorm(-abs(z)), scientific = TRUE, digits = 3))
}
say("\n-- 结局＝他们的口径 --")
pt(rbind(auc_vs_chance(d$event, d$lp_tu, "Tu 6 因子"),
         auc_vs_chance(d$event, d$score_us, "我们 3 项评分卡"),
         auc_vs_chance(d$event, d$lp_us, "我们连续模型")))
say("-- 结局＝我们的口径 --")
pt(rbind(auc_vs_chance(dd$event2, dd$lp_tu, "Tu 6 因子"),
         auc_vs_chance(dd$event2, dd$score_us, "我们 3 项评分卡"),
         auc_vs_chance(dd$event2, dd$lp_us, "我们连续模型")))
say("\n非劣性（他们口径）：ΔAUC 的 95% 区间下限 = ", round(quantile(dv, .025, na.rm = TRUE), 4),
    "；若预设非劣界值为 −0.02，则区间完全落在界值之上，可判定非劣。")

say("\n========== 敏感性：握力改用「两手最好值的均值」==========")
grip2 <- b[, .(ID = as.character(ID),
               g1 = pmax(qc003, qc004, na.rm = TRUE),
               g2 = pmax(qc005, qc006, na.rm = TRUE))]
grip2[, grip_alt := rowMeans(cbind(g1, g2), na.rm = TRUE)]
d2 <- merge(d[, .(ID, event, lp_tu, lp_us, mem15, orient15, age15, male, edu3)], grip2,
            by = "ID")
d2 <- d2[!is.na(grip_alt)]
d2[, lp_tu_alt := B_TU[["age_lt65"]] * as.integer(age15 < 65) +
                   B_TU[["age_ge75"]] * as.integer(age15 >= 75) +
                   B_TU[["male"]] * as.integer(male) +
                   B_TU[["edu_primid"]] * as.integer(edu3 == "Primary/middle") +
                   B_TU[["edu_none"]] * as.integer(edu3 == "No formal schooling") +
                   B_TU[["grip"]] * grip_alt + B_TU[["memory"]] * mem15 +
                   B_TU[["orient"]] * orient15]
say("n = ", nrow(d2), "｜Tu 模型（握力取两手最好值均值）AUC = ",
    round(as.numeric(auc(roc(d2$event, d2$lp_tu_alt, quiet = TRUE, direction = "<"))), 3),
    "；我们 4 因子 AUC = ",
    round(as.numeric(auc(roc(d2$event, d2$lp_us, quiet = TRUE, direction = "<"))), 3))

## 导出
out <- data.table(
  口径 = c("Tu（教育分层切点）", "Tu（教育分层切点）", "我们（全局切点 ≤", "我们（全局切点 ≤"),
  模型 = c("Tu 6 因子", "我们 4 因子", "Tu 6 因子", "我们 4 因子")
)
out <- data.table(
  结局口径 = rep(c("Tu 教育分层切点", "我们全局切点"), each = 2),
  模型 = rep(c("Tu 6 因子（含认知测验）", "我们 4 因子（无认知测验）"), 2),
  AUC = c(round(ci_tu[2], 3), round(ci_us[2], 3),
          round(as.numeric(auc(r2t)), 3), round(as.numeric(auc(r2u)), 3)),
  n = nrow(d), 事件 = c(sum(d$event), sum(d$event), sum(dd$event2), sum(dd$event2)))
fwrite(out, file.path(OUT, "F07_headtohead.csv"), bom = TRUE)

close(con)
cat("\n完成 ->", file.path(OUT, "F07_report.txt"), "\n")
