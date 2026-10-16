# ============================================================================
# F05：结局的测量问题——错分有多少、换口径后模型表现如何、模型究竟在预测什么
# 数据：waves_cognition（2018 与 2020 记忆分，即时试次1 + 延迟，0–20）
# 输出：C:/y03tmp/out/F05_report.txt、F05_*.csv、F05_*.png
# 说明：全部为事后分析；口径预先指定，不做事后挑选
# ============================================================================
invisible(try(Sys.setlocale("LC_ALL", "zh_CN.UTF-8"), silent = TRUE))
suppressPackageStartupMessages({ library(data.table); library(pROC) })
if (dir.exists("D:/R/library")) .libPaths(c("D:/R/library", .libPaths()))
TMP <- "C:/y03tmp/tmp"; OUT <- "C:/y03tmp/out"
con <- file(file.path(OUT, "F05_report.txt"), open = "wt", encoding = "UTF-8")
say <- function(...) { s <- paste0(..., collapse = ""); cat(s, "\n"); cat(s, "\n", file = con) }
pt  <- function(x) say(paste(capture.output(print(x)), collapse = "\n"))

BETA_OLD <- c("(Intercept)" = -4.19524508778947, age = 0.0492039250378111,
              eduP = -0.290105771000328, eduM = -0.994711800202235, eduH = -1.43588388035287,
              dep = 0.412166484954086, adl = -0.0110119344848014)
LV <- c("Less than primary", "Primary school", "Middle school", "High school or above")

## ---------------- 数据 ----------------
w <- as.data.table(readRDS(file.path(TMP, "waves_cognition.rds")))
a18 <- w[year == 2018, .(ID = IDn, mem18 = as.numeric(imrc) + as.numeric(dlrc),
                         age18 = as.numeric(age), edu18 = edu_group,
                         cesd18 = as.numeric(cesd10), adl18 = as.numeric(adl),
                         wt18 = as.numeric(weight))]
a20 <- w[year == 2020, .(ID = IDn, mem20 = as.numeric(imrc) + as.numeric(dlrc))]
m <- merge(a18, a20, by = "ID")
pt_paired <- m[!is.na(mem18) & !is.na(mem20)]
say("配对样本 n = ", nrow(pt_paired))

## 波内标准化（用配对样本的波内均值/SD）
mu18 <- mean(pt_paired$mem18); sd18 <- sd(pt_paired$mem18)
mu20 <- mean(pt_paired$mem20); sd20 <- sd(pt_paired$mem20)
pt_paired[, z18 := (mem18 - mu18)/sd18][, z20 := (mem20 - mu20)/sd20]
pt_paired[, dz := z20 - z18][, draw := mem20 - mem18]
say("2018：mean=", round(mu18, 2), " SD=", round(sd18, 2),
    "｜2020：mean=", round(mu20, 2), " SD=", round(sd20, 2))
say("变化量（0–20 分）：mean=", round(mean(pt_paired$draw), 2),
    " SD=", round(sd(pt_paired$draw), 2))
say("标准化变化 dz：mean=", round(mean(pt_paired$dz), 3), " SD=", round(sd(pt_paired$dz), 3))

## ---------------- 1. 两波相关、信度与最小可测变化 ----------------
say("\n========== 1. 测量信度与最小可测变化 ==========")
r <- cor(pt_paired$mem18, pt_paired$mem20)
sem <- sd18 * sqrt(1 - r)
mdc <- 1.96 * sqrt(2) * sem
say("两波相关 r = ", round(r, 3), "（注意：2018/2020 换过词表，此值会被压低，")
say("      按 ρ≈r 处理即高估测量误差，属保守方向）")
say("估计 SEM = ", round(sem, 2), " 分；最小可测变化 MDC95 = ", round(mdc, 2), " 分")

## ---------------- 2. 转归矩阵与切点附近的往返 ----------------
say("\n========== 2. 转归矩阵（配对样本）==========")
cut_pop <- floor(mu18 - sd18)
pt_paired[, g18 := cut(mem18, c(-Inf, cut_pop, 6, 9, Inf),
                       labels = c(paste0("≤", cut_pop), "4-6", "7-9", "≥10"))]
pt_paired[, g20 := cut(mem20, c(-Inf, cut_pop, 6, 9, Inf),
                       labels = c(paste0("≤", cut_pop), "4-6", "7-9", "≥10"))]
tt <- dcast(pt_paired[, .(n = .N), by = .(g18, g20)], g18 ~ g20, value.var = "n", fill = 0)
pt(tt)
say("\n按 2018 分数的分组，2020 跨阈值（≤", cut_pop, "）的比例：")
say("（只列 2018 未受损者）")
sub <- pt_paired[mem18 > cut_pop]
sub[, g := cut(mem18, c(cut_pop, 6, 9, Inf), labels = c("4-6（紧贴切点上方）", "7-9", "≥10"))]
pt(sub[, .(n = .N, 事件 = sum(mem20 <= cut_pop),
           跨线率 = paste0(round(100*mean(mem20 <= cut_pop), 1), "%")), by = g][order(g)])
say("→ 紧贴切点者的跨线率是高分组的 ",
    round(sub[mem18 %in% 4:6, mean(mem20 <= cut_pop)] / sub[mem18 >= 10, mean(mem20 <= cut_pop)], 1),
    " 倍")

## 事件中的下降幅度：有多少超过最小可测变化
ev <- pt_paired[mem18 > cut_pop & mem20 <= cut_pop]
say("\n跨阈值事件 n = ", nrow(ev),
    "；其中下降幅度 ≥ MDC95（", round(mdc, 2), " 分）者 ",
    sum(ev$draw <= -mdc), " 例（", round(100*mean(ev$draw <= -mdc), 1), "%）")
say("   下降幅度中位数 = ", round(median(ev$draw), 1), " 分；≥5 分者 ",
    sum(ev$draw <= -5), " 例（", round(100*mean(ev$draw <= -5), 1), "%）")
say("→ 其余跨线者的下降幅度在测量波动可解释范围内（含回归均值）")
say("\n对照：2018 未受损者中，非事件者的分数变化中位数 = ",
    round(median(sub[mem20 > cut_pop]$draw), 1), " 分")

## ---------------- 3. 不同结局口径下的模型表现 ----------------
say("\n========== 3. 不同结局口径下的模型表现（冻结 vs 同批人重估）==========")
## m 本身已含协变量；直接在其上补 z 分与变化量，避免 merge 产生 .x/.y 后缀
mm <- m[!is.na(mem18) & !is.na(mem20)]
mm[, `:=`(z18 = (mem18 - mu18)/sd18, z20 = (mem20 - mu20)/sd20)]
mm[, `:=`(dz = z20 - z18, draw = mem20 - mem18)]
mm_full <- copy(mm)
mm <- mm[!is.na(age18) & !is.na(edu18) & !is.na(cesd18) & !is.na(adl18) & age18 >= 45]
mm_full <- mm_full[!is.na(age18) & !is.na(edu18) & !is.na(cesd18) & !is.na(adl18) & age18 >= 45]
mm[, dep := as.integer(cesd18 >= 10)]
mm[, edu_group := factor(edu18, levels = LV)]
mm[, p_frozen := plogis(BETA_OLD[["(Intercept)"]] + BETA_OLD[["age"]]*age18 +
       BETA_OLD[["dep"]]*dep + BETA_OLD[["adl"]]*adl18 +
       BETA_OLD[["eduP"]]*as.integer(edu_group == LV[2]) +
       BETA_OLD[["eduM"]]*as.integer(edu_group == LV[3]) +
       BETA_OLD[["eduH"]]*as.integer(edu_group == LV[4]))]

## 等百分位锚定的 2020 阈值
p_anchor <- mean(pt_paired$mem18 <= cut_pop)
cut20p <- as.numeric(quantile(pt_paired$mem20, p_anchor, type = 1))
say("O2 等百分位锚定：2018 ≤", cut_pop, " 对应百分位 ", round(100*p_anchor, 1),
    "% → 2020 切点 ≤", cut20p)

mm[, o1 := as.integer(mem20 <= cut_pop)]
mm[, o2 := as.integer(mem20 <= cut20p)]
mm[, o3 := as.integer(dz <= -1.0)]
mm[, o4 := as.integer(dz <= -0.5)]
mm[, o6 := as.integer(mem20 <= cut_pop & draw <= -5)]
mm_sub <- mm[mem18 > cut_pop]                 # 人群 A：2018 未受损者
mm_risk <- copy(mm_sub)                       # 跨线类结局的评估人群
mm_full2 <- mm                               # 人群 B：全部配对者（用于下降类结局，避免截断）

aucv <- function(y, p) if (length(unique(y)) < 2) NA_real_ else
  as.numeric(auc(roc(y, p, quiet = TRUE, direction = "<")))
eval_def <- function(col, tag, data = mm_risk) {
  d <- copy(data)[, event := get(col)]
  fit <- glm(event ~ age18 + edu_group + dep + adl18, data = d, family = binomial())
  d[, p_new := predict(fit, newdata = d, type = "response")]
  cal <- glm(event ~ qlogis(pmin(pmax(p_frozen, 1e-9), 1 - 1e-9)), data = d, family = binomial())
  a_old <- aucv(d$event, d$p_frozen); a_new <- aucv(d$event, d$p_new)
  data.table(口径 = tag, 人群 = if (nrow(d) == nrow(mm_risk)) "A 未受损者" else "B 全部配对者",
             n = nrow(d), 事件 = sum(d$event),
             事件率 = round(100*mean(d$event), 2),
             冻结AUC = round(a_old, 4), 重估AUC = round(a_new, 4),
             ΔAUC = round(a_new - a_old, 4),
             冻结校准斜率 = round(coef(cal)[2], 3), 冻结校准截距 = round(coef(cal)[1], 3),
             Brier冻结 = round(mean((d$p_frozen - d$event)^2), 4),
             Brier重估 = round(mean((d$p_new - d$event)^2), 4))
}
defs <- rbind(
  eval_def("o1", paste0("O1 绝对切点 ≤", cut_pop, "（现状）"), mm_risk),
  eval_def("o2", "O2 等百分位锚定", mm_risk),
  eval_def("o6", "O6 跨阈值且下降 ≥5 分", mm_risk),
  eval_def("o3", "O3 波内标准化下降 ≥1 SD", mm_full2),
  eval_def("o4", "O4 波内标准化下降 ≥0.5 SD", mm_full2),
  eval_def("o3", "O3 波内标准化下降 ≥1 SD（限未受损者）", mm_risk),
  eval_def("o4", "O4 波内标准化下降 ≥0.5 SD（限未受损者）", mm_risk))
pt(defs)
fwrite(defs, file.path(OUT, "F05_outcome_definitions.csv"), bom = TRUE)

## ---------------- 4. 模型究竟在预测什么（R²） ----------------
say("\n========== 4. 四变量能解释什么：水平 vs 变化 ==========")
r2 <- function(form, data) summary(lm(form, data = data))$r.squared
mm_full[, `:=`(dep = as.integer(cesd18 >= 10),
               edu_group = factor(edu18, levels = LV))]
mm_sub <- copy(mm_full)[mem18 > cut_pop]
r2tab <- data.table(
  样本 = rep(c("2018 未受损者", "全部配对者"), each = 3),
  结局 = rep(c("2020 记忆分（水平）", "2 年变化量（0–20 分）", "2 年标准化变化 dz"), 2),
  R2 = c(r2(mem20 ~ age18 + edu_group + dep + adl18, mm_sub),
         r2(draw ~ age18 + edu_group + dep + adl18, mm_sub),
         r2(dz ~ age18 + edu_group + dep + adl18, mm_sub),
         r2(mem20 ~ age18 + edu_group + dep + adl18, mm_full),
         r2(draw ~ age18 + edu_group + dep + adl18, mm_full),
         r2(dz ~ age18 + edu_group + dep + adl18, mm_full)))
r2tab[, R2 := round(R2, 4)]
pt(r2tab)
fwrite(r2tab, file.path(OUT, "F05_R2.csv"), bom = TRUE)

say("\n四变量对「2020 水平」与「2 年变化量」的线性系数（同一量纲 0–20 分，全部配对者）：")
lv <- summary(lm(mem20 ~ age18 + edu_group + dep + adl18, data = mm_full2))$coefficients
ch <- summary(lm(draw ~ age18 + edu_group + dep + adl18, data = mm_full2))$coefficients
cm <- data.table(项 = rownames(lv), 对2020水平 = round(lv[, 1], 3), 对2年变化量 = round(ch[, 1], 3))
pt(cm)
fwrite(cm, file.path(OUT, "F05_coefficients_level_vs_change.csv"), bom = TRUE)

## ---------------- 图 ----------------
png(file.path(OUT, "F05_fig_change_dist.png"), width = 1150, height = 560, res = 120)
par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 1.5))
hist(pt_paired$draw, breaks = seq(-20.5, 20.5, 1), col = "grey80", border = "white",
     xlab = "Change in memory score 2018 to 2020", ylab = "Count",
     main = "Distribution of 2-year change (n = 12,571)")
abline(v = 0, lty = 3); abline(v = -mdc, col = "red", lwd = 2)
legend("topleft", c("MDC95", "no change"), col = c("red", "black"), lwd = c(2, 1),
       lty = c(1, 3), bty = "n", cex = 0.8)
lvl <- c(paste0("≤", cut_pop), "4-6", "7-9", "≥10")
bandmat <- cbind(w2018 = as.integer(table(factor(pt_paired$g18, levels = lvl))),
                 w2020 = as.integer(table(factor(pt_paired$g20, levels = lvl))))
barplot(t(bandmat), beside = TRUE, col = c("grey60", "grey25"), xlab = "Memory score band",
        ylab = "Count", main = "Band distribution: 2018 vs 2020", las = 1)
legend("topright", c("2018", "2020"), fill = c("grey60", "grey25"), bty = "n", cex = 0.8)
dev.off()

close(con)
cat("\n完成 ->", file.path(OUT, "F05_report.txt"), "\n")
