# ============================================================================
# F02 哑火检查：能否用 CHARLS 2018→2020 建立/验证模型
# 输入：waves_cognition_2026-09-13.rds（五波逐波数据）、2020 原始 HSF（词表版本）
# 输出：C:/y03tmp/out/F02_report.txt、F02_*.csv
# 要点：冻结模型 = A14（2011→2015，4 年记忆跨阈值）
#       结局口径：2018 记忆分（即时试次1 + 延迟，0–20）→ 2020 同分跨阈值
#       切点：按原设计逻辑用 2018 自身 mean−1SD 取整；另报"照搬 2011 切点 ≤4"作对照
# ============================================================================
invisible(try(Sys.setlocale("LC_ALL", "zh_CN.UTF-8"), silent = TRUE))
suppressPackageStartupMessages({ library(data.table); library(pROC); library(haven) })
if (dir.exists("D:/R/library")) .libPaths(c("D:/R/library", .libPaths()))

TMP <- "C:/y03tmp/tmp"; OUT <- "C:/y03tmp/out"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
con <- file(file.path(OUT, "F02_report.txt"), open = "wt", encoding = "UTF-8")
say <- function(...) { s <- paste0(..., collapse = ""); cat(s, "\n"); cat(s, "\n", file = con) }
pt <- function(x) say(paste(capture.output(print(x)), collapse = "\n"))

## ---------- 冻结模型系数（A14_最终模型系数.csv） ----------
beta <- c("(Intercept)" = -4.19524508778947,
          age = 0.0492039250378111,
          eduPrimary = -0.290105771000328,
          eduMiddle = -0.994711800202235,
          eduHigh = -1.43588388035287,
          dep = 0.412166484954086,
          adl = -0.0110119344848014)

## ---------- 1. 数据 ----------
w <- as.data.table(readRDS(file.path(TMP, "waves_cognition.rds")))
say("五波数据集行数：", nrow(w), "｜年份：", paste(sort(unique(w$year)), collapse = ", "))

d18 <- w[year == 2018]; d20 <- w[year == 2020]
say("2018 行数：", nrow(d18), "｜2020 行数：", nrow(d20))

a18 <- d18[, .(ID = IDn, imrc18 = as.numeric(imrc), dlrc18 = as.numeric(dlrc),
               age18 = as.numeric(age), edu18 = edu_group, cesd18 = as.numeric(cesd10),
               adl18 = as.numeric(adl), wt18 = as.numeric(weight))]
a18[, mem18 := imrc18 + dlrc18]
a20 <- d20[, .(ID = IDn, imrc20 = as.numeric(imrc), dlrc20 = as.numeric(dlrc))]
a20[, mem20 := imrc20 + dlrc20]
m <- merge(a18, a20, by = "ID", all = FALSE)
say("2018 与 2020 双波齐全（按 ID）：", nrow(m), " 人")

## ---------- 2. 逐步样本流程 ----------
say("\n========== 2. 样本流程 ==========")
say("① 2018 记忆分可计分：", sum(!is.na(m$mem18)))
say("② 2020 记忆分可计分：", sum(!is.na(m$mem20)))
say("③ 两波均可计分：", sum(!is.na(m$mem18) & !is.na(m$mem20)))
say("④ ③ 且 2018 四个协变量齐全：",
    sum(!is.na(m$mem18) & !is.na(m$mem20) & !is.na(m$age18) & !is.na(m$edu18) &
          !is.na(m$cesd18) & !is.na(m$adl18)))
say("   2018 教育缺失：", sum(is.na(m$edu18)), "｜CES-D 缺失：", sum(is.na(m$cesd18)),
    "｜ADL 缺失：", sum(is.na(m$adl18)), "｜年龄缺失：", sum(is.na(m$age18)))

## ---------- 3. 两波记忆分分布对照 ----------
say("\n========== 3. 记忆分分布对照（同一批人）==========")
pp <- m[!is.na(mem18) & !is.na(mem20)]
say("配对人数：", nrow(pp))
tab <- data.table(
  指标 = c("记忆分(0-20)", "即时回忆(0-10)", "延迟回忆(0-10)"),
  m2018 = round(c(mean(pp$mem18), mean(pp$imrc18), mean(pp$dlrc18)), 2),
  sd2018 = round(c(sd(pp$mem18), sd(pp$imrc18), sd(pp$dlrc18)), 2),
  m2020 = round(c(mean(pp$mem20), mean(pp$imrc20), mean(pp$dlrc20)), 2),
  sd2020 = round(c(sd(pp$mem20), sd(pp$imrc20), sd(pp$dlrc20)), 2),
  diff = round(c(mean(pp$mem20 - pp$mem18), mean(pp$imrc20 - pp$imrc18),
                 mean(pp$dlrc20 - pp$dlrc18)), 2),
  diff_sd = round(c(sd(pp$mem20 - pp$mem18), sd(pp$imrc20 - pp$imrc18),
                    sd(pp$dlrc20 - pp$dlrc18)), 2))
pt(tab)
say("配对差值 t 检验（记忆分）：p = ",
    format(t.test(pp$mem20, pp$mem18, paired = TRUE)$p.value, scientific = TRUE, digits = 3))
fwrite(tab, file.path(OUT, "F02_memory_distribution.csv"), bom = TRUE)

## ---------- 4. 切点与事件 + 冻结模型搬运 ----------
say("\n========== 4. 切点、事件与冻结模型搬运 ==========")
cut18 <- floor(mean(pp$mem18) - sd(pp$mem18))
say("2018 记忆分：mean=", round(mean(pp$mem18), 2), " SD=", round(sd(pp$mem18), 2),
    " → mean−1SD 取整 = ≤", cut18)
say("（原开发口径 2011 切点 = ≤4）")
say("2018 受损比例（≤", cut18, "）：", round(100 * mean(pp$mem18 <= cut18), 1), "%")

res <- data.table()
d18_keep <- NULL
for (k in c(cut18, 4)) {
  nm <- if (k == cut18) "2018自身切点" else "照搬2011切点"
  d <- m[!is.na(mem18) & !is.na(mem20) & !is.na(age18) & !is.na(edu18) &
           !is.na(cesd18) & !is.na(adl18)]
  d <- d[mem18 > k]
  d[, event := as.integer(mem20 <= k)]
  d[, dep := as.integer(cesd18 >= 10)]
  d[, edu_group := factor(edu18, levels = c("Less than primary", "Primary school",
                                            "Middle school", "High school or above"))]
  d[, lp := beta[["(Intercept)"]] + beta[["age"]] * age18 + beta[["dep"]] * dep +
        beta[["adl"]] * adl18 +
        beta[["eduPrimary"]] * as.integer(edu_group == "Primary school") +
        beta[["eduMiddle"]] * as.integer(edu_group == "Middle school") +
        beta[["eduHigh"]] * as.integer(edu_group == "High school or above")]
  d[, phat := plogis(lp)]
  r <- roc(d$event, d$phat, quiet = TRUE, direction = "<")
  ci <- as.numeric(ci.auc(r))
  cal <- glm(event ~ lp, data = d, family = binomial())
  brier <- mean((d$phat - d$event)^2)
  say("\n-- ", nm, "（≤", k, "）--")
  say("   风险人群 n=", nrow(d), "（2018 未受损）｜事件=", sum(d$event),
      "｜事件率=", round(100 * mean(d$event), 2), "%")
  say("   冻结模型：AUC=", round(ci[2], 3), " (", round(ci[1], 3), "–", round(ci[3], 3),
      ")｜校准斜率=", round(coef(cal)[2], 3), " 截距=", round(coef(cal)[1], 3),
      "｜Brier=", round(brier, 4))
  res <- rbind(res, data.table(口径 = nm, 切点 = k, n = nrow(d), 事件 = sum(d$event),
                               事件率 = round(100 * mean(d$event), 2),
                               AUC = round(ci[2], 3), AUC低 = round(ci[1], 3), AUC高 = round(ci[3], 3),
                               校准斜率 = round(coef(cal)[2], 3), 校准截距 = round(coef(cal)[1], 3),
                               Brier = round(brier, 4)))
  if (k == cut18) d18_keep <- copy(d)
}
fwrite(res, file.path(OUT, "F02_frozen_model_transport.csv"), bom = TRUE)

## ---------- 5. 同批人重估权重（预览）----------
say("\n========== 5. 同批人重估：2018 权重 vs 2011 冻结权重（预览）==========")
d <- copy(d18_keep)
fit <- glm(event ~ age18 + edu_group + dep + adl18, data = d, family = binomial())
say("重估系数："); pt(round(coef(fit), 4))
d[, phat_new := predict(fit, type = "response")]
r_old <- roc(d$event, d$phat, quiet = TRUE, direction = "<")
r_new <- roc(d$event, d$phat_new, quiet = TRUE, direction = "<")
ci_old <- as.numeric(ci.auc(r_old)); ci_new <- as.numeric(ci.auc(r_new))
say("旧权重 AUC=", round(ci_old[2], 3), " (", round(ci_old[1], 3), "–", round(ci_old[3], 3), ")")
say("新权重 AUC=", round(ci_new[2], 3), " (", round(ci_new[1], 3), "–", round(ci_new[3], 3), ")")
say("配对 DeLong：p = ", format(roc.test(r_old, r_new, paired = TRUE)$p.value, digits = 3))
best <- function(r) {
  co <- coords(r, "best", best.method = "youden", ret = c("threshold", "sensitivity", "specificity"))
  data.table(阈值 = round(co$threshold, 4), 敏感度 = round(co$sensitivity, 3), 特异度 = round(co$specificity, 3))
}
say("\n旧权重（约登切点）："); pt(best(r_old))
say("新权重（约登切点）："); pt(best(r_new))

## ---------- 6. 2020 词表版本 A/B 与记忆分 ----------
say("\n========== 6. 2020 词表版本 A/B 与记忆分 ==========")
x5 <- as.data.table(read_dta(file.path(TMP, "w5_hsf.dta"),
                             col_select = all_of(c("ID", paste0("xwordlist_", 1:10, "_")))))
say("2020 原始文件 ID 样例：", paste(head(as.character(x5$ID), 3), collapse = ", "))
x5[, ver := fifelse(as.character(xwordlist_1_) == "\u6cb3\u6d41", "A",
                    fifelse(as.character(xwordlist_1_) == "\u5929\u7a7a", "B", "未施测/缺失"))]
x5 <- x5[, .(ID = as.character(ID), ver)]
mm <- merge(m[!is.na(mem18) & !is.na(mem20)], x5, by = "ID", all.x = TRUE)
pt(mm[, .(n = .N, mem2018 = round(mean(mem18), 2), mem2020 = round(mean(mem20), 2),
          即时2020 = round(mean(imrc20), 2), 延迟2020 = round(mean(dlrc20), 2)), by = ver][order(ver)])

close(con)
cat("\n完成 ->", file.path(OUT, "F02_report.txt"), "\n")
