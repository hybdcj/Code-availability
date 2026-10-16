# ============================================================================
# F16 分析：外部队列（HRS / ELSA / SHARE）全人群口径
# 与 A11/A12 唯一的差别：不再要求基线记忆可得、不再要求基线未受损
# 切点定义与 A11 v2 主口径（baseline_all_mean_sd）一致，便于直接比较
# ============================================================================
invisible(try(Sys.setlocale("LC_ALL", "zh_CN.UTF-8"), silent = TRUE))
suppressPackageStartupMessages({ library(data.table); library(pROC) })
if (dir.exists("D:/R/library")) .libPaths(c("D:/R/library", .libPaths()))
TMP <- "C:/y03tmp/tmp"; OUT <- "C:/y03tmp/out"
con <- file(file.path(OUT, "F16_report.txt"), open = "wt", encoding = "UTF-8")
say <- function(...) { s <- paste0(..., collapse = ""); cat(s, "\n"); cat(s, "\n", file = con) }
pt <- function(x) say(paste(capture.output(print(x)), collapse = "\n"))

B <- c(intercept = -4.19524508778947, age = 0.0492039250378111,
       eduPri = -0.290105771000328, eduMid = -0.994711800202235,
       eduHigh = -1.43588388035287, dep = 0.412166484954086, adl = -0.0110119344848014)
CO <- c("HRS", "ELSA", "SHARE")
LB <- c(U = "U_全人口", N = "N_基线未受损", I = "I_基线已受损", M = "M_基线记忆缺失")

prep <- function(co) {
  d <- as.data.table(readRDS(file.path(TMP, paste0("F16_base_", co, ".rds"))))
  d <- d[!is.na(age) & age >= 45]
  d[, score := 2L * pmax(floor((age - 45) / 5), 0) + edu_pts + 3L * dep]
  d[, lp := B[["intercept"]] + B[["age"]] * age + B[["dep"]] * dep + B[["adl"]] * adl +
      B[["eduPri"]] * as.integer(edu_pts == -2) +
      B[["eduMid"]] * as.integer(edu_pts == -7) +
      B[["eduHigh"]] * as.integer(edu_pts == -10)]
  d[, phat := plogis(lp)]
  cs <- d[!is.na(mem_base)]
  d[, cut := mean(cs$mem_base) - sd(cs$mem_base)]
  d[, cut_int := floor(cut)]
  d[, event := as.integer(mem_fu <= cut_int)]
  d[, low0 := as.integer(mem_base <= cut_int)]
  d[, cohort := co]
  d[]
}
D <- lapply(CO, prep); names(D) <- CO

say("切点核对（期望 HRS 6.244/6、ELSA 6.77/6、SHARE 5.28/5）：")
pt(rbindlist(lapply(CO, function(co) data.table(
  队列 = co, 切点原始 = round(unique(D[[co]]$cut), 3), 整数阈值 = unique(D[[co]]$cut_int),
  mem_base可得 = sum(!is.na(D[[co]]$mem_base)), 年龄教育抑郁齐全 = nrow(D[[co]])))))

say("\n========== 1. 管线校验（复现 A11 主风险集）==========")
say("期望：HRS 14051 / 事件 1842；ELSA 6725 / 554；SHARE 26879 / 2930")
ver <- rbindlist(lapply(CO, function(co) {
  d <- D[[co]][!is.na(mem_base) & mem_base > cut_int & !is.na(mem_fu)]
  data.table(队列 = co, 风险集n = nrow(d), 事件 = sum(d$event),
             AUC_评分卡 = round(as.numeric(auc(roc(d$event, d$score, quiet = TRUE, direction = "<"))), 4))
}))
pt(ver); fwrite(ver, file.path(OUT, "F16_verify_A11.csv"), bom = TRUE)

fl <- rbindlist(lapply(CO, function(co) {
  d <- D[[co]]; fu <- !is.na(d$mem_fu)
  data.table(队列 = co, 年龄教育抑郁齐全 = nrow(d), 随访记忆可计分 = sum(fu),
             U_全人口 = sum(fu),
             U_基线记忆缺失 = sum(fu & is.na(d$mem_base)),
             N_基线未受损 = sum(fu & !is.na(d$mem_base) & d$mem_base > d$cut_int),
             I_基线已受损 = sum(fu & !is.na(d$mem_base) & d$mem_base <= d$cut_int))
}))
say("\n========== 2. 样本流程 ==========")
pt(fl); fwrite(fl, file.path(OUT, "F16_flow.csv"), bom = TRUE)

aucrow <- function(d, v, co, arm, vlab) {
  ok <- !is.na(d[[v]]) & !is.na(d$event)
  if (sum(ok) < 30 || length(unique(d$event[ok])) < 2) return(NULL)
  r <- roc(d$event[ok], d[[v]][ok], quiet = TRUE, direction = "<")
  ci <- as.numeric(ci.auc(r))
  data.table(队列 = co, 口径 = arm, 预测子 = vlab, n = sum(ok), 事件 = sum(d$event[ok]),
             事件率 = round(100 * mean(d$event[ok]), 2),
             AUC = round(ci[2], 4), AUC低 = round(ci[1], 4), AUC高 = round(ci[3], 4))
}
mk_sets <- function(d) {
  list(U = d,
       N = d[!is.na(mem_base) & mem_base > cut_int],
       I = d[!is.na(mem_base) & mem_base <= cut_int],
       M = d[is.na(mem_base)])
}
aucT <- rbindlist(lapply(CO, function(co) {
  st <- mk_sets(D[[co]][!is.na(mem_fu)])
  rbindlist(lapply(names(st), function(s) {
    if (nrow(st[[s]]) == 0) return(NULL)
    rbind(aucrow(st[[s]], "score", co, LB[[s]], "三项评分卡"),
          aucrow(st[[s]], "phat", co, LB[[s]], "连续模型"))
  }))
}))
say("\n========== 3. 各口径 AUC ==========")
pt(aucT); fwrite(aucT, file.path(OUT, "F16_auc.csv"), bom = TRUE)

dec <- rbindlist(lapply(CO, function(co) {
  d <- D[[co]][!is.na(mem_fu) & !is.na(mem_base)]
  r_low <- roc(d$event, d$low0, quiet = TRUE, direction = "<")
  r_scr <- roc(d$event, d$score, quiet = TRUE, direction = "<")
  fit <- glm(event ~ score + low0, data = d, family = binomial())
  r_cmb <- roc(d$event, predict(fit, type = "response"), quiet = TRUE, direction = "<")
  data.table(队列 = co, n = nrow(d), 事件 = sum(d$event),
             事件中基线已低分占比 = round(100 * mean(d$low0[d$event == 1]), 1),
             AUC_仅基线已低分 = round(as.numeric(auc(r_low)), 4),
             AUC_仅评分卡 = round(as.numeric(auc(r_scr)), 4),
             AUC_合并 = round(as.numeric(auc(r_cmb)), 4),
             P_合并优于评分卡 = signif(roc.test(r_cmb, r_scr, paired = TRUE)$p.value, 3))
}))
say("\n========== 4. 分解：基线已低分的预测力 ==========")
pt(dec); fwrite(dec, file.path(OUT, "F16_decomposition.csv"), bom = TRUE)

scr <- rbindlist(lapply(CO, function(co) {
  d <- D[[co]][!is.na(mem_fu)]
  s <- d$score; ev <- d$event; n <- length(ev)
  out <- rbindlist(lapply(c(0.10, 0.20, 0.30), function(fr) {
    k <- floor(fr * n); th <- sort(s, decreasing = TRUE)[k]; pos <- s >= th
    sub <- d[pos]; yld <- sum(sub$event == 1 & sub$low0 == 0, na.rm = TRUE)
    data.table(队列 = co, 目标筛查比例 = paste0(fr * 100, "%"), 阈值 = th, n = length(ev), 事件 = sum(ev), 检出人数 = sum(pos), 实际筛查比例 = round(100 * mean(pos), 2), TP = sum(pos & ev == 1), FP = sum(pos & ev == 0), TN = sum(!pos & ev == 0), FN = sum(!pos & ev == 1),
               灵敏度 = round(100 * sum(ev[pos]) / sum(ev), 1),
               特异度 = round(100 * sum(!pos & ev == 0) / sum(ev == 0), 1),
               阳性预测值 = round(100 * mean(ev[pos]), 1),
               阴性预测值 = round(100 * mean(ev[!pos] == 0), 1),
               检出者中基线已低分占比 = round(100 * mean(sub$low0 == 1, na.rm = TRUE), 1),
               新发真事件数 = yld, 每检出1例新发需筛人数 = if (yld > 0) round(n / yld, 1) else NA_real_)
  }))
  c2 <- coords(roc(ev, s, quiet = TRUE, direction = "<"), "best", best.method = "youden",
               ret = c("threshold", "sensitivity", "specificity"))
  c2 <- c2[which.max(c2$specificity), , drop = FALSE]
  pos <- s >= c2$threshold; sub <- d[pos]
  yld <- sum(sub$event == 1 & sub$low0 == 0, na.rm = TRUE)
  rbind(out, data.table(队列 = co, 目标筛查比例 = "约登最优", 阈值 = c2$threshold, n = length(ev), 事件 = sum(ev), 检出人数 = sum(pos), 实际筛查比例 = round(100 * mean(pos), 2), TP = sum(pos & ev == 1), FP = sum(pos & ev == 0), TN = sum(!pos & ev == 0), FN = sum(!pos & ev == 1),
                        灵敏度 = round(100 * c2$sensitivity, 1), 特异度 = round(100 * c2$specificity, 1),
                        阳性预测值 = round(100 * mean(ev[pos]), 1),
                        阴性预测值 = round(100 * mean(ev[!pos] == 0), 1),
                        检出者中基线已低分占比 = round(100 * mean(sub$low0 == 1, na.rm = TRUE), 1),
                        新发真事件数 = yld, 每检出1例新发需筛人数 = if (yld > 0) round(n / yld, 1) else NA_real_))
}))
say("\n========== 5. 全人群筛查模拟 ==========")
pt(scr); fwrite(scr, file.path(OUT, "F16_screening.csv"), bom = TRUE)

strat_auc <- function(s, co, kind, lab) {
  if (nrow(s) < 30 || length(unique(s$event)) < 2) return(NULL)
  data.table(队列 = co, 分层 = kind, 层 = lab, n = nrow(s),
             事件率 = round(100 * mean(s$event), 1),
             AUC = round(as.numeric(auc(roc(s$event, s$score, quiet = TRUE, direction = "<"))), 4))
}
st <- rbindlist(lapply(CO, function(co) {
  d <- D[[co]][!is.na(mem_fu)]
  d[, ag := fcase(age < 60, "45-59", age < 70, "60-69", age < 80, "70-79", default = "80+")]
  out <- rbindlist(lapply(c("45-59", "60-69", "70-79", "80+"),
                          function(a) strat_auc(d[ag == a], co, "年龄组", a)))
  if (!all(is.na(d$sex))) {
    sx <- d[!is.na(sex)]
    sx[, sexl := paste0("取值", as.character(as.integer(sex)))]
    out <- rbind(out, rbindlist(lapply(unique(sx$sexl),
                                       function(v) strat_auc(sx[sexl == v], co, "性别", v))))
    say("   [", co, "] 性别原值频数：",
        paste(names(table(as.integer(sx$sex))), as.integer(table(as.integer(sx$sex))),
              sep = "=", collapse = "  "))
  }
  eg <- d[!is.na(mem_base)]
  eg[, elab := fcase(edu_pts == 0, "未上学", edu_pts == -2, "小学", edu_pts == -7, "初中",
                     edu_pts == -10, "高中及以上", default = "其他")]
  out <- rbind(out, rbindlist(lapply(unique(eg$elab),
                                     function(v) strat_auc(eg[elab == v], co, "教育", v))))
  out
}))
say("\n========== 6. 分层稳定性 ==========")
pt(st); fwrite(st, file.path(OUT, "F16_strata.csv"), bom = TRUE)

cal <- rbindlist(lapply(CO, function(co) {
  d0 <- D[[co]][!is.na(mem_fu) & !is.na(phat)]
  rbindlist(lapply(c("U", "N"), function(s) {
    d <- if (s == "U") d0 else d0[!is.na(mem_base) & mem_base > cut_int]
    fit <- glm(event ~ lp, data = d, family = binomial())
    data.table(队列 = co, 口径 = LB[[s]], n = nrow(d),
               期望事件率 = round(100 * mean(d$phat), 1), 实际事件率 = round(100 * mean(d$event), 1),
               校准斜率 = round(coef(fit)[2], 3), 校准截距 = round(coef(fit)[1], 3),
               Brier = round(mean((d$phat - d$event)^2), 4))
  }))
}))
say("\n========== 7. 校准 ==========")
pt(cal); fwrite(cal, file.path(OUT, "F16_calibration.csv"), bom = TRUE)

sh <- D[["SHARE"]][!is.na(mem_fu) & !is.na(country_label)]
cn <- sh[, .(n = .N, 事件率 = round(100 * mean(event), 1),
             全人群AUC = round(as.numeric(auc(roc(event, score, quiet = TRUE, direction = "<"))), 3),
             基线已低分占比 = round(100 * mean(low0 == 1, na.rm = TRUE), 1)),
         by = country_label][order(-n)]
cn <- cn[n >= 200]
say("\n========== 8. SHARE 分国家（n>=200，全人群口径）==========")
pt(cn); fwrite(cn, file.path(OUT, "F16_share_country.csv"), bom = TRUE)

close(con)
cat("\n完成 ->", file.path(OUT, "F16_report.txt"), "\n")
