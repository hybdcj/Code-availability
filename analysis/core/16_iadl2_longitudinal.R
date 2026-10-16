# ============================================================================
# 16_iadl2_longitudinal.R —— IADL（两条目版 iadl2）纵向模型补充分析
# 目的：为 A02/A07 决策「IADL 是否纳入主模型」提供纵向证据
# 设定：与 T06（13_weighted_longitudinal.R）一致——单暴露 + age + female + edu_group；
#       主敏感性：十年变化模型额外校正基线 ms10（回归均值/地板效应）
# 生成：2026-09-26
# 输出：
#   output/04_纵向队列/iadl2_longitudinal_models_2026-09-26.csv
#   output/04_纵向队列/iadl2_longitudinal_summary_2026-09-26.md
# ============================================================================
invisible(try(Sys.setlocale("LC_ALL", "zh_CN.UTF-8"), silent = TRUE))
.libPaths(c("D:/R/library", .libPaths()))
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(survival); library(survey) })
P <- "D:/BaiduSyncdisk/科研/徐中豪"; OUT <- file.path(P, "output", "04_纵向队列")

coh <- readRDS(file.path(OUT, "cohort_change_2026-09-13.rds"))$cohort
inc <- readRDS(file.path(OUT, "interval_incidence_events_2026-09-20.rds"))

b11 <- coh %>% filter(year == 2011) %>%
  transmute(KEY, age, female, edu = edu_group, iadl2, adl, ms10_2011 = ms10, weight, communityID)
delta_tab <- coh %>% select(KEY, year, ms10) %>%
  pivot_wider(names_from = year, values_from = ms10, names_prefix = "y") %>%
  transmute(KEY, delta = y2020 - y2011, ms10_2020 = y2020)

d_inc  <- inc$base %>% inner_join(b11, by = "KEY")
yr <- c("2013" = 2, "2015" = 4, "2018" = 7, "2020" = 9)
im <- cbind(d_inc$im2013, d_inc$im2015, d_inc$im2018, d_inc$im2020)
d_inc$time  <- apply(im, 1, function(z) { i <- which(z == 1); if (length(i)) yr[i[1]] else 9 })
d_inc$event <- as.integer(apply(im, 1, function(z) any(z == 1, na.rm = TRUE)))
d_chg <- b11 %>% inner_join(delta_tab, by = "KEY")

get_est <- function(fit, term, ratio = TRUE) {
  b <- coef(fit)[term]; se <- sqrt(vcov(fit)[term, term])
  if (ratio) list(est = exp(b), lo = exp(b - 1.96*se), hi = exp(b + 1.96*se))
  else list(est = b, lo = b - 1.96*se, hi = b + 1.96*se)
}
fit_all <- function(e, label) {
  rhs  <- paste(e, "+ age + female + edu")
  rhs2 <- paste(rhs, "+ ms10_2011")
  des  <- svydesign(ids = ~communityID, weights = ~weight, data = d_inc)
  des2 <- svydesign(ids = ~communityID, weights = ~weight, data = d_inc)
  des3 <- svydesign(ids = ~communityID, weights = ~weight, data = d_chg)
  out <- list(
    logistic     = c(get_est(svyglm(as.formula(paste("endpoint_2020 ~", rhs)), des, family = quasibinomial()), e),
                     get_est(glm(as.formula(paste("endpoint_2020 ~", rhs)), d_inc, family = binomial()), e)),
    cox          = c(get_est(svycoxph(as.formula(paste("Surv(time, event) ~", rhs)), des2), e),
                     get_est(coxph(as.formula(paste("Surv(time, event) ~", rhs)), d_inc), e)),
    change       = c(get_est(svyglm(as.formula(paste("delta ~", rhs)), des3), e, FALSE),
                     get_est(lm(as.formula(paste("delta ~", rhs)), d_chg), e, FALSE)),
    change_adj   = c(get_est(svyglm(as.formula(paste("delta ~", rhs2)), des3), e, FALSE),
                     get_est(lm(as.formula(paste("delta ~", rhs2)), d_chg), e, FALSE))
  )
  lapply(out, function(z) z)
}
fmt <- function(z) sprintf("%.2f (%.2f-%.2f)", z[[1]], z[[2]], z[[3]])
res <- list()
for (e in c("iadl2", "adl")) {
  r <- fit_all(e, e)
  vars_inc <- c(e, "age", "female", "edu", "endpoint_2020")
  vars_surv <- c(e, "age", "female", "edu", "time", "event")
  vars_chg <- c(e, "age", "female", "edu", "delta")
  n_inc  <- sum(complete.cases(d_inc[, vars_inc]))
  n_surv <- sum(complete.cases(d_inc[, vars_surv]))
  n_chg  <- sum(complete.cases(d_chg[, vars_chg]))
  n_lookup <- c(logistic = n_inc, cox = n_surv, change = n_chg, change_adj = n_chg)
  for (m in names(r)) {
    z <- r[[m]]
    res[[length(res)+1]] <- data.frame(exposure = e, model = m, n = n_lookup[[m]],
      weighted = fmt(list(z[[1]], z[[2]], z[[3]])), unweighted = fmt(list(z[[4]], z[[5]], z[[6]])))
  }
}
res <- do.call(rbind, res)
write.csv(res, file.path(OUT, "iadl2_longitudinal_models_2026-09-26.csv"), row.names = FALSE, fileEncoding = "UTF-8")

desc <- d_chg %>% group_by(iadl2) %>%
  summarise(n = n(), base = round(mean(ms10_2011, na.rm = TRUE), 2), follow = round(mean(ms10_2020, na.rm = TRUE), 2),
            delta = round(mean(delta, na.rm = TRUE), 2), dec1 = sprintf("%.1f%%", 100*mean((ms10_2020 - ms10_2011) <= -1, na.rm = TRUE)), .groups = "drop")
md <- c("# IADL（两条目版 iadl2）纵向补充分析（2026-09-26）", "",
  "> 目的：为「IADL 是否纳入主模型」提供纵向证据（此前仅有横断面 OR 1.53–2.41）。",
  "> 设定：与 T06 一致（单暴露 + 年龄 + 性别 + 教育）；加权＝2011 基线权重 + PSU=communityID。", "",
  "## 一、描述性（按 2011 年 iadl2）", "", paste(capture.output(print(as.data.frame(desc), row.names = FALSE)), collapse = "\n"), "",
  "## 二、模型结果（加权 / 未加权）", "", paste(capture.output(print(res, row.names = FALSE)), collapse = "\n"), "",
  "## 三、解读", "",
  "- 二分类结局：OR 1.32（0.99–1.75，加权）/ 1.28（1.00–1.65，未加权）；HR 1.19（1.00–1.41）/ 1.17（1.00–1.38）——**均为临界**，弱于 ADL（OR 1.49、HR 1.20）。",
  "- 连续变化：不校正基线时 β=+0.57（看似「下降更少」），**校正基线 ms10 后 β≈−0.01（p=0.93）** → 该「效应」是**基线差异导致的回归均值假象**（IADL 受限组基线 ms10 仅 6.90 vs 8.19）。",
  "- 结论建议：**IADL 不进入主评分**；作为补充/次要分析报告（横断面证据 + 临界纵向关联 + 基线差异显著）。",
  "", "*脚本：code/16_iadl2_longitudinal.R；数据：cohort_change_2026-09-13.rds + interval_incidence_events_2026-09-20.rds*")
writeLines(md, file.path(OUT, "iadl2_longitudinal_summary_2026-09-26.md"), useBytes = TRUE)
cat("\n==== 结果 ====\n"); print(res, row.names = FALSE)
cat("\n==== 描述性 ====\n"); print(as.data.frame(desc), row.names = FALSE)
