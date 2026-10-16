# ============================================================================
# 05_cohort_risk_factors.R —— 固定队列：基线危险因素 → 十年认知下降
# 版本：v1.0（2026-09-13）
# 输入：output/04_纵向队列/cohort_change_2026-09-13.rds（由 code/04_cohort_change.R 生成）
# 输出：output/04_纵向队列/cohort_risk_2026-09-13.csv（四种模型结果汇总）
# 四个分析（均以 2011 基线暴露为自变量，校正年龄、性别、教育）：
#   ① 新发受损（Logistic，OR）：基线正常者 2011→2020 变为受损
#   ② 首次受损时间（Cox，HR）：2013/2015/2018/2020 四波随访
#   ③ 认知下降速度（线性混合模型，nlme）：时间×暴露 交互项
#   ④ 十年变化量（线性回归，β）：2020 − 2011 分差
# 说明：2011 期没有 IADL，故基线暴露不含 IADL；睡眠按小时，疼痛/体力活动为二分类。
# ============================================================================
invisible(try(Sys.setlocale("LC_ALL", "zh_CN.UTF-8"), silent = TRUE))
if (dir.exists("D:/R/library")) .libPaths(c("D:/R/library", .libPaths()))
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(survival); library(nlme) })

P <- Sys.getenv("CHARLS_PROJ", "D:/BaiduSyncdisk/科研/徐中豪")
OUT <- file.path(P, "output", "04_纵向队列")
obj <- readRDS(file.path(OUT, "cohort_change_2026-09-13.rds"))
coh <- obj$cohort

## ---- 基线暴露（2011） ----
base <- coh %>% filter(year == 2011) %>%
  transmute(KEY, age, female, edu = edu_group,
            dep = as.integer(cesd10 >= 10), adl, sleep_h,
            pain = pain_any, pa = pa_any, hibp, diab,
            ms10_2011 = ms10, impaired_2011 = impaired_fixed)

## ---- 随访结局 ----
# ① 新发受损（2011 正常 → 2020 受损）
inc <- base %>% filter(impaired_2011 == 0) %>%
  left_join(coh %>% filter(year == 2020) %>% select(KEY, impaired_2020 = impaired_fixed), by = "KEY")
# ② 首次受损时间（年）：随访时点 2013=2, 2015=4, 2018=7, 2020=9
years <- c("2013" = 2, "2015" = 4, "2018" = 7, "2020" = 9)
fu <- coh %>% filter(year != 2011) %>% select(KEY, year, impaired_fixed) %>%
  mutate(t = years[as.character(year)]) %>% arrange(KEY, t)
surv_df <- fu %>% group_by(KEY) %>% summarise(
  time = if (any(impaired_fixed == 1, na.rm = TRUE)) min(t[impaired_fixed == 1], na.rm = TRUE) else 9,
  event = as.integer(any(impaired_fixed == 1, na.rm = TRUE)), .groups = "drop")
surv_df <- base %>% filter(impaired_2011 == 0) %>% left_join(surv_df, by = "KEY")
# ④ 十年变化量
delta <- coh %>% select(KEY, year, ms10) %>%
  pivot_wider(names_from = year, values_from = ms10, names_prefix = "y") %>%
  transmute(KEY, delta = y2020 - y2011)
base <- base %>% left_join(delta, by = "KEY")

exposures <- c("dep", "adl", "sleep_h", "pain", "pa", "hibp", "diab")
labels <- c(dep = "抑郁（CES-D-10 ≥10）", adl = "ADL 受限", sleep_h = "睡眠时长（每 1 小时）",
            pain = "疼痛", pa = "规律体力活动", hibp = "高血压", diab = "糖尿病")
adj <- " + age + female + edu"
fmt <- function(est, lo, hi, p) sprintf("%.2f (%.2f-%.2f), p=%s", est, lo, hi,
                                        ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))

res <- do.call(rbind, lapply(exposures, function(e) {
  # ① Logistic
  m1 <- glm(as.formula(paste("impaired_2020 ~", e, adj)), data = inc, family = binomial)
  s1 <- summary(m1)$coefficients[2, ]
  or <- exp(s1[1]); ci1 <- exp(s1[1] - 1.96*s1[2]); ci2 <- exp(s1[1] + 1.96*s1[2])
  # ② Cox
  m2 <- coxph(as.formula(paste("Surv(time, event) ~", e, adj)), data = surv_df)
  s2 <- summary(m2)$coefficients[e, ]; s2c <- summary(m2)$conf.int[e, ]
  # ③ 混合模型：时间 × 暴露
  d <- coh %>% select(KEY, year, ms10) %>% mutate(t = c("2011"=0,"2013"=2,"2015"=4,"2018"=7,"2020"=9)[as.character(year)]) %>%
    left_join(base %>% select(KEY, all_of(exposures), age, female, edu), by = "KEY") %>% filter(!is.na(ms10))
  m3 <- try(lme(as.formula(paste("ms10 ~ t *", e, adj)), random = ~ 1 + t | KEY, data = d, na.action = na.omit), silent = TRUE)
  if (inherits(m3, "try-error")) { beta <- NA; blo <- NA; bhi <- NA; bp <- NA } else {
    s3 <- summary(m3)$tTable; rn <- grep(paste0("^t:", e), rownames(s3))
    if (length(rn)) { beta <- s3[rn[1],1]; blo <- beta - 1.96*s3[rn[1],2]; bhi <- beta + 1.96*s3[rn[1],2]; bp <- s3[rn[1],5] }
    else { beta <- NA; blo <- NA; bhi <- NA; bp <- NA } }
  # ④ 十年变化量
  m4 <- lm(as.formula(paste("delta ~", e, adj)), data = base)
  s4 <- summary(m4)$coefficients[2, ]
  data.frame(暴露 = labels[[e]],
             新发受损_OR = fmt(or, ci1, ci2, s1[4]),
             首次受损_HR = fmt(s2c[1], s2c[3], s2c[4], s2[5]),
             下降速度_beta = ifelse(is.na(beta), "—", sprintf("%.3f (%.3f-%.3f), p=%s", beta, blo, bhi,
                                                        ifelse(bp < 0.001, "<0.001", sprintf("%.3f", bp)))),
             十年变化_beta = sprintf("%.3f (%.3f-%.3f), p=%s", s4[1], s4[1]-1.96*s4[2], s4[1]+1.96*s4[2],
                                     ifelse(s4[4] < 0.001, "<0.001", sprintf("%.3f", s4[4]))),
             row.names = NULL)
}))

cat("=== 固定队列危险因素分析（n=", nrow(base), "；新发分析 n=", nrow(inc), "）===\n", sep = "")
print(res, row.names = FALSE)
write.csv(res, file.path(OUT, "cohort_risk_2026-09-13.csv"), row.names = FALSE, fileEncoding = "UTF-8")
cat("\n事件数：随访期间发生认知受损 ", sum(surv_df$event, na.rm = TRUE), " 例（",
    sprintf("%.1f%%", 100*mean(surv_df$event, na.rm = TRUE)), "）\n", sep = "")
cat("已存：output/04_纵向队列/cohort_risk_2026-09-13.csv\n")
