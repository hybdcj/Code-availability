# ============================================================================
# 04_cohort_change.R —— 同一批人十年（五波）认知变化（固定队列）
# 版本：v1.0（2026-09-13）
# 输入：output/02_五波数据与质控/waves_cognition_2026-09-13.rds（由 code/03_waves_build.R 生成）
# 输出：output/04_纵向队列/cohort_change_table_2026-09-13.csv、output/04_纵向队列/cohort_change_2026-09-13.md
# 口径要点：
#   ① 队列定义：五波均有"定向 + 连续减 7"（ms10，0–10）的人（这是五波唯一完全可比的认知指标）
#   ② 跨波比较用波内 z（避免 2018 量表差异）；同时报绝对分与"2011 固定切点"下的受损率
#   ③ 转移分析：2011 正常/受损 → 2020 正常/受损
# ============================================================================
invisible(try(Sys.setlocale("LC_ALL", "zh_CN.UTF-8"), silent = TRUE))
if (dir.exists("D:/R/library")) .libPaths(c("D:/R/library", .libPaths()))
suppressPackageStartupMessages({ library(dplyr); library(tidyr) })

P <- Sys.getenv("CHARLS_PROJ", "D:/BaiduSyncdisk/科研/徐中豪")
OUT <- file.path(P, "output", "04_纵向队列")
w <- readRDS(file.path(P, "output", "02_五波数据与质控", "waves_cognition_2026-09-13.rds"))

## ---- 1. 可比指标 ms10 = 定向 + 连续减 7（0–10），五波都有 ----
w <- w %>% mutate(ms10 = orient + ser7)

## ---- 2. 固定队列：2011 年≥45岁 + 五波 ms10 均可用 ----
# 先建立统一 KEY，再做年龄筛选；年龄口径与 2020 期主分析一致。
w$KEY <- ifelse(is.na(w$IDn), w$ID, w$IDn)
candidate_ids <- w %>% filter(!is.na(ms10)) %>% count(KEY) %>%
  filter(n == 5) %>% pull(KEY)
baseline_eligible_ids <- w %>% filter(year == 2011, !is.na(age), age >= 45) %>% pull(KEY) %>% unique()
cohort_ids <- intersect(candidate_ids, baseline_eligible_ids)
coh <- w %>% filter(KEY %in% cohort_ids) %>% arrange(KEY, year)
cat("五波 ms10 均可用的候选者：", length(candidate_ids), " 人；2011≥45岁合格者：",
    length(baseline_eligible_ids), " 人；最终固定队列：", length(cohort_ids), " 人\n")

## ---- 3. 每波：均值（原始 + 波内 z）与受损率 ----
base2011 <- coh %>% filter(year == 2011) %>%
  transmute(KEY, age2011 = age, female, edu_group, ms10_2011 = ms10,
            anchor_thr = ave(ms10, edu_group, FUN = function(x) mean(x, na.rm = TRUE) - sd(x, na.rm = TRUE)))
coh <- coh %>% left_join(base2011 %>% select(KEY, anchor_thr), by = "KEY") %>%
  mutate(impaired_fixed = as.integer(ms10 < anchor_thr))          # 2011 固定切点（教育分层）

by_wave <- coh %>% group_by(year) %>%
  summarise(n = n(),
            原始分均值 = round(mean(ms10, na.rm = TRUE), 2),
            原始分SD = round(sd(ms10, na.rm = TRUE), 2),
            波内z均值 = round(mean(scale(ms10)[, 1], na.rm = TRUE), 3),
            受损率_固定切点 = sprintf("%.1f%%", 100*mean(impaired_fixed, na.rm = TRUE)),
            抑郁率 = sprintf("%.1f%%", 100*mean(cesd10 >= 10, na.rm = TRUE)),
            .groups = "drop")
cat("\n=== 固定队列：各波认知（定向+连续减7，0–10）===\n"); print(as.data.frame(by_wave), row.names = FALSE)

## ---- 4. 变化幅度（2011 → 2020）与分层 ----
ms_wide <- coh %>% select(KEY, year, ms10) %>%
  pivot_wider(names_from = year, values_from = ms10, names_prefix = "y")
base <- coh %>% filter(year == 2011) %>% transmute(KEY, age2011 = age, female, edu_group)
wide <- ms_wide %>% left_join(base, by = "KEY") %>% mutate(delta = y2020 - y2011)
cat("\n=== 十年变化（2020 − 2011，原始分）===\n")
cat(sprintf("均值变化 %.2f 分（SD %.2f）｜下降≥2 分者 %.1f%%｜上升≥2 分者 %.1f%%\n",
            mean(wide$delta, na.rm = TRUE), sd(wide$delta, na.rm = TRUE),
            100*mean(wide$delta <= -2, na.rm = TRUE), 100*mean(wide$delta >= 2, na.rm = TRUE)))

cat("\n按基线年龄组：\n")
print(wide %>% mutate(ageg = cut(age2011, c(0, 54, 64, 74, Inf), labels = c("<55", "55-64", "65-74", "75+"))) %>%
        group_by(ageg) %>% summarise(n = n(), 变化均值 = round(mean(delta, na.rm = TRUE), 2),
                                     decline2 = sprintf("%.1f%%", 100*mean(delta <= -2, na.rm = TRUE)), .groups = "drop") %>%
        as.data.frame(), row.names = FALSE)
cat("\n按教育：\n")
print(wide %>% group_by(edu_group) %>% summarise(n = n(), 变化均值 = round(mean(delta, na.rm = TRUE), 2),
                                                  decline2 = sprintf("%.1f%%", 100*mean(delta <= -2, na.rm = TRUE)), .groups = "drop") %>%
        as.data.frame(), row.names = FALSE)
cat("\n按性别：\n")
print(wide %>% mutate(性别 = ifelse(female == 1, "女", "男")) %>% group_by(性别) %>%
        summarise(n = n(), 变化均值 = round(mean(delta, na.rm = TRUE), 2), .groups = "drop") %>% as.data.frame(), row.names = FALSE)

## ---- 5. 认知状态转移（2011 → 2020，固定切点） ----
tr <- coh %>% filter(year %in% c(2011, 2020)) %>% select(KEY, year, impaired_fixed) %>%
  pivot_wider(names_from = year, values_from = impaired_fixed, names_prefix = "im")
tab <- table(im2011 = tr$im2011, im2020 = tr$im2020, useNA = "no")
cat("\n=== 认知状态转移（行=2011，列=2020；1=受损）===\n"); print(tab)
cat(sprintf("持续正常 %.1f%%｜2011 正常→2020 受损 %.1f%%｜持续受损 %.1f%%｜受损→正常 %.1f%%\n",
            100*tab[1,1]/sum(tab), 100*tab[1,2]/sum(tab), 100*tab[2,2]/sum(tab), 100*tab[2,1]/sum(tab)))

## ---- 6. 基线特征（2011） ----
b <- coh %>% filter(year == 2011)
cat("\n=== 固定队列基线特征（2011）===\n")
cat(sprintf("n=%d｜年龄 %.1f±%.1f 岁｜女性 %.1f%%｜农村 %.1f%%\n", nrow(b),
            mean(b$age, na.rm = TRUE), sd(b$age, na.rm = TRUE),
            100*mean(b$female == 1, na.rm = TRUE), 100*mean(b$rural == 1, na.rm = TRUE)))
cat(sprintf("教育：小学以下 %.1f%%｜小学 %.1f%%｜初中 %.1f%%｜高中及以上 %.1f%%\n",
            100*mean(b$edu_group == "Less than primary", na.rm = TRUE),
            100*mean(b$edu_group == "Primary school", na.rm = TRUE),
            100*mean(b$edu_group == "Middle school", na.rm = TRUE),
            100*mean(b$edu_group == "High school or above", na.rm = TRUE)))
cat(sprintf("ADL 受限 %.1f%%｜抑郁 %.1f%%｜高血压 %.1f%%｜糖尿病 %.1f%%｜睡眠 %.2f h｜疼痛 %.1f%%｜体力活动 %.1f%%\n",
            100*mean(b$adl == 1, na.rm = TRUE), 100*mean(b$cesd10 >= 10, na.rm = TRUE),
            100*mean(b$hibp == 1, na.rm = TRUE), 100*mean(b$diab == 1, na.rm = TRUE),
            mean(b$sleep_h, na.rm = TRUE), 100*mean(b$pain_any == 1, na.rm = TRUE),
            100*mean(b$pa_any == 1, na.rm = TRUE)))

## ---- 7. 输出 ----
dir.create(OUT, showWarnings = FALSE)
write.csv(by_wave, file.path(OUT, "cohort_change_table_2026-09-13.csv"), row.names = FALSE, fileEncoding = "UTF-8")
saveRDS(list(cohort = coh, wide = wide, by_wave = by_wave, transition = tab),
        file.path(OUT, "cohort_change_2026-09-13.rds"))
cat("\n已存：output/04_纵向队列/cohort_change_table_2026-09-13.csv、output/04_纵向队列/cohort_change_2026-09-13.rds\n")
