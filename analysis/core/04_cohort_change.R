# ============================================================================
# 04_cohort_change.R —— same participants across ten years(five waves)cognitive change(fixed cohort)
# version: v1.0(2026-09-13)
# input: output/02_five_wave_data_and_quality_control/waves_cognition_2026-09-13.rds(from code/03_waves_build.R generate)
# output: output/04_longitudinal_cohort/cohort_change_table_2026-09-13.csv, output/04_longitudinal_cohort/cohort_change_2026-09-13.md
# key definitions: 
#   ① cohort definition: available at all five waves"orientation + serial subtraction of 7"(ms10, 0–10)participants(only fully comparable cognitive measure across five waves)
#   ② cross wave comparison uses within wave z(avoid 2018 scale differences); report absolute scores alongside"2011 fixed cutoff"impairment rate under
#   ③ transition analysis: 2011 normal/impaired → 2020 normal/impaired
# ============================================================================
invisible(try(Sys.setlocale("LC_ALL", "zh_CN.UTF-8"), silent = TRUE))
if (dir.exists("D:/R/library")) .libPaths(c("D:/R/library", .libPaths()))
suppressPackageStartupMessages({ library(dplyr); library(tidyr) })

P <- Sys.getenv("CHARLS_PROJ", "D:/BaiduSyncdisk/research/Xu_Zhonghao")
OUT <- file.path(P, "output", "04_longitudinal_cohort")
w <- readRDS(file.path(P, "output", "02_five_wave_data_and_quality_control", "waves_cognition_2026-09-13.rds"))

## ---- 1. comparable measure ms10 = orientation + serial subtraction of 7(0–10), present at all five waves ----
w <- w %>% mutate(ms10 = orient + ser7)

## ---- 2. fixed cohort: 2011 year≥45 years old + five waves ms10 all available ----
# first construct common KEY, then apply age eligibility; age criteria and 2020 consistent with wave primary analysis.
w$KEY <- ifelse(is.na(w$IDn), w$ID, w$IDn)
candidate_ids <- w %>% filter(!is.na(ms10)) %>% count(KEY) %>%
  filter(n == 5) %>% pull(KEY)
baseline_eligible_ids <- w %>% filter(year == 2011, !is.na(age), age >= 45) %>% pull(KEY) %>% unique()
cohort_ids <- intersect(candidate_ids, baseline_eligible_ids)
coh <- w %>% filter(KEY %in% cohort_ids) %>% arrange(KEY, year)
cat("five waves ms10 candidates with all values available: ", length(candidate_ids), " participants; 2011≥45 age eligible participants: ",
    length(baseline_eligible_ids), " participants; final fixed cohort: ", length(cohort_ids), " participants\n")

## ---- 3. each wave: mean(original + within wave z)and impairment rate ----
base2011 <- coh %>% filter(year == 2011) %>%
  transmute(KEY, age2011 = age, female, edu_group, ms10_2011 = ms10,
            anchor_thr = ave(ms10, edu_group, FUN = function(x) mean(x, na.rm = TRUE) - sd(x, na.rm = TRUE)))
coh <- coh %>% left_join(base2011 %>% select(KEY, anchor_thr), by = "KEY") %>%
  mutate(impaired_fixed = as.integer(ms10 < anchor_thr))          # 2011 fixed cutoff(education strata)

by_wave <- coh %>% group_by(year) %>%
  summarise(n = n(),
            mean_raw_score = round(mean(ms10, na.rm = TRUE), 2),
            raw_score_SD = round(sd(ms10, na.rm = TRUE), 2),
            within_wave_z_mean = round(mean(scale(ms10)[, 1], na.rm = TRUE), 3),
            impairment_rate_fixed_cutoff = sprintf("%.1f%%", 100*mean(impaired_fixed, na.rm = TRUE)),
            depression_rate = sprintf("%.1f%%", 100*mean(cesd10 >= 10, na.rm = TRUE)),
            .groups = "drop")
cat("\n=== fixed cohort: cognition at each wave(orientation+serial subtraction of 7, 0–10)===\n"); print(as.data.frame(by_wave), row.names = FALSE)

## ---- 4. magnitude of change(2011 → 2020)and strata ----
ms_wide <- coh %>% select(KEY, year, ms10) %>%
  pivot_wider(names_from = year, values_from = ms10, names_prefix = "y")
base <- coh %>% filter(year == 2011) %>% transmute(KEY, age2011 = age, female, edu_group)
wide <- ms_wide %>% left_join(base, by = "KEY") %>% mutate(delta = y2020 - y2011)
cat("\n=== ten year change(2020 − 2011, raw score)===\n")
cat(sprintf("mean change %.2f points(SD %.2f)｜decline≥2 scoring participants %.1f%%｜increase≥2 scoring participants %.1f%%\n",
            mean(wide$delta, na.rm = TRUE), sd(wide$delta, na.rm = TRUE),
            100*mean(wide$delta <= -2, na.rm = TRUE), 100*mean(wide$delta >= 2, na.rm = TRUE)))

cat("\n by baseline age group: \n")
print(wide %>% mutate(ageg = cut(age2011, c(0, 54, 64, 74, Inf), labels = c("<55", "55-64", "65-74", "75+"))) %>%
        group_by(ageg) %>% summarise(n = n(), mean_change = round(mean(delta, na.rm = TRUE), 2),
                                     decline2 = sprintf("%.1f%%", 100*mean(delta <= -2, na.rm = TRUE)), .groups = "drop") %>%
        as.data.frame(), row.names = FALSE)
cat("\n by education: \n")
print(wide %>% group_by(edu_group) %>% summarise(n = n(), mean_change = round(mean(delta, na.rm = TRUE), 2),
                                                  decline2 = sprintf("%.1f%%", 100*mean(delta <= -2, na.rm = TRUE)), .groups = "drop") %>%
        as.data.frame(), row.names = FALSE)
cat("\n by sex: \n")
print(wide %>% mutate(sex = ifelse(female == 1, "female", "male")) %>% group_by(sex) %>%
        summarise(n = n(), mean_change = round(mean(delta, na.rm = TRUE), 2), .groups = "drop") %>% as.data.frame(), row.names = FALSE)

## ---- 5. cognitive status transitions(2011 → 2020, fixed cutoff) ----
tr <- coh %>% filter(year %in% c(2011, 2020)) %>% select(KEY, year, impaired_fixed) %>%
  pivot_wider(names_from = year, values_from = impaired_fixed, names_prefix = "im")
tab <- table(im2011 = tr$im2011, im2020 = tr$im2020, useNA = "no")
cat("\n=== cognitive status transitions(rows=2011, columns=2020; 1=impaired)===\n"); print(tab)
cat(sprintf("persistently unimpaired %.1f%%｜2011 normal→2020 impaired %.1f%%｜persistently impaired %.1f%%｜impaired→normal %.1f%%\n",
            100*tab[1,1]/sum(tab), 100*tab[1,2]/sum(tab), 100*tab[2,2]/sum(tab), 100*tab[2,1]/sum(tab)))

## ---- 6. baseline characteristics(2011) ----
b <- coh %>% filter(year == 2011)
cat("\n=== fixed cohort baseline characteristics(2011)===\n")
cat(sprintf("n=%d｜age %.1f±%.1f years old｜female %.1f%%｜rural %.1f%%\n", nrow(b),
            mean(b$age, na.rm = TRUE), sd(b$age, na.rm = TRUE),
            100*mean(b$female == 1, na.rm = TRUE), 100*mean(b$rural == 1, na.rm = TRUE)))
cat(sprintf("education: below primary school %.1f%%｜primary school %.1f%%｜middle school %.1f%%｜high school or above %.1f%%\n",
            100*mean(b$edu_group == "Less than primary", na.rm = TRUE),
            100*mean(b$edu_group == "Primary school", na.rm = TRUE),
            100*mean(b$edu_group == "Middle school", na.rm = TRUE),
            100*mean(b$edu_group == "High school or above", na.rm = TRUE)))
cat(sprintf("ADL limited %.1f%%｜depression %.1f%%｜hypertension %.1f%%｜diabetes %.1f%%｜sleep %.2f h｜pain %.1f%%｜physical activity %.1f%%\n",
            100*mean(b$adl == 1, na.rm = TRUE), 100*mean(b$cesd10 >= 10, na.rm = TRUE),
            100*mean(b$hibp == 1, na.rm = TRUE), 100*mean(b$diab == 1, na.rm = TRUE),
            mean(b$sleep_h, na.rm = TRUE), 100*mean(b$pain_any == 1, na.rm = TRUE),
            100*mean(b$pa_any == 1, na.rm = TRUE)))

## ---- 7. output ----
dir.create(OUT, showWarnings = FALSE)
write.csv(by_wave, file.path(OUT, "cohort_change_table_2026-09-13.csv"), row.names = FALSE, fileEncoding = "UTF-8")
saveRDS(list(cohort = coh, wide = wide, by_wave = by_wave, transition = tab),
        file.path(OUT, "cohort_change_2026-09-13.rds"))
cat("\n saved: output/04_longitudinal_cohort/cohort_change_table_2026-09-13.csv, output/04_longitudinal_cohort/cohort_change_2026-09-13.rds\n")
