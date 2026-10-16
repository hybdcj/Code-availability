# ============================================================
# 血液统计第二步：2011 与 2015 两波血液数据统一、固定队列对接
# 输出：
#   - 两波血液变量统一后的固定队列分析数据
#   - 两波重叠、缺失、分布和变化摘要
# ============================================================

invisible(try(Sys.setlocale("LC_ALL", "zh_CN.UTF-8"), silent = TRUE))
.libPaths(c("D:/R/library", .libPaths()))
suppressPackageStartupMessages({
  library(haven)
  library(dplyr)
  library(tidyr)
})

in_2011_zip <- "data/raw/wave1_2011/blood/Blood_20140429.zip"
in_2015_zip <- "data/raw/wave3_2015/blood/Blood.zip"
in_cohort   <- "output/04_纵向队列/cohort_change_2026-09-13.rds"
out_dir     <- "output/07_血液指标"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

normalize_id <- function(x) {
  x <- as.character(x)
  ifelse(
    nchar(x) == 11,
    paste0(substr(x, 1, 9), sprintf("%03d", as.integer(substr(x, 10, 11)))),
    x
  )
}

read_first_dta <- function(zip_path) {
  td <- tempfile()
  dir.create(td)
  unzip(zip_path, exdir = td)
  f <- list.files(td, pattern = "\\.dta$", full.names = TRUE)
  stopifnot(length(f) >= 1)
  read_dta(f[1])
}

b11 <- read_first_dta(in_2011_zip) %>%
  transmute(
    KEY = normalize_id(ID),
    bloodweight_2011 = as.numeric(bloodweight),
    plt_2011 = as.numeric(qc1_vb009),
    wbc_2011 = as.numeric(qc1_vb002),
    crp_2011 = as.numeric(newcrp),
    hgb_2011 = as.numeric(qc1_vb004),
    hct_2011 = as.numeric(qc1_vb005),
    mcv_2011 = as.numeric(qc1_vb006)
  ) %>%
  filter(!is.na(KEY)) %>%
  distinct(KEY, .keep_all = TRUE)

b15 <- read_first_dta(in_2015_zip) %>%
  transmute(
    KEY = normalize_id(ID),
    bloodweight_2015 = as.numeric(Blood_weight),
    communityID_2015 = as.character(communityID),
    fasting_2015 = as.numeric(bl_fasting),
    plt_2015 = as.numeric(bl_plt),
    wbc_2015 = as.numeric(bl_wbc),
    crp_2015 = as.numeric(bl_crp),
    hgb_2015 = as.numeric(bl_hgb),
    hct_2015 = as.numeric(bl_hct),
    mcv_2015 = as.numeric(bl_mcv)
  ) %>%
  filter(!is.na(KEY)) %>%
  distinct(KEY, .keep_all = TRUE)

cohort_obj <- readRDS(in_cohort)
cohort_long <- cohort_obj$cohort %>%
  mutate(KEY = normalize_id(ifelse(!is.na(IDn), IDn, ID))) %>%
  group_by(KEY) %>%
  summarise(
    n_waves = n(),
    n_ms10 = sum(!is.na(ms10)),
    y2011 = ms10[year == 2011][1],
    y2013 = ms10[year == 2013][1],
    y2015 = ms10[year == 2015][1],
    y2018 = ms10[year == 2018][1],
    y2020 = ms10[year == 2020][1],
    impaired_2011 = impaired_fixed[year == 2011][1],
    impaired_2020 = impaired_fixed[year == 2020][1],
    .groups = "drop"
  ) %>%
  mutate(
    fixed_cohort = n_waves == 5 & n_ms10 == 5,
    incident_2020 = fixed_cohort & impaired_2011 == 0 & impaired_2020 == 1,
    decline_2plus = fixed_cohort & (y2020 - y2011) <= -2
  )

analysis <- cohort_long %>%
  filter(fixed_cohort) %>%
  left_join(cohort_obj$wide %>% select(KEY, age2011, female, edu_group, delta), by = "KEY") %>%
  left_join(b11, by = "KEY") %>%
  left_join(b15, by = "KEY") %>%
  mutate(
    any_2011 = !is.na(plt_2011) | !is.na(wbc_2011) | !is.na(crp_2011),
    any_2015 = !is.na(plt_2015) | !is.na(wbc_2015) | !is.na(crp_2015),
    core3_2011 = !is.na(plt_2011) & !is.na(wbc_2011) & !is.na(crp_2011),
    core3_2015 = !is.na(plt_2015) & !is.na(wbc_2015) & !is.na(crp_2015),
    core3_both = core3_2011 & core3_2015,
    d_plt = plt_2015 - plt_2011,
    d_wbc = wbc_2015 - wbc_2011,
    d_crp = crp_2015 - crp_2011,
    logcrp_2011 = log1p(pmax(crp_2011, 0)),
    logcrp_2015 = log1p(pmax(crp_2015, 0)),
    d_logcrp = logcrp_2015 - logcrp_2011
  )

# 两波重叠统计
overlap <- tibble(
  fixed_cohort = nrow(analysis),
  any_2011 = sum(analysis$any_2011),
  any_2015 = sum(analysis$any_2015),
  both_any = sum(analysis$any_2011 & analysis$any_2015),
  core3_2011 = sum(analysis$core3_2011),
  core3_2015 = sum(analysis$core3_2015),
  core3_both = sum(analysis$core3_both),
  incident_both_any = sum(analysis$incident_2020 & analysis$any_2011 & analysis$any_2015, na.rm = TRUE),
  decline2_both_any = sum(analysis$decline_2plus & analysis$any_2011 & analysis$any_2015, na.rm = TRUE),
  incident_core3_both = sum(analysis$incident_2020 & analysis$core3_both, na.rm = TRUE),
  decline2_core3_both = sum(analysis$decline_2plus & analysis$core3_both, na.rm = TRUE)
)

# 两波逐变量摘要
summ_one <- function(x, wave, variable) {
  tibble(
    wave = wave, variable = variable,
    n = sum(!is.na(x)), missing = sum(is.na(x)),
    mean = mean(x, na.rm = TRUE), sd = sd(x, na.rm = TRUE),
    p1 = quantile(x, .01, na.rm = TRUE, names = FALSE),
    p25 = quantile(x, .25, na.rm = TRUE, names = FALSE),
    median = median(x, na.rm = TRUE),
    p75 = quantile(x, .75, na.rm = TRUE, names = FALSE),
    p99 = quantile(x, .99, na.rm = TRUE, names = FALSE),
    min = min(x, na.rm = TRUE), max = max(x, na.rm = TRUE)
  )
}
summary_two <- bind_rows(
  summ_one(analysis$plt_2011, 2011, "PLT"),
  summ_one(analysis$wbc_2011, 2011, "WBC"),
  summ_one(analysis$crp_2011, 2011, "CRP"),
  summ_one(analysis$hgb_2011, 2011, "HGB"),
  summ_one(analysis$hct_2011, 2011, "HCT"),
  summ_one(analysis$mcv_2011, 2011, "MCV"),
  summ_one(analysis$plt_2015, 2015, "PLT"),
  summ_one(analysis$wbc_2015, 2015, "WBC"),
  summ_one(analysis$crp_2015, 2015, "CRP"),
  summ_one(analysis$hgb_2015, 2015, "HGB"),
  summ_one(analysis$hct_2015, 2015, "HCT"),
  summ_one(analysis$mcv_2015, 2015, "MCV")
)

summary_change <- bind_rows(
  summ_one(analysis$d_plt, "2015-2011", "ΔPLT"),
  summ_one(analysis$d_wbc, "2015-2011", "ΔWBC"),
  summ_one(analysis$d_crp, "2015-2011", "ΔCRP"),
  summ_one(analysis$d_logcrp, "2015-2011", "Δlog(CRP+1)")
)

write.csv(overlap, file.path(out_dir, "blood_two_wave_overlap_2026-09-19.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
write.csv(summary_two, file.path(out_dir, "blood_two_wave_summary_2026-09-19.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
write.csv(summary_change, file.path(out_dir, "blood_two_wave_change_summary_2026-09-19.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
saveRDS(analysis, file.path(out_dir, "analysis_blood_two_waves_fixedcohort_2026-09-19.rds"))

fmt <- function(x, d = 3) formatC(x, format = "f", digits = d)
lines <- c(
  "# 血液统计第二步 · 两波统一与固定队列对接",
  "",
  "> 生成日期：2026-09-20（年龄口径修正后重算）",
  "",
  "## 一、结论",
  "",
  paste0("- 2011 血液：n=", nrow(b11), "；2015 血液：n=", nrow(b15), "。"),
  paste0("- 五波固定队列：n=", overlap$fixed_cohort, "。"),
  paste0("- 两波至少一项血液指标均有值：**", overlap$both_any, " 人**。"),
  paste0("- 两波 PLT+WBC+CRP 均完整：**", overlap$core3_both, " 人**。"),
  paste0("- 两波血液子集中，十年新发认知受损：", overlap$incident_both_any, " 人；下降≥2分：", overlap$decline2_both_any, " 人。"),
  paste0("- 三项核心指标两波均完整者中，新发受损：", overlap$incident_core3_both, " 人；下降≥2分：", overlap$decline2_core3_both, " 人。"),
  "- 2015 数据同样没有中性粒细胞/淋巴细胞分类，Neu、NLR、PLR 仍无法计算。",
  "",
  "## 二、两波变量摘要",
  "",
  "| 波次 | 指标 | n | 缺失 | 均值 | SD | 中位数 | P1 | P99 | 最小值 | 最大值 |",
  "|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|"
)
for (i in seq_len(nrow(summary_two))) {
  r <- summary_two[i, ]
  lines <- c(lines, paste0("| ", r$wave, " | ", r$variable, " | ", r$n, " | ", r$missing,
                           " | ", fmt(r$mean), " | ", fmt(r$sd), " | ", fmt(r$median),
                           " | ", fmt(r$p1), " | ", fmt(r$p99), " | ", fmt(r$min), " | ", fmt(r$max), " |"))
}
lines <- c(lines,
  "",
  "## 三、2015−2011 变化摘要",
  "",
  "| 指标 | n | 均值变化 | SD | 中位数变化 | P1 | P99 |",
  "|---|---:|---:|---:|---:|---:|---:|"
)
for (i in seq_len(nrow(summary_change))) {
  r <- summary_change[i, ]
  lines <- c(lines, paste0("| ", r$variable, " | ", r$n, " | ", fmt(r$mean), " | ", fmt(r$sd),
                           " | ", fmt(r$median), " | ", fmt(r$p1), " | ", fmt(r$p99), " |"))
}
lines <- c(lines,
  "",
  "## 四、后续分析数据",
  "",
  "可直接读取：",
  "",
  "`output/07_血液指标/analysis_blood_two_waves_fixedcohort_2026-09-19.rds`",
  "",
  "下一阶段可做：",
  "1. PLT、WBC、CRP 的连续值、四分位和异常值处理；",
  "2. 基线 2011 指标预测 2020 新发受损和十年认知变化；",
  "3. 2011→2015 变化量预测后续认知变化；",
  "4. 校正年龄、性别、教育、城乡、基线认知和体重/炎症相关协变量；",
  "5. 完成敏感性分析后再进入主论文 Table 5/Figure 6。"
)
writeLines(lines, file.path(out_dir, "blood_two_wave_report_2026-09-19.md"), useBytes = TRUE)

cat("完成：", file.path(out_dir, "blood_two_wave_report_2026-09-19.md"), "\n")
print(overlap)
print(summary_change)
