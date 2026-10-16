invisible(try(Sys.setlocale("LC_ALL", "zh_CN.UTF-8"), silent = TRUE))
# ============================================================
# 血液统计第一步：2011 血样变量清查与固定队列对接
# 目标：
#   1) 核对 2011 血液数据是否含 PLT/WBC/CRP 等指标
#   2) 与2011年≥45岁、五波ms10可计分固定队列 对接
#   3) 统计缺失、分布范围和可分析样本量
#   4) 输出可供后续回归使用的 2011 血液子集
# 说明：
#   - 原始数据只读，不修改；
#   - CHARLS 已释放血液数据不含中性粒细胞/淋巴细胞分类，故不能计算 Neu/NLR/PLR；
#   - 2015 血液数据尚未下载，待取得后再扩展第二时点。
# ============================================================

.libPaths(c("D:/R/library", .libPaths()))
suppressPackageStartupMessages({
  library(haven)
  library(dplyr)
  library(tidyr)
})

in_blood_zip <- "data/raw/wave1_2011/blood/Blood_20140429.zip"
in_cohort    <- "output/04_纵向队列/cohort_change_2026-09-13.rds"
out_dir      <- "output/07_血液指标"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# 统一 CHARLS 个体 ID：
# 2011 为 11 位（户号9位 + 个人序号2位），其余为 12 位。
normalize_id <- function(x) {
  x <- as.character(x)
  ifelse(
    nchar(x) == 11,
    paste0(substr(x, 1, 9), sprintf("%03d", as.integer(substr(x, 10, 11)))),
    x
  )
}

# 1. 读取 2011 血液数据
tmp <- tempfile()
dir.create(tmp)
unzip(in_blood_zip, exdir = tmp)
blood_file <- list.files(tmp, pattern = "\\.dta$", full.names = TRUE)[1]
stopifnot(length(blood_file) == 1)
blood <- read_dta(blood_file)

blood <- blood %>%
  mutate(
    KEY = normalize_id(ID),
    PLT = as.numeric(qc1_vb009),       # Platelets, 10^9/L
    WBC = as.numeric(qc1_vb002),       # White blood cell, 10^3/uL
    CRP = as.numeric(newcrp),          # C-reactive protein, mg/L
    HGB = as.numeric(qc1_vb004),       # Hemoglobin, g/dL
    HCT = as.numeric(qc1_vb005),       # Hematocrit, %
    MCV = as.numeric(qc1_vb006)        # Mean corpuscular volume, fL
  )

# 2. 读取五波固定队列
cohort_obj <- readRDS(in_cohort)
cohort_long <- cohort_obj$cohort %>%
  mutate(KEY = normalize_id(ifelse(!is.na(IDn), IDn, ID)))

cohort_wide <- cohort_obj$wide
key_year <- cohort_long %>%
  group_by(KEY) %>%
  summarise(
    n_waves = n(),
    n_ms10 = sum(!is.na(ms10)),
    y2011 = ms10[year == 2011][1],
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

# 3. 血液与固定队列对接
dat <- key_year %>%
  left_join(blood %>% select(KEY, PLT, WBC, CRP, HGB, HCT, MCV, bloodweight = bloodweight),
            by = "KEY") %>%
  left_join(cohort_obj$wide %>% select(KEY, age2011), by = "KEY")

# 4. 变量清单
var_inventory <- tibble(
  variable = c("PLT", "WBC", "CRP", "HGB", "HCT", "MCV"),
  raw_variable = c("qc1_vb009", "qc1_vb002", "newcrp", "qc1_vb004", "qc1_vb005", "qc1_vb006"),
  label = c("Platelets", "White blood cell", "C-reactive protein",
            "Hemoglobin", "Hematocrit", "Mean corpuscular volume"),
  unit = c("10^9/L", "10^3/uL", "mg/L", "g/dL", "%", "fL"),
  available_2011 = TRUE,
  available_2015_expected = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE),
  available_neutrophil = FALSE
)
write.csv(var_inventory, file.path(out_dir, "blood_variable_inventory_2026-09-19.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

# 5. 对接与缺失统计
blood_flags <- dat %>%
  summarise(
    fixed_cohort = sum(fixed_cohort, na.rm = TRUE),
    blood_any = sum(!is.na(PLT) | !is.na(WBC) | !is.na(CRP)),
    PLT_nonmissing = sum(!is.na(PLT)),
    WBC_nonmissing = sum(!is.na(WBC)),
    CRP_nonmissing = sum(!is.na(CRP)),
    all3_nonmissing = sum(!is.na(PLT) & !is.na(WBC) & !is.na(CRP)),
    incident_available = sum(!is.na(incident_2020)),
    incident_among_blood = sum(incident_2020 & (!is.na(PLT) | !is.na(WBC) | !is.na(CRP)), na.rm = TRUE),
    decline2_among_blood = sum(decline_2plus & (!is.na(PLT) | !is.na(WBC) | !is.na(CRP)), na.rm = TRUE)
  )
write.csv(blood_flags, file.path(out_dir, "blood_overlap_2026-09-19.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

# 6. 分布摘要（固定队列中血液非缺失者）
summary_rows <- lapply(c("PLT", "WBC", "CRP", "HGB", "HCT", "MCV"), function(v) {
  x <- dat[[v]][dat$fixed_cohort]
  tibble(
    variable = v,
    n = sum(!is.na(x)),
    missing = sum(is.na(x)),
    mean = mean(x, na.rm = TRUE),
    sd = sd(x, na.rm = TRUE),
    p1 = quantile(x, 0.01, na.rm = TRUE, names = FALSE),
    p25 = quantile(x, 0.25, na.rm = TRUE, names = FALSE),
    median = median(x, na.rm = TRUE),
    p75 = quantile(x, 0.75, na.rm = TRUE, names = FALSE),
    p99 = quantile(x, 0.99, na.rm = TRUE, names = FALSE),
    min = min(x, na.rm = TRUE),
    max = max(x, na.rm = TRUE)
  )
}) %>% bind_rows()
write.csv(summary_rows, file.path(out_dir, "blood_summary_fixedcohort_2026-09-19.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

# 7. 仅保留固定队列，供后续回归直接使用
analysis <- dat %>% filter(fixed_cohort)
saveRDS(analysis, file.path(out_dir, "analysis_blood_2011_fixedcohort_2026-09-19.rds"))

# 8. 生成 Markdown 清查报告
fmt <- function(x, digits = 3) formatC(x, format = "f", digits = digits)
lines <- c(
  "# 血液统计第一步 · 数据可用性与变量清查",
  "",
  "> 生成日期：2026-09-20（年龄口径修正后重算）",
  "",
  "## 一、结论",
  "",
  paste0("- 2011 年血样数据：**已下载并成功读取**，n=", nrow(blood), "。"),
  paste0("- 与五波固定队列可对接：**", blood_flags$fixed_cohort, " 人**。"),
  paste0("- 三项核心指标均非缺失：PLT=", blood_flags$PLT_nonmissing,
         "；WBC=", blood_flags$WBC_nonmissing,
         "；CRP=", blood_flags$CRP_nonmissing, "。"),
  paste0("- PLT、WBC、CRP 三项同时非缺失：**", blood_flags$all3_nonmissing, " 人**。"),
  "- CHARLS 已释放血液数据中没有中性粒细胞或淋巴细胞分类，**Neu、NLR、PLR 无法计算**。",
  "- 2015 年血液数据已于 2026-09-19 下载（Blood.zip），并已在第二步完成两波对接。",
  "",
  "## 二、2011 年可用变量",
  "",
  "| 指标 | 原始变量 | 单位 |",
  "|---|---|---|",
  "| PLT | qc1_vb009 | 10^9/L |",
  "| WBC | qc1_vb002 | 10^3/uL |",
  "| CRP | newcrp | mg/L |",
  "| HGB | qc1_vb004 | g/dL |",
  "| HCT | qc1_vb005 | % |",
  "| MCV | qc1_vb006 | fL |",
  "",
  "## 三、固定队列分布摘要",
  "",
  "| 指标 | n | 缺失 | 均值 | SD | 中位数 | P1 | P99 |",
  "|---|---:|---:|---:|---:|---:|---:|---:|"
)
for (i in seq_len(nrow(summary_rows))) {
  r <- summary_rows[i, ]
  lines <- c(lines, paste0("| ", r$variable, " | ", r$n, " | ", r$missing, " | ",
                           fmt(r$mean), " | ", fmt(r$sd), " | ", fmt(r$median), " | ",
                           fmt(r$p1), " | ", fmt(r$p99), " |"))
}
lines <- c(lines,
  "",
  "## 四、后续分析入口",
  "",
  paste0("可直接读取：`", file.path(out_dir, "analysis_blood_2011_fixedcohort_2026-09-19.rds"), "`"),
  "",
  "建议模型顺序：",
  "1. PLT、WBC、CRP 分别按连续值和三分位/四分位建模；",
  "2. 结局先做十年新发受损 Logistic，再做认知变化量线性回归；",
  "3. 校正年龄、性别、教育、城乡及基线认知；",
  "4. 对极端值做截尾或排除敏感性分析；",
  "5. 2015 血样下载后再做第二时点及重复测量分析。",
  "",
  "## 五、当前数据状态",
  "",
  "- 2011 与 2015 血液数据均已完成变量清查；",
  "- 两波对接及后续模型见 `output/07_血液指标/blood_two_wave_report_2026-09-19.md`；",
  "- Neu/NLR/PLR 不可计算，后续以 PLT、WBC、CRP 作为血液炎症相关指标。"
)
writeLines(lines, file.path(out_dir, "blood_inventory_report_2026-09-19.md"), useBytes = TRUE)

cat("完成：", file.path(out_dir, "blood_inventory_report_2026-09-19.md"), "\n")
print(blood_flags)
print(summary_rows)
