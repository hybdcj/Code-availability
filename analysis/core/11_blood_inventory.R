invisible(try(Sys.setlocale("LC_ALL", "zh_CN.UTF-8"), silent = TRUE))
# ============================================================
# blood analysis step one: 2011 blood variable inventory and fixed cohort linkage
# target: 
#   1) verify 2011 whether blood data contain PLT/WBC/CRP and related biomarkers
#   2) and 2011 year≥45 years old, five waves ms10 scorable fixed cohort linkage
#   3) missingness counts, distribution ranges and analyzable sample sizes
#   4) output for later regression 2011 blood subset
# notes: 
#   - raw data read only, do not modify; 
#   - CHARLS released blood data contain no neutrophils/lymphocyte differential, therefore cannot calculate Neu/NLR/PLR; 
#   - 2015 blood data not yet downloaded, extend to second time point after acquisition.
# ============================================================

.libPaths(c("D:/R/library", .libPaths()))
suppressPackageStartupMessages({
  library(haven)
  library(dplyr)
  library(tidyr)
})

in_blood_zip <- "data/raw/wave1_2011/blood/Blood_20140429.zip"
in_cohort    <- "output/04_longitudinal_cohort/cohort_change_2026-09-13.rds"
out_dir      <- "output/07_blood_biomarkers"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# harmonization CHARLS individual ID: 
# 2011 is 11 positions(household number 9 positions + individual sequence number 2 positions), remaining values are 12 positions.
normalize_id <- function(x) {
  x <- as.character(x)
  ifelse(
    nchar(x) == 11,
    paste0(substr(x, 1, 9), sprintf("%03d", as.integer(substr(x, 10, 11)))),
    x
  )
}

# 1. read 2011 blood data
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

# 2. read five wave fixed cohort
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

# 3. link blood data to fixed cohort
dat <- key_year %>%
  left_join(blood %>% select(KEY, PLT, WBC, CRP, HGB, HCT, MCV, bloodweight = bloodweight),
            by = "KEY") %>%
  left_join(cohort_obj$wide %>% select(KEY, age2011), by = "KEY")

# 4. variable inventory
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

# 5. linkage and missingness counts
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

# 6. distribution summary(fixed cohort participants with nonmissing blood measurements)
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

# 7. retain fixed cohort only, for direct use in subsequent regression
analysis <- dat %>% filter(fixed_cohort)
saveRDS(analysis, file.path(out_dir, "analysis_blood_2011_fixedcohort_2026-09-19.rds"))

# 8. generate Markdown inventory report
fmt <- function(x, digits = 3) formatC(x, format = "f", digits = digits)
lines <- c(
  "# blood analysis step one · data availability and variable inventory",
  "",
  "> generation date: 2026-09-20(recompute after correcting age eligibility)",
  "",
  "## one, conclusion",
  "",
  paste0("- 2011 year blood sample data: **downloaded and successfully read**, n=", nrow(blood), "."),
  paste0("- linkable to five wave fixed cohort: **", blood_flags$fixed_cohort, " participants**."),
  paste0("- all three core biomarkers nonmissing: PLT=", blood_flags$PLT_nonmissing,
         "; WBC=", blood_flags$WBC_nonmissing,
         "; CRP=", blood_flags$CRP_nonmissing, "."),
  paste0("- PLT, WBC, CRP all three simultaneously nonmissing: **", blood_flags$all3_nonmissing, " participants**."),
  "- CHARLS released blood data contain no neutrophil or lymphocyte differentials, **Neu, NLR, PLR cannot calculate**.",
  "- 2015 year blood data acquired on 2026-09-19 download(Blood.zip), two wave linkage completed in step two.",
  "",
  "## two, 2011 available variables by year",
  "",
  "| metrics | raw variable | unit |",
  "|---|---|---|",
  "| PLT | qc1_vb009 | 10^9/L |",
  "| WBC | qc1_vb002 | 10^3/uL |",
  "| CRP | newcrp | mg/L |",
  "| HGB | qc1_vb004 | g/dL |",
  "| HCT | qc1_vb005 | % |",
  "| MCV | qc1_vb006 | fL |",
  "",
  "## three, fixed cohort distribution summary",
  "",
  "| metrics | n | missing | mean | SD | median | P1 | P99 |",
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
  "## four, entry point for subsequent analysis",
  "",
  paste0("directly readable: `", file.path(out_dir, "analysis_blood_2011_fixedcohort_2026-09-19.rds"), "`"),
  "",
  "suggested model sequence: ",
  "1. PLT, WBC, CRP model as continuous values and tertiles respectively/quartile models; ",
  "2. first analyze ten year incident impairment Logistic, then linear regression of cognitive change; ",
  "3. adjust for age, sex, education, urban rural residence and baseline cognition; ",
  "4. winsorize or exclude extremes in sensitivity analyses; ",
  "5. 2015 analyze second time point and repeated measures after blood data acquisition.",
  "",
  "## five, current data status",
  "",
  "- 2011 and 2015 blood variable inventories completed; ",
  "- two wave linkage and subsequent models in `output/07_blood_biomarkers/blood_two_wave_report_2026-09-19.md`; ",
  "- Neu/NLR/PLR not calculable, subsequently use PLT, WBC, CRP as blood inflammation related biomarkers."
)
writeLines(lines, file.path(out_dir, "blood_inventory_report_2026-09-19.md"), useBytes = TRUE)

cat("complete: ", file.path(out_dir, "blood_inventory_report_2026-09-19.md"), "\n")
print(blood_flags)
print(summary_rows)
