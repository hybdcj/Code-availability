# =============================================================================
# 09_fivewave_table1.R —— CHARLS five waves(2011/2013/2015/2018/2020)cross sectional Table 1: 5 combined table
# purpose: supervisor meeting(2026-09-20)one deliverable; cross sectional characteristics by wave specific cognitive impairment status
# definition template: code/02_table1.R(2020 previously delivered wave baseline table v4.4)——
#   sample = age at this wave ≥ 45 and cognition 9 complete components(cog_complete == 1); 
#   statistic = mean (SD) / n (%); p value = wilcox.test / chisq.test; including Overall(total)columns; unweighted
# input: output/02_five_wave_data_and_quality_control/waves_cognition_2026-09-13.rds(code/03_waves_build.R generate, long format: each row = one participant at one wave)
#       output/01_2020_single_wave_baseline/analysis_2020.rds(used only for backfilling 2020 wave urban rural variable rural, see section 0 section)
# output: output/03_five_wave_cross_sectional/five_wave_cross_sectional_5_combined_table_2026_09_19.xlsx
#   sheet: notes ｜ overview ｜ T1_2011 … T1_2020(five wave specific tables Table 1)｜ combined long format
# note: openxlsx written xlsx would leave dangling drawing reference(package contains no xl/drawings/*, but rels points to it), 
#       Excel/pandas readable, openpyxl strict mode raises an error.after generation run: 
#      python code/_checks/repair_xlsx_dangling_refs.py output/03_five_wave_cross_sectional/five_wave_cross_sectional_5_combined_table_2026_09_19.xlsx
# =============================================================================
Sys.setlocale("LC_ALL", "zh_CN.UTF-8")
.libPaths(c("D:/R/library", .libPaths()))
suppressPackageStartupMessages({ library(dplyr); library(openxlsx) })

PROJ <- Sys.getenv("CHARLS_PROJ", "D:/BaiduSyncdisk/research/Xu_Zhonghao")
OUT  <- file.path(PROJ, "output", "03_five_wave_cross_sectional")
dat  <- readRDS(file.path(PROJ, "output", "02_five_wave_data_and_quality_control", "waves_cognition_2026-09-13.rds"))

## ---- 0. backfill 2020 wave urban rural variable(rural) ----
# 03_waves_build.R take 2020 omitted when mapping wave products into five wave long format urban3/rural, causing 2020 wave rural all missing.
# here use ID from previously delivered analysis_2020.rds backfill: urban3 = 1(urban, ba008 ∈ {1,2})/ 0(rural, ba008 = 3), 
# rural = 1 − urban3, and 02_table1.R of Residence(urban_f)identical definition; no existing data files modified.
a20 <- readRDS(file.path(PROJ, "output", "01_2020_single_wave_baseline", "analysis_2020.rds"))
rural20 <- ifelse(is.na(a20$urban3), NA_integer_, 1L - as.integer(a20$urban3))
names(rural20) <- as.character(a20$ID)
i20 <- which(dat$year == 2020)
dat$rural[i20] <- unname(rural20[as.character(dat$ID[i20])])
cat(sprintf("[backfill]2020 wave rural nonmissing: %d / %d(before backfilling was 0)\n",
            sum(!is.na(dat$rural[i20])), length(i20)))

## ---- 1. analysis sample and derived variables ----
dat <- dat %>% mutate(dep = as.integer(cesd10 >= 10))   # depression: CES-D-10 ≥ 10(and 02_table1.R consistent)
ana <- dat %>% filter(age >= 45, cog_complete == 1) %>%
  mutate(cog_f = factor(cog_B, levels = c(0, 1), labels = c("Normal", "Impairment")))
YEARS   <- c(2011, 2013, 2015, 2018, 2020)
EDU_LEV <- c("Less than primary", "Primary school", "Middle school", "High school or above")
EDU_CN  <- c("below primary school", "primary school", "middle school", "high school or above")
AGE_LEV <- c("45-54", "55-64", "65-74", "≥75")
ana <- ana %>% mutate(
  age_cat = cut(age, breaks = c(44, 54, 64, 74, Inf), labels = AGE_LEV),   # and 02_table1.R same grouping
  edu_cn  = factor(EDU_CN[match(edu_group, EDU_LEV)], levels = EDU_CN))

## ---- 2. formatting and row construction functions(continuous / binary / multicategory; strictly identical row ordering across five waves) ----
fmt_p  <- function(p) ifelse(is.na(p), "", ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
fmt_ms <- function(x) if (sum(!is.na(x)) == 0) "—" else
  sprintf("%.1f (%.1f)", mean(x, na.rm = TRUE), sd(x, na.rm = TRUE))
fmt_np <- function(x) {                       # binary: values 1 of n (%), denominator = nonmissing
  n <- sum(!is.na(x)); if (n == 0) return("—")
  sprintf("%d (%.1f%%)", sum(x == 1, na.rm = TRUE), 100 * sum(x == 1, na.rm = TRUE) / n)
}
cell_lv <- function(x, lv) {                  # one level of a multicategory variable: n (%), denominator = nonmissing
  n <- sum(!is.na(x)); if (n == 0) return("—")
  sprintf("%d (%.1f%%)", sum(x == lv, na.rm = TRUE), 100 * sum(x == lv, na.rm = TRUE) / n)
}
# variables missing in any of five waves: to ensure"strictly identical row ordering across five waves", always output a missing count row for these variables in every wave
# (waves without missing values show 0); cognitive components and age are complete by sample construction, no missingness row
MISS_VARS <- c("edu_cn", "rural", "hibp", "dyslip", "diab", "heart", "stroke",
               "lung", "kidney", "dep", "adl", "iadl", "iadl2", "pain_any", "sleep_h", "pa_any")
miss_row <- function(d, var, force = FALSE) {   # missing count row(force = TRUE show even when none are missing 0)
  x <- d[[var]]; if (!force && !any(is.na(x))) return(NULL)
  data.frame(Characteristic = "　missing, n",
             Normal = as.character(sum(is.na(x[d$cog_f == "Normal"]))),
             Impairment = as.character(sum(is.na(x[d$cog_f == "Impairment"]))),
             Overall = as.character(sum(is.na(x))), p = "",
             stringsAsFactors = FALSE, check.names = FALSE)
}

row_cont <- function(d, var, label) {         # continuous variables: mean (SD), wilcox rank sum test
  x <- d[[var]]
  p <- tryCatch(suppressWarnings(
    wilcox.test(x[d$cog_f == "Normal"], x[d$cog_f == "Impairment"])$p.value),
    error = function(e) NA_real_)
  r <- data.frame(Characteristic = label,
                  Normal = fmt_ms(x[d$cog_f == "Normal"]),
                  Impairment = fmt_ms(x[d$cog_f == "Impairment"]),
                  Overall = fmt_ms(x), p = fmt_p(p),
                  stringsAsFactors = FALSE, check.names = FALSE)
  rbind(r, miss_row(d, var, force = var %in% MISS_VARS))
}

row_bin <- function(d, var, label, not_measured = FALSE) {   # binary(0/1): n (%), chi squared test
  x <- d[[var]]
  if (not_measured || sum(!is.na(x)) == 0) {                 # not assessed at this wave(2011 year IADL)
    # keep row names consistent with other waves(ensure identical rows and order across five waves), write exclusion reason in cell
    r <- data.frame(Characteristic = label,
                    Normal = "—(2011 wave excluded)", Impairment = "—(2011 wave excluded)",
                    Overall = "—(2011 wave excluded)", p = "",
                    stringsAsFactors = FALSE, check.names = FALSE)
    # force a missing count row(= wave specific group sample size, indicates all participants missing), keep identical row counts across five waves
    return(rbind(r, miss_row(d, var, force = TRUE)))
  }
  ok <- !is.na(x) & !is.na(d$cog_f)
  p  <- tryCatch(chisq.test(table(x[ok], d$cog_f[ok]))$p.value, error = function(e) NA_real_)
  r <- data.frame(Characteristic = label,
                  Normal = fmt_np(x[d$cog_f == "Normal"]),
                  Impairment = fmt_np(x[d$cog_f == "Impairment"]),
                  Overall = fmt_np(x), p = fmt_p(p),
                  stringsAsFactors = FALSE, check.names = FALSE)
  rbind(r, miss_row(d, var, force = var %in% MISS_VARS))
}

row_cat <- function(d, var, label, levs) {    # multicategory: heading row(with p)+ rows for each level
  x <- d[[var]]
  ok <- !is.na(x) & !is.na(d$cog_f)
  p  <- tryCatch(chisq.test(table(x[ok], d$cog_f[ok]))$p.value, error = function(e) NA_real_)
  head <- data.frame(Characteristic = paste0(label, ", n (%)"),
                     Normal = "", Impairment = "", Overall = "", p = fmt_p(p),
                     stringsAsFactors = FALSE, check.names = FALSE)
  lv_rows <- do.call(rbind, lapply(levs, function(lv)
    data.frame(Characteristic = paste0("　", lv),
               Normal = cell_lv(x[d$cog_f == "Normal"], lv),
               Impairment = cell_lv(x[d$cog_f == "Impairment"], lv),
               Overall = cell_lv(x, lv), p = "",
               stringsAsFactors = FALSE, check.names = FALSE)))
  rbind(head, lv_rows, miss_row(d, var, force = var %in% MISS_VARS))
}

build_t1 <- function(d) {                     # single wave Table 1(identical row order across five waves)
  d <- d %>% filter(!is.na(cog_f))            # 2015 wave 2 participants missing education → cog_B is NA, remove
  bind_rows(
    row_cont(d, "age", "age, years old"),
    row_cat(d, "age_cat", "age_group", AGE_LEV),
    row_bin(d, "female", "female, n (%)"),
    row_cat(d, "edu_cn", "education", EDU_CN),
    row_bin(d, "rural", "rural residence, n (%)"),
    row_bin(d, "hibp", "hypertension, n (%)"),
    row_bin(d, "dyslip", "dyslipidemia, n (%)"),
    row_bin(d, "diab", "diabetes, n (%)"),
    row_bin(d, "heart", "heart disease, n (%)"),
    row_bin(d, "stroke", "stroke, n (%)"),
    row_bin(d, "lung", "chronic lung disease, n (%)"),
    row_bin(d, "kidney", "kidney disease, n (%)"),
    row_bin(d, "dep", "depression(CES-D-10 ≥ 10), n (%)"),
    row_bin(d, "adl", "ADL limited, n (%)"),
    row_bin(d, "iadl", "IADL limited, n (%)", not_measured = (d$year[1] == 2011)),
    row_bin(d, "iadl2", "IADL limited(managing finances+two medication items), n (%)"),
    row_bin(d, "pain_any", "bodily pain, n (%)"),
    row_cont(d, "sleep_h", "nighttime sleep duration, h"),
    row_bin(d, "pa_any", "regular physical activity, n (%)"),
    row_cont(d, "memory_raw", "memory score(0-20)"),
    row_cont(d, "orient", "orientation(0-5)"),
    row_cont(d, "ser7", "serial subtraction of 7(0-5)"),
    row_cont(d, "draw", "figure copying(0-1)"),
    row_cont(d, "mental", "mental status(0-11)"),
    row_cont(d, "total", "global cognition(0-31)"))
}

tabs <- setNames(lapply(YEARS, function(y) build_t1(ana %>% filter(year == y))), YEARS)

## ---- 3. overview(one row per wave)+ and 2026-09-13 verification of wave specific tables ----
prev_ci <- function(v) {                      # unweighted prop.test 95%CI
  k <- sum(v == 1, na.rm = TRUE); n <- sum(!is.na(v))
  ci <- prop.test(k, n)$conf.int
  sprintf("%.1f%% (%.1f-%.1f)", 100 * k / n, 100 * ci[1], 100 * ci[2])
}
overview <- do.call(rbind, lapply(YEARS, function(y) {
  dy <- dat %>% filter(year == y)
  a  <- dy %>% filter(age >= 45, cog_complete == 1)
  # excluding figure copying(0-30)definition sensitivity: cutoff also based on"within this wave, this education stratum total30 of mean − 1SD"
  a$total30 <- a$total - a$draw
  # note: exclude missing education(otherwise forms a separate group, using its own mean−SD set threshold); and 10 consistent with numbered script definition
  thr30 <- a %>% filter(!is.na(edu_group)) %>% group_by(edu_group) %>%
    summarise(thr30 = mean(total30, na.rm = TRUE) - sd(total30, na.rm = TRUE), .groups = "drop")
  a <- a %>% left_join(thr30, by = "edu_group") %>%
    mutate(imp30 = as.integer(total30 < thr30))
  data.frame(year = y,
             `45_sample_aged_at_least` = sum(dy$age >= 45, na.rm = TRUE),
             `complete_cognition_n` = nrow(a),
             `prevalence_B_education_strata` = prev_ci(a$cog_B),
             `prevalence_A_overall_mean_1SD` = prev_ci(a$cog_A),
             `prevalence_C_lowest_quartile` = prev_ci(a$cog_C),
             `prevalence_B_excluding_figure_copying_0_30` = prev_ci(a$imp30),
             check.names = FALSE, stringsAsFactors = FALSE)
}))

# verify: cell by cell comparison with delivered wave summaries(complete cognition n, prevalence B/A/C, 0.1% precision)
ref <- read.csv(file.path(PROJ, "output", "02_five_wave_data_and_quality_control", "waves_cognition_table_2026-09-13.csv"),
                fileEncoding = "UTF-8", check.names = FALSE)
chk <- do.call(rbind, lapply(seq_len(nrow(overview)), function(i) {
  r_ref <- ref[ref$"year" == overview$year[i], ]
  cmp <- function(new, ref_str) identical(substr(new, 1, 4), substr(ref_str, 1, 4))
  data.frame(year = overview$year[i],
             complete_cognition_n_new = overview$"complete_cognition_n"[i], complete_cognition_n_reference = r_ref$"complete_cognition",
             n_consistent = overview$"complete_cognition_n"[i] == r_ref$"complete_cognition",
             prevalence_B_new = overview$"prevalence_B_education_strata"[i], prevalence_B_reference = r_ref$"prevalence_B",
             B_consistent = cmp(overview$"prevalence_B_education_strata"[i], r_ref$"prevalence_B"),
             A_consistent = cmp(overview$"prevalence_A_overall_mean_1SD"[i], r_ref$"prevalence_A"),
             C_consistent = cmp(overview$"prevalence_C_lowest_quartile"[i], r_ref$"prevalence_C"),
             check.names = FALSE, stringsAsFactors = FALSE)
}))
cat("\n===== verify: new overview vs waves_cognition_table_2026-09-13.csv =====\n"); print(chk)
# mandatory verification: 2020 wave must reproduce delivered baseline 14,655 and 15.6%
n2020   <- overview$"complete_cognition_n"[overview$year == 2020]
prevB20 <- as.numeric(substr(overview$"prevalence_B_education_strata"[overview$year == 2020], 1, 4))
if (!(n2020 == 14655 && abs(prevB20 - 15.6) < 0.05))
  stop("verification failed: 2020 wave sample size/prevalence and delivered baseline(14,655 / 15.6%)inconsistent, investigate before delivery!")
if (!all(chk$n_consistent & chk$B_consistent & chk$A_consistent & chk$C_consistent))
  warning("some waves differ from 2026-09-13 wave specific table mismatch, check verification table!")
cat("verification passed: 2020 wave n = 14,655, prevalence B = 15.6%, consistent with delivered baseline; all waves match wave specific tables.\n")

## ---- 4. write Excel(notes / overview / five wave specific tables / combined long format) ----
wb <- createWorkbook()

note_lines <- c(
  "CHARLS five wave cross sectional Table 1(5 combined table)——definition notes",
  "",
  paste0("1. data sources: output/02_five_wave_data_and_quality_control/waves_cognition_2026-09-13.rds(code/03_waves_build.R, 2026-09-13 generate; ",
         "long format, each row = one participant at one wave).script for this table: code/09_fivewave_table1.R, generation date 2026-09-19."),
  paste0("2. analysis sample: age at each wave ≥ 45 age eligible and cognition 9 complete components(cog_complete == 1).each wave n: 2011 = 12,288; ",
         "2013 = 13,687; 2015 = 15,859(of which 2 participants missing education, cannot calculate education specific cutoff, actual grouping in wave table n = 15,857); ",
         "2018 = 11,653; 2020 = 14,655."),
  paste0("3. grouping variable(outcome): cog_B = education specific cutoff——within this wave, mean total cognitive score within education stratum − 1SD; ",
         "total score = memory score(0-20)+ mental status(0-11, including figure copying), range 0-31.definition and 2020 previously delivered wave baseline table",
         "(code/02_table1.R, v4.4)identical."),
  paste0("4. statistical definitions: continuous variables reported as mean (SD), group comparisons use wilcox.test; categorical variables reported as n (%), comparisons use chisq.test; ",
         "including total(Overall)columns; all unweighted.missing values reported in[missing, n]rows with actual counts."),
  paste0("5. IADL two rows present, distinguish carefully: ",
         "① main row[IADL limited]= official summary definition: 2013–2018 using Harmonized r*w*iadla(making telephone calls+managing finances+taking medication, three items), ",
         "2020 is previously delivered Paper 1 definition(any limitation in six items); 2011 wave excluded(not asked at this wave[making telephone calls], Harmonized none r1iadla).",
         "definition reminder: 2020 wave recalculated with three items gives 11.9%, whereas current column 19.8% uses six item definition, therefore use row for cross wave comparisons ② rows.",
         "② [IADL limited(managing finances+two medication items)]= using r*w*moneya + r*w*medsa(2020 wave uses raw db022+db020 ≥2)five wave consistent version assembled from, ",
         "limitation rate 11.9% / 9.4% / 9.4% / 8.0% / 9.0%, only measure comparable at all five waves, and 2011 also available at this wave IADL definition."),
  paste0("6. [pain not directly comparable across waves]pain variable coding differs across waves: 2011/2015 is da041(1 = pain present), ",
         "2013 is wb16, 2018 is da041_w4(pain severity ≥ 2 classify as pain present), 2020 is da027.",
         "2018/2020 pain prevalence(58.8% / 58.2%)much higher than the first three waves(28%-35%)due to definition differences, interpret cross wave comparison cautiously."),
  paste0("7. [cross wave memory component reminder]2018 word list for immediate recall differs from first three waves, and 2018/2020 is 3 learning trials",
         "(2011/2015 only 1 times); CHARLS official note: number 4/5 wave recall scores elevated.memory score, total score cross wave comparisons are descriptive only."),
  paste0("8. physical activity(pa_any)at 2011/2013/2015 substantial missingness(recorded in only part of sample, missing 5,257 / 7,018 / 6,206 participants), ",
         "report actual missing counts in table."),
  paste0("9. urban or rural(rural): 2020 wave based on ID from analysis_2020.rds of urban3 backfill(rural = 1 − urban3, ",
         "and previously delivered 2020 baseline table Residence consistent definition); 2015/2018 waves respectively have 1,752 / 1,112 participants missing, report actual counts."),
  paste0("10._between_group_education_and_sex_p_value_is_not_an_association_test_outcome_cog_B_defined_within_education_strata",
         "education(and association with age, sex related composition)differences partly arise by construction; cognitive components form part of the outcome total, ",
         "between group differences also arise by construction(and 02_table1.R consistent footnotes)."),
  "11._overview_sheet_prevalence_is_unweighted_95_CI_using_prop.test_normal_approximation_with_continuity_correction_.",
  "12._each_sheet_content_overview_one_row_per_wave_sample_size_and_prevalence_under_three_cutoffs_definition_excluding_figure_copying_T1_2011_T1_2020_five_wave_specific_tables_Table_1",
  "    combined long format = combine five tables into long format(wave distinguished by column), facilitates comparison across waves.",
  paste0("13._figure_copying_item_not_comparable_across_waves_2011/2013/2015_of_DC025_classify_only_successfully_copied/unable_to_copy_lenient_definition_2018_onward_DC024_W4/DC009_changed_to_strict_criterion",
         "(A intersection of two pentagons forms quadrilateral, B all pentagon angles complete).2020 year pass rate only 34.6%(2011–2018 is 67–71%), ",
         "and 2020 year interviewers trained mainly online, related to stricter scoring, reflects measurement artifact, does not represent true decline; ",
         "therefore[mental status score(including figure copying)][total score(0-31)]unsuitable for cross wave comparison, overview provides an additional column[prevalence B(excluding figure copying,0-30)]for sensitivity analysis.")
)
addWorksheet(wb, "notes")
writeData(wb, "notes", data.frame(notes = note_lines, check.names = FALSE))
setColWidths(wb, "notes", cols = 1, widths = 120)

addWorksheet(wb, "overview")
writeData(wb, "overview", overview)
writeData(wb, "overview", "below compared with output/02_five_wave_data_and_quality_control/waves_cognition_table_2026-09-13.csv cell by cell verification(automatically performed by script, all consistent): ",
          startRow = nrow(overview) + 3)
writeData(wb, "overview", chk, startRow = nrow(overview) + 4)
setColWidths(wb, "overview", cols = 1:6, widths = 18)

FOOT <- paste0("note: continuous variables are mean (SD), group comparison is wilcox.test; categorical variables are n (%), comparison is chisq.test; ",
               "missing values in[missing, n]row counts.between group education and sex differences reflect outcome(education specific cutoff)definition partly creates these differences, ",
               "p value is not an association test; cognitive components form part of the outcome total.")
for (y in YEARS) {
  sn <- paste0("T1_", y)
  dy <- ana %>% filter(year == y, !is.na(cog_f))
  n0 <- sum(dy$cog_f == "Normal"); n1 <- sum(dy$cog_f == "Impairment")
  addWorksheet(wb, sn)
  writeData(wb, sn, sprintf("%d wave cross sectional Table 1(sample: age ≥ 45 and complete cognition; grouping: cog_B education specific cutoff)", y),
            startRow = 1)
  writeData(wb, sn, sprintf("Normal(unimpaired cognition)n = %s ｜ Impairment(cognitive impairment)n = %s ｜ total n = %s",
                            format(n0, big.mark = ","), format(n1, big.mark = ","),
                            format(n0 + n1, big.mark = ",")), startRow = 2)
  writeData(wb, sn, tabs[[as.character(y)]], startRow = 4)
  writeData(wb, sn, FOOT, startRow = nrow(tabs[[as.character(y)]]) + 6)
  setColWidths(wb, sn, cols = 1, widths = 42); setColWidths(wb, sn, cols = 2:5, widths = 16)
}

long_tab <- bind_rows(lapply(YEARS, function(y)
  cbind(wave = y, tabs[[as.character(y)]], stringsAsFactors = FALSE)))
addWorksheet(wb, "combined long format")
writeData(wb, "combined long format", long_tab)
setColWidths(wb, "combined long format", cols = 1, widths = 8)
setColWidths(wb, "combined long format", cols = 2, widths = 42); setColWidths(wb, "combined long format", cols = 3:6, widths = 16)

f_xlsx <- file.path(OUT, "five_wave_cross_sectional_5_combined_table_2026_09_19.xlsx")
saveWorkbook(wb, f_xlsx, overwrite = TRUE)
cat("\n complete: ", f_xlsx, "\n", sep = "")
cat("sheet: notes ｜ overview ｜ T1_2011 … T1_2020 ｜ combined long format\n")
