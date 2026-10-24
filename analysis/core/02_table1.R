# =============================================================================
# 02_table1.R  --  CHARLS 2020 baseline characteristics table(Table 1)+ prevalence(v4.4, 2026-09-12)
# definition: primary analysis = cognition 9 complete components(education completed, no missing values); cutoff = education strata(<within stratum mean−1SD)
# version: v4.4(2026-09-12)——① TableS1 restore edu_f(education gradient in selection bias, no missing education in either group after backfilling); 
#                      ② S2/S4/S5 of CSV quote final footnote row as a valid field(unquoted commas previously produced extra columns); 
#                      ③ education strata(B/A two definitions)harmonize to weighted estimates + 95%CI, also report unweighted point estimates; 
#                      ④ TableS2/S3/S4 additional outputs n columns(actual sample size of each model)
# input: output/01_2020_single_wave_baseline/analysis_2020.rds
# output: output/01_2020_single_wave_baseline/Table1_baseline.docx / .html, Table1_and_prevalence.xlsx, 
#       TableS1_complete_vs_incomplete.csv, TableS2_pain_sleep_varform.csv(weighted), 
#       TableS3_comorbidity_adj.csv(chronic conditions/function/weighted adjusted associations of behavioral variables), 
#       TableS4_sensitivityA_adj.csv(sensitivity A weighted adjusted associations under this definition, and S3 symmetric)
# =============================================================================
Sys.setlocale("LC_ALL", "zh_CN.UTF-8")
.libPaths(c("D:/R/library", .libPaths()))
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(gtsummary); library(survey)
  library(openxlsx); library(flextable); library(srvyr)
})
proj <- "D:/BaiduSyncdisk/research/Xu_Zhonghao"; out <- file.path(proj, "output", "01_2020_single_wave_baseline")
dat_all <- readRDS(file.path(out, "analysis_2020.rds"))

## ---- 1. sample partition ----
dat    <- dat_all %>% filter(primary == 1, !is.na(cog_impair))   # primary analysis sample
# note: noedu(complete cognition but missing education)at 2026-09-12 always empty after completing education; 
#     retained only for backward compatibility, below S1 group construction automatically skips empty groups.
noedu  <- dat_all %>% filter(primary == 1,  is.na(cog_impair))
incomp <- dat_all %>% filter(primary == 0)                       # incomplete assessment(sensitivity)
cat(sprintf("primary analysis n=%d | complete cognition but missing education n=%d | incomplete assessment n=%d\n",
            nrow(dat), nrow(noedu), nrow(incomp)))

## ---- 2. variable preparation and labels ----
prep <- function(d) d %>% mutate(
  age_cat   = cut(age, breaks = c(44, 54, 64, 74, Inf),
                  labels = c("45-54", "55-64", "65-74", "≥75")),
  female_f  = factor(female, c(0,1), c("Male","Female")),
  edu_f     = factor(edu_group,
                     levels = c("Less than primary","Primary school","Middle school","High school or above")),
  married_f = factor(ifelse(marital %in% c(1,2), "Married/partnered", "Not partnered"),
                     levels = c("Married/partnered","Not partnered")),
  urban_f   = factor(urban3, c(0,1), c("Rural","Urban")),
  hukou_f   = factor(case_when(hukou == 1 ~ "Agricultural",
                               hukou %in% c(2,3) ~ "Non-agricultural/Unified"),
                     levels = c("Agricultural","Non-agricultural/Unified")),
  pain_f    = factor(pain_any, c(0,1), c("No pain","Any pain")),
  pain3_f   = factor(pain_3cat, c(0,1,2), c("No pain","Mild","Severe")),
  sleep_f   = factor(sleep_cat, levels = c("<6h","6-8h",">8h")),
  sleep_ref = factor(sleep_cat, levels = c("6-8h","<6h",">8h")),   # model reference is 6-8h as reference
  dep_f     = factor(depression, c(0,1), c("No","Yes")),
  adl_f     = factor(adl_any,  c(0,1), c("No","Yes")),
  iadl_f    = factor(iadl_any, c(0,1), c("No","Yes")),
  multi_f   = factor(multimorbidity, c(0,1), c("No","Yes")),
  cog_f     = factor(cog_impair, levels = c(0,1), labels = c("Normal","Impairment")),
  hypertension_f   = factor(hypertension,   c(0,1), c("No","Yes")),
  dyslipidemia_f   = factor(dyslipidemia,   c(0,1), c("No","Yes")),
  diabetes_f       = factor(diabetes,       c(0,1), c("No","Yes")),
  heart_disease_f  = factor(heart_disease,  c(0,1), c("No","Yes")),
  stroke_f         = factor(stroke,         c(0,1), c("No","Yes")),
  lung_disease_f   = factor(lung_disease,   c(0,1), c("No","Yes")),
  arthritis_f      = factor(arthritis,      c(0,1), c("No","Yes")),
  kidney_disease_f = factor(kidney_disease, c(0,1), c("No","Yes")),
  smoke_f   = factor(smoke, c(0,1,2), c("Never","Former","Current")),
  drink_f   = factor(drink_freq, c(3,2,1),
                     c("None in past year","Less than monthly","At least monthly")),
  pa_any_f  = factor(pa_any, c(0,1), c("No","Yes")),
  sleep_restless_f = factor(sleep_restless, c(0,1), c("No","Yes"))
)
dat <- prep(dat); noedu <- prep(noedu); incomp <- prep(incomp)

## ---- 3. Table 1(primary analysis sample, grouped by cognitive status) ----
t1_vars <- c("age", "age_cat", "female_f", "edu_f", "married_f", "urban_f", "hukou_f",
             "hypertension_f","dyslipidemia_f","diabetes_f","heart_disease_f","stroke_f",
             "lung_disease_f","arthritis_f","kidney_disease_f","multi_f","n_chronic",
             "pain3_f","sleep_h","sleep_f","sleep_restless_f",
             "smoke_f","drink_f","pa_any_f",
             "dep_f","cesd10","adl_f","iadl_f",
             "imrc","dlrc","memory","orient","ser7","draw","mental","global")
labels1 <- list(
  age ~ "Age, years", age_cat ~ "Age group, n (%)", female_f ~ "Sex", edu_f ~ "Education",
  married_f ~ "Marital status", urban_f ~ "Residence", hukou_f ~ "Hukou type",
  hypertension_f ~ "Hypertension", dyslipidemia_f ~ "Dyslipidemia", diabetes_f ~ "Diabetes",
  heart_disease_f ~ "Heart disease", stroke_f ~ "Stroke", lung_disease_f ~ "Chronic lung disease",
  arthritis_f ~ "Arthritis/rheumatism", kidney_disease_f ~ "Kidney disease",
  multi_f ~ "Multimorbidity (≥2)", n_chronic ~ "Number of chronic diseases",
  pain_f ~ "Body pain", pain3_f ~ "Body pain severity", sleep_h ~ "Night sleep duration, h", sleep_f ~ "Sleep duration category",
  sleep_restless_f ~ "Restless sleep (≥3-4 days/week)",
  smoke_f ~ "Smoking status", drink_f ~ "Alcohol use, past year",
  pa_any_f ~ "Any regular physical activity (≥10 min)",
  dep_f ~ "Depression (CES-D-10 ≥ 10)", cesd10 ~ "CES-D-10 score",
  adl_f ~ "ADL limitation", iadl_f ~ "IADL limitation",
  imrc ~ "Immediate word recall (0-10)", dlrc ~ "Delayed word recall (0-10)",
  memory ~ "Episodic memory (0-20)", orient ~ "Orientation (0-5)", ser7 ~ "Serial 7s (0-5)",
  draw ~ "Figure drawing (0-1)", mental ~ "Mental status (0-11)", global ~ "Global cognition (0-31)")
# select labels dynamically for variables actually tabulated(avoid label reference nonexistent columns)
labels_for <- function(vars) labels1[vapply(labels1, function(f) as.character(f[[2]]), character(1)) %in% vars]

tab1 <- dat %>% select(all_of(c("cog_f", t1_vars))) %>%
  tbl_summary(by = cog_f,
    statistic = list(all_continuous() ~ "{mean} ({sd})", all_categorical() ~ "{n} ({p}%)"),
    digits = list(all_continuous() ~ 1, all_categorical() ~ c(0, 1)),
    missing = "ifany", label = labels_for(t1_vars)) %>%
  add_p(test = list(all_continuous() ~ "wilcox.test", all_categorical() ~ "chisq.test")) %>%
  add_overall(last = TRUE) %>% bold_labels() %>% modify_header(label ~ "**Characteristic**") %>%
  modify_table_styling(columns = label,
    footnote = "Cognitive impairment was defined with education-stratified cutoffs; between-group comparisons for education and sex are by design and not interpretable as tests of association.")

## ---- 4. Table S1: primary analysis vs incomplete assessment(missingness bias assessment) ----
# 2026-09-12(v4.2): edu_f restored s1_vars.
#   previous rationale("missing education defines group 2 group definition")corresponds to groups before backfilling: at that time group 2 groups = complete cognition but missing education.
#   after backfilling group 2 group is now"incomplete cognitive assessment", and neither group has missing education, education gradient(below primary school 71.7% vs 36.8%)
#   is a key selection bias indicator, must accompany the table, otherwise manuscript numbers cannot be found in the supplement.
s1_vars <- c("age", "age_cat", "female_f", "edu_f", "urban_f", "married_f",
             "hypertension_f","diabetes_f","stroke_f","arthritis_f",
             "pain3_f","dep_f","adl_f","iadl_f","sleep_h")
# 2026-09-12: education via zredu+ba010 no missing values remain after combining, empty groups removed automatically
grp_list <- list(
  "Analytic sample (complete cognition)"       = dat,
  "Complete cognition, education missing"      = noedu,
  "Incomplete cognition (sensitivity only)"    = incomp
)
grp_list <- grp_list[vapply(grp_list, nrow, integer(1)) > 0]
grp <- bind_rows(lapply(names(grp_list), function(g)
  grp_list[[g]] %>% select(all_of(s1_vars)) %>% mutate(grp = g))) %>%
  mutate(grp = factor(grp, levels = names(grp_list)))
tabS1 <- grp %>% tbl_summary(by = grp,
    statistic = list(all_continuous() ~ "{mean} ({sd})", all_categorical() ~ "{n} ({p}%)"),
    digits = list(all_categorical() ~ c(0,1)), missing = "no",
    label = labels_for(s1_vars)) %>%
  add_p(test = list(all_continuous() ~ "kruskal.test", all_categorical() ~ "chisq.test")) %>%
  bold_labels()

## ---- 4b. Table S2: pain/sleep"variable specification"sensitivity analysis(crude model + adjusted model OR) ----
# purpose: show both binary and ordinal specifications, avoid"reporting only significant findings"concern about
# definition: survey weighted logistic(INDV_weight, community clustering); adjustment age+sex+education (four groups)+urban or rural
des_main <- dat %>% filter(!is.na(INDV_weight), !is.na(communityID)) %>%
  as_survey_design(weights = INDV_weight, ids = communityID)
# note: report exposure coefficients only; covariates(including edu_f)coefficients not reported——education defines the outcome specific cutoffs, 
#     its coefficient is not interpretable, reporting could be misread as"higher education, higher cognitive impairment risk".
or_block <- function(expo, label) {
  adj <- " + age + female_f + edu_f + urban_f"
  grab <- function(f, tag) {
    m <- svyglm(f, design = des_main, family = quasibinomial)
    s <- summary(m)$coefficients
    s <- s[grepl(paste0("^", expo), rownames(s)), , drop = FALSE]   # exposure terms only
    data.frame(exposure = label, term = rownames(s), model = tag,
               OR_95CI = sprintf("%.2f (%.2f-%.2f)", exp(s[,1]),
                                 exp(s[,1] - 1.96*s[,2]), exp(s[,1] + 1.96*s[,2])),
               p = sprintf("%.3g", s[,4]),
               # actual sample size entering this model(2026-09-12 supplement: missingness differs across models, n different)
               n = nrow(model.frame(m)), row.names = NULL, stringsAsFactors = FALSE)
  }
  rbind(grab(as.formula(paste0("cog_impair ~ ", expo)), "Crude"),
        grab(as.formula(paste0("cog_impair ~ ", expo, adj)), "Adjusted"))
}
NOTE_S2 <- paste0("Note: all models are survey-weighted logistic regressions (individual weight INDV_weight, ",
                  "community as cluster), adjusted for age, sex, education (4 groups) and residence. ",
                  "Only exposure coefficients are reported; covariate coefficients (including education) are ",
                  "deliberately omitted because education defines the stratified cut-off used for the outcome ",
                  "and is therefore not interpretable.")
tabS2 <- bind_rows(
  or_block("pain_f",        "Body pain, binary (any vs none)"),
  or_block("pain3_f",       "Body pain severity (3 levels, ref = no pain)"),
  or_block("sleep_h",       "Sleep duration per 1 h (continuous)"),
  or_block("sleep_ref",     "Sleep duration grouped (ref = 6-8 h)"),
  or_block("sleep_restless_f", "Restless sleep")
)

## ---- 4c. Table S3: chronic conditions and function/weighted adjusted associations of behavioral variables ----
# and Table S2 same definition: svyglm weighted, adjustment age+sex+education (four groups)+urban or rural
run_adj <- function(v, label){
  f <- as.formula(paste0("cog_impair ~ ", v, " + age + female_f + edu_f + urban_f"))
  m <- svyglm(f, design = des_main, family = quasibinomial)
  s <- summary(m)$coefficients
  i <- grep(paste0("^", v), rownames(s))
  data.frame(exposure = label, term = rownames(s)[i],
             OR_95CI = sprintf("%.2f (%.2f-%.2f)", exp(s[i,1]),
                               exp(s[i,1] - 1.96*s[i,2]), exp(s[i,1] + 1.96*s[i,2])),
             p = sprintf("%.3g", s[i,4]),
             n = nrow(model.frame(m)), row.names = NULL, stringsAsFactors = FALSE)
}
tabS3 <- bind_rows(
  run_adj("hypertension",   "Hypertension"),
  run_adj("stroke",         "Stroke"),
  run_adj("arthritis",      "Arthritis/rheumatism"),
  run_adj("diabetes",       "Diabetes"),
  run_adj("lung_disease",   "Chronic lung disease"),
  run_adj("heart_disease",  "Heart disease"),
  run_adj("kidney_disease", "Kidney disease"),
  run_adj("n_chronic",      "Number of chronic diseases"),
  run_adj("multimorbidity", "Multimorbidity (>=2)"),
  run_adj("depression",     "Depression (CES-D-10 >= 10)"),
  run_adj("adl_any",        "ADL limitation"),
  run_adj("iadl_any",       "IADL limitation"),
  run_adj("pa_any",         "Regular physical activity"),
  run_adj("dyslipidemia",   "Dyslipidemia"),
  run_adj("drink_f",        "Alcohol use, past year (ref = none)")
)

## ---- 4d. Table S4: sensitivity A definition(overall mean−1SD cutoff)weighted adjusted associations under ----
# and S3 symmetric: same exposures, same weights and adjustment set, only outcome changes to cog_low_1sd.
# purpose: stroke/multimorbidity/alcohol and related findings depend on cutoff definitions, report alongside primary definition.
run_adj_A <- function(v, label){
  f <- as.formula(paste0("cog_low_1sd ~ ", v, " + age + female_f + edu_f + urban_f"))
  m <- svyglm(f, design = des_main, family = quasibinomial)
  s <- summary(m)$coefficients
  i <- grep(paste0("^", v), rownames(s))
  data.frame(exposure = label, term = rownames(s)[i],
             OR_95CI = sprintf("%.2f (%.2f-%.2f)", exp(s[i,1]),
                               exp(s[i,1] - 1.96*s[i,2]), exp(s[i,1] + 1.96*s[i,2])),
             p = sprintf("%.3g", s[i,4]),
             n = nrow(model.frame(m)), row.names = NULL, stringsAsFactors = FALSE)
}
tabS4 <- bind_rows(
  run_adj_A("hypertension",   "Hypertension"),
  run_adj_A("stroke",         "Stroke"),
  run_adj_A("arthritis",      "Arthritis/rheumatism"),
  run_adj_A("diabetes",       "Diabetes"),
  run_adj_A("lung_disease",   "Chronic lung disease"),
  run_adj_A("heart_disease",  "Heart disease"),
  run_adj_A("kidney_disease", "Kidney disease"),
  run_adj_A("n_chronic",      "Number of chronic diseases"),
  run_adj_A("multimorbidity", "Multimorbidity (>=2)"),
  run_adj_A("depression",     "Depression (CES-D-10 >= 10)"),
  run_adj_A("adl_any",        "ADL limitation"),
  run_adj_A("iadl_any",       "IADL limitation"),
  run_adj_A("pa_any",         "Regular physical activity"),
  run_adj_A("dyslipidemia",   "Dyslipidemia"),
  run_adj_A("drink_f",        "Alcohol use, past year (ref = none)")
)

## ---- 4e. Table S5: three cutoff definitions(B / A / C)adjusted associations of core exposures under ----
# purpose: SAP §6 of S1(definition C)and S2(definition A)sensitivity analysis requires core exposure models under both definitions; 
#       previously only chronic conditions/function/behavioral variables analyzed A definition, **sleep categories and three pain levels not modeled under A/C modeled under this definition**.
# definition: survey weighted logistic; adjustment age+sex+education (four groups)+urban or rural; outcomes respectively use B/A/C three cutoff definitions.
core <- list(
  list(v = "depression", lab = "Depression (CES-D-10 >= 10)"),
  list(v = "adl_any",    lab = "ADL limitation"),
  list(v = "iadl_any",   lab = "IADL limitation"),
  list(v = "sleep_ref",  lab = "Sleep duration grouped (ref = 6-8 h)"),
  list(v = "pa_any",     lab = "Regular physical activity"),
  list(v = "pain3_f",    lab = "Body pain severity (ref = no pain)")
)
outc <- c(B = "cog_impair", A = "cog_low_1sd", C = "cog_low_q1")
one_oc <- function(v, oc) {
  f <- as.formula(paste0(oc, " ~ ", v, " + age + female_f + edu_f + urban_f"))
  m <- svyglm(f, design = des_main, family = quasibinomial)
  s <- summary(m)$coefficients
  i <- grep(paste0("^", v), rownames(s))
  data.frame(term = rownames(s)[i],
             OR = sprintf("%.2f (%.2f-%.2f)", exp(s[i,1]),
                          exp(s[i,1] - 1.96*s[i,2]), exp(s[i,1] + 1.96*s[i,2])),
             p  = sprintf("%.3g", s[i,4]), row.names = NULL, stringsAsFactors = FALSE)
}
tabS5 <- do.call(rbind, lapply(core, function(cx) {
  res <- lapply(names(outc), function(k) one_oc(cx$v, outc[[k]]))
  base <- res[[1]][, "term", drop = FALSE]
  for (j in seq_along(outc)) {
    base[[paste0(names(outc)[j], "_OR_95CI")]] <- res[[j]]$OR
    base[[paste0(names(outc)[j], "_p")]]       <- res[[j]]$p
  }
  cbind(exposure = cx$lab, base, row.names = NULL)
}))
NOTE_S5 <- paste0("Note: B = education-stratified cut-off (primary), A = overall mean-1SD, ",
                  "C = lowest quartile. All models are survey-weighted logistic regressions ",
                  "(individual weight INDV_weight, community as cluster), adjusted for age, sex, ",
                  "education (4 groups) and residence. Only exposure coefficients are reported.")

## ---- 5. prevalence(primary definition + sensitivity) ----
# des_main already in 4b created in section(of the primary analysis sample survey design)
mk_design <- function(d) d %>% filter(!is.na(INDV_weight), !is.na(communityID)) %>%
  as_survey_design(weights = INDV_weight, ids = communityID)
prev_row <- function(var, label, d, design = NULL) {
  y <- d[[var]]; n <- sum(!is.na(y)); k <- sum(y == 1, na.rm = TRUE)
  ci_u <- prop.test(k, n)$conf.int
  if (!is.null(design)) {
    f <- svymean(as.formula(paste0("~", var)), design, na.rm = TRUE); ci_w <- confint(f)
    w <- sprintf("%.1f%% (%.1f-%.1f)", 100*coef(f)[1], 100*ci_w[1], 100*ci_w[2])
    nw <- nrow(design$variables)
    yw <- design$variables[[var]]; nwv <- sum(!is.na(yw)); kw <- sum(yw == 1, na.rm = TRUE)
    ci_uw <- prop.test(kw, nwv)$conf.int
    uw <- sprintf("%.1f%% (%.1f-%.1f)", 100*kw/nwv, 100*ci_uw[1], 100*ci_uw[2])
  } else { w <- NA_character_; uw <- NA_character_; nw <- NA_integer_ }
  data.frame(definition = label, n = n, unweighted = sprintf("%.1f%% (%.1f-%.1f)", 100*k/n, 100*ci_u[1], 100*ci_u[2]),
             unweighted_wsub = uw, weighted = w, n_weighted = nw, stringsAsFactors = FALSE)
}
des_all  <- mk_design(dat_all %>% filter(!is.na(cog_impair)))
prev_tab <- bind_rows(
  prev_row("cog_impair",   "[primary definition]education specific cutoff, participants with complete cognition",        dat, des_main),
  prev_row("cog_low_1sd",  "sensitivity A: overall mean−1SD, participants with complete cognition",        dat, des_main),
  prev_row("cog_low_q1",   "sensitivity C: lowest quartile, participants with complete cognition",          dat, des_main),
  prev_row("cog_impair",   "sensitivity D: education specific cutoff, include incomplete assessments",      dat_all %>% filter(!is.na(cog_impair)), des_all)
)

## ---- 6. stratified prevalence ----
strat <- function(d, byvar, label, design, def = "cog_impair", tag = "primary") {
  dd <- d %>% filter(!is.na(.data[[byvar]]))
  ds <- dd %>% filter(!is.na(INDV_weight), !is.na(communityID)) %>%
    as_survey_design(weights = INDV_weight, ids = communityID)
  ds %>% group_by(.data[[byvar]]) %>%
    summarise(p = survey_mean(.data[[def]], na.rm = TRUE, vartype = "ci"), n = unweighted(n())) %>%
    mutate(stratum = paste0(label, " [", tag, "]"), level = as.character(.data[[byvar]]),
           prevalence = sprintf("%.1f%% (%.1f-%.1f)", 100*p, 100*p_low, 100*p_upp)) %>%
    select(stratum, level, prevalence, n) %>% as.data.frame()
}
strat_tab <- bind_rows(
  strat(dat, "age_cat",  "Age group",  des_main, "cog_impair",  "primary definition B"),
  strat(dat, "female_f", "Sex",        des_main, "cog_impair",  "primary definition B"),
  strat(dat, "urban_f",  "Residence",  des_main, "cog_impair",  "primary definition B"),
  strat(dat, "age_cat",  "Age group",  des_main, "cog_low_1sd", "sensitivity A"),
  strat(dat, "female_f", "Sex",        des_main, "cog_low_1sd", "sensitivity A"),
  strat(dat, "urban_f",  "Residence",  des_main, "cog_low_1sd", "sensitivity A")
) %>% mutate(unweighted = NA_character_)
# education gradient(v4.3, 2026-09-12): harmonize with other strata as[weighted + 95%CI], listed separately[unweighted]point estimate(full sample n=14,655).
#   primary definition B cutoffs defined within education strata, B education gradient is compressed by construction(weighted approximately 13-17%, between groups CI overlap); 
#   sensitivity A(overall cutoff)education gradient under this definition reflects the underlying pattern(weighted 34.8% -> 1.8%).
edu_unw <- bind_rows(lapply(c("primary definition B", "sensitivity A"), function(tag) {
  v <- if (tag == "primary definition B") "cog_low_edu" else "cog_low_1sd"
  dat %>% filter(!is.na(edu_f)) %>% group_by(edu_f) %>%
    summarise(p = mean(.data[[v]], na.rm = TRUE), .groups = "drop") %>%
    transmute(tag = tag, level = as.character(edu_f), unweighted = sprintf("%.1f%%", 100*p))
}))
strat_edu <- bind_rows(lapply(c("primary definition B", "sensitivity A"), function(tag) {
  v <- if (tag == "primary definition B") "cog_low_edu" else "cog_low_1sd"
  ds <- dat %>% filter(!is.na(edu_f), !is.na(INDV_weight), !is.na(communityID)) %>%
    as_survey_design(weights = INDV_weight, ids = communityID)
  ds %>% group_by(edu_f) %>%
    summarise(p = survey_mean(.data[[v]], na.rm = TRUE, vartype = "ci"), n = unweighted(n())) %>%
    mutate(stratum = paste0("Education [", tag, ", weighted]"), level = as.character(edu_f),
           prevalence = sprintf("%.1f%% (%.1f-%.1f)", 100*p, 100*p_low, 100*p_upp), tag = tag) %>%
    select(tag, stratum, level, prevalence, n) %>% as.data.frame()
})) %>% left_join(edu_unw, by = c("tag", "level")) %>%
  select(stratum, level, prevalence, n, unweighted) %>% as.data.frame()

strat_all <- bind_rows(strat_tab, strat_edu)
NOTE_STRAT <- paste0("Education rows: prevalence = survey-weighted estimate (95% CI) on the weighted-eligible ",
  "subsample; unweighted = point estimate on the full analytic sample (n=14,655). Other strata: weighted ",
  "estimates only; n = weighted-eligible respondents.")

## ---- 7. output ----
dir.create(out, showWarnings = FALSE)
print(tab1)
as_gt(tab1) %>% gt::gtsave(file.path(out, "Table1_baseline.html"))
as_gt(tabS1) %>% gt::gtsave(file.path(out, "TableS1_complete_vs_incomplete.html"))
write.csv(as.data.frame(tabS1$table_body), file.path(out, "TableS1_complete_vs_incomplete.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
write.csv(tabS2, file.path(out, "TableS2_pain_sleep_varform.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
cat(sprintf('"# %s"\n', NOTE_S2), file = file.path(out, "TableS2_pain_sleep_varform.csv"), append = TRUE)
write.csv(tabS3, file.path(out, "TableS3_comorbidity_adj.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
write.csv(tabS4, file.path(out, "TableS4_sensitivityA_adj.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
cat(sprintf('"# %s"\n', NOTE_S2), file = file.path(out, "TableS4_sensitivityA_adj.csv"), append = TRUE)
write.csv(tabS5, file.path(out, "TableS5_three_cutpoints_core.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
cat(sprintf('"# %s"\n', NOTE_S5), file = file.path(out, "TableS5_three_cutpoints_core.csv"), append = TRUE)

w1 <- createWorkbook()
addWorksheet(w1, "Table1_raw");            writeData(w1, "Table1_raw", as_tibble(tab1$table_body))
addWorksheet(w1, "TableS1_raw");           writeData(w1, "TableS1_raw", as_tibble(tabS1$table_body))
addWorksheet(w1, "TableS2_varform");       writeData(w1, "TableS2_varform", tabS2)
writeData(w1, "TableS2_varform", NOTE_S2, startRow = nrow(tabS2) + 3, startCol = 1)
addWorksheet(w1, "TableS3_comorbidity");   writeData(w1, "TableS3_comorbidity", tabS3)
writeData(w1, "TableS3_comorbidity", NOTE_S2, startRow = nrow(tabS3) + 3, startCol = 1)
addWorksheet(w1, "TableS4_sensitivityA");  writeData(w1, "TableS4_sensitivityA", tabS4)
writeData(w1, "TableS4_sensitivityA", NOTE_S2, startRow = nrow(tabS4) + 3, startCol = 1)
addWorksheet(w1, "TableS5_three_cutpoints"); writeData(w1, "TableS5_three_cutpoints", tabS5)
writeData(w1, "TableS5_three_cutpoints", NOTE_S5, startRow = nrow(tabS5) + 3, startCol = 1)
addWorksheet(w1, "Prevalence");            writeData(w1, "Prevalence", prev_tab)
addWorksheet(w1, "Prevalence_stratified"); writeData(w1, "Prevalence_stratified", strat_all)
writeData(w1, "Prevalence_stratified", NOTE_STRAT, startRow = nrow(strat_all) + 3, startCol = 1)
saveWorkbook(w1, file.path(out, "Table1_and_prevalence.xlsx"), overwrite = TRUE)

# docx write last: if locked by Word in use, does not affect html/xlsx generation
save_as_docx(as_flex_table(tab1), path = file.path(out, "Table1_baseline.docx"))

cat("\n=== cognitive impairment prevalence ===\n"); print(prev_tab)
cat("\n=== stratification ===\n"); print(strat_all)
cat("\n=== Table S2: variable specification sensitivity analysis ===\n"); print(tabS2)
cat("\n=== Table S3: chronic conditions and function/weighted adjusted associations of behavioral variables ===\n"); print(tabS3)
cat("\n=== Table S4: sensitivity A weighted adjusted associations by definition ===\n"); print(tabS4)
cat("\n=== Table S5: core exposures under three definitions ===\n"); print(tabS5)
cat("\n complete: Table1_baseline.docx / .html, TableS1_*.html/.csv, TableS2/S3/S4/S5_*.csv, Table1_and_prevalence.xlsx\n")
