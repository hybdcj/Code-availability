# ============================================================================
# 08_measurement_equivalence.R
# comparability of cognitive measurement across five waves/measurement equivalence sensitivity analysis
# version: v1.1(2026-09-13, with 03_waves_build.R of 2018 updated alongside wave variable correction)
# input: output/02_five_wave_data_and_quality_control/waves_cognition_2026-09-13.rds + raw cognitive modules across five waves
# output: output/05_measurement_equivalence_ratio/measurement_*_2026-09-13.csv / .rds
# purpose: 
#   ① decompose“orientation 5 questions + serial subtraction of 7(0-10)”item missingness and distribution changes in common item set; 
#   ② using multigroup ordinal CFA test 2011-2020 five wave measurement invariance, focus on adjacent wave comparisons; 
#   ③ compare memory scores A(prorated), B(within wave z), C(exclude 2018)three schemes; 
#   ④ at 2018 aged at least 60 years old, and completed 2018 rerun ten year change and risk factors among participants with complete cognitive modules
#      (sensitivity scenarios, for testing"2018 year 60 all items administered only above age"whether conclusions are affected).
# note: 2018 wave memory scores now use full sample assessments 10 word list(DC028_W4 immediate recall trial 1 + DC047_W4 delayed, 0-20), 
#     no longer use items administered only to 60 administered above age 3 word list(DC013_W4_1 / DC015_W4).
# ============================================================================
invisible(try(Sys.setlocale("LC_ALL", "zh_CN.UTF-8"), silent = TRUE))
if (dir.exists("D:/R/library")) .libPaths(c("D:/R/library", .libPaths()))
suppressPackageStartupMessages({
  library(haven); library(dplyr); library(tidyr); library(nlme)
  library(survival); library(MASS); library(lavaan)
})
select <- dplyr::select  # MASS would mask dplyr::select

PROJ <- Sys.getenv("CHARLS_PROJ", getwd())
RAW <- file.path(PROJ, "data", "raw")
OUT <- file.path(PROJ, "output", "05_measurement_equivalence_ratio")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
num <- function(x) as.numeric(x)
norm_id <- function(x) {
  x <- as.character(x); n <- nchar(x)
  ifelse(n == 11, paste0(substr(x, 1, 9), sprintf("%03d", as.integer(substr(x, 10, 11)))), x)
}
fmt_p <- function(p) ifelse(is.na(p), "NA", ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
fmt_est <- function(est, lo, hi, p) sprintf("%.3f (%.3f-%.3f), p=%s", est, lo, hi, fmt_p(p))
paired_ci <- function(x) {
  x <- x[!is.na(x)]; n <- length(x)
  if (n < 2) return(c(n = n, mean = NA, lo = NA, hi = NA, sd = NA))
  m <- mean(x); se <- sd(x) / sqrt(n)
  c(n = n, mean = m, lo = m - 1.96 * se, hi = m + 1.96 * se, sd = sd(x))
}

## ---------------------------------------------------------------------------
## 1. read constructed five wave data and fixed common item cohort
## ---------------------------------------------------------------------------
w <- readRDS(file.path(PROJ, "output", "02_five_wave_data_and_quality_control", "waves_cognition_2026-09-13.rds"))
w <- w %>% mutate(KEY = ifelse(is.na(IDn), ID, IDn), ms10 = orient + ser7)
candidate_ids <- w %>% filter(!is.na(ms10)) %>% count(KEY) %>% filter(n == 5) %>% pull(KEY)
baseline_age_ids <- w %>% filter(year == 2011, !is.na(age), age >= 45) %>% pull(KEY) %>% unique()
cohort_ids <- intersect(candidate_ids, baseline_age_ids)
coh <- w %>% filter(KEY %in% cohort_ids) %>% arrange(KEY, year)
cat("2011 year≥45 age eligible fixed cohort with scorable common items across five waves n =", length(cohort_ids), "\n")

## within wave overall z used only for B scheme description; memory_z retain if already present, recompute to ensure consistent cohort definition.
w <- w %>% group_by(year) %>% mutate(memory_z_pop = as.numeric(scale(memory_prop))) %>% ungroup()

## ---------------------------------------------------------------------------
## 2. extract from raw modules“orientation 5 questions + serial subtraction of 7 five steps”item level responses
## ---------------------------------------------------------------------------
read_orient_items <- function(year, file, vars, ok) {
  d <- read_dta(file, col_select = all_of(c("ID", vars)))
  out <- data.frame(KEY = norm_id(d$ID), year = year, check.names = FALSE)
  item_names <- c("ori_year","ori_month","ori_date","ori_day","ori_season")
  for (i in seq_along(vars)) {
    x <- num(d[[vars[i]]])
    x[x %in% c(999, -1)] <- NA
    out[[item_names[i]]] <- as.integer(!is.na(x) & x == ok[i])  # consistent with original scoring: not answered/do not know=incorrect
  }
  out
}
read_s7_items <- function(year, file, ans_vars, flag_vars = NULL) {
  cols <- c("ID", ans_vars, flag_vars)
  d <- read_dta(file, col_select = all_of(cols))
  a <- sapply(ans_vars, function(v) num(d[[v]]))
  a <- matrix(a, nrow = nrow(d), ncol = length(ans_vars))
  a[!is.na(a) & (a < 0 | a > 1000)] <- NA
  eng <- !is.na(a)
  if (!is.null(flag_vars)) {
    f <- sapply(flag_vars, function(v) num(d[[v]]))
    f <- matrix(f, nrow = nrow(d), ncol = length(flag_vars))
    eng <- eng | (!is.na(f) & f == 997)
  }
  ok <- matrix(FALSE, nrow = nrow(a), ncol = ncol(a))
  ok[, 1] <- !is.na(a[, 1]) & a[, 1] == 93
  if (ncol(a) >= 2) for (j in 2:ncol(a)) {
    ok[, j] <- !is.na(a[, j - 1]) & !is.na(a[, j]) & a[, j] == a[, j - 1] - 7
  }
  ok <- ok & eng
  out <- data.frame(KEY = norm_id(d$ID), year = year, check.names = FALSE)
  for (j in seq_len(ncol(ok))) out[[paste0("s7_", j)]] <- ifelse(rowSums(!is.na(a)) == 0, NA_integer_, as.integer(ok[, j]))
  out
}

ori_list <- list(
  read_orient_items(2011, file.path(RAW, "wave1_2011", "health_status_and_functioning.dta"),
                    c("dc001s1","dc001s2","dc001s3","dc002","dc003"), c(1,2,3,1,1)),
  read_orient_items(2013, file.path(RAW, "wave2_2013", "Health_Status_and_Functioning.dta"),
                    c("dc001s1","dc001s2","dc001s3","dc002","dc003"), c(1,2,3,1,1)),
  read_orient_items(2015, file.path(RAW, "wave3_2015", "Health_Status_and_Functioning.dta"),
                    c("dc001s1","dc001s2","dc001s3","dc002","dc003"), c(1,2,3,1,1)),
  read_orient_items(2018, file.path(RAW, "wave4_2018", "Cognition.dta"),
                    c("dc001_w4","dc006_w4","dc003_w4","dc005_w4","dc002_w4"), rep(1, 5)),
  read_orient_items(2020, file.path(RAW, "wave5_2020", "Health_Status_and_Functioning.dta"),
                    c("dc001","dc005","dc003","dc004","dc002"), rep(1, 5))
)
s7_list <- list(
  read_s7_items(2011, file.path(RAW, "wave1_2011", "health_status_and_functioning.dta"),
                paste0("dc0", 19:23)),
  read_s7_items(2013, file.path(RAW, "wave2_2013", "Health_Status_and_Functioning.dta"),
                paste0("dc0", 19:23)),
  read_s7_items(2015, file.path(RAW, "wave3_2015", "Health_Status_and_Functioning.dta"),
                paste0("dc0", 19:23)),
  read_s7_items(2018, file.path(RAW, "wave4_2018", "Cognition.dta"),
                paste0("dc014_w4_", 1:5, "_1")),
  read_s7_items(2020, file.path(RAW, "wave5_2020", "Health_Status_and_Functioning.dta"),
                paste0("dc007_", 1:5, "_1"), paste0("dc007_", 1:5))
)
items_raw <- full_join(bind_rows(ori_list), bind_rows(s7_list), by = c("KEY", "year")) %>%
  arrange(KEY, year)
common_items <- c("ori_year","ori_month","ori_date","ori_day","ori_season","s7_1","s7_2","s7_3","s7_4","s7_5")
coh_i <- coh %>% select(KEY, year, age, female, ms10, orient, ser7) %>%
  left_join(items_raw, by = c("KEY", "year")) %>%
  mutate(item_n = rowSums(!is.na(across(all_of(common_items)))),
         item_complete = as.integer(item_n == length(common_items)))

item_by_wave <- do.call(rbind, lapply(sort(unique(coh_i$year)), function(y) {
  d <- coh_i %>% filter(year == y)
  ori_mat <- as.matrix(d[, c("ori_year","ori_month","ori_date","ori_day","ori_season")])
  s7_mat <- as.matrix(d[, paste0("s7_", 1:5)])
  data.frame(
    year = y, n = nrow(d), item_complete = sum(d$item_complete == 1, na.rm = TRUE),
    complete_pct = sprintf("%.1f%%", 100 * mean(d$item_complete == 1, na.rm = TRUE)),
    orient_total_available = sprintf("%.1f%%", 100 * mean(!is.na(d$orient))),
    ser7_total_available = sprintf("%.1f%%", 100 * mean(!is.na(d$ser7))),
    mean_ms10 = round(mean(d$ms10, na.rm = TRUE), 3),
    orient_item_sum_mean = round(mean(rowSums(ori_mat, na.rm = TRUE), na.rm = TRUE), 3),
    s7_item_sum_mean = round(mean(rowSums(s7_mat, na.rm = TRUE), na.rm = TRUE), 3)
  )
}))
write.csv(item_by_wave, file.path(OUT, "measurement_item_completeness_2026-09-13.csv"), row.names = FALSE, fileEncoding = "UTF-8")

## item correct response rate(fixed cohort)
item_prev <- do.call(rbind, lapply(sort(unique(coh_i$year)), function(y) {
  d <- coh_i %>% filter(year == y)
  data.frame(year = y, item = common_items,
             p_correct = round(colMeans(d[, common_items], na.rm = TRUE), 3),
             n_available = colSums(!is.na(d[, common_items])), row.names = NULL)
}))
write.csv(item_prev, file.path(OUT, "measurement_item_prevalence_2026-09-13.csv"), row.names = FALSE, fileEncoding = "UTF-8")

## item level 2015 -> 2018 paired change(participants responding at both waves only)
item_change <- do.call(rbind, lapply(common_items, function(v) {
  d <- coh_i %>% filter(year %in% c(2015, 2018), !is.na(.data[[v]])) %>%
    select(KEY, year, val = all_of(v)) %>% pivot_wider(names_from = year, values_from = val, names_prefix = "y")
  if (nrow(d) < 2) return(NULL)
  delta <- d$y2018 - d$y2015
  data.frame(item = v, n_pair = nrow(d), change = round(mean(delta), 4),
             lo = round(mean(delta) - 1.96 * sd(delta) / sqrt(nrow(d)), 4),
             hi = round(mean(delta) + 1.96 * sd(delta) / sqrt(nrow(d)), 4),
             p = fmt_p(t.test(delta)$p.value))
}))
write.csv(item_change, file.path(OUT, "measurement_item_change_2015_2018_2026-09-13.csv"), row.names = FALSE, fileEncoding = "UTF-8")

## ---------------------------------------------------------------------------
## 3. multigroup ordinal CFA: measurement invariance of common items
## ---------------------------------------------------------------------------
all_complete_ids <- coh_i %>% group_by(KEY) %>% summarise(k = sum(item_complete == 1), .groups = "drop") %>%
  filter(k == 5) %>% pull(KEY)
moddat <- coh_i %>% filter(KEY %in% all_complete_ids, item_complete == 1) %>%
  mutate(wave = factor(year, levels = c(2011, 2013, 2015, 2018, 2020))) %>%
  droplevels()
cat("complete item responses across five waves, available for all wave CFA fixed cohort n =", length(all_complete_ids), "\n")

cfa_model <- paste("F1 =~", paste(common_items[1:5], collapse = " + "), "\nF2 =~", paste(common_items[6:10], collapse = " + "))
fit_measure <- function(fit) {
  m <- tryCatch(lavaan::fitMeasures(fit), error = function(e) NULL)
  if (is.null(m)) return(NULL)
  nm <- function(x) if (x %in% names(m)) as.numeric(m[[x]]) else NA_real_
  data.frame(chisq = nm("chisq.scaled"), df = nm("df.scaled"), p = nm("pvalue.scaled"),
             cfi = nm("cfi.scaled"), rmsea = nm("rmsea.scaled"), srmr = nm("srmr"))
}
fit_pair <- function(y1, y2) {
  d <- coh_i %>% filter(year %in% c(y1, y2), item_complete == 1) %>%
    mutate(wave = factor(year, levels = c(y1, y2)))
  if (nrow(d) < 100 || length(unique(d$wave)) < 2) return(NULL)
  fit0 <- tryCatch(cfa(cfa_model, data = d, group = "wave", ordered = common_items,
                       estimator = "WLSMV", std.lv = TRUE, missing = "listwise", meanstructure = TRUE), error = function(e) e)
  fit1 <- tryCatch(cfa(cfa_model, data = d, group = "wave", ordered = common_items,
                       estimator = "WLSMV", std.lv = TRUE, missing = "listwise", meanstructure = TRUE,
                       group.equal = c("loadings", "thresholds")), error = function(e) e)
  if (inherits(fit0, "error") || inherits(fit1, "error")) return(NULL)
  f0 <- fit_measure(fit0); f1 <- fit_measure(fit1)
  data.frame(pair = paste(y1, y2, sep = "-"), n = nrow(d),
             config_cfi = f0$cfi, scalar_cfi = f1$cfi, delta_cfi = f1$cfi - f0$cfi,
             config_rmsea = f0$rmsea, scalar_rmsea = f1$rmsea, delta_rmsea = f1$rmsea - f0$rmsea,
             scalar_chisq = f1$chisq, scalar_df = f1$df, scalar_p = f1$p)
}
pairs <- list(c(2011, 2013), c(2013, 2015), c(2015, 2018), c(2018, 2020), c(2011, 2020))
cfa_pairs <- bind_rows(lapply(pairs, function(z) fit_pair(z[1], z[2])))
if (nrow(cfa_pairs)) {
  cfa_pairs$delta_cfi <- round(cfa_pairs$delta_cfi, 4)
  cfa_pairs$delta_rmsea <- round(cfa_pairs$delta_rmsea, 4)
  write.csv(cfa_pairs, file.path(OUT, "measurement_cfa_pairs_2026-09-13.csv"), row.names = FALSE, fileEncoding = "UTF-8")
}

fit0_all <- tryCatch(cfa(cfa_model, data = moddat, group = "wave", ordered = common_items,
                         estimator = "WLSMV", std.lv = TRUE, missing = "listwise", meanstructure = TRUE), error = function(e) e)
fit1_all <- tryCatch(cfa(cfa_model, data = moddat, group = "wave", ordered = common_items,
                         estimator = "WLSMV", std.lv = TRUE, missing = "listwise", meanstructure = TRUE,
                         group.equal = c("loadings", "thresholds")), error = function(e) e)
if (!inherits(fit0_all, "error") && !inherits(fit1_all, "error")) {
  cfa_all <- rbind(cbind(model = "Configural", fit_measure(fit0_all)),
                   cbind(model = "Scalar", fit_measure(fit1_all)))
  write.csv(cfa_all, file.path(OUT, "measurement_cfa_allwaves_2026-09-13.csv"), row.names = FALSE, fileEncoding = "UTF-8")
}

## ---------------------------------------------------------------------------
## 4. 2018 year 60 age restriction / complete 2018 complete cognitive module sensitivity analysis
## ---------------------------------------------------------------------------
ids60 <- coh %>% filter(year == 2018, age >= 60) %>% pull(KEY)
ids60_mem <- coh %>% filter(year == 2018, age >= 60, !is.na(imrc), !is.na(dlrc)) %>% pull(KEY)
ids60_full <- coh_i %>% filter(year == 2018, age >= 60, item_complete == 1) %>% pull(KEY)
scenarios <- list(
  all = cohort_ids,
  age60_2018 = ids60,
  age60_2018_memory = ids60_mem,
  age60_2018_common_complete = ids60_full
)
scen_desc <- c(all = "entire five wave common item cohort", age60_2018 = "2018 aged at least 60 years old",
               age60_2018_memory = "2018 aged at least 60 years old with complete memory items",
               age60_2018_common_complete = "2018 aged at least 60 age eligible with complete common items")

scen_tab <- do.call(rbind, lapply(names(scenarios), function(s) {
  d <- coh %>% filter(KEY %in% scenarios[[s]])
  wide <- d %>% select(KEY, year, ms10) %>% pivot_wider(names_from = year, values_from = ms10, names_prefix = "y")
  ch1518 <- paired_ci(wide$y2018 - wide$y2015)
  ch1120 <- paired_ci(wide$y2020 - wide$y2011)
  means <- d %>% group_by(year) %>% summarise(m = mean(ms10, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = year, values_from = m, names_prefix = "m")
  data.frame(
    scenario = s, label = scen_desc[[s]], n = length(scenarios[[s]]),
    m2011 = round(means$m2011, 3), m2013 = round(means$m2013, 3), m2015 = round(means$m2015, 3),
    m2018 = round(means$m2018, 3), m2020 = round(means$m2020, 3),
    change_2015_2018 = round(ch1518["mean"], 3),
    ci_2015_2018 = sprintf("%.3f to %.3f", ch1518["lo"], ch1518["hi"]),
    change_2011_2020 = round(ch1120["mean"], 3),
    ci_2011_2020 = sprintf("%.3f to %.3f", ch1120["lo"], ch1120["hi"])
  )
}))
write.csv(scen_tab, file.path(OUT, "measurement_common_score_sensitivity_2026-09-13.csv"), row.names = FALSE, fileEncoding = "UTF-8")

## risk factors: rerun depression separately under four sample definitions, ADL -> incident impairment / first impairment
run_risk <- function(d, label) {
  d <- d %>% mutate(dep = as.integer(cesd10 >= 10), edu = factor(edu_group))
  base0 <- d %>% filter(year == 2011)
  thr <- base0 %>% group_by(edu) %>%
    summarise(cut_sub = mean(ms10, na.rm = TRUE) - sd(ms10, na.rm = TRUE), .groups = "drop")
  cutoff <- setNames(thr$cut_sub, as.character(thr$edu))
  d <- d %>% mutate(impaired_sub = as.integer(ms10 < cutoff[as.character(edu)]))
  base <- d %>% filter(year == 2011) %>%
    select(KEY, age, female, edu, dep, adl, ms10, impaired_2011 = impaired_sub)
  inc <- base %>% filter(impaired_2011 == 0) %>%
    left_join(d %>% filter(year == 2020) %>% select(KEY, impaired_2020 = impaired_sub), by = "KEY")
  fu <- d %>% filter(year != 2011) %>% select(KEY, year, impaired_sub) %>%
    mutate(t = c("2013"=2, "2015"=4, "2018"=7, "2020"=9)[as.character(year)]) %>% arrange(KEY, t)
  surv_events <- fu %>% group_by(KEY) %>% summarise(
    time = if (any(impaired_sub == 1, na.rm = TRUE)) min(t[impaired_sub == 1], na.rm = TRUE) else 9,
    event = as.integer(any(impaired_sub == 1, na.rm = TRUE)), .groups = "drop")
  surv_df <- base %>% filter(impaired_2011 == 0) %>% select(KEY, age, female, edu, dep, adl) %>%
    left_join(surv_events, by = "KEY") %>%
    mutate(time = ifelse(is.na(time), 9, time), event = ifelse(is.na(event), 0, event))
  delta <- d %>% select(KEY, year, ms10) %>% pivot_wider(names_from = year, values_from = ms10, names_prefix = "y") %>%
    transmute(KEY, delta = y2020 - y2011)
  base <- base %>% left_join(delta, by = "KEY")
  out <- list()
  for (e in c("dep","adl")) {
    f <- as.formula(paste("impaired_2020 ~", e, "+ age + female + edu"))
    m1 <- tryCatch(glm(f, data = inc, family = binomial), error = function(e) NULL)
    m2 <- tryCatch(coxph(as.formula(paste("Surv(time, event) ~", e, "+ age + female + edu")), data = surv_df), error = function(e) NULL)
    m3 <- tryCatch(lm(as.formula(paste("delta ~", e, "+ age + female + edu")), data = base), error = function(e) NULL)
    if (!is.null(m1)) {
      s <- summary(m1)$coefficients[e, ]; or <- exp(s[1]); lo <- exp(s[1]-1.96*s[2]); hi <- exp(s[1]+1.96*s[2]); op <- s[4]
    } else or <- lo <- hi <- op <- NA
    if (!is.null(m2)) {
      s2 <- summary(m2)$coefficients[e, ]; c2 <- summary(m2)$conf.int[e, ]; hr <- c2[1]; hlo <- c2[3]; hhi <- c2[4]; hp <- s2[5]
    } else hr <- hlo <- hhi <- hp <- NA
    if (!is.null(m3)) {
      s3 <- summary(m3)$coefficients[e, ]; b <- s3[1]; blo <- b-1.96*s3[2]; bhi <- b+1.96*s3[2]; bp <- s3[4]
    } else b <- blo <- bhi <- bp <- NA
    out[[e]] <- data.frame(scenario = label, exposure = e,
                           n_inc = nrow(inc), events = sum(inc$impaired_2020, na.rm = TRUE),
                           OR = fmt_est(or, lo, hi, op), HR = fmt_est(hr, hlo, hhi, hp),
                           delta_beta = fmt_est(b, blo, bhi, bp), row.names = NULL)
  }
  bind_rows(out)
}

risk_sens <- bind_rows(lapply(names(scenarios), function(s) run_risk(coh %>% filter(KEY %in% scenarios[[s]]), scen_desc[[s]])))
write.csv(risk_sens, file.path(OUT, "measurement_common_score_risk_sensitivity_2026-09-13.csv"), row.names = FALSE, fileEncoding = "UTF-8")
## ---------------------------------------------------------------------------
## 5. three memory score definitions: A prorated, B within wave z, C exclude 2018
## ---------------------------------------------------------------------------
mem_ids <- w %>% filter(!is.na(memory_prop)) %>% count(KEY) %>% filter(n == 5) %>% pull(KEY)
mem <- w %>% filter(KEY %in% mem_ids) %>% arrange(KEY, year)
cat("fixed cohort with scorable memory at all five waves n =", length(mem_ids), "\n")

mem_desc <- mem %>% group_by(year) %>% summarise(
  n = n(), list_len = first(list_len_imm + list_len_dlr),
  raw_mean = mean(memory_raw, na.rm = TRUE), prop_mean = mean(memory_prop, na.rm = TRUE),
  prop_sd = sd(memory_prop, na.rm = TRUE), floor_pct = 100*mean(memory_prop == 0, na.rm = TRUE),
  ceiling_pct = 100*mean(memory_prop == 1, na.rm = TRUE),
  imm_delay_r = cor(imrc/list_len_imm, dlrc/list_len_dlr, method = "spearman", use = "complete.obs"),
  mem_ms10_r = cor(memory_prop, ms10, method = "spearman", use = "complete.obs"),
  .groups = "drop")
mem_desc$alpha_2item <- 2 * mem_desc$imm_delay_r / (1 + mem_desc$imm_delay_r)
write.csv(mem_desc, file.path(OUT, "measurement_memory_wave_descriptives_2026-09-13.csv"), row.names = FALSE, fileEncoding = "UTF-8")

base_mem <- mem %>% filter(year == 2011) %>%
  transmute(KEY, age2011 = age, female, edu_group, dep = as.integer(cesd10 >= 10), adl)
mem_long <- mem %>% transmute(KEY, year, prop = memory_prop, zp = memory_z_pop) %>%
  left_join(base_mem, by = "KEY") %>%
  mutate(t = c("2011"=0, "2013"=2, "2015"=4, "2018"=7, "2020"=9)[as.character(year)]) %>%
  filter(!is.na(edu_group))

fit_mem_model <- function(dat, score, scheme, excl2018 = FALSE) {
  dd <- dat
  if (excl2018) dd <- dd %>% filter(year != 2018)
  m <- tryCatch(lme(as.formula(paste(score, "~ t + age2011 + female + edu_group")),
                    random = ~ 1 + t | KEY, data = dd, na.action = na.omit), error = function(e) e)
  if (inherits(m, "error")) return(data.frame(scheme = scheme, n = length(unique(dd$KEY)), slope_year = NA, slope_decade = NA, p = NA))
  s <- summary(m)$tTable; b <- s["t",1]; se <- s["t",2]
  data.frame(scheme = scheme, n = length(unique(dd$KEY)), slope_year = b, slope_decade = 10*b,
             ci_decade = sprintf("%.3f to %.3f", 10*(b-1.96*se), 10*(b+1.96*se)), p = fmt_p(s["t",5]))
}
mem_schemes <- bind_rows(
  fit_mem_model(mem_long, "prop", "A: prorated(including 2018)", FALSE),
  fit_mem_model(mem_long, "zp", "B: within wave z(including 2018)", FALSE),
  fit_mem_model(mem_long, "prop", "C: prorated(exclude 2018)", TRUE)
)
write.csv(mem_schemes, file.path(OUT, "measurement_memory_schemes_2026-09-13.csv"), row.names = FALSE, fileEncoding = "UTF-8")

## time under two prespecified exposures and three memory definitions×exposure interaction
mem_risk <- do.call(rbind, lapply(c("dep","adl"), function(e) {
  do.call(rbind, lapply(list(
    list(score="prop", scheme="A: prorated(including 2018)", excl=FALSE),
    list(score="zp", scheme="B: within wave z(including 2018)", excl=FALSE),
    list(score="prop", scheme="C: prorated(exclude 2018)", excl=TRUE)
  ), function(z) {
    dd <- mem_long; if (z$excl) dd <- dd %>% filter(year != 2018)
    m <- tryCatch(lme(as.formula(paste(z$score, "~ t *", e, "+ age2011 + female + edu_group")),
                      random = ~ 1 + t | KEY, data = dd, na.action = na.omit), error = function(e) e)
    if (inherits(m, "error")) return(data.frame(exposure=e, scheme=z$scheme, beta_decade=NA, ci="", p=NA))
    s <- summary(m)$tTable; rn <- grep(paste0("^t:", e), rownames(s), value = TRUE)[1]
    b <- s[rn,1]; se <- s[rn,2]
    data.frame(exposure=e, scheme=z$scheme, beta_decade=10*b,
               ci=sprintf("%.3f to %.3f", 10*(b-1.96*se), 10*(b+1.96*se)), p=fmt_p(s[rn,5]))
  }))
}))
write.csv(mem_risk, file.path(OUT, "measurement_memory_risk_schemes_2026-09-13.csv"), row.names = FALSE, fileEncoding = "UTF-8")

saveRDS(list(coh=coh, coh_item=coh_i, item_by_wave=item_by_wave, item_prev=item_prev,
             item_change=item_change, cfa_pairs=cfa_pairs, cfa_all=if (exists("cfa_all")) cfa_all else NULL,
             scenarios=scen_tab, risk_sensitivity=risk_sens, memory_desc=mem_desc,
             memory_schemes=mem_schemes, memory_risk=mem_risk),
        file.path(OUT, "measurement_equivalence_2026-09-13.rds"))

cat("\n===== common item completeness =====\n"); print(item_by_wave, row.names = FALSE)
cat("\n===== 2015->2018 item specific change =====\n"); print(item_change, row.names = FALSE)
cat("\n===== sample sensitivity =====\n"); print(scen_tab, row.names = FALSE)
cat("\n===== CFA adjacent waves =====\n"); print(cfa_pairs, row.names = FALSE)
cat("\n===== memory score description =====\n"); print(mem_desc, row.names = FALSE)
cat("\n===== memory score time slope =====\n"); print(mem_schemes, row.names = FALSE)
cat("\n outputs generated measurement_*_2026-09-13.csv/.rds\n")
