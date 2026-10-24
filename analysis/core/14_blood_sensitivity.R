# ============================================================
# blood analysis sensitivity: blood sample weights + blood/follow up selection IPW
# ============================================================

invisible(try(Sys.setlocale("LC_ALL", "zh_CN.UTF-8"), silent = TRUE))
.libPaths(c("D:/R/library", .libPaths()))
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(survey)
})

in_analysis <- "output/07_blood_biomarkers/analysis_blood_two_waves_fixedcohort_2026-09-19.rds"
in_cohort   <- "output/04_longitudinal_cohort/cohort_change_2026-09-13.rds"
in_waves    <- "output/02_five_wave_data_and_quality_control/waves_cognition_2026-09-13.rds"
out_dir     <- "output/07_blood_biomarkers"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

norm_id <- function(x) {
  x <- as.character(x)
  ifelse(nchar(x) == 11, paste0(substr(x, 1, 9), sprintf("%03d", as.integer(substr(x, 10, 11)))), x)
}

a <- readRDS(in_analysis)
cohort_obj <- readRDS(in_cohort)
cohort_long <- cohort_obj$cohort %>% mutate(KEY = norm_id(ifelse(!is.na(IDn), IDn, ID)))
base_cov <- cohort_long %>%
  filter(year == 2011) %>%
  transmute(
    KEY,
    rural = as.numeric(rural),
    communityID = as.character(communityID),
    adl = as.numeric(adl),
    cesd10 = as.numeric(cesd10),
    hibp = as.numeric(hibp),
    diab = as.numeric(diab),
    heart = as.numeric(heart),
    stroke = as.numeric(stroke),
    lung = as.numeric(lung),
    kidney = as.numeric(kidney),
    dyslip = as.numeric(dyslip),
    sleep_h = as.numeric(sleep_h),
    pain_any = as.numeric(pain_any)
  )

d <- a %>% left_join(base_cov, by = "KEY") %>%
  mutate(
    edu_group = factor(edu_group),
    edu_group = relevel(edu_group, ref = "Less than primary"),
    y_change_10 = y2020 - y2011,
    y_change_5 = y2020 - y2015,
    logcrp_2011 = log1p(pmax(crp_2011, 0)),
    logcrp_2015 = log1p(pmax(crp_2015, 0)),
    z_plt_2011 = as.numeric(scale(plt_2011)),
    z_wbc_2011 = as.numeric(scale(wbc_2011)),
    z_logcrp_2011 = as.numeric(scale(logcrp_2011)),
    z_dplt = as.numeric(scale(d_plt)),
    z_dwbc = as.numeric(scale(d_wbc)),
    z_dlogcrp = as.numeric(scale(d_logcrp))
  )

# ---------------- 1. selection model: fixed cohort and 2011 blood data available / all three core biomarkers complete at both waves ----------------
w2011 <- readRDS(in_waves) %>%
  filter(year == 2011, age >= 45, !is.na(mental)) %>%
  mutate(KEY = norm_id(IDn)) %>%
  distinct(KEY, .keep_all = TRUE) %>%
  mutate(
    edu_group = factor(edu_group),
    edu_group = relevel(edu_group, ref = "Less than primary")
  )

flags <- d %>%
  transmute(KEY, fixed_cohort, any_2011, core3_2011, core3_both,
            selected_2011 = fixed_cohort & any_2011,
            selected_both = fixed_cohort & core3_both)

seldat <- w2011 %>%
  left_join(flags, by = "KEY") %>%
  mutate(
    selected_2011 = coalesce(selected_2011, FALSE),
    selected_both = coalesce(selected_both, FALSE)
  ) %>%
  filter(!is.na(age), !is.na(female), !is.na(edu_group), !is.na(mental)) 

fit_selection <- function(outcome) {
  if (outcome == "selected_2011") {
    f <- selected_2011 ~ age + female + edu_group + rural + mental
  } else {
    f <- selected_both ~ age + female + edu_group + rural + mental
  }
  dd <- seldat %>% select(all_of(all.vars(f)), rural) %>% drop_na()
  glm(f, data = dd, family = binomial())
}
sel1 <- fit_selection("selected_2011")
sel2 <- fit_selection("selected_both")

seldat <- seldat %>%
  mutate(
    p_sel_2011 = predict(sel1, newdata = seldat, type = "response"),
    p_sel_both = predict(sel2, newdata = seldat, type = "response"),
    p_sel_2011 = pmin(pmax(p_sel_2011, quantile(p_sel_2011[selected_2011], .01, na.rm = TRUE)), quantile(p_sel_2011[selected_2011], .99, na.rm = TRUE)),
    p_sel_both = pmin(pmax(p_sel_both, quantile(p_sel_both[selected_both], .01, na.rm = TRUE)), quantile(p_sel_both[selected_both], .99, na.rm = TRUE)),
    ipw_2011 = ifelse(selected_2011, 1 / p_sel_2011, NA_real_),
    ipw_both = ifelse(selected_both, 1 / p_sel_both, NA_real_)
  )

d <- d %>% left_join(seldat %>% select(KEY, p_sel_2011, p_sel_both, ipw_2011, ipw_both), by = "KEY")

# selection bias description: SMD
smd_fun <- function(x, g) {
  ok <- !is.na(x) & !is.na(g)
  x <- x[ok]; g <- g[ok]
  m1 <- mean(x[g]); m0 <- mean(x[!g])
  s1 <- sd(x[g]); s0 <- sd(x[!g])
  (m1 - m0) / sqrt((s1^2 + s0^2) / 2)
}
sel_smd <- tibble(
  variable = c("age", "female", "rural", "mental", "depression", "ADL", "pain", "hypertension"),
  label = c("age", "female", "rural", "baseline cognition", "depression", "ADL_limited", "pain", "hypertension"),
  smd_2011_blood = c(
    smd_fun(w2011$age, w2011$KEY %in% d$KEY[d$any_2011]),
    smd_fun(w2011$female, w2011$KEY %in% d$KEY[d$any_2011]),
    smd_fun(w2011$rural, w2011$KEY %in% d$KEY[d$any_2011]),
    smd_fun(w2011$mental, w2011$KEY %in% d$KEY[d$any_2011]),
    smd_fun(w2011$cesd10, w2011$KEY %in% d$KEY[d$any_2011]),
    smd_fun(w2011$adl, w2011$KEY %in% d$KEY[d$any_2011]),
    smd_fun(w2011$pain_any, w2011$KEY %in% d$KEY[d$any_2011]),
    smd_fun(w2011$hibp, w2011$KEY %in% d$KEY[d$any_2011])
  ),
  smd_both_core = c(
    smd_fun(w2011$age, w2011$KEY %in% d$KEY[d$core3_both]),
    smd_fun(w2011$female, w2011$KEY %in% d$KEY[d$core3_both]),
    smd_fun(w2011$rural, w2011$KEY %in% d$KEY[d$core3_both]),
    smd_fun(w2011$mental, w2011$KEY %in% d$KEY[d$core3_both]),
    smd_fun(w2011$cesd10, w2011$KEY %in% d$KEY[d$core3_both]),
    smd_fun(w2011$adl, w2011$KEY %in% d$KEY[d$core3_both]),
    smd_fun(w2011$pain_any, w2011$KEY %in% d$KEY[d$core3_both]),
    smd_fun(w2011$hibp, w2011$KEY %in% d$KEY[d$core3_both])
  )
)

# ---------------- 2. weight sensitivity model ----------------
base_covars <- c("age2011", "female", "edu_group", "rural", "y2011")

extract_svy_glm <- function(fit, outcome_label, model_label, exposure_label) {
  co <- coef(summary(fit))
  ci <- confint(fit)
  term <- attr(terms(fit), "term.labels")[1]
  keep <- grepl(paste0("^", term), rownames(co))
  data.frame(
    outcome = outcome_label,
    model = model_label,
    exposure = exposure_label,
    term = rownames(co)[keep],
    n = nrow(fit$model),
    estimate = co[keep, 1],
    se = co[keep, 2],
    p = co[keep, 4],
    OR = if (fit$family$family %in% c("quasibinomial", "binomial")) exp(co[keep, 1]) else NA_real_,
    lower = if (fit$family$family %in% c("quasibinomial", "binomial")) exp(ci[keep, 1]) else ci[keep, 1],
    upper = if (fit$family$family %in% c("quasibinomial", "binomial")) exp(ci[keep, 2]) else ci[keep, 2],
    beta = if (fit$family$family %in% c("quasibinomial", "binomial")) NA_real_ else co[keep, 1],
    stringsAsFactors = FALSE
  )
}

fit_weighted <- function(data, exposure, covars, idvar, weightvar, outcome, family = c("binomial", "gaussian")) {
  family <- match.arg(family)
  cadjust <- switch(family, binomial = quasibinomial(), gaussian = gaussian())
  vars <- unique(c(outcome, exposure, covars, idvar, weightvar, "impaired_2011"))
  dd <- data %>% select(all_of(vars)) %>% drop_na()
  if (family == "binomial") dd <- dd %>% filter(impaired_2011 == 0)
  f <- as.formula(paste(outcome, "~", paste(c(exposure, covars), collapse = " + ")))
  des <- svydesign(ids = as.formula(paste0("~", idvar)), weights = as.formula(paste0("~", weightvar)),
                   data = dd, nest = TRUE)
  svyglm(f, design = des, family = cadjust)
}

weighted_rows <- list()
idx <- 0
for (outcome_info in list(
  list(outcome = "incident_2020", label = "Incident cognitive impairment 2020", family = "binomial"),
  list(outcome = "y_change_10", label = "Cognitive change 2011-2020", family = "gaussian")
)) {
  for (i in seq_along(c("z_plt_2011", "z_wbc_2011", "z_logcrp_2011"))) {
    ex <- c("z_plt_2011", "z_wbc_2011", "z_logcrp_2011")[i]
    lab <- c("PLT per 1 SD", "WBC per 1 SD", "log(CRP+1) per 1 SD")[i]
    for (w in c("bloodweight_2011", "ipw_2011", "w_comb_2011")) {
      if (w == "w_comb_2011") {
        d$w_comb_2011 <- d$bloodweight_2011 * d$ipw_2011
        d$w_comb_2011 <- pmin(pmax(d$w_comb_2011, quantile(d$w_comb_2011, .01, na.rm = TRUE)), quantile(d$w_comb_2011, .99, na.rm = TRUE))
      }
      m <- tryCatch(fit_weighted(d, ex, base_covars, "communityID", w, outcome_info$outcome, outcome_info$family), error = function(e) NULL)
      if (!is.null(m)) {
        idx <- idx + 1
        weighted_rows[[idx]] <- extract_svy_glm(m, outcome_info$label, paste0(w, " base-adjusted"), lab)
      }
    }
  }
}

# blood biomarker change model: 2015 weight/selection IPW
d <- d %>% mutate(w_comb_2015 = pmin(pmax(bloodweight_2015 * ipw_both, quantile(bloodweight_2015 * ipw_both, .01, na.rm = TRUE)), quantile(bloodweight_2015 * ipw_both, .99, na.rm = TRUE)))
for (i in seq_along(c("z_dplt", "z_dwbc", "z_dlogcrp"))) {
  ex <- c("z_dplt", "z_dwbc", "z_dlogcrp")[i]
  lab <- c("delta PLT per 1 SD", "delta WBC per 1 SD", "delta log(CRP+1) per 1 SD")[i]
  for (w in c("bloodweight_2015", "ipw_both", "w_comb_2015")) {
    m <- tryCatch(fit_weighted(d, ex, c("age2011", "female", "edu_group", "rural", "y2015"), "communityID_2015", w, "y_change_5", "gaussian"), error = function(e) NULL)
    if (!is.null(m)) {
      idx <- idx + 1
      weighted_rows[[idx]] <- extract_svy_glm(m, "Cognitive change 2015-2020", paste0(w, " base-adjusted"), lab)
    }
  }
}
weighted <- bind_rows(weighted_rows)
weighted$q <- NA_real_
for (key in unique(weighted[, c("outcome", "model")])) {}
# only after complete filling apply outcome/model apply to three continuous exposures BH-FDR
out <- weighted; out$q <- NA_real_
keys <- unique(out[, c("outcome", "model")])
for (i in seq_len(nrow(keys))) {
  idx2 <- which(out$outcome == keys$outcome[i] & out$model == keys$model[i])
  if (length(idx2)) out$q[idx2] <- p.adjust(out$p[idx2], method = "BH")
}
weighted <- out

write.csv(sel_smd, file.path(out_dir, "blood_selection_smd_2026-09-19.csv"), row.names = FALSE)
write.csv(weighted, file.path(out_dir, "blood_weighted_sensitivity_2026-09-19.csv"), row.names = FALSE)

# supplementary table: retain exposure rows only
supp <- bind_rows(
  read.csv(file.path(out_dir, "blood_main_continuous_2026-09-19.csv")) %>% mutate(Table="Main continuous"),
  read.csv(file.path(out_dir, "blood_main_quartile_2026-09-19.csv")) %>% mutate(Table="Main quartile"),
  read.csv(file.path(out_dir, "blood_change_models_2026-09-19.csv")) %>% mutate(Table="Change exposure"),
  read.csv(file.path(out_dir, "blood_winsor_sensitivity_2026-09-19.csv")) %>% mutate(Table="Winsorized"),
  weighted %>% mutate(Table="Weighted/IPW")
)
write.csv(supp, file.path(out_dir, "blood_supplementary_all_models_2026-09-19.csv"), row.names = FALSE)

cat("complete: blood_weighted_sensitivity_2026-09-19.csv\n")
print(weighted %>% select(outcome, model, exposure, OR, beta, lower, upper, p, q, n))
