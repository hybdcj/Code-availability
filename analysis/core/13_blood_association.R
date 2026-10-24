# ============================================================
# blood analysis step three: PLT/WBC/CRP and ten year cognitive impairment/association with change
# primary analysis: 
#   1) 2011 PLT/WBC/logCRP prediction 2020 incident cognitive impairment
#   2) 2011 PLT/WBC/logCRP and 2011→2020 cognitive change
#   3) 2011→2015 blood biomarker change predicts 2015→2020 cognitive change(exploration)
# notes: 
#   - continuous exposures scaled by 1 SD standardization; 
#   - CRP use log(CRP+1), quartiles still ranked by raw values; 
#   - primarily unweighted, blood weights reserved for later sensitivity analyses; 
#   - extreme values handled using 1%/99% winsorization sensitivity analysis.
# ============================================================

invisible(try(Sys.setlocale("LC_ALL", "zh_CN.UTF-8"), silent = TRUE))
.libPaths(c("D:/R/library", .libPaths()))
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(openxlsx)
})

in_analysis <- "output/07_blood_biomarkers/analysis_blood_two_waves_fixedcohort_2026-09-19.rds"
in_cohort   <- "output/04_longitudinal_cohort/cohort_change_2026-09-13.rds"
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
    pain_any = as.numeric(pain_any),
    communityID = as.character(communityID)
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
    z_dlogcrp = as.numeric(scale(d_logcrp)),
    q_plt_2011 = ntile(plt_2011, 4),
    q_wbc_2011 = ntile(wbc_2011, 4),
    q_crp_2011 = ntile(crp_2011, 4),
    q_dplt = ntile(d_plt, 4),
    q_dwbc = ntile(d_wbc, 4),
    q_dcrp = ntile(d_logcrp, 4),
    q_plt_2011 = factor(q_plt_2011, levels = 1:4, labels = paste0("Q", 1:4)),
    q_wbc_2011 = factor(q_wbc_2011, levels = 1:4, labels = paste0("Q", 1:4)),
    q_crp_2011 = factor(q_crp_2011, levels = 1:4, labels = paste0("Q", 1:4)),
    q_dplt = factor(q_dplt, levels = 1:4, labels = paste0("Q", 1:4)),
    q_dwbc = factor(q_dwbc, levels = 1:4, labels = paste0("Q", 1:4)),
    q_dcrp = factor(q_dcrp, levels = 1:4, labels = paste0("Q", 1:4))
  )

# 1%/99% winsorization sensitivity variables
wins <- function(x) {
  lo <- quantile(x, 0.01, na.rm = TRUE, names = FALSE)
  hi <- quantile(x, 0.99, na.rm = TRUE, names = FALSE)
  pmin(pmax(x, lo), hi)
}
d <- d %>%
  mutate(
    plt_2011_w = wins(plt_2011),
    wbc_2011_w = wins(wbc_2011),
    crp_2011_w = wins(crp_2011),
    logcrp_2011_w = log1p(pmax(crp_2011_w, 0)),
    z_plt_2011_w = as.numeric(scale(plt_2011_w)),
    z_wbc_2011_w = as.numeric(scale(wbc_2011_w)),
    z_logcrp_2011_w = as.numeric(scale(logcrp_2011_w))
  )

# basic model and Fully-adjusted model
base_covars <- c("age2011", "female", "edu_group", "rural", "y2011")
full_covars <- c(base_covars, "cesd10", "adl", "hibp", "diab", "heart", "stroke", "lung", "kidney", "dyslip", "sleep_h", "pain_any")

extract_glm <- function(fit, exposure_label, model_label, outcome_label) {
  co <- coef(summary(fit))
  ci <- confint.default(fit)
  keep <- grepl(paste0("^", attr(terms(fit), "term.labels")[1]), rownames(co))
  data.frame(
    outcome = outcome_label,
    model = model_label,
    exposure = exposure_label,
    term = rownames(co)[keep],
    n = length(fit$fitted.values),
    events = if (!is.null(fit$y)) sum(fit$y == 1, na.rm = TRUE) else NA_integer_,
    estimate = co[keep, 1],
    se = co[keep, 2],
    p = co[keep, 4],
    OR = exp(co[keep, 1]),
    lower = exp(ci[keep, 1]),
    upper = exp(ci[keep, 2]),
    stringsAsFactors = FALSE
  )
}

extract_lm <- function(fit, exposure_label, model_label, outcome_label) {
  co <- coef(summary(fit))
  ci <- confint.default(fit)
  keep <- grepl(paste0("^", attr(terms(fit), "term.labels")[1]), rownames(co))
  data.frame(
    outcome = outcome_label,
    model = model_label,
    exposure = exposure_label,
    term = rownames(co)[keep],
    n = length(fit$fitted.values),
    events = NA_integer_,
    estimate = co[keep, 1],
    se = co[keep, 2],
    p = co[keep, 4],
    beta = co[keep, 1],
    lower = ci[keep, 1],
    upper = ci[keep, 2],
    stringsAsFactors = FALSE
  )
}

fit_logit <- function(data, exposure, covars) {
  vars <- c("incident_2020", exposure, covars)
  dd <- data %>% select(all_of(vars)) %>% drop_na()
  # incidence analysis restricted to baseline unimpaired
  dd <- dd[data$impaired_2011[match(rownames(dd), rownames(data))] == 0, , drop = FALSE]
  f <- as.formula(paste("incident_2020 ~", paste(c(exposure, covars), collapse = " + ")))
  list(fit = glm(f, data = dd, family = binomial()), data = dd)
}

# incidence model: handle baseline impairment restriction separately
make_incident_data <- function(data, exposure, covars) {
  vars <- unique(c("incident_2020", "impaired_2011", exposure, covars))
  data %>% select(all_of(vars)) %>% filter(impaired_2011 == 0) %>% drop_na()
}
fit_logit2 <- function(data, exposure, covars) {
  dd <- make_incident_data(data, exposure, covars)
  f <- as.formula(paste("incident_2020 ~", paste(c(exposure, covars), collapse = " + ")))
  list(fit = glm(f, data = dd, family = binomial()), data = dd)
}
fit_lm2 <- function(data, outcome, exposure, covars) {
  vars <- unique(c(outcome, exposure, covars))
  dd <- data %>% select(all_of(vars)) %>% drop_na()
  f <- as.formula(paste(outcome, "~", paste(c(exposure, covars), collapse = " + ")))
  list(fit = lm(f, data = dd), data = dd)
}

cont_exposures <- c("z_plt_2011", "z_wbc_2011", "z_logcrp_2011")
cont_labels <- c("PLT per 1 SD", "WBC per 1 SD", "log(CRP+1) per 1 SD")
quart_exposures <- c("q_plt_2011", "q_wbc_2011", "q_crp_2011")
quart_labels <- c("PLT quartile", "WBC quartile", "CRP quartile")

main_cont <- bind_rows(lapply(seq_along(cont_exposures), function(i) {
  ex <- cont_exposures[i]
  lab <- cont_labels[i]
  m1 <- fit_logit2(d, ex, base_covars)$fit
  m2 <- fit_logit2(d, ex, full_covars)$fit
  bind_rows(
    extract_glm(m1, lab, "Base-adjusted", "Incident cognitive impairment 2020"),
    extract_glm(m2, lab, "Fully-adjusted", "Incident cognitive impairment 2020")
  )
}))
main_cont <- bind_rows(
  main_cont,
  bind_rows(lapply(seq_along(cont_exposures), function(i) {
    ex <- cont_exposures[i]; lab <- cont_labels[i]
    m1 <- fit_lm2(d, "y_change_10", ex, base_covars)$fit
    m2 <- fit_lm2(d, "y_change_10", ex, full_covars)$fit
    bind_rows(
      extract_lm(m1, lab, "Base-adjusted", "Cognitive change 2011-2020"),
      extract_lm(m2, lab, "Fully-adjusted", "Cognitive change 2011-2020")
    )
  }))
)

main_quart <- bind_rows(lapply(seq_along(quart_exposures), function(i) {
  ex <- quart_exposures[i]; lab <- quart_labels[i]
  m1 <- fit_logit2(d, ex, base_covars)$fit
  m2 <- fit_logit2(d, ex, full_covars)$fit
  bind_rows(
    extract_glm(m1, lab, "Base-adjusted", "Incident cognitive impairment 2020"),
    extract_glm(m2, lab, "Fully-adjusted", "Incident cognitive impairment 2020")
  )
}))
main_quart <- bind_rows(
  main_quart,
  bind_rows(lapply(seq_along(quart_exposures), function(i) {
    ex <- quart_exposures[i]; lab <- quart_labels[i]
    m1 <- fit_lm2(d, "y_change_10", ex, base_covars)$fit
    m2 <- fit_lm2(d, "y_change_10", ex, full_covars)$fit
    bind_rows(
      extract_lm(m1, lab, "Base-adjusted", "Cognitive change 2011-2020"),
      extract_lm(m2, lab, "Fully-adjusted", "Cognitive change 2011-2020")
    )
  }))
)

# change exposures: using 2015 cognition as baseline, test 2011→2015 blood biomarker change
change_exposures <- c("z_dplt", "z_dwbc", "z_dlogcrp")
change_labels <- c("delta PLT per 1 SD", "delta WBC per 1 SD", "delta log(CRP+1) per 1 SD")
change_covars <- c("age2011", "female", "edu_group", "rural", "y2015")
change_models <- bind_rows(lapply(seq_along(change_exposures), function(i) {
  ex <- change_exposures[i]; lab <- change_labels[i]
  m1 <- fit_lm2(d, "y_change_5", ex, change_covars)$fit
  m2 <- fit_lm2(d, "y_change_5", ex, c(change_covars, "cesd10", "adl", "hibp", "diab", "heart", "stroke", "lung", "kidney", "dyslip", "sleep_h", "pain_any"))$fit
  bind_rows(
    extract_lm(m1, lab, "Base-adjusted", "Cognitive change 2015-2020"),
    extract_lm(m2, lab, "Fully-adjusted", "Cognitive change 2015-2020")
  )
}))

# outlier sensitivity: 1%/99% winsorization
sens_cont <- bind_rows(lapply(seq_along(cont_exposures), function(i) {
  ex <- c("z_plt_2011_w", "z_wbc_2011_w", "z_logcrp_2011_w")[i]
  lab <- cont_labels[i]
  m1 <- fit_logit2(d, ex, base_covars)$fit
  m2 <- fit_lm2(d, "y_change_10", ex, base_covars)$fit
  bind_rows(
    extract_glm(m1, lab, "Winsorized base-adjusted", "Incident cognitive impairment 2020"),
    extract_lm(m2, lab, "Winsorized base-adjusted", "Cognitive change 2011-2020")
  )
}))

# FDR: primary continuous exposure, by outcome/adjust separately within model families
add_fdr <- function(df) {
  out <- df
  out$q <- NA_real_
  keys <- unique(out[, c("outcome", "model")])
  for (i in seq_len(nrow(keys))) {
    idx <- which(out$outcome == keys$outcome[i] & out$model == keys$model[i] & grepl("per 1 SD", out$exposure))
    if (length(idx)) out$q[idx] <- p.adjust(out$p[idx], method = "BH")
  }
  out
}
main_cont <- add_fdr(main_cont)
change_models <- add_fdr(change_models)
sens_cont <- add_fdr(sens_cont)

# output CSV(all data labels use ASCII, avoid Chinese encoding problems)
write.csv(main_cont, file.path(out_dir, "blood_main_continuous_2026-09-19.csv"), row.names = FALSE)
write.csv(main_quart, file.path(out_dir, "blood_main_quartile_2026-09-19.csv"), row.names = FALSE)
write.csv(change_models, file.path(out_dir, "blood_change_models_2026-09-19.csv"), row.names = FALSE)
write.csv(sens_cont, file.path(out_dir, "blood_winsor_sensitivity_2026-09-19.csv"), row.names = FALSE)

fmt <- function(x, d = 3) formatC(x, format = "f", digits = d)
ci_or <- function(r) paste0(fmt(r$OR), " (", fmt(r$lower), "-", fmt(r$upper), ")")
ci_b <- function(r) paste0(fmt(r$beta), " (", fmt(r$lower, 3), "-", fmt(r$upper, 3), ")")
exposure_cn <- c(
  "PLT per 1 SD" = "PLT(per 1 SD)",
  "WBC per 1 SD" = "WBC(per 1 SD)",
  "log(CRP+1) per 1 SD" = "log(CRP+1)(per 1 SD)",
  "PLT quartile" = "PLT quartiles",
  "WBC quartile" = "WBC quartiles",
  "CRP quartile" = "CRP quartiles",
  "delta PLT per 1 SD" = "ΔPLT(per 1 SD)",
  "delta WBC per 1 SD" = "ΔWBC(per 1 SD)",
  "delta log(CRP+1) per 1 SD" = "Δlog(CRP+1)(per 1 SD)"
)
tr <- function(x) ifelse(x %in% names(exposure_cn), unname(exposure_cn[x]), x)

lines <- c(
  "# blood analysis step three · PLT/WBC/CRP and cognitive outcomes",
  "",
  "> generation date: 2026-09-20(recompute after correcting age eligibility)",
  "",
  "## one, primary model results summary",
  "",
  "### 2020 incident cognitive impairment(Base-adjusted)",
  "",
  "| exposure | OR (95%CI) | p | q | n |",
  "|---|---:|---:|---:|---:|"
)
tmp <- main_cont %>% filter(outcome == "Incident cognitive impairment 2020", model == "Base-adjusted")
for (i in seq_len(nrow(tmp))) lines <- c(lines, paste0("| ", tr(tmp$exposure[i]), " | ", ci_or(tmp[i, ]), " | ", fmt(tmp$p[i], 4), " | ", fmt(tmp$q[i], 4), " | ", tmp$n[i], " |"))
lines <- c(lines, "", "### 2011→2020 cognitive change(Base-adjusted)", "", "| exposure | β (95%CI) | p | q | n |", "|---|---:|---:|---:|---:|")
tmp <- main_cont %>% filter(outcome == "Cognitive change 2011-2020", model == "Base-adjusted")
for (i in seq_len(nrow(tmp))) lines <- c(lines, paste0("| ", tr(tmp$exposure[i]), " | ", ci_b(tmp[i, ]), " | ", fmt(tmp$p[i], 4), " | ", fmt(tmp$q[i], 4), " | ", tmp$n[i], " |"))
lines <- c(lines, "", "### 2011→2015 blood biomarker change → 2015→2020 cognitive change(Base-adjusted)", "", "| exposure | β (95%CI) | p | q | n |", "|---|---:|---:|---:|---:|")
tmp <- change_models %>% filter(model == "Base-adjusted")
for (i in seq_len(nrow(tmp))) lines <- c(lines, paste0("| ", tr(tmp$exposure[i]), " | ", ci_b(tmp[i, ]), " | ", fmt(tmp$p[i], 4), " | ", fmt(tmp$q[i], 4), " | ", tmp$n[i], " |"))
lines <- c(lines, "", "## two, quartile results", "", "quartile models in Excel of `Main_quartile`; Q1 as reference.", "", "## three, sensitivity analysis", "", "1%/99% winsorization results in `Winsorized`.", "", "## four, interpretation boundaries", "", "- blood biomarkers represent observational associations, not causal effects.", "- blood subsample subject to healthy participant/survivor selection, disclose in manuscript.", "- Neu/NLR/PLR cannot derive from CHARLS released data.", "- next add: blood sample weights, longitudinal missingness, extreme value exclusion and main manuscript Table 5/Figure 6 integration with.")
writeLines(lines, file.path(out_dir, "blood_analysis_report_2026-09-19.md"), useBytes = TRUE)

cat("complete: blood model CSV outputs\\n")
cat("primary continuous model: \n"); print(main_cont %>% filter(model == "Base-adjusted") %>% select(outcome, exposure, OR, lower, upper, p, q, n))
cat("change model: \n"); print(change_models %>% filter(model == "Base-adjusted") %>% select(exposure, beta, lower, upper, p, q, n))
