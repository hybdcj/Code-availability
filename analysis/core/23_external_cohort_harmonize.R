# ============================================================================
# 23_external_cohort_harmonize.R —— A11 v2(third round audited correction)
# external cohort harmonization and fixed score mapping
# key corrections(2026-09-26 third round audit): 
#   1) cutoff and CHARLS consistent with development definition: by“all baseline eligible participants”calculate mean−1SD, 
#      follow up memory availability not required; integer threshold = floor(mean−1SD).
#   2) SHARE ISCED-97 of 95(still in education)/97(other)treat as missing, no longer misclassified as high school or above.
#   3) SHARE retain country_mod(14 countries)and country_label, for A12 country clustering/stratification.
#   4) also output cutoff definition sensitivity: baseline_all, complete_fu, fixed CHARLS ≤4.
# dependencies: D:/codex_tmp/{hrs,elsa,share} extracted under DTA(raw archives in external_data/<cohort>/raw|harmonized/)
# output: external_data/<cohort>/analysis/<cohort>_analysis_2026-09-26.rds(primary analysis)
#       external_data/mapping/A11_three_cohort_score_summary_2026_09_26.csv
#       external_data/mapping/A11_cutoff_sensitivity_2026_09_26.csv
#       external_data/mapping/A11_cutoff_sensitivity_stratification_2026_09_26.csv
#       external_data/mapping/A11_cohort_flowchart_2026_09_26.csv
# ============================================================================
invisible(try(Sys.setlocale("LC_ALL", "zh_CN.UTF-8"), silent = TRUE))
.libPaths(c("D:/R/library", .libPaths()))
suppressPackageStartupMessages({ library(haven); library(dplyr); library(tidyr); library(pROC) })

BASE <- "D:/BaiduSyncdisk/research/Xu_Zhonghao/external_data"
OUT  <- Sys.getenv("A11_OUT", unset = BASE)
dir.create(file.path(OUT, "mapping"), showWarnings = FALSE, recursive = TRUE)
for (co in c("HRS","ELSA","SHARE")) dir.create(file.path(OUT, co, "analysis"), showWarnings = FALSE, recursive = TRUE)

score_age <- function(age) 2 * pmax(floor((age - 45)/5), 0)
clean_neg <- function(x) { x <- as.numeric(x); x[x < 0] <- NA; x }

derive_edu_group <- function(pts) {
  factor(dplyr::case_when(
    is.na(pts) ~ NA_character_,
    pts ==   0 ~ "Less than primary",
    pts ==  -2 ~ "Primary school",
    pts ==  -7 ~ "Middle school",
    pts == -10 ~ "High school or above",
    TRUE       ~ NA_character_
  ), levels = c("Less than primary","Primary school","Middle school","High school or above"))
}

add_score <- function(d) {
  d$score     <- score_age(d$age) + d$edu_pts     + 3 * d$dep
  d$score_alt <- score_age(d$age) + d$edu_pts_alt + 3 * d$dep
  d$edu_group <- derive_edu_group(d$edu_pts)
  d
}

build_hrs <- function() {
  rd <- as.data.frame(read_dta("D:/codex_tmp/hrs/randhrs1992_2022v1.dta",
    col_select = any_of(c("hhidpn","raedyrs","r10agey_b","r10cesd","r10adl5a",
                          "r10imrc","r10dlrc","r12imrc","r12dlrc"))))
  raw_n <- nrow(rd)
  rd <- as.data.frame(lapply(rd, function(x) { x <- as.numeric(x); x[x < 0] <- NA; x }))
  d <- rd %>% transmute(
    ID = hhidpn, age = r10agey_b,
    edu_pts = case_when(is.na(raedyrs) ~ NA_real_, raedyrs == 0 ~ 0, raedyrs <= 6 ~ -2,
                        raedyrs <= 11 ~ -7, raedyrs >= 12 ~ -10),
    edu_pts_alt = case_when(is.na(raedyrs) ~ NA_real_, raedyrs <= 5 ~ 0, raedyrs <= 11 ~ -2,
                            raedyrs == 12 ~ -7, raedyrs >= 13 ~ -10),
    dep = as.integer(r10cesd >= 4), dep_raw = r10cesd,
    adl = as.integer(r10adl5a >= 1),
    mem_base = r10imrc + r10dlrc,
    mem_fu   = r12imrc + r12dlrc)
  age45_n <- sum(!is.na(d$age) & d$age >= 45)
  d <- d %>% filter(!is.na(age), age >= 45) %>% filter(!is.na(edu_pts), !is.na(dep), !is.na(mem_base))
  d <- add_score(as.data.frame(d))
  list(d = d, raw_n = raw_n, age45_n = age45_n, base_n = nrow(d))
}

build_elsa <- function() {
  el <- as.data.frame(read_dta("D:/codex_tmp/elsa/gh_elsa_h.dta",
    col_select = any_of(c("idauniq","raeducl","r4agey","r4cesd","r4adlfive",
                          "r4imrc","r4dlrc","r6imrc","r6dlrc"))))
  raw_n <- nrow(el)
  el <- as.data.frame(lapply(el, function(x) { x <- as.numeric(x); x[x < 0] <- NA; x }))
  d <- el %>% transmute(
    ID = idauniq, age = r4agey,
    edu_pts = case_when(is.na(raeducl) ~ NA_real_, raeducl == 1 ~ -7, raeducl >= 2 ~ -10),
    edu_pts_alt = case_when(is.na(raeducl) ~ NA_real_, raeducl == 1 ~ -2,
                            raeducl == 2 ~ -7, raeducl == 3 ~ -10),
    dep = as.integer(r4cesd >= 4), dep_raw = r4cesd,
    adl = as.integer(r4adlfive >= 1),
    mem_base = r4imrc + r4dlrc,
    mem_fu   = r6imrc + r6dlrc)
  age45_n <- sum(!is.na(d$age) & d$age >= 45)
  d <- d %>% filter(!is.na(age), age >= 45) %>% filter(!is.na(edu_pts), !is.na(dep), !is.na(mem_base))
  d <- add_score(as.data.frame(d))
  list(d = d, raw_n = raw_n, age45_n = age45_n, base_n = nrow(d))
}

build_share <- function() {
  sh <- as.data.frame(read_dta("D:/codex_tmp/share/easySHARE_rel9-0-0.dta",
    col_select = c(mergeid, wave, country_mod, age, isced1997_r, eurod, adla, recall_1, recall_2)))
  raw_n <- nrow(sh)
  lab <- attr(sh$country_mod, "labels")
  lab_map <- setNames(sub("^[0-9]+\\.\\s*", "", names(lab)), as.integer(lab))
  sh$mergeid <- as.character(sh$mergeid)
  for (v in setdiff(names(sh), "mergeid")) sh[[v]] <- clean_neg(sh[[v]])
  # ISCED-97 valid values 0–6; 95=still in education, 97=otherwise treat as missing
  sh$isced1997_r[!is.na(sh$isced1997_r) & sh$isced1997_r > 6] <- NA

  a <- sh[sh$wave == 4, c("mergeid","country_mod","age","isced1997_r","eurod","adla","recall_1","recall_2")]
  names(a) <- c("ID","country_mod","age","isced","eurod","adla","rec1","rec2")
  b <- sh[sh$wave == 6, c("mergeid","recall_1","recall_2")]
  names(b) <- c("ID","rec1_fu","rec2_fu")
  d <- merge(a, b, by = "ID", all.x = TRUE)  # retain all wave4 baseline eligible participants; wave6 participants with missing values mem_fu=NA
  d <- d %>% transmute(
    ID, country_mod = as.numeric(country_mod),
    country_label = unname(lab_map[as.character(country_mod)]),
    age,
    edu_pts = case_when(is.na(isced) ~ NA_real_, isced <= 1 ~ 0, isced == 2 ~ -7, isced >= 3 ~ -10),
    edu_pts_alt = case_when(is.na(isced) ~ NA_real_, isced <= 1 ~ 0, isced == 2 ~ -2,
                            isced == 3 ~ -7, isced >= 4 ~ -10),
    dep = as.integer(eurod >= 4), dep_raw = eurod,
    adl = as.integer(adla >= 1),
    mem_base = rec1 + rec2,
    mem_fu   = rec1_fu + rec2_fu)
  age45_n <- sum(!is.na(d$age) & d$age >= 45)
  d <- d %>% filter(!is.na(age), age >= 45) %>% filter(!is.na(edu_pts), !is.na(dep), !is.na(mem_base))
  d <- add_score(as.data.frame(d))
  wave4 <- sh[sh$wave == 4, ]
  wave4_n <- nrow(wave4)
  wave4_age45_n <- sum(!is.na(wave4$age) & wave4$age >= 45)
  wave4_base_n <- sum(!is.na(wave4$age) & wave4$age >= 45 & !is.na(wave4$isced1997_r) &
                      !is.na(wave4$eurod) & !is.na(wave4$recall_1 + wave4$recall_2))
  list(d = d, raw_n = raw_n, age45_n = age45_n, base_n = nrow(d),
       wave4_n = wave4_n, wave4_age45_n = wave4_age45_n, wave4_base_n = wave4_base_n)
}

analyse_cut <- function(cohort, d, cut_def) {
  if (cut_def == "baseline_all_mean_sd") {
    cut_raw <- mean(d$mem_base) - sd(d$mem_base)
    cut_int <- floor(cut_raw)
    cut_sample_n <- nrow(d)
  } else if (cut_def == "complete_fu_mean_sd") {
    dc <- d[!is.na(d$mem_fu), ]
    cut_raw <- mean(dc$mem_base) - sd(dc$mem_base)
    cut_int <- floor(cut_raw)
    cut_sample_n <- nrow(dc)
  } else if (cut_def == "fixed_charls_le4") {
    cut_raw <- 4
    cut_int <- 4
    cut_sample_n <- NA_integer_
  } else stop("unknown cut_def")

  risk <- d %>% filter(mem_base > cut_int, !is.na(mem_fu)) %>%
    mutate(event = as.integer(mem_fu <= cut_int),
           cut = cut_raw, cut_int = cut_int, cohort = cohort)
  strata <- risk %>% mutate(stratum = cut(score, c(-Inf, -3, 4, Inf), labels = c("low","in","high"))) %>%
    group_by(stratum) %>% summarise(n = n(), events = sum(event), rate = mean(event), .groups = "drop") %>%
    mutate(cohort = cohort, cut_definition = cut_def) %>% relocate(cohort, cut_definition)
  r1 <- roc(risk$event, risk$score, direction = "<", quiet = TRUE)
  r2 <- roc(risk$event, risk$score_alt, direction = "<", quiet = TRUE)
  overall <- data.frame(
    cohort = cohort, cut_definition = cut_def,
    baseline_n = nrow(d), followup_n = sum(!is.na(d$mem_fu)), cut_sample_n = cut_sample_n,
    cut_raw = cut_raw, cut_int = cut_int, risk_n = nrow(risk), events = sum(risk$event),
    event_rate = mean(risk$event), score_mean = mean(risk$score), score_sd = sd(risk$score),
    auc_score = as.numeric(auc(r1)), auc_score_low = ci.auc(r1)[1], auc_score_high = ci.auc(r1)[3],
    auc_score_alt = as.numeric(auc(r2)), auc_score_alt_low = ci.auc(r2)[1], auc_score_alt_high = ci.auc(r2)[3])
  list(overall = overall, strata = strata, risk = risk)
}

hrs  <- build_hrs();  elsa <- build_elsa();  share <- build_share()
cohorts <- list(HRS = hrs$d, ELSA = elsa$d, SHARE = share$d)
cut_defs <- c("baseline_all_mean_sd","complete_fu_mean_sd","fixed_charls_le4")
all_overall <- list(); all_strata <- list(); primary <- list(); sens_risk <- list()
for (nm in names(cohorts)) {
  for (cd in cut_defs) {
    res <- analyse_cut(nm, cohorts[[nm]], cd)
    all_overall[[paste(nm, cd)]] <- res$overall
    all_strata[[paste(nm, cd)]] <- res$strata
    if (cd == "baseline_all_mean_sd") primary[[nm]] <- res$risk
    sens_risk[[paste(nm, cd, sep = "__")]] <- res$risk
  }
}
overall <- bind_rows(all_overall)
strata  <- bind_rows(all_strata)
primary_overall <- overall %>% filter(cut_definition == "baseline_all_mean_sd")

flow <- bind_rows(
  data.frame(cohort = "HRS",
    step = c("raw rows","age>=45","baseline eligible (age/edu/dep/mem_base)",
             "follow-up memory available","primary risk set (baseline>cut_int & follow-up available)"),
    n = c(hrs$raw_n, hrs$age45_n, hrs$base_n, sum(!is.na(hrs$d$mem_fu)), nrow(primary$HRS))),
  data.frame(cohort = "ELSA",
    step = c("raw rows","age>=45","baseline eligible (age/edu/dep/mem_base)",
             "follow-up memory available","primary risk set (baseline>cut_int & follow-up available)"),
    n = c(elsa$raw_n, elsa$age45_n, elsa$base_n, sum(!is.na(elsa$d$mem_fu)), nrow(primary$ELSA))),
  data.frame(cohort = "SHARE",
    step = c("raw person-wave rows","wave 4 records","wave 4 age>=45",
             "wave 4 baseline eligible (age/edu/dep/mem_base)",
             "follow-up memory available","primary risk set (baseline>cut_int & follow-up available)"),
    n = c(share$raw_n, share$wave4_n, share$wave4_age45_n, share$base_n,
          sum(!is.na(share$d$mem_fu)), nrow(primary$SHARE)))
)

write.csv(primary_overall, file.path(OUT, "mapping/A11_three_cohort_score_summary_2026_09_26.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(overall, file.path(OUT, "mapping/A11_cutoff_sensitivity_2026_09_26.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(strata, file.path(OUT, "mapping/A11_cutoff_sensitivity_stratification_2026_09_26.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(flow, file.path(OUT, "mapping/A11_cohort_flowchart_2026_09_26.csv"), row.names = FALSE, fileEncoding = "UTF-8")

saveRDS(primary$HRS,  file.path(OUT, "HRS/analysis/HRS_analysis_2026-09-26.rds"))
saveRDS(primary$ELSA, file.path(OUT, "ELSA/analysis/ELSA_analysis_2026-09-26.rds"))
saveRDS(primary$SHARE,file.path(OUT, "SHARE/analysis/SHARE_analysis_2026-09-26.rds"))

for (nm in names(cohorts)) {
  keys <- grep(paste0("^", nm, "__"), names(sens_risk), value = TRUE)
  saveRDS(sens_risk[keys], file.path(OUT, nm, "analysis", paste0(nm, "_sensitivity_risk_sets_2026-09-26.rds")))
}

cat("\n================ A11 v2 primary analysis(baseline_all_mean_sd) ================\n")
print(primary_overall %>% select(cohort, baseline_n, followup_n, cut_raw, cut_int, risk_n, events,
                                 event_rate, score_mean, score_sd, auc_score, auc_score_low, auc_score_high) %>% mutate(across(where(is.numeric), ~round(.x, 4))))
cat("\n================ stratum specific primary analysis event rates ================\n")
print(strata %>% filter(cut_definition == "baseline_all_mean_sd") %>% mutate(rate = round(rate,4)) %>% as.data.frame())
cat("\n================ cutoff definition sensitivity ================\n")
print(overall %>% select(cohort, cut_definition, cut_raw, cut_int, risk_n, events, event_rate, auc_score) %>%
      mutate(across(where(is.numeric), ~round(.x, 4))) %>% arrange(cohort, cut_definition) %>% as.data.frame())
cat("\n written: ", OUT, "\n")