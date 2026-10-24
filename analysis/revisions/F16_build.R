# ============================================================================
# F16 build: reconstruct external cohorts, **no baseline memory restriction**
# usage: Rscript F16_build.R {HRS|ELSA|SHARE}
# variable mapping and missing codes copied exactly code/23_external_cohort_harmonize_v2.R(A11 v2)
# only change: no longer filter(!is.na(mem_base)); no longer requires mem_base > cut
# output: C:/y03tmp/tmp/F16_base_<cohort>.rds
# ============================================================================
args <- commandArgs(trailingOnly = TRUE)
co <- args[1]
suppressPackageStartupMessages({ library(haven); library(dplyr) })
if (dir.exists("D:/R/library")) .libPaths(c("D:/R/library", .libPaths()))
TMP <- "C:/y03tmp/tmp"
neg2na <- function(x) { x <- as.numeric(x); x[x < 0] <- NA; x }
say <- function(...) cat(paste0(..., collapse = ""), "\n")

if (co == "HRS") {
  v <- c("hhidpn","raedyrs","ragender","r10agey_b","r10cesd","r10adl5a",
         "r10imrc","r10dlrc","r12imrc","r12dlrc")
  rd <- as.data.frame(read_dta("D:/codex_tmp/hrs/randhrs1992_2022v1.dta", col_select = any_of(v)))
  raw_n <- nrow(rd)
  say("HRS original row count: ", raw_n, "｜column read: ", paste(names(rd), collapse = ", "))
  rd <- as.data.frame(lapply(rd, neg2na))
  sx <- if ("ragender" %in% names(rd)) rd$ragender else rep(NA_real_, nrow(rd))
  d <- data.frame(
    ID = rd$hhidpn, age = rd$r10agey_b, sex = sx,
    edu_pts = dplyr::case_when(is.na(rd$raedyrs) ~ NA_real_, rd$raedyrs == 0 ~ 0,
                               rd$raedyrs <= 6 ~ -2, rd$raedyrs <= 11 ~ -7,
                               rd$raedyrs >= 12 ~ -10),
    dep = as.integer(rd$r10cesd >= 4),
    dep_raw = rd$r10cesd,
    adl = as.integer(rd$r10adl5a >= 1),
    mem_base = rd$r10imrc + rd$r10dlrc,
    mem_fu   = rd$r12imrc + rd$r12dlrc,
    cohort = "HRS")
}
if (co == "ELSA") {
  v <- c("idauniq","raeducl","ragender","r4sex","r4agey","r4cesd","r4adlfive",
         "r4imrc","r4dlrc","r6imrc","r6dlrc")
  rd <- as.data.frame(read_dta("D:/codex_tmp/elsa/gh_elsa_h.dta", col_select = any_of(v)))
  raw_n <- nrow(rd)
  say("ELSA original row count: ", raw_n, "｜column read: ", paste(names(rd), collapse = ", "))
  rd <- as.data.frame(lapply(rd, neg2na))
  sx <- if ("ragender" %in% names(rd)) rd$ragender else
        if ("r4sex" %in% names(rd)) rd$r4sex else rep(NA_real_, nrow(rd))
  d <- data.frame(
    ID = rd$idauniq, age = rd$r4agey, sex = sx,
    edu_pts = dplyr::case_when(is.na(rd$raeducl) ~ NA_real_, rd$raeducl == 1 ~ -7,
                               rd$raeducl >= 2 ~ -10),
    dep = as.integer(rd$r4cesd >= 4),
    dep_raw = rd$r4cesd,
    adl = as.integer(rd$r4adlfive >= 1),
    mem_base = rd$r4imrc + rd$r4dlrc,
    mem_fu   = rd$r6imrc + rd$r6dlrc,
    cohort = "ELSA")
}
if (co == "SHARE") {
  v <- c("mergeid","wave","country_mod","age","isced1997_r","eurod","adla",
         "gender","female","recall_1","recall_2")
  sh <- as.data.frame(read_dta("D:/codex_tmp/share/easySHARE_rel9-0-0.dta", col_select = any_of(v)))
  raw_n <- nrow(sh)
  say("SHARE original row count: ", raw_n, "｜column read: ", paste(names(sh), collapse = ", "))
  lab <- attr(sh$country_mod, "labels")
  lab_map <- setNames(sub("^[0-9]+\\.\\s*", "", names(lab)), as.integer(lab))
  sh$mergeid <- as.character(sh$mergeid)
  for (nm in setdiff(names(sh), "mergeid")) sh[[nm]] <- neg2na(sh[[nm]])
  sh$isced1997_r[!is.na(sh$isced1997_r) & sh$isced1997_r > 6] <- NA
  w4 <- sh[sh$wave == 4, ]
  w6 <- sh[sh$wave == 6, c("mergeid","recall_1","recall_2")]
  names(w6) <- c("ID","rec1_fu","rec2_fu")
  sx <- if ("gender" %in% names(w4)) w4$gender else
        if ("female" %in% names(w4)) w4$female else rep(NA_real_, nrow(w4))
  a <- data.frame(ID = w4$mergeid, country_mod = w4$country_mod, age = w4$age,
                  isced = w4$isced1997_r, eurod = w4$eurod, adla = w4$adla,
                  rec1 = w4$recall_1, rec2 = w4$recall_2, sex = sx,
                  country_label = unname(lab_map[as.character(w4$country_mod)]))
  # and repository version code/23_external_cohort_harmonize.R consistent: left join, retain all wave 4 respondents
  m <- merge(a, w6, by = "ID", all.x = TRUE)
  d <- data.frame(
    ID = m$ID, age = m$age, sex = m$sex, country_mod = m$country_mod,
    country_label = m$country_label,
    edu_pts = dplyr::case_when(is.na(m$isced) ~ NA_real_, m$isced <= 1 ~ 0,
                               m$isced == 2 ~ -7, m$isced >= 3 ~ -10),
    dep = as.integer(m$eurod >= 4), dep_raw = m$eurod,
    adl = as.integer(m$adla >= 1),
    mem_base = m$rec1 + m$rec2, mem_fu = m$rec1_fu + m$rec2_fu,
    cohort = "SHARE")
}

# harmonization: apply only necessary age and three score item filters, **do not change mem_base and mem_fu**
d <- d[!is.na(d$age) & !is.na(d$edu_pts) & !is.na(d$dep), ]
say("age/education/complete depression: ", nrow(d), "｜of which mem_base available ", sum(!is.na(d$mem_base)),
    "｜mem_fu available ", sum(!is.na(d$mem_fu)))
saveRDS(d, file.path(TMP, paste0("F16_base_", co, ".rds")))
say("written F16_base_", co, ".rds; raw_n = ", raw_n)
