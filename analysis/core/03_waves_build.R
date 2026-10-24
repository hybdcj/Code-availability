# ============================================================================
# 03_waves_build.R —— five wave cognitive dataset construction(route B: reconstruct each wave from raw modules using 2020 reconstruct definitions)
# version: v1.0(2026-09-13)
# input: data/raw/wave1_2011, wave2_2013, wave3_2015, wave4_2018 raw module of
#       (2020 reuse previously constructed wave data directly output/01_2020_single_wave_baseline/analysis_2020.rds)
# output: output/02_five_wave_data_and_quality_control/waves_cognition_2026-09-13.rds(combined longitudinal cognition data across five waves)
#       output/02_five_wave_data_and_quality_control/waves_cognition_table_2026-09-13.csv(wave specific tables)
# key definitions: 
#   ① 9 components: orientation 5 items, serial subtraction of 7(5 step), figure copying, immediate recall, delayed recall
#   ② three memory score versions: raw score(raw), proportional score(prorated, comparable across waves), within wave z points
#   ③ prevalence: each wave uses"its own education specific cutoff"(primary definition B)and overall/quartile cutoff(A/C)
#   ④ 2018 wave definition(2026-09-13 correction): primary definition uses DC028_W4(10 word list·immediate recall trial 1)+ DC047_W4(delayed)
#      and DC024_W4(copy overlapping pentagons); these three items follow the questionnaire PROCEDURE yes 45 assessed in full sample above age.
#      60 participants older than MMSE subset additionally includes 3 word list(DC013_W4_1 / DC015_W4)and place orientation, excluded from primary longitudinal definition.
#      additionally: number 4 round word list content and administration(3 learning trials)and 2011–2015 different, official recommendation is longitudinal equipercentile equating.
# ============================================================================
invisible(try(Sys.setlocale("LC_ALL", "zh_CN.UTF-8"), silent = TRUE))
if (dir.exists("D:/R/library")) .libPaths(c("D:/R/library", .libPaths()))
suppressPackageStartupMessages({ library(haven); library(dplyr) })

PROJ <- Sys.getenv("CHARLS_PROJ", "D:/BaiduSyncdisk/research/Xu_Zhonghao")
RAW  <- file.path(PROJ, "data/raw")
OUT  <- file.path(PROJ, "output", "02_five_wave_data_and_quality_control")
if (!dir.exists(OUT)) dir.create(OUT, recursive = TRUE)

## ---------- general utilities ----------
num <- function(x) as.numeric(x)
# word list recall: before ns values in word positions 1..ns treat as"word recalled"; handle indicator flags separately
#   none_var   : "no word recalled" flag(11 / 4), used to determine"actually administered"
#   refuse_var : "refused recall" flag(12)→ code as missing, and 2020 wave 01_build_2020.R consistent definition
recall_cnt <- function(d, vars, ns, none_var = NULL, refuse_var = NULL) {
  m <- sapply(vars, function(v) num(d[[v]]))
  m <- matrix(m, ncol = length(vars))
  cnt <- rowSums(m >= 1 & m <= ns, na.rm = TRUE)
  any_asked <- rowSums(!is.na(m)) > 0
  if (!is.null(none_var) && none_var %in% names(d)) {
    none <- num(d[[none_var]]); any_asked <- any_asked | (!is.na(none) & none != 0)
  }
  if (!is.null(refuse_var) && refuse_var %in% names(d)) {
    rv <- num(d[[refuse_var]]); cnt[!is.na(rv) & rv == 12] <- NA
  }
  cnt[!any_asked] <- NA
  as.numeric(cnt)
}
# serial subtraction of 7: answers assessed using"previous step − 7"classify as correct; first step uses 93(100−7)as correct
ser7_score <- function(d, ans_vars) {
  a <- sapply(ans_vars, function(v) { x <- num(d[[v]]); x[!is.na(x) & (x < 0 | x > 100)] <- NA; x })
  a <- matrix(a, ncol = length(ans_vars))
  n <- ncol(a); ok <- matrix(FALSE, nrow = nrow(a), ncol = n)
  ok[, 1] <- !is.na(a[, 1]) & a[, 1] == 93
  for (i in 2:n) ok[, i] <- !is.na(a[, i-1]) & !is.na(a[, i]) & a[, i] == a[, i-1] - 7
  ans <- rowSums(!is.na(a)) > 0
  s <- rowSums(ok, na.rm = TRUE); s[!ans] <- NA
  as.numeric(s)
}
# four education groups(and 2020 consistent definition: <primary school / primary school / middle school / high school or above)
edu4_raw <- function(x) {
  x <- num(x)
  dplyr::case_when(x %in% 1:3 ~ "Less than primary", x == 4 ~ "Primary school",
                   x == 5 ~ "Middle school", x >= 6 ~ "High school or above", TRUE ~ NA_character_)
}

## ---------- wave configuration ----------
## demographics consistently from Harmonized D(harmonized names, codes verified); cognition always from raw modules
cfg_year <- c(2011, 2013, 2015, 2018)     # and cfg consistent ordering, used to locate Harmonized of r1–r4 field
norm_id <- function(x) {
  x <- as.character(x); n <- nchar(x)
  ifelse(n == 11, paste0(substr(x, 1, 9), sprintf("%03d", as.integer(substr(x, 10, 11)))), x)
}
hd <- read_dta(file.path(RAW, "harmonized_d/H_CHARLS_D_Data.dta"),
               col_select = any_of(c("ID", "communityID",
                                     paste0("r", 1:4, "agey"), "ragender", "raeduc_c", "raeducl",
                                     paste0("r", 1:4, "hukou"), paste0("r", 1:4, "rural2"),
                                     paste0("r", 1:4, "adlab_c"), paste0("r", 2:4, "iadla"),
                                     paste0("r", 1:4, "cesd10"),
                                     paste0("r", 1:4, "hibpe"), paste0("r", 1:4, "diabe"),
                                     paste0("r", 1:4, "hearte"), paste0("r", 1:4, "stroke"),
                                     paste0("r", 1:4, "lunge"), paste0("r", 1:4, "kidneye"),
                                     paste0("r", 1:4, "dyslipe"),
                                     paste0("r", 1:4, "moneya"), paste0("r", 1:4, "medsa"),
                                     paste0("r", 1:4, "wtresp"))))
hd$IDn <- norm_id(hd$ID)

cfg <- list(
  list(year = 2011, dir = "wave1_2011", hsf = "health_status_and_functioning.dta",
       demo = "demographic_background.dta", wt = "weight.dta", wt_var = "ind_weight",
       orient = c("dc001s1","dc001s2","dc001s3","dc002","dc003"),
       orient_ok = c(1, 2, 3, 1, 1),
       ser7 = paste0("dc0", 19:23), draw = "dc025",
       pain = "da041", pain_ge2 = FALSE, pa_prefix = "da051_",
       imm = paste0("dc006s", 1:10), imm_none = "dc006s11",
       dlr = paste0("dc027s", 1:10), dlr_none = "dc027s11", draw_ok = 1),
  list(year = 2013, dir = "wave2_2013", hsf = "Health_Status_and_Functioning.dta",
       pain = "wb16", pain_ge2 = TRUE, pa_prefix = "da051_",
       demo = "Demographic_Background.dta", wt = "Weights.dta", wt_var = "INDV_weight",
       orient = c("dc001s1","dc001s2","dc001s3","dc002","dc003"),
       orient_ok = c(1, 2, 3, 1, 1),
       ser7 = paste0("dc0", 19:23), draw = "dc025",
       imm = paste0("dc006_1_s", 1:10), imm_none = "dc006_1_s11",
       dlr = paste0("dc027s", 1:10), dlr_none = "dc027s11", draw_ok = 1),
   list(year = 2015, dir = "wave3_2015", hsf = "Health_Status_and_Functioning.dta",
        pain = "da041", pain_ge2 = FALSE, pa_prefix = "da051_",
        demo = "Demographic_Background.dta", wt = "Weights.dta", wt_var = "INDV_weight",
        orient = c("dc001s1","dc001s2","dc001s3","dc002","dc003"),
        orient_ok = c(1, 2, 3, 1, 1),
        ser7 = paste0("dc0", 19:23), draw = "dc025",
        imm = paste0("dc006s", 1:10), imm_none = "dc006s11", imm_refuse = "dc006s12",
        dlr = paste0("dc027s", 1:10), dlr_none = "dc027s11", dlr_refuse = "dc027s12", draw_ok = 1),
  list(year = 2018, dir = "wave4_2018", hsf = "Cognition.dta",
       extra = "Health_Status_and_Functioning.dta",     # 2018 sleep from/chronic conditions in health module
       demo = "Demographic_Background.dta", wt = "Weights.dta", wt_var = "INDV_weight",
       orient = c("dc001_w4","dc006_w4","dc005_w4","dc003_w4","dc002_w4"),
       orient_ok = c(1, 1, 1, 1, 1),                       # 1=for / 5=incorrect / 97=do not know
       ser7 = paste0("dc014_w4_", 1:5, "_1"),
       draw = "dc024_w4", draw_err = 5,                     # copy overlapping pentagons: 1=for / 5=incorrect / 97=not assessed
       pain = "da041_w4", pain_ge2 = TRUE, pa_prefix = "da051_",
       imm = paste0("dc028_w4_s", 1:10), imm_none = "dc028_w4_s11", imm_refuse = "dc028_w4_s12",
       dlr = paste0("dc047_w4_s", 1:10), dlr_none = "dc047_w4_s11", dlr_refuse = "dc047_w4_s12",
       draw_ok = 1)
)

build_wave <- function(cf) {
  hs  <- read_dta(file.path(RAW, cf$dir, cf$hsf))
  dm  <- read_dta(file.path(RAW, cf$dir, cf$demo))
  wt  <- read_dta(file.path(RAW, cf$dir, cf$wt))
  cat("\n== ", cf$year, " ==\n  HSF row count ", nrow(hs), " ｜ demographics ", nrow(dm), " ｜ weights ", nrow(wt), "\n", sep = "")

  # orientation: value equals"correct code for the item"code as 1, otherwise code 0(NA treat as incorrect = 0; and 2020 take"do not know"code as 0 consistent)
  orient_m <- sapply(seq_along(cf$orient), function(i) {
    x <- num(hs[[cf$orient[i]]]); as.integer(!is.na(x) & x == cf$orient_ok[i])
  })
  orient_m <- matrix(orient_m, ncol = length(cf$orient))
  orient <- rowSums(orient_m)

  ser7 <- ser7_score(hs, cf$ser7)

  draw <- if (is.na(cf$draw)) rep(NA_integer_, nrow(hs)) else {
    x <- num(hs[[cf$draw]])
    err_code <- if (!is.null(cf$draw_err)) cf$draw_err else 2      # each wave"incorrect figure copy"codes differ: 2011/2013/2015=2, 2018=5
    x[!is.na(x) & !(x %in% c(1, err_code))] <- NA                 # 97/3/997/999 etc"not assessed/do not know"code as missing
    ifelse(is.na(x), NA_integer_, as.integer(x == cf$draw_ok))
  }
  imm <- recall_cnt(hs, cf$imm, length(cf$imm), cf$imm_none, cf$imm_refuse)
  dlr <- recall_cnt(hs, cf$dlr, length(cf$dlr), cf$dlr_none, cf$dlr_refuse)

  df <- data.frame(ID = as.character(hs$ID), year = cf$year,
                   orient = orient, ser7 = ser7, draw = draw, imrc = imm, dlrc = dlr)
  df$list_len_imm <- length(cf$imm); df$list_len_dlr <- length(cf$dlr)
  df$memory_raw <- df$imrc + df$dlrc
  df$memory_prop <- (df$imrc + df$dlrc) / (df$list_len_imm + df$list_len_dlr)   # A: prorated
  # mental status: when figure copying is available 3 item composite(0–11), when figure copying unavailable 2 items(0–10)
  df$mental <- if (all(is.na(draw))) rowSums(cbind(df$orient, df$ser7), na.rm = FALSE) else
    rowSums(cbind(df$orient, df$ser7, df$draw), na.rm = FALSE)
  df$draw_available <- !all(is.na(draw))
  df$cog_complete <- as.integer(!is.na(df$orient) & !is.na(df$ser7) &
                                  (!df$draw_available | !is.na(df$draw)) &
                                  !is.na(df$imrc) & !is.na(df$dlrc))

  # demographics: from Harmonized D corresponding wave field uses(w = 1..4)
  w <- which(sapply(cfg_year, function(y) y == cf$year))[1]
  dm2 <- data.frame(
    IDn = hd$IDn,
    age = num(hd[[paste0("r", w, "agey")]]),
    female = as.integer(num(hd$ragender) == 2),
    edu = num(hd$raeduc_c),
    rural = num(hd[[paste0("r", w, "rural2")]]),
    communityID = hd$communityID,
    weight = num(hd[[paste0("r", w, "wtresp")]]),
    adl    = as.integer(num(hd[[paste0("r", w, "adlab_c")]]) >= 1),
    iadl   = if (paste0("r", w, "iadla") %in% names(hd)) as.integer(num(hd[[paste0("r", w, "iadla")]]) >= 1) else NA_integer_,
    # IADL two item version(managing finances + taking medication): official summary r*iadla including"making telephone calls", whereas 2011 telephone item not asked in this wave, 
    # therefore r1iadla does not exist; these two items 2011 also available in this wave(r1moneya / r1medsa), supports a comparable two item definition across five waves
    iadl2  = {
      mo <- num(hd[[paste0("r", w, "moneya")]]); me <- num(hd[[paste0("r", w, "medsa")]])
      ifelse(is.na(mo) & is.na(me), NA_integer_,
             as.integer((!is.na(mo) & mo == 1) | (!is.na(me) & me == 1)))
    },
    cesd10 = num(hd[[paste0("r", w, "cesd10")]]),
    hibp   = num(hd[[paste0("r", w, "hibpe")]]),
    diab   = num(hd[[paste0("r", w, "diabe")]]),
    heart  = num(hd[[paste0("r", w, "hearte")]]),
    stroke = num(hd[[paste0("r", w, "stroke")]]),
    lung   = num(hd[[paste0("r", w, "lunge")]]),
    kidney = num(hd[[paste0("r", w, "kidneye")]]),
    dyslip = num(hd[[paste0("r", w, "dyslipe")]]))
  dm2$edu_group <- edu4_raw(dm2$edu)
  df$IDn <- norm_id(df$ID)
  slp_var <- if ("da049" %in% names(hs)) "da049" else if ("da030" %in% names(hs)) "da030" else NA
  if (is.na(slp_var) && !is.null(cf$extra)) {         # 2018: sleep in another module
    ex <- read_dta(file.path(RAW, cf$dir, cf$extra))
    if ("da049" %in% names(ex)) {
      x <- num(ex$da049); x[!is.na(x) & (x < 0 | x > 24)] <- NA
      df$sleep_h <- x[match(df$ID, as.character(ex$ID))]
      slp_var <- "da049"
    }
  }
  if (!"sleep_h" %in% names(df)) df$sleep_h <- if (is.na(slp_var)) NA_real_ else {
    x <- num(hs[[slp_var]]); x[!is.na(x) & (x < 0 | x > 24)] <- NA; x }
  # pain(coding differs by wave: da041 is 1=pain present/2=no pain; da041_w4, wb16, da027 is 1=no pain, 2-5=pain severity)
  px <- hs                                   # pain/module containing physical activity(2018 at extra within module)
  if (!is.null(cf$extra)) px <- read_dta(file.path(RAW, cf$dir, cf$extra))
  idx <- match(df$ID, as.character(px$ID))
  if (!is.null(cf$pain) && cf$pain %in% names(px)) {
    x <- num(px[[cf$pain]]); x[!is.na(x) & x %in% c(997, 999, -1)] <- NA
    v <- if (isTRUE(cf$pain_ge2)) as.integer(x >= 2) else as.integer(x == 1)
    df$pain_any <- v[idx]
  } else df$pain_any <- NA_integer_
  pa_vars <- intersect(paste0(cf$pa_prefix, 1:3, "_"), names(px))
  if (length(pa_vars)) {
    m <- sapply(pa_vars, function(v) num(px[[v]])); m <- matrix(m, ncol = length(pa_vars))
    v <- as.integer(rowSums(m == 1, na.rm = TRUE) > 0); v[rowSums(!is.na(m)) == 0] <- NA
    df$pa_any <- v[idx]
  } else df$pa_any <- NA_integer_
  df <- left_join(df, dm2, by = "IDn")

  # cutoff(calculate separately by wave): education strata B, overall −1SD of A, lowest quartile C
  prim <- df %>% filter(cog_complete == 1, !is.na(memory_raw), !is.na(mental))
  gm <- mean(prim$memory_raw + prim$mental, na.rm = TRUE); gs <- sd(prim$memory_raw + prim$mental, na.rm = TRUE)
  thrA <- gm - gs; thrC <- quantile(prim$memory_raw + prim$mental, 0.25, na.rm = TRUE)
  edu_thr <- prim %>% filter(!is.na(edu_group)) %>% group_by(edu_group) %>%
    summarise(thr = mean(memory_raw + mental, na.rm = TRUE) - sd(memory_raw + mental, na.rm = TRUE), .groups = "drop")
  df <- df %>% left_join(edu_thr, by = "edu_group") %>%
    mutate(total = memory_raw + mental,
           cog_B = as.integer(total < thr), cog_A = as.integer(total < thrA),
           cog_C = as.integer(total <= thrC))
  # within wave z(B scheme: for cross wave comparisons and trajectory models)
  df$total_z <- (df$total - gm) / gs
  df$memory_z <- (df$memory_raw - mean(prim$memory_raw, na.rm = TRUE)) / sd(prim$memory_raw, na.rm = TRUE)
  df$wave <- as.character(cf$year)
  df
}

waves <- lapply(cfg, build_wave)

## 2020 wave: reuse previously constructed dataset, map to common fields
a20 <- readRDS(file.path(PROJ, "output", "01_2020_single_wave_baseline", "analysis_2020.rds"))
w2020 <- data.frame(
  ID = as.character(a20$ID), year = 2020, wave = "2020",
  # [2026-10-11 repair]previously omitted IDn, bind_rows would convert 2020 rows with IDn into NA, 
  # prevents longitudinal linkage to other waves(see details in code/33_fix_waves_cognition_idn.R).
  # 2020 of ID already is 12 digit standardization ID, norm_id() returns unchanged.
  IDn = norm_id(as.character(a20$ID)),
  orient = a20$orient, ser7 = a20$ser7, draw = a20$draw, imrc = a20$imrc, dlrc = a20$dlrc,
  list_len_imm = 10, list_len_dlr = 10,
  memory_raw = a20$memory, memory_prop = a20$memory / 20,
  mental = a20$mental, draw_available = TRUE, cog_complete = a20$cog_complete,
  age = a20$age, female = a20$female, edu_group = a20$edu_group,
  communityID = a20$communityID, weight = a20$INDV_weight,
  total = a20$global, cog_B = a20$cog_low_edu, cog_A = a20$cog_low_1sd,
  cog_C = a20$cog_low_q1,
  adl = a20$adl_any, iadl = a20$iadl_any, cesd10 = a20$cesd10, sleep_h = a20$sleep_h,
  hibp = a20$hypertension, diab = a20$diabetes, heart = a20$heart_disease,
  stroke = a20$stroke, lung = a20$lung_disease, kidney = a20$kidney_disease,
  dyslip = a20$dyslipidemia, pain_any = a20$pain_any, pa_any = a20$pa_any)
# 2020 wave IADL two items(managing finances db022 + taking medication db020): and 2011-2018 of Harmonized consistent single item definitions(>=2 classify as difficulty)
hs20 <- read_dta(file.path(RAW, "wave5_2020/Health_Status_and_Functioning.dta"),
                 col_select = any_of(c("ID", "db020", "db022")))
mo20 <- num(hs20$db022); me20 <- num(hs20$db020)
mo20 <- ifelse(!is.na(mo20), as.integer(mo20 >= 2), NA_integer_)
me20 <- ifelse(!is.na(me20), as.integer(me20 >= 2), NA_integer_)
iadl2_20 <- ifelse(is.na(mo20) & is.na(me20), NA_integer_, as.integer(mo20 == 1 | me20 == 1))
w2020$iadl2 <- iadl2_20[match(as.character(a20$ID), as.character(hs20$ID))]
prim20 <- w2020 %>% filter(cog_complete == 1, !is.na(total))
w2020$total_z <- (w2020$total - mean(prim20$total, na.rm = TRUE)) / sd(prim20$total, na.rm = TRUE)
w2020$memory_z <- (w2020$memory_raw - mean(prim20$memory_raw, na.rm = TRUE)) / sd(prim20$memory_raw, na.rm = TRUE)
waves <- c(waves, list(w2020))

allw <- bind_rows(lapply(waves, function(d) d %>% mutate(across(where(is.numeric), ~ .))))
saveRDS(allw, file.path(OUT, "waves_cognition_2026-09-13.rds"))

## ---------- wave specific tables ----------
tab <- do.call(rbind, lapply(waves, function(d) {
  inw <- d[!is.na(d$age) & d$age >= 45, ]
  comp <- inw[inw$cog_complete == 1, ]
  data.frame(
    year = d$year[1],
    available_sample = nrow(inw),
    complete_cognition = nrow(comp),
    complete_memory = sum(!is.na(comp$memory_raw)),
    mean_age = round(mean(inw$age, na.rm = TRUE), 1),
    female = sprintf("%.1f%%", 100*mean(inw$female, na.rm = TRUE)),
    mean_raw_memory_score = round(mean(comp$memory_raw, na.rm = TRUE), 2),
    mean_proportional_memory_score = round(mean(comp$memory_prop, na.rm = TRUE), 3),
    mean_mental_status_score = round(mean(comp$mental, na.rm = TRUE), 2),
    mean_global_cognition_score = round(mean(comp$total, na.rm = TRUE), 2),
    prevalence_B = sprintf("%.1f%%", 100*mean(comp$cog_B, na.rm = TRUE)),
    prevalence_A = sprintf("%.1f%%", 100*mean(comp$cog_A, na.rm = TRUE)),
    prevalence_C = sprintf("%.1f%%", 100*mean(comp$cog_C, na.rm = TRUE)),
    ADL_limited = sprintf("%.1f%%", 100*mean(comp$adl, na.rm = TRUE)),
    IADL_limited = if (all(is.na(comp$iadl))) "-" else sprintf("%.1f%%", 100*mean(comp$iadl, na.rm = TRUE)),
    depression = sprintf("%.1f%%", 100*mean(comp$cesd10 >= 10, na.rm = TRUE)),
    hypertension = sprintf("%.1f%%", 100*mean(comp$hibp == 1, na.rm = TRUE)),
    diabetes = sprintf("%.1f%%", 100*mean(comp$diab == 1, na.rm = TRUE)),
    sleep_hours = round(mean(comp$sleep_h, na.rm = TRUE), 2),
    pain = sprintf("%.1f%%", 100*mean(comp$pain_any, na.rm = TRUE)),
    physical_activity = sprintf("%.1f%%", 100*mean(comp$pa_any, na.rm = TRUE)),
    including_figure_copying_item = unique(comp$draw_available)[1])
}))
write.csv(tab, file.path(OUT, "waves_cognition_table_2026-09-13.csv"), row.names = FALSE, fileEncoding = "UTF-8")
cat("\n================ wave specific tables(cognition) ================\n"); print(tab, row.names = FALSE)
cat("\n saved: output/02_five_wave_data_and_quality_control/waves_cognition_2026-09-13.rds, output/02_five_wave_data_and_quality_control/waves_cognition_table_2026-09-13.csv\n")
