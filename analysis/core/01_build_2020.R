# =============================================================================
# 01_build_2020.R  --  CHARLS 2020 (Wave 5) analysis dataset construction
# project: Prevalence of cognitive impairment and multidimensional preoperative vulnerability factors in middle-aged and older Chinese adults
# input: data/raw/wave5_2020/*.dta
# output: output/01_2020_single_wave_baseline/analysis_2020.rds, output/01_2020_single_wave_baseline/analysis_2020.csv, output/01_2020_single_wave_baseline/build_log.txt
# notes: raw/ read only; this script does not modify any raw files
# version: v3.2(2026-09-12)——education changed to zredu + ba010 combine(2020 complementary within wave, no missing values after combining)
# =============================================================================

## ---- 0. environment ----
Sys.setlocale("LC_ALL", "zh_CN.UTF-8")          # support Chinese paths
.libPaths(c("D:/R/library", .libPaths()))
suppressPackageStartupMessages({
  library(haven); library(dplyr); library(tidyr)
})

proj <- "D:/BaiduSyncdisk/research/Xu_Zhonghao"
raw_dir <- file.path(proj, "data/raw/wave5_2020")
out_dir <- file.path(proj, "output", "01_2020_single_wave_baseline")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

log_con <- file(file.path(out_dir, "build_log.txt"), open = "wt", encoding = "UTF-8")
say <- function(...) { msg <- paste0(...); cat(msg, "\n"); writeLines(msg, log_con) }

say("== CHARLS 2020 data construction ", format(Sys.time()), " ==")
say("definition notes(v3.2, 2026-09-12): orientation / serial subtraction of 7 / in figure copying[do not know (997)]always score incorrect as 0(counts as attempted); ")
say("         [refused (999)/not applicable (-1)]code as missing; figure copying[not assessed (3)]code as missing; ")
say("         in word recall[no words recalled]score 0 points, [refused recall/unable to test]code as missing; ")
say("         primary analysis restricted to cognition 9 complete components and known education; cutoff stratified by education(<within stratum mean−1SD).")

## ---- 1. read raw data(select required columns only) ----
hs_file   <- file.path(raw_dir, "Health_Status_and_Functioning.dta")
demo_file <- file.path(raw_dir, "Demographic_Background.dta")
wt_file   <- file.path(raw_dir, "Weights.dta")

hs_cols <- c(
  "ID",
  paste0("dc00", 1:5),                                    # orientation: year/season/day/day of week/month
  paste0("dc007_", 1:5), paste0("dc007_", 1:5, "_1"),     # serial subtraction of 7
  "dc009",                                                # figure copying
  paste0("dc012_s", 1:13), paste0("dc013_s", 1:12),       # immediate recall trial 1, 2(including s11/s12/s13 indicator flags)
  paste0("dc014_s", 1:12),                                # immediate recall trial 3
  paste0("dc028_s", 1:12),                                # delayed recall
  paste0("dc0", 16:25),                                   # CES-D-10
  paste0("da002_", 1:15, "_"),                            # chronic conditions(present in prior wave, comparison)
  paste0("da003_", 1:15, "_"),                            # chronic conditions(ever diagnosed by a doctor)
  paste0("db", sprintf("%03d", c(1,3,5,7,9,11,12,14,16,18,20,22))),  # ADL/IADL
  "da027", "da030", "da031",                              # pain, sleep
  "da046", "da047", "zsmoke", "da051",                    # smoking, alcohol use(da051 is past year drinking frequency)
  "da032_1_", "da032_2_", "da032_3_"                      # physical activity(high/in/low intensity ≥10 minutes)
)
hs <- read_dta(hs_file, col_select = all_of(hs_cols))

demo_cols <- c("ID","householdID","communityID","proxy_2","ba008","ba009","ba010","ba011",
               "xrage","xrgender","zrgender","zredu","xrpartner","zrbirthyear")
demo <- read_dta(demo_file, col_select = all_of(demo_cols))
wt   <- read_dta(wt_file, col_select = all_of(c("ID","HH_weight","HH_weight_ad1","INDV_weight","INDV_weight_ad2")))

say("HS row count: ", nrow(hs), " | Demo row count: ", nrow(demo), " | weight row count: ", nrow(wt))
say("note: Demographic module ", nrow(demo), " participants, of which ", nrow(demo) - nrow(hs),
    " participants without Health Status module records; analysis starts from HS module as starting point.")
say("education source: zredu nonmissing ", sum(!is.na(demo$zredu)), " ｜ ba010 nonmissing ", sum(!is.na(demo$ba010)),
    " ｜ missing after combining ", sum(is.na(dplyr::coalesce(as.numeric(demo$zredu), as.numeric(demo$ba010)))), "")

## ---- 2. helper functions ----
# special codes -> NA(997 do not know / 999 refused / -1 not applicable)
to_na <- function(x, codes = c(997, 999, -1)) { x <- as.numeric(x); x[x %in% codes] <- NA; x }

# word list recall scoring: s1..s10 value is the word index in the list (1-10), not recalled=0
# "no words recalled"(s11) -> 0 points; "refused recall"(s12) / "unable to test"(s13) -> missing
recall_score <- function(df, prefix, ns = 10, refuse_cols = character(0)) {
  m <- as.matrix(df[paste0(prefix, "_s", seq_len(ns))])
  m <- apply(m, 2, as.numeric)
  cnt <- rowSums(m >= 1 & m <= 10, na.rm = TRUE)
  asked <- rowSums(!is.na(m)) > 0
  cnt[!asked] <- NA
  for (fc in refuse_cols) {
    if (fc %in% names(df)) {
      v <- as.numeric(df[[fc]])
      cnt[!is.na(v) & v > 0] <- NA            # refused / unable to test -> missing
    }
  }
  as.numeric(cnt)
}

## ---- 3. cognitive domains ----
# 3.1 orientation dc001-dc005(1=correct; 2 incorrect and 997 do not know all scored as 0; 999 refused / -1 not applicable = missing)
orient_items <- lapply(paste0("dc00", 1:5), function(v) {
  x <- as.numeric(hs[[v]])
  x[x %in% c(999, -1)] <- NA                          # only refusal/not applicable coded missing
  ifelse(is.na(x), NA, ifelse(x == 1, 1, 0))          # incorrect and do not know both scored 0(HRS convention)
})
orient_df <- as.data.frame(orient_items); names(orient_df) <- paste0("orient", 1:5)

# 3.2 serial subtraction of 7 dc007_1..5(answers in _1_ variables; following official Harmonized logic relative to the previous answer)
# [v3 correction]response determination = numeric value present or explicit response"do not know (997)"; 997 scored incorrect(no points added).
#            999=refused / -1=not applicable / not assessed -> missing.
s7_num <- lapply(1:5, function(i) {
  x <- as.numeric(hs[[paste0("dc007_", i, "_1")]])
  x[!is.na(x) & (x < 0 | x > 1000)] <- NA       # obvious entry error / sentinel value
  x
})
s7_flg  <- lapply(1:5, function(i) as.numeric(hs[[paste0("dc007_", i)]]))
s7_eng  <- Map(function(n, f) (!is.na(n)) | (!is.na(f) & f == 997), s7_num, s7_flg)
s7 <- as.data.frame(s7_num); names(s7) <- paste0("s7_", 1:5)
ser7 <- with(s7, {
  ok1 <- !is.na(s7_1) & s7_1 == 93
  ok2 <- !is.na(s7_1) & !is.na(s7_2) & s7_2 == s7_1 - 7
  ok3 <- !is.na(s7_2) & !is.na(s7_3) & s7_3 == s7_2 - 7
  ok4 <- !is.na(s7_3) & !is.na(s7_4) & s7_4 == s7_3 - 7
  ok5 <- !is.na(s7_4) & !is.na(s7_5) & s7_5 == s7_4 - 7
  ok <- cbind(ok1, ok2, ok3, ok4, ok5) & do.call(cbind, s7_eng)   # unanswered steps receive no points
  score   <- rowSums(ok, na.rm = TRUE)
  any_ans <- Reduce(`|`, s7_eng)
  score[!any_ans] <- NA
  as.numeric(score)
})

# 3.3 figure copying dc009(1 correct→1; 2 incorrect, 997 do not know→0(incorrect answer); 3 not assessed / 999 refused / not assessed→missing)
draw_x <- as.numeric(hs$dc009); draw_x[draw_x %in% c(999, -1)] <- NA
draw <- ifelse(is.na(draw_x), NA,
               ifelse(draw_x == 1, 1, ifelse(draw_x %in% c(2, 997), 0, NA)))

# 3.4 word recall(v3: refused recall / unable to test -> missing; "no words recalled" -> 0 points)
imrc1 <- recall_score(hs, "dc012", refuse_cols = c("dc012_s12","dc012_s13"))  # immediate recall trial 1
imrc2 <- recall_score(hs, "dc013", refuse_cols = c("dc013_s12"))
imrc3 <- recall_score(hs, "dc014", refuse_cols = c("dc014_s12"))
imrc_mean <- rowMeans(cbind(imrc1, imrc2, imrc3), na.rm = TRUE)   # mean of three learning trials
dlrc  <- recall_score(hs, "dc028", refuse_cols = c("dc028_s12"))   # delayed recall

# 3.5 composite score(follow Harmonized CHARLS convention: rowtotal ignore missing items; code only when all items missing NA)
sum_avail <- function(mat) {
  mat <- as.matrix(mat)
  out <- rowSums(mat, na.rm = TRUE)
  out[rowSums(!is.na(mat)) == 0] <- NA
  as.numeric(out)
}
orient <- sum_avail(orient_df)                                   # orientation 0-5
memory <- sum_avail(cbind(imrc1, dlrc))                          # episodic memory 0-20
mental <- sum_avail(cbind(orient_df, ser7, draw))                # mental status 0-11
global <- ifelse(is.na(memory) | is.na(mental), NA, memory + mental)   # global cognition 0-31
# complete case indicator(9 cognitive components: immediate recall, delayed recall, orientation 5 items, serial subtraction of 7, figure copying, all nonmissing)
cog_mat <- cbind(imrc1, dlrc, orient_df, ser7, draw)
cog_items_avail <- rowSums(!is.na(cog_mat))
cog_complete <- as.integer(cog_items_avail == ncol(cog_mat))

## ---- 4. depression CES-D-10(dc016-dc025; 1-4 -> 0-3; dc020, dc023 reverse) ----
cesd_items <- paste0("dc0", 16:25)
cesd_raw <- as.data.frame(lapply(cesd_items, function(v) to_na(hs[[v]])))
names(cesd_raw) <- cesd_items
cesd <- as.data.frame(lapply(names(cesd_raw), function(v) {
  x <- cesd_raw[[v]]
  if (v %in% c("dc020","dc023")) { 4 - x } else { x - 1 }   # reverse scored items
}))
names(cesd) <- cesd_items
cesd10 <- rowSums(cesd)
cesd10[rowSums(!is.na(cesd)) < 10] <- NA                     # require 10 complete questions
cesd_avail <- rowSums(!is.na(cesd))
dep <- ifelse(is.na(cesd10), NA, as.integer(cesd10 >= 10))    # conventional threshold ≥10

## ---- 5. chronic conditions(da002 and da003 joint classification) ----
dis_labels <- c("hypertension","dyslipidemia","diabetes","cancer","lung_disease",
                "liver_disease","heart_disease","stroke","kidney_disease","digestive_disease",
                "psychiatric","memory_disease","parkinson","arthritis","asthma")
dis <- data.frame(ID = hs$ID)
for (i in 1:15) {
  prev <- to_na(hs[[paste0("da002_", i, "_")]])      # 1/2/3 = condition present in prior wave; 99 = absent at preceding wave
  ever <- to_na(hs[[paste0("da003_", i, "_")]])      # 1 = ever diagnosed by a doctor; 2 = no
  has <- ifelse(!is.na(prev) & prev %in% c(1,2,3), 1L,
                ifelse(!is.na(ever) & ever == 1, 1L,
                       ifelse(!is.na(ever) & ever == 2, 0L, NA_integer_)))
  dis[[dis_labels[i]]] <- has
}
dis$n_chronic <- rowSums(dis[dis_labels], na.rm = TRUE)
dis$multimorbidity <- ifelse(rowSums(!is.na(dis[dis_labels])) == 0, NA,
                             as.integer(dis$n_chronic >= 2))

## ---- 6. ADL / IADL ----
adl_items <- paste0("db", sprintf("%03d", c(1,3,5,7,9,11)))
iadl_items <- paste0("db", sprintf("%03d", c(12,14,16,18,20,22)))
adl_x <- as.data.frame(lapply(adl_items, function(v) {
  x <- to_na(hs[[v]]); ifelse(is.na(x), NA, as.integer(x >= 2))   # ≥2 indicates difficulty
}))
iadl_x <- as.data.frame(lapply(iadl_items, function(v) {
  x <- to_na(hs[[v]]); ifelse(is.na(x), NA, as.integer(x >= 2))
}))
names(adl_x) <- adl_items; names(iadl_x) <- iadl_items
adl_count  <- rowSums(adl_x,  na.rm = TRUE);  adl_count[rowSums(!is.na(adl_x))  == 0] <- NA
iadl_count <- rowSums(iadl_x, na.rm = TRUE); iadl_count[rowSums(!is.na(iadl_x)) == 0] <- NA
adl_any  <- ifelse(is.na(adl_count),  NA, as.integer(adl_count  >= 1))
iadl_any <- ifelse(is.na(iadl_count), NA, as.integer(iadl_count >= 1))

## ---- 7. pain / sleep / lifestyle ----
pain <- to_na(hs$da027)
pain_any <- ifelse(is.na(pain), NA, as.integer(pain >= 2))
# three pain categories(2026-09-10 decision: dichotomization masks J shaped dose response): 1 none / 2-3 mild / 4-5 severe
pain_3cat <- ifelse(is.na(pain), NA,
                    ifelse(pain == 1, 0L, ifelse(pain %in% c(2, 3), 1L, 2L)))
sleep_h  <- to_na(hs$da030); sleep_h[!is.na(sleep_h) & (sleep_h < 0 | sleep_h > 24)] <- NA
nap_min  <- to_na(hs$da031); nap_min[!is.na(nap_min) & (nap_min < 0 | nap_min > 600)] <- NA
sleep_cat <- cut(sleep_h, breaks = c(-Inf, 5.999, 8.0001, Inf),
                 labels = c("<6h","6-8h",">8h"), right = FALSE)
# proxy for sleep quality: CES-D items dc022"restless sleep"(1-4), take 3-4(per week ≥3-4 days)indicates restless sleep
# note: dc022 also include CES-D total score, when jointly included in regression explain in SAP explain in
slp_restless_x <- to_na(hs$dc022)
sleep_restless <- ifelse(is.na(slp_restless_x), NA, as.integer(slp_restless_x >= 3))

## ---- 7b. smoking / alcohol use / physical activity ----
# smoking: zsmoke==1(smoked in preceding wave)participants skip this question at current wave da046, direct response da047(1 current smoker / 2 former smoker / 3 never); 
#       others first answer da046(1 ever smoked / 2 never), ever smokers then answer da047.using da047 is final status.
da046v <- to_na(hs$da046); da047v <- to_na(hs$da047)
smoke <- rep(NA_integer_, nrow(hs))
smoke[!is.na(da047v) & da047v == 3] <- 0L          # never
smoke[!is.na(da046v) & da046v == 2] <- 0L          # never
smoke[!is.na(da047v) & da047v == 2] <- 1L          # former smoker
smoke[!is.na(da047v) & da047v == 1] <- 2L          # current smoking
# alcohol use(past year, da051): 1 ≥once per month / 2 less than monthly / 3 no alcohol use
drink_freq <- to_na(hs$da051)
drink_any  <- ifelse(is.na(drink_freq), NA, as.integer(drink_freq %in% c(1, 2)))
# physical activity(da032, per occasion ≥10 minutes): 1 yes / 2 no
pa_vigorous <- { x <- to_na(hs$da032_1_); ifelse(is.na(x), NA, as.integer(x == 1)) }
pa_moderate <- { x <- to_na(hs$da032_2_); ifelse(is.na(x), NA, as.integer(x == 1)) }
pa_light    <- { x <- to_na(hs$da032_3_); ifelse(is.na(x), NA, as.integer(x == 1)) }
pa_mat <- cbind(pa_vigorous, pa_moderate, pa_light)
pa_any <- ifelse(rowSums(!is.na(pa_mat)) == 0, NA,
                 as.integer(rowSums(pa_mat == 1, na.rm = TRUE) > 0))

## ---- 7c. BMI: exclude ----
# 2020 core module lacks height and weight measurements; Harmonized D most recent measured values from 2015 year(r3mbmi), 
# across 5 combining years introduces temporal mismatch and nonrandom missingness.as decided(2026-09-10): BMI not included in this study, 
# document data unavailability as a limitation.

## ---- 8. demographics ----
urba <- to_na(demo$ba008)      # 1 city center/town center; 2 urban rural fringe; 3 village; 4 special area
hukou <- to_na(demo$ba009)     # 1 agricultural; 2 nonagricultural; 3 resident household registration; 4 no household registration
df <- data.frame(
  ID          = hs$ID,
  # cognition
  imrc        = imrc1, imrc2 = imrc2, imrc3 = imrc3, imrc_mean = imrc_mean,
  dlrc        = dlrc, memory = memory, orient = orient,
  ser7        = ser7, draw = draw, mental = mental, global = global,
  cog_items_avail = cog_items_avail, cog_complete = cog_complete,
  # psychological/function
  cesd10 = cesd10, cesd_avail = cesd_avail, depression = dep,
  adl_count = adl_count, adl_any = adl_any, iadl_count = iadl_count, iadl_any = iadl_any,
  pain = pain, pain_any = pain_any, pain_3cat = pain_3cat,
  sleep_h = sleep_h, nap_min = nap_min, sleep_cat = sleep_cat,
  sleep_restless = sleep_restless,
  # lifestyle
  smoke = smoke, drink_freq = drink_freq, drink_any = drink_any,
  pa_vigorous = pa_vigorous, pa_moderate = pa_moderate, pa_light = pa_light, pa_any = pa_any,
  stringsAsFactors = FALSE
)
df <- cbind(df, dis[dis_labels], n_chronic = dis$n_chronic, multimorbidity = dis$multimorbidity)

demo_small <- data.frame(
  ID = demo$ID, age = demo$xrage,
  female = ifelse(is.na(demo$xrgender), NA, as.integer(demo$xrgender == 2)),
  # education(2026-09-12 revision): zredu and ba010 at 2020 fully complementary within wave——
  #   zredu nonmissing 17,687(carried forward from previous wave), ba010 nonmissing 1,708(newly asked this wave), no overlap, no gaps, total 19,395.
  #   education after combining 0 missing, therefore no cross wave backfilling required; both variables use the same 11 level coding(verified against birth cohort distributions).
  edu = dplyr::coalesce(to_na(demo$zredu), to_na(demo$ba010)),
  marital = to_na(demo$ba011),
  partner = to_na(demo$xrpartner),
  urban = ifelse(is.na(urba), NA, ifelse(urba == 1, 1L, ifelse(urba == 3, 0L, NA_integer_))),
  urban3 = ifelse(is.na(urba), NA, ifelse(urba %in% c(1,2), 1L, ifelse(urba == 3, 0L, NA_integer_))),
  hukou = hukou,
  proxy = ifelse(is.na(demo$proxy_2), NA, as.integer(demo$proxy_2 == 1)),
  communityID = demo$communityID, householdID = demo$householdID,
  stringsAsFactors = FALSE
)

dat <- df %>% left_join(demo_small, by = "ID") %>% left_join(wt, by = "ID")

## ---- 9. analysis sample definition(v3: primary analysis = cognition 9 complete components) ----
n0 <- nrow(dat)
dat <- dat %>% mutate(
  edu_group = case_when(
    edu %in% 1:3 ~ "Less than primary",
    edu == 4     ~ "Primary school",
    edu == 5     ~ "Middle school",
    edu >= 6     ~ "High school or above",
    TRUE         ~ NA_character_),
  eligible = !is.na(age) & age >= 45,
  nonproxy = is.na(proxy) | proxy == 0
) %>% filter(eligible, nonproxy)
n1 <- nrow(dat)
dat <- dat %>% mutate(primary = as.integer(cog_complete == 1))
n2 <- sum(dat$primary, na.rm = TRUE)
n2b <- sum(dat$primary == 1 & !is.na(dat$edu_group), na.rm = TRUE)
say("sample flow: all ", n0, " -> age≥45 and nonproxy ", n1,
    " -> cognition 9 complete components ", n2, "(of whom education is known, enter primary analysis ", n2b,
    "; missing education ", n2 - n2b, "), incomplete cognitive assessment(sensitivity only) ", n1 - n2)

## ---- 10. cognitive outcome definition(cutoff and z scores always reference the primary analysis sample) ----
gm <- mean(dat$global[dat$primary == 1], na.rm = TRUE); gs <- sd(dat$global[dat$primary == 1], na.rm = TRUE)
mm <- mean(dat$memory[dat$primary == 1], na.rm = TRUE); ms <- sd(dat$memory[dat$primary == 1], na.rm = TRUE)
em <- mean(dat$mental[dat$primary == 1], na.rm = TRUE); es <- sd(dat$mental[dat$primary == 1], na.rm = TRUE)
thr_all <- gm - gs
q1_all  <- quantile(dat$global[dat$primary == 1], 0.25, na.rm = TRUE)
thr_edu_tab <- dat %>% filter(primary == 1, !is.na(edu_group)) %>% group_by(edu_group) %>%
  summarise(thr_edu = mean(global, na.rm = TRUE) - sd(global, na.rm = TRUE), .groups = "drop")
dat <- dat %>% left_join(thr_edu_tab, by = "edu_group") %>% mutate(
  cog_z     = (global - gm) / gs,
  cog_z_mem = (memory - mm) / ms,
  cog_z_men = (mental - em) / es,
  # primary definition(scheme B, 2026-09-10 set): education specific cutoff = below mean in this education stratum − 1SD
  cog_low_edu = ifelse(is.na(global) | is.na(thr_edu), NA, as.integer(global < thr_edu)),
  # sensitivity definition A: overall primary sample mean − 1SD
  cog_low_1sd = ifelse(is.na(global), NA, as.integer(global < thr_all)),
  # sensitivity definition C: lowest quartile
  cog_low_q1  = ifelse(is.na(global), NA, as.integer(global <= q1_all))
)
dat$cog_impair <- dat$cog_low_edu          # primary outcome alias
say("cutoff: all participants with complete cognition mean−1SD = ", round(thr_all, 2), "; lowest quartile = ", round(q1_all, 2))
say("education specific cutoff: ", paste(sprintf("%s %.2f", thr_edu_tab$edu_group, thr_edu_tab$thr_edu), collapse = " | "))

## ---- 11. save ----
saveRDS(dat, file.path(out_dir, "analysis_2020.rds"))
write.csv(dat, file.path(out_dir, "analysis_2020.csv"), row.names = FALSE, fileEncoding = "UTF-8")

## ---- 12. build log: key variable distributions ----
say("\n---- key variable distributions(after construction) ----")
say("nonmissing count by component: imrc1=", sum(!is.na(dat$imrc)), " dlrc=", sum(!is.na(dat$dlrc)),
    " orient=", sum(!is.na(dat$orient)), " ser7=", sum(!is.na(dat$ser7)),
    " draw=", sum(!is.na(dat$draw)))
d1 <- dat[dat$primary == 1, ]
say("participants with complete cognition n=", nrow(d1), "(of whom education is known, enter primary analysis ", sum(!is.na(d1$edu_group)),
    "; missing education ", sum(is.na(d1$edu_group)), ").distributions below use participants with complete cognition(n=", nrow(d1), ")")
say("global cognition global(participants with complete cognition): mean=", round(mean(d1$global, na.rm=TRUE),2),
    " sd=", round(sd(d1$global, na.rm=TRUE),2))
say("episodic memory memory: mean=", round(mean(d1$memory, na.rm=TRUE),2),
    " | immediate (trial 1) mean=", round(mean(d1$imrc, na.rm=TRUE),2),
    " | delayed mean=", round(mean(d1$dlrc, na.rm=TRUE),2))
say("mental status mental: mean=", round(mean(d1$mental, na.rm=TRUE),2),
    " | orientation mean=", round(mean(d1$orient, na.rm=TRUE),2),
    " | serial subtraction of 7 mean=", round(mean(d1$ser7, na.rm=TRUE),2),
    " | figure copying mean=", round(mean(d1$draw, na.rm=TRUE),2))
say("CES-D-10: n=", sum(!is.na(d1$cesd10)), " mean=", round(mean(d1$cesd10, na.rm=TRUE),2),
    " | depression (≥10) prevalence=", round(100*mean(d1$depression, na.rm=TRUE),1), "%")
say("pain (any): ", round(100*mean(d1$pain_any, na.rm=TRUE),1), "%",
    " | three categories [none/mild/severe]: ", sum(d1$pain_3cat==0, na.rm=TRUE), " / ",
    sum(d1$pain_3cat==1, na.rm=TRUE), " / ", sum(d1$pain_3cat==2, na.rm=TRUE),
    " | sleep duration mean=", round(mean(d1$sleep_h, na.rm=TRUE),2), "h")
say("ADL limited: ", round(100*mean(d1$adl_any, na.rm=TRUE),1), "%",
    " | IADL limited: ", round(100*mean(d1$iadl_any, na.rm=TRUE),1), "%")
say("smoking: never ", sum(d1$smoke==0, na.rm=TRUE), " / former smoker ", sum(d1$smoke==1, na.rm=TRUE),
    " / current ", sum(d1$smoke==2, na.rm=TRUE), " / missing ", sum(is.na(d1$smoke)))
say("alcohol use (past year): ≥monthly ", sum(d1$drink_freq==1, na.rm=TRUE), " / <monthly ", sum(d1$drink_freq==2, na.rm=TRUE),
    " / nondrinker ", sum(d1$drink_freq==3, na.rm=TRUE), " / missing ", sum(is.na(d1$drink_freq)))
say("physical activity (any,≥10 minutes/times): ", round(100*mean(d1$pa_any, na.rm=TRUE),1), "%",
    " | high intensity: ", round(100*mean(d1$pa_vigorous, na.rm=TRUE),1), "%")
say("restless sleep (≥3-4 days/week): ", round(100*mean(d1$sleep_restless, na.rm=TRUE),1), "%")
say("cognitive impairment(participants with complete cognition; primary definition among participants with known education): primary definition [education specific cutoff] ", round(100*mean(d1$cog_low_edu, na.rm=TRUE),1), "%",
    " | sensitivity A[overall mean−1SD] ", round(100*mean(d1$cog_low_1sd, na.rm=TRUE),1), "%",
    " | sensitivity C[lowest quartile] ", round(100*mean(d1$cog_low_q1, na.rm=TRUE),1), "%")
say("participants with incomplete assessment(n=", nrow(dat) - nrow(d1), ")if including: prevalence under primary definition ",
    round(100*mean(dat$cog_low_edu, na.rm=TRUE),1), "%(for sensitivity analyses)")
say("weights available (INDV_weight nonmissing): ", sum(!is.na(dat$INDV_weight)), " / ", nrow(dat))

close(log_con)
cat("\n complete.output: output/01_2020_single_wave_baseline/analysis_2020.rds\n")
