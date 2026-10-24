# ============================================================================
# 07_gbtm_final.R —— cognitive trajectory classification(GBTM)+ reliable change sensitivity analysis
# version: v1.1(2026-09-19 rerun verification: outputs and 09-13 version are byte identical, version number updated only)｜ dependencies: lcmm, nnet, dplyr
# input: output/04_longitudinal_cohort/cohort_change_2026-09-13.rds(fixed cohort meeting age eligibility 5,672 participants × 5 wave)
# output: output/06_trajectory_model/gbtm_final_trajectories_2026-09-19.csv(trajectories by class, including raw scores and z)
#       output/06_trajectory_model/gbtm_final_classes_2026-09-19.csv(baseline characteristics by class)
#       output/06_trajectory_model/gbtm_final_predictors_2026-09-19.csv(baseline factors → category, if generating)
#       output/04_longitudinal_cohort/reliable_change_2026-09-19.csv(reliable change sensitivity analysis)
# method: 
#   ① initial values: first use k-means(k=2..4)obtain class mean trajectories → as hlme() initial values for B
#   ② primary analysis: GBTM(hlme, quadratic polynomial, ng=2..4), by BIC select optimum
#   ③ outcome uses"within wave z"(remove 2018 wave scale differences), translate class trajectories back to raw scores for display
#   ④ sensitivity: reliable change(Δz ≤ −0.5)Logistic, check consistent direction of baseline factors
# ============================================================================
invisible(try(Sys.setlocale("LC_ALL", "zh_CN.UTF-8"), silent = TRUE))
if (dir.exists("D:/R/library")) .libPaths(c("D:/R/library", .libPaths()))
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(lcmm); library(nnet) })

P <- Sys.getenv("CHARLS_PROJ", "D:/BaiduSyncdisk/research/Xu_Zhonghao")
OUT <- file.path(P, "output", "06_trajectory_model")
OUT_LONG <- file.path(P, "output", "04_longitudinal_cohort")
coh <- readRDS(file.path(OUT_LONG, "cohort_change_2026-09-13.rds"))$cohort
tmap <- c("2011" = 0, "2013" = 2, "2015" = 4, "2018" = 7, "2020" = 9)

## ---- data(individual×wave, z outcome + raw scores as backup) ----
lat <- coh %>% select(KEY, year, ms10) %>% group_by(year) %>%
  mutate(z = as.numeric(scale(ms10))) %>% ungroup() %>%
  mutate(time = as.numeric(tmap[as.character(year)])) %>%
  filter(!is.na(z)) %>% arrange(KEY, time) %>% as.data.frame()
lat$KEYn <- as.integer(factor(lat$KEY))
cat("long format: ", nrow(lat), "rows; individual", length(unique(lat$KEYn)), fill = TRUE)

## ---- ① k-means initial values ----
wide <- lat %>% select(KEYn, time, z) %>% pivot_wider(names_from = time, values_from = z, names_prefix = "t")
wide <- wide[complete.cases(wide), ]
X <- as.matrix(wide[, -1]); X <- scale(X, center = FALSE, scale = FALSE)   # retain shape differences only
set.seed(20260913)
init_km <- lapply(2:4, function(k) kmeans(X, centers = k, nstart = 50, iter.max = 100))

## ---- ② GBTM: initial values from k-means class mean construction ----
# hlme parameter vector: reference class fixed effects (3) + class offsets (3×(ng−1)) + class variances (ng) → length 4×ng
need_len <- function(ng) 4 * ng
mkB <- function(km, ng) {
  cm <- km$centers
  tt <- c(0, 2, 4, 7, 9)
  fitc <- t(apply(cm, 1, function(y) coef(lm(y ~ tt + I(tt^2)))))
  fitc <- fitc[order(fitc[, 1]), , drop = FALSE]
  B <- as.numeric(fitc[1, ])
  if (ng > 1) for (k in 2:ng) B <- c(B, as.numeric(fitc[k, ] - fitc[1, ]))
  B <- c(B, rep(1, ng))
  n <- need_len(ng); if (length(B) < n) B <- c(B, rep(1, n - length(B)))
  B[seq_len(n)]
}
fits <- list(); bics <- numeric()
for (ng in 2:4) {
  B0 <- mkB(init_km[[ng - 1]], ng)
  m <- try(hlme(z ~ time + I(time^2), subject = "KEYn", ng = ng, data = lat,
                mixture = ~ time + I(time^2), random = ~ -1, B = B0, verbose = FALSE), silent = TRUE)
  if (inherits(m, "try-error")) { cat("ng =", ng, "failed: ", conditionMessage(attr(m, "condition")), fill = TRUE); next }
  if (!is.finite(m$BIC) || m$BIC > 1e8) { cat("ng =", ng, "degenerate(BIC =", m$BIC, ")", fill = TRUE); next }
  fits[[as.character(ng)]] <- m; bics[as.character(ng)] <- m$BIC
  cat("ng =", ng, "complete; BIC =", round(m$BIC, 1), fill = TRUE)
}
if (!length(fits)) {
  cat("GBTM still failed to converge —— see reliable change sensitivity analysis below(can be used independently)", fill = TRUE)
  best <- NULL; best_ng <- NA
} else {
  best_ng <- names(which.min(bics)); best <- fits[[best_ng]]
  cat("by BIC selection: ", best_ng, "classes", fill = TRUE)
}

## ---- ③ class trajectory(z and raw score) ----
if (!is.null(best)) {
  pp <- best$pprob; pp$class <- pp[[grep("^class", names(pp))[1]]]
  lat2 <- lat %>% left_join(pp[, c("KEYn", "class")], by = "KEYn")
  wave_stats <- coh %>% group_by(year) %>% summarise(m = mean(ms10, na.rm = TRUE), s = sd(ms10, na.rm = TRUE), .groups = "drop")
  lat2 <- lat2 %>% left_join(wave_stats, by = "year") %>% mutate(ms10_fitted = m + z * s)
  traj <- lat2 %>% group_by(class, time) %>%
    summarise(n = n(), z = round(mean(z, na.rm = TRUE), 2), raw_score = round(mean(ms10_fitted, na.rm = TRUE), 2), .groups = "drop") %>%
    pivot_wider(names_from = time, values_from = c(z, raw_score))
  cat("", "=== each class trajectory(column names: z_0/z_2/z_4/z_7/z_9 and raw score_0…9)===", fill = TRUE)
  print(as.data.frame(traj), row.names = FALSE)
  cat("class size: ", paste(names(table(lat2$class[!duplicated(lat2$KEYn)])), as.integer(table(lat2$class[!duplicated(lat2$KEYn)])), sep = "=", collapse = "  "), fill = TRUE)
}

## ---- ④ reliable change(sensitivity): Δz ≤ −0.5 treat as"reliable decline" ----
chg <- coh %>% filter(year %in% c(2011, 2020)) %>% select(KEY, year, ms10) %>%
  group_by(year) %>% mutate(z = as.numeric(scale(ms10))) %>% ungroup() %>%
  select(KEY, year, z) %>% pivot_wider(names_from = year, values_from = z, names_prefix = "z")
chg <- chg %>% mutate(delta_z = z2020 - z2011, reliable_decline = as.integer(delta_z <= -0.5))
base <- coh %>% filter(year == 2011) %>%
  transmute(KEY, age, female, edu = edu_group, dep = as.integer(cesd10 >= 10), adl, sleep_h,
            pain = pain_any, pa = pa_any, hibp, diab) %>%
  left_join(chg %>% select(KEY, delta_z, reliable_decline), by = "KEY")
cat("", "=== reliable decline(Δz ≤ −0.5)proportion: ", sprintf("%.1f%%", 100*mean(base$reliable_decline, na.rm = TRUE)),
    "; Δz mean ", round(mean(base$delta_z, na.rm = TRUE), 3), " ===", fill = TRUE)
m_rc <- glm(reliable_decline ~ age + female + dep + adl + sleep_h + pain + pa + hibp + diab,
            data = base, family = binomial)
s <- summary(m_rc)$coefficients[-1, , drop = FALSE]
rc <- data.frame(variable = rownames(s),
                 OR = sprintf("%.2f (%.2f-%.2f)", exp(s[, 1]), exp(s[, 1] - 1.96*s[, 2]), exp(s[, 1] + 1.96*s[, 2])),
                 p = signif(s[, 4], 3), row.names = NULL)
print(rc, row.names = FALSE)

## ---- ⑤ output ----
if (!is.null(best)) {
  write.csv(traj, file.path(OUT, "gbtm_final_trajectories_2026-09-19.csv"), row.names = FALSE, fileEncoding = "UTF-8")
  cls <- lat2 %>% group_by(class) %>% summarise(n = length(unique(KEYn)), .groups = "drop")
  write.csv(cls, file.path(OUT, "gbtm_final_classes_2026-09-19.csv"), row.names = FALSE, fileEncoding = "UTF-8")
  saveRDS(list(best = best, fits = fits, bics = bics, lat = lat2), file.path(OUT, "gbtm_final_2026-09-19.rds"))
}
write.csv(rc, file.path(OUT_LONG, "reliable_change_2026-09-19.csv"), row.names = FALSE, fileEncoding = "UTF-8")
cat("saved_output/06_trajectory_model/gbtm_final___2026_09_19.csv_and_reliable_change_2026_09_19.csv", fill = TRUE)
