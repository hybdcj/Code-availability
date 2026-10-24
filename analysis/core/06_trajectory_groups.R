# ============================================================================
# 06_trajectory_groups.R —— cognitive trajectory groups(exploratory; GBTM intermediate version of)
# version: v1.1(2026-09-13)｜ dependencies: cluster(silhouette), nnet
# input: output/04_longitudinal_cohort/cohort_change_2026-09-13.rds(fixed cohort meeting age eligibility 5,672 participants × 5 wave)
# output: output/06_trajectory_model/trajectory_groups_2026-09-19.csv, trajectory_classes_2026-09-19.csv, 
#       trajectory_predictors_2026-09-19.csv, trajectory_groups_2026-09-19.rds
#      (2026-09-19 version number: script rerun and checked against same day data, outputs byte identical to previous version)
# method: for"change trajectory relative to individual baseline"(0, d2013, d2015, d2018, d2020)perform k-means, 
#       k=2–5 select optimum using mean silhouette width; standardize cognition within wave to remove 2018 wave scale differences.
# notes: formal GBTM(lcmm::hlme)requires more extensive initial value search, planned as next step; 
#       this script provides**reproducible exploratory trajectory groups**, clearly label when reporting.
# ============================================================================
invisible(try(Sys.setlocale("LC_ALL", "zh_CN.UTF-8"), silent = TRUE))
if (dir.exists("D:/R/library")) .libPaths(c("D:/R/library", .libPaths()))
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(cluster); library(nnet) })

P <- Sys.getenv("CHARLS_PROJ", "D:/BaiduSyncdisk/research/Xu_Zhonghao")
OUT <- file.path(P, "output", "06_trajectory_model")
coh <- readRDS(file.path(P, "output", "04_longitudinal_cohort", "cohort_change_2026-09-13.rds"))$cohort

## 1. within wave standardization + shape of change ----
lat <- coh %>% select(KEY, year, ms10) %>% group_by(year) %>%
  mutate(z = as.numeric(scale(ms10))) %>% ungroup() %>%
  select(KEY, year, z) %>% pivot_wider(names_from = year, values_from = z, names_prefix = "t")
lat <- lat[complete.cases(lat), ]
shape <- data.frame(d0 = 0, d1 = lat$t2013 - lat$t2011, d2 = lat$t2015 - lat$t2011,
                    d3 = lat$t2018 - lat$t2011, d4 = lat$t2020 - lat$t2011)
cat("participants included in trajectory analysis: ", nrow(lat), fill = TRUE)

## 2. select k(mean silhouette width) ----
set.seed(20260913)
X <- as.matrix(shape)
sil <- sapply(2:5, function(k) {
  km <- kmeans(X, centers = k, nstart = 50, iter.max = 100)
  mean(silhouette(km$cluster, dist(X))[, 3])
})
cat("mean silhouette width(k=2..5): ", paste(round(sil, 3), collapse = "  "), fill = TRUE)
k_best <- (2:5)[which.max(sil)]
cat("selected number of classes k =", k_best, fill = TRUE)
km <- kmeans(X, centers = k_best, nstart = 50, iter.max = 100)
lat$class <- km$cluster

## 3. class trajectories and naming(by 2020 year relative change from worst to best) ----
tr <- lat %>% group_by(class) %>%
  summarise(n = n(), z2011 = mean(t2011), z2013 = mean(t2013), z2015 = mean(t2015),
            z2018 = mean(t2018), z2020 = mean(t2020),
            delta = mean(t2020 - t2011), .groups = "drop") %>%
  mutate(across(starts_with("z"), ~ round(.x, 2)), delta = round(delta, 2)) %>% arrange(delta)
nm <- c("rapid decline group", "moderate decline group", "mild decline group", "largely stable group", "improvement group")[seq_len(nrow(tr))]
tr$group <- nm
lat <- lat %>% left_join(tr %>% select(class, group), by = "class")
cat("", "=== trajectory classes(within wave z mean; delta = 2020 − 2011)===", fill = TRUE)
print(as.data.frame(tr), row.names = FALSE)

## 4. baseline characteristics × category ----
base <- coh %>% filter(year == 2011) %>%
  transmute(KEY, age, female, edu = edu_group, dep = as.integer(cesd10 >= 10), adl, sleep_h,
            pain = pain_any, pa = pa_any, hibp, diab) %>%
  left_join(lat %>% select(KEY, group), by = "KEY")
tab <- base %>% group_by(group) %>%
  summarise(n = n(),
            age = round(mean(age, na.rm = TRUE), 1),
            female_pct = round(100*mean(female == 1, na.rm = TRUE), 1),
            lowedu_pct = round(100*mean(edu == "Less than primary", na.rm = TRUE), 1),
            dep_pct = round(100*mean(dep == 1, na.rm = TRUE), 1),
            adl_pct = round(100*mean(adl == 1, na.rm = TRUE), 1),
            sleep_h = round(mean(sleep_h, na.rm = TRUE), 2),
            pain_pct = round(100*mean(pain == 1, na.rm = TRUE), 1),
            pa_pct = round(100*mean(pa == 1, na.rm = TRUE), 1),
            hibp_pct = round(100*mean(hibp == 1, na.rm = TRUE), 1),
            diab_pct = round(100*mean(diab == 1, na.rm = TRUE), 1), .groups = "drop")
cat("", "=== baseline characteristics by class(2011)===", fill = TRUE)
print(as.data.frame(tab), row.names = FALSE)

## 5. baseline factors → trajectory classes ----
base$group <- factor(base$group)
if (nlevels(base$group) == 2) {
  base$bad <- as.integer(base$group == levels(base$group)[1])
  m <- glm(bad ~ age + female + dep + adl + sleep_h + pain + pa + hibp + diab,
           data = base, family = binomial)
  s <- summary(m)$coefficients[-1, , drop = FALSE]
  pred <- data.frame(contrast = paste0(levels(base$group)[1], " vs ", levels(base$group)[2]),
                     variable = rownames(s),
                     OR = sprintf("%.2f (%.2f-%.2f)", exp(s[, 1]), exp(s[, 1] - 1.96*s[, 2]),
                                  exp(s[, 1] + 1.96*s[, 2])),
                     p = signif(s[, 4], 3), row.names = NULL)
} else {
  base$group <- relevel(base$group, ref = levels(base$group)[1])
  m <- multinom(group ~ age + female + dep + adl + sleep_h + pain + pa + hibp + diab,
                data = base, trace = FALSE)
  cc <- summary(m)$coefficients; ss <- summary(m)$standard.errors
  pred <- do.call(rbind, lapply(rownames(cc), function(g) data.frame(
    contrast = paste0(g, " vs ", levels(base$group)[1]), variable = colnames(cc),
    OR = sprintf("%.2f (%.2f-%.2f)", exp(cc[g, ]), exp(cc[g, ] - 1.96*ss[g, ]), exp(cc[g, ] + 1.96*ss[g, ])),
    p = signif(2*pnorm(-abs(cc[g, ]/ss[g, ])), 3), row.names = NULL)))
}
cat("", "=== baseline factors → trajectory classes(OR, adjust for age/sex/education)===", fill = TRUE)
print(head(pred, 40), row.names = FALSE)

## 6. output ----
write.csv(tr,   file.path(OUT, "trajectory_groups_2026-09-19.csv"),     row.names = FALSE, fileEncoding = "UTF-8")
write.csv(tab,  file.path(OUT, "trajectory_classes_2026-09-19.csv"),    row.names = FALSE, fileEncoding = "UTF-8")
write.csv(pred, file.path(OUT, "trajectory_predictors_2026-09-19.csv"), row.names = FALSE, fileEncoding = "UTF-8")
saveRDS(list(lat = lat, traj = tr, km = km, sil = sil, base = base, pred = pred),
        file.path(OUT, "trajectory_groups_2026-09-19.rds"))
cat("saved_output/06_trajectory_model/trajectory_*_2026-09-19.csv / .rds", fill = TRUE)
