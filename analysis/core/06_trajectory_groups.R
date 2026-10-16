# ============================================================================
# 06_trajectory_groups.R —— 认知轨迹分组（探索性；GBTM 的过渡版本）
# 版本：v1.1（2026-09-13）｜ 依赖：cluster（silhouette）、nnet
# 输入：output/04_纵向队列/cohort_change_2026-09-13.rds（年龄合格固定队列 5,672 人 × 5 波）
# 输出：output/06_轨迹模型/trajectory_groups_2026-09-19.csv、trajectory_classes_2026-09-19.csv、
#       trajectory_predictors_2026-09-19.csv、trajectory_groups_2026-09-19.rds
#      （2026-09-19 版本号：脚本已在当日数据上复跑核对，输出与旧版逐字节一致）
# 方法：对"相对自身基线的变化轨迹"（0, d2013, d2015, d2018, d2020）做 k-means，
#       k=2–5 按平均轮廓宽度选最优；认知分先做波内标准化以消除 2018 期量表差异。
# 说明：正式的 GBTM（lcmm::hlme）需要更充分的初值搜索，列为下一步；
#       本脚本给的是**可复现的探索性轨迹分组**，报告时须明确标注。
# ============================================================================
invisible(try(Sys.setlocale("LC_ALL", "zh_CN.UTF-8"), silent = TRUE))
if (dir.exists("D:/R/library")) .libPaths(c("D:/R/library", .libPaths()))
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(cluster); library(nnet) })

P <- Sys.getenv("CHARLS_PROJ", "D:/BaiduSyncdisk/科研/徐中豪")
OUT <- file.path(P, "output", "06_轨迹模型")
coh <- readRDS(file.path(P, "output", "04_纵向队列", "cohort_change_2026-09-13.rds"))$cohort

## 1. 波内标准化 + 变化形状 ----
lat <- coh %>% select(KEY, year, ms10) %>% group_by(year) %>%
  mutate(z = as.numeric(scale(ms10))) %>% ungroup() %>%
  select(KEY, year, z) %>% pivot_wider(names_from = year, values_from = z, names_prefix = "t")
lat <- lat[complete.cases(lat), ]
shape <- data.frame(d0 = 0, d1 = lat$t2013 - lat$t2011, d2 = lat$t2015 - lat$t2011,
                    d3 = lat$t2018 - lat$t2011, d4 = lat$t2020 - lat$t2011)
cat("纳入轨迹分析个体：", nrow(lat), fill = TRUE)

## 2. 选 k（平均轮廓宽度） ----
set.seed(20260913)
X <- as.matrix(shape)
sil <- sapply(2:5, function(k) {
  km <- kmeans(X, centers = k, nstart = 50, iter.max = 100)
  mean(silhouette(km$cluster, dist(X))[, 3])
})
cat("平均轮廓宽度（k=2..5）：", paste(round(sil, 3), collapse = "  "), fill = TRUE)
k_best <- (2:5)[which.max(sil)]
cat("选定类别数 k =", k_best, fill = TRUE)
km <- kmeans(X, centers = k_best, nstart = 50, iter.max = 100)
lat$class <- km$cluster

## 3. 类别轨迹与命名（按 2020 年相对变化由差到好） ----
tr <- lat %>% group_by(class) %>%
  summarise(n = n(), z2011 = mean(t2011), z2013 = mean(t2013), z2015 = mean(t2015),
            z2018 = mean(t2018), z2020 = mean(t2020),
            delta = mean(t2020 - t2011), .groups = "drop") %>%
  mutate(across(starts_with("z"), ~ round(.x, 2)), delta = round(delta, 2)) %>% arrange(delta)
nm <- c("快速下降组", "中度下降组", "轻度下降组", "基本稳定组", "改善组")[seq_len(nrow(tr))]
tr$group <- nm
lat <- lat %>% left_join(tr %>% select(class, group), by = "class")
cat("", "=== 轨迹类别（波内 z 均值；delta = 2020 − 2011）===", fill = TRUE)
print(as.data.frame(tr), row.names = FALSE)

## 4. 基线特征 × 类别 ----
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
cat("", "=== 各类别基线特征（2011）===", fill = TRUE)
print(as.data.frame(tab), row.names = FALSE)

## 5. 基线因素 → 轨迹类别 ----
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
cat("", "=== 基线因素 → 轨迹类别（OR，校正年龄/性别/教育）===", fill = TRUE)
print(head(pred, 40), row.names = FALSE)

## 6. 输出 ----
write.csv(tr,   file.path(OUT, "trajectory_groups_2026-09-19.csv"),     row.names = FALSE, fileEncoding = "UTF-8")
write.csv(tab,  file.path(OUT, "trajectory_classes_2026-09-19.csv"),    row.names = FALSE, fileEncoding = "UTF-8")
write.csv(pred, file.path(OUT, "trajectory_predictors_2026-09-19.csv"), row.names = FALSE, fileEncoding = "UTF-8")
saveRDS(list(lat = lat, traj = tr, km = km, sil = sil, base = base, pred = pred),
        file.path(OUT, "trajectory_groups_2026-09-19.rds"))
cat("已存 output/06_轨迹模型/trajectory_*_2026-09-19.csv / .rds", fill = TRUE)
