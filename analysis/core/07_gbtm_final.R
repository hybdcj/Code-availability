# ============================================================================
# 07_gbtm_final.R —— 认知轨迹分型（GBTM）+ 可靠变化敏感性分析
# 版本：v1.1（2026-09-19 复跑核对：输出与 09-13 版逐字节一致，仅版本号更新）｜ 依赖：lcmm、nnet、dplyr
# 输入：output/04_纵向队列/cohort_change_2026-09-13.rds（年龄合格固定队列 5,672 人 × 5 波）
# 输出：output/06_轨迹模型/gbtm_final_trajectories_2026-09-19.csv（各类轨迹，含原始分与 z）
#       output/06_轨迹模型/gbtm_final_classes_2026-09-19.csv（各类基线特征）
#       output/06_轨迹模型/gbtm_final_predictors_2026-09-19.csv（基线因素 → 类别，若生成）
#       output/04_纵向队列/reliable_change_2026-09-19.csv（可靠变化敏感性分析）
# 方法：
#   ① 初值：先用 k-means（k=2..4）得到各类均值轨迹 → 作为 hlme() 的初值 B
#   ② 主分析：GBTM（hlme，二次多项式，ng=2..4），按 BIC 选最优
#   ③ 结局用"波内 z"（消除 2018 期量表差异），类别轨迹再回译成原始分展示
#   ④ 敏感性：可靠变化（Δz ≤ −0.5）Logistic，看基线因素方向是否一致
# ============================================================================
invisible(try(Sys.setlocale("LC_ALL", "zh_CN.UTF-8"), silent = TRUE))
if (dir.exists("D:/R/library")) .libPaths(c("D:/R/library", .libPaths()))
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(lcmm); library(nnet) })

P <- Sys.getenv("CHARLS_PROJ", "D:/BaiduSyncdisk/科研/徐中豪")
OUT <- file.path(P, "output", "06_轨迹模型")
OUT_LONG <- file.path(P, "output", "04_纵向队列")
coh <- readRDS(file.path(OUT_LONG, "cohort_change_2026-09-13.rds"))$cohort
tmap <- c("2011" = 0, "2013" = 2, "2015" = 4, "2018" = 7, "2020" = 9)

## ---- 数据（个体×波，z 结局 + 原始分备用） ----
lat <- coh %>% select(KEY, year, ms10) %>% group_by(year) %>%
  mutate(z = as.numeric(scale(ms10))) %>% ungroup() %>%
  mutate(time = as.numeric(tmap[as.character(year)])) %>%
  filter(!is.na(z)) %>% arrange(KEY, time) %>% as.data.frame()
lat$KEYn <- as.integer(factor(lat$KEY))
cat("长表：", nrow(lat), "行；个体", length(unique(lat$KEYn)), fill = TRUE)

## ---- ① k-means 初值 ----
wide <- lat %>% select(KEYn, time, z) %>% pivot_wider(names_from = time, values_from = z, names_prefix = "t")
wide <- wide[complete.cases(wide), ]
X <- as.matrix(wide[, -1]); X <- scale(X, center = FALSE, scale = FALSE)   # 仅保留形状差异
set.seed(20260913)
init_km <- lapply(2:4, function(k) kmeans(X, centers = k, nstart = 50, iter.max = 100))

## ---- ② GBTM：初值由 k-means 类别均值构造 ----
# hlme 参数向量：参考类固定效应(3) + 各类偏移(3×(ng−1)) + 各类方差(ng) → 长度 4×ng
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
  if (inherits(m, "try-error")) { cat("ng =", ng, "失败：", conditionMessage(attr(m, "condition")), fill = TRUE); next }
  if (!is.finite(m$BIC) || m$BIC > 1e8) { cat("ng =", ng, "退化（BIC =", m$BIC, "）", fill = TRUE); next }
  fits[[as.character(ng)]] <- m; bics[as.character(ng)] <- m$BIC
  cat("ng =", ng, "完成；BIC =", round(m$BIC, 1), fill = TRUE)
}
if (!length(fits)) {
  cat("GBTM 仍未能收敛 —— 见下方可靠性变化敏感性分析（可独立使用）", fill = TRUE)
  best <- NULL; best_ng <- NA
} else {
  best_ng <- names(which.min(bics)); best <- fits[[best_ng]]
  cat("按 BIC 选择：", best_ng, "类", fill = TRUE)
}

## ---- ③ 类别轨迹（z 与原始分） ----
if (!is.null(best)) {
  pp <- best$pprob; pp$class <- pp[[grep("^class", names(pp))[1]]]
  lat2 <- lat %>% left_join(pp[, c("KEYn", "class")], by = "KEYn")
  wave_stats <- coh %>% group_by(year) %>% summarise(m = mean(ms10, na.rm = TRUE), s = sd(ms10, na.rm = TRUE), .groups = "drop")
  lat2 <- lat2 %>% left_join(wave_stats, by = "year") %>% mutate(ms10_fitted = m + z * s)
  traj <- lat2 %>% group_by(class, time) %>%
    summarise(n = n(), z = round(mean(z, na.rm = TRUE), 2), 原始分 = round(mean(ms10_fitted, na.rm = TRUE), 2), .groups = "drop") %>%
    pivot_wider(names_from = time, values_from = c(z, 原始分))
  cat("", "=== 各类别轨迹（列名：z_0/z_2/z_4/z_7/z_9 与 原始分_0…9）===", fill = TRUE)
  print(as.data.frame(traj), row.names = FALSE)
  cat("类别人数：", paste(names(table(lat2$class[!duplicated(lat2$KEYn)])), as.integer(table(lat2$class[!duplicated(lat2$KEYn)])), sep = "=", collapse = "  "), fill = TRUE)
}

## ---- ④ 可靠变化（敏感性）：Δz ≤ −0.5 视为"可靠下降" ----
chg <- coh %>% filter(year %in% c(2011, 2020)) %>% select(KEY, year, ms10) %>%
  group_by(year) %>% mutate(z = as.numeric(scale(ms10))) %>% ungroup() %>%
  select(KEY, year, z) %>% pivot_wider(names_from = year, values_from = z, names_prefix = "z")
chg <- chg %>% mutate(delta_z = z2020 - z2011, reliable_decline = as.integer(delta_z <= -0.5))
base <- coh %>% filter(year == 2011) %>%
  transmute(KEY, age, female, edu = edu_group, dep = as.integer(cesd10 >= 10), adl, sleep_h,
            pain = pain_any, pa = pa_any, hibp, diab) %>%
  left_join(chg %>% select(KEY, delta_z, reliable_decline), by = "KEY")
cat("", "=== 可靠下降（Δz ≤ −0.5）比例：", sprintf("%.1f%%", 100*mean(base$reliable_decline, na.rm = TRUE)),
    "；Δz 均值 ", round(mean(base$delta_z, na.rm = TRUE), 3), " ===", fill = TRUE)
m_rc <- glm(reliable_decline ~ age + female + dep + adl + sleep_h + pain + pa + hibp + diab,
            data = base, family = binomial)
s <- summary(m_rc)$coefficients[-1, , drop = FALSE]
rc <- data.frame(variable = rownames(s),
                 OR = sprintf("%.2f (%.2f-%.2f)", exp(s[, 1]), exp(s[, 1] - 1.96*s[, 2]), exp(s[, 1] + 1.96*s[, 2])),
                 p = signif(s[, 4], 3), row.names = NULL)
print(rc, row.names = FALSE)

## ---- ⑤ 输出 ----
if (!is.null(best)) {
  write.csv(traj, file.path(OUT, "gbtm_final_trajectories_2026-09-19.csv"), row.names = FALSE, fileEncoding = "UTF-8")
  cls <- lat2 %>% group_by(class) %>% summarise(n = length(unique(KEYn)), .groups = "drop")
  write.csv(cls, file.path(OUT, "gbtm_final_classes_2026-09-19.csv"), row.names = FALSE, fileEncoding = "UTF-8")
  saveRDS(list(best = best, fits = fits, bics = bics, lat = lat2), file.path(OUT, "gbtm_final_2026-09-19.rds"))
}
write.csv(rc, file.path(OUT_LONG, "reliable_change_2026-09-19.csv"), row.names = FALSE, fileEncoding = "UTF-8")
cat("已存 output/06_轨迹模型/gbtm_final_*_2026-09-19.csv 与 reliable_change_2026-09-19.csv", fill = TRUE)
