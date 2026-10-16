# ============================================================================
# 18_screening_index_v3.R —— A07：构建认知脆弱性筛查指数（v3 冻结版）
# 依据：SAP_v3_2026-09-26.md §7（主模型 4 项 + 整数化规则 + 风险分层）
# 变量：年龄、教育（低/中/高）、抑郁（CES-D-10 ≥10）、ADL（任一受限），均为 2011 基线
# 结局：impaired_2020（2011 正常者的 2020 终点评定）；次要：首次受损 Cox
# 主分析不加权；加权（2011 基线权重 + PSU=communityID）作敏感性
# 生成：2026-09-26
# ============================================================================
invisible(try(Sys.setlocale("LC_ALL", "zh_CN.UTF-8"), silent = TRUE))
.libPaths(c("D:/R/library", .libPaths()))
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(survival); library(survey); library(pROC) })
P <- "D:/BaiduSyncdisk/科研/徐中豪"; OUT <- file.path(P, "output", "04_纵向队列")
DIR <- file.path(P, "output", "10_筛查指数"); if (!dir.exists(DIR)) dir.create(DIR, recursive = TRUE)

coh <- readRDS(file.path(OUT, "cohort_change_2026-09-13.rds"))$cohort
inc <- readRDS(file.path(OUT, "interval_incidence_events_2026-09-20.rds"))

b11 <- coh %>% filter(year == 2011) %>%
  transmute(KEY, age, female, edu_group, dep = as.integer(cesd10 >= 10), adl,
            ms10_2011 = ms10, weight, communityID)
d <- inc$base %>% inner_join(b11, by = "KEY")
yr <- c("2013" = 2, "2015" = 4, "2018" = 7, "2020" = 9)
im <- cbind(d$im2013, d$im2015, d$im2018, d$im2020)
d$time  <- apply(im, 1, function(z) { i <- which(z == 1); if (length(i)) yr[i[1]] else 9 })
d$event <- as.integer(apply(im, 1, function(z) any(z == 1, na.rm = TRUE)))
d$edu_group <- factor(d$edu_group, levels = c("Less than primary","Primary school","Middle school","High school or above"))

cat("分析样本:", nrow(d), "｜2020 受损事件:", sum(d$endpoint_2020), "｜首次受损事件:", sum(d$event), "\n")
cat("完整病例（4 变量 + 结局）:", sum(complete.cases(d[, c("age","edu_group","dep","adl","endpoint_2020")])), "\n\n")

## ---------- 1) 主模型（未加权 = 主分析；加权 = 敏感性） ----------
f <- endpoint_2020 ~ age + edu_group + dep + adl
m_u <- glm(f, data = d, family = binomial())
m_w <- svyglm(f, design = svydesign(ids = ~communityID, weights = ~weight, data = d), family = quasibinomial())
cf <- data.frame(
  term = names(coef(m_u)),
  beta_unweighted = as.numeric(coef(m_u)),
  se_unweighted = as.numeric(sqrt(diag(vcov(m_u)))),
  beta_weighted = as.numeric(coef(m_w)[names(coef(m_u))]),
  se_weighted = as.numeric(sqrt(diag(vcov(m_w)))[names(coef(m_u))])
)
cf$OR_unweighted <- exp(cf$beta_unweighted); cf$OR_weighted <- exp(cf$beta_weighted)
cf$ci_u <- sprintf("%.2f (%.2f-%.2f)", cf$OR_unweighted, exp(cf$beta_unweighted - 1.96*cf$se_unweighted), exp(cf$beta_unweighted + 1.96*cf$se_unweighted))
cf$ci_w <- sprintf("%.2f (%.2f-%.2f)", cf$OR_weighted, exp(cf$beta_weighted - 1.96*cf$se_weighted), exp(cf$beta_weighted + 1.96*cf$se_weighted))
print(cf[, c("term","ci_u","ci_w")], row.names = FALSE)
write.csv(cf, file.path(DIR, "screening_index_v3_models_2026-09-26.csv"), row.names = FALSE, fileEncoding = "UTF-8")

## ---------- 2) 整数化评分（SAP §7.3；基于未加权主模型） ----------
b <- coef(m_u)
b_age5 <- 5 * b[["age"]]                       # 年龄按每 5 岁折算
terms_beta <- c(age5 = b_age5,
                edu_primary = b[["edu_groupPrimary school"]],
                edu_middle  = b[["edu_groupMiddle school"]],
                edu_high    = b[["edu_groupHigh school or above"]],
                dep = b[["dep"]], adl = b[["adl"]])
names(terms_beta) <- c("年龄（每 5 岁）","教育：小学（vs 未上学）","教育：初中（vs 未上学）","教育：高中及以上（vs 未上学）","抑郁","ADL 受限")
bmin <- min(abs(terms_beta))
pts <- round(terms_beta / bmin)
score_tab <- data.frame(项目 = names(terms_beta), beta_以5岁折算 = round(terms_beta, 4),
                        每单位分值 = as.integer(pts), 说明 = c(
  "45–49=0分；每增 5 岁 +1 个单位分值", "相对 Less than primary", "相对 Less than primary", "相对 Less than primary",
  "CES-D-10 ≥10", "任一 ADL 受限"))
print(score_tab, row.names = FALSE)
write.csv(score_tab, file.path(DIR, "screening_index_v3_score_table_2026-09-26.csv"), row.names = FALSE, fileEncoding = "UTF-8")

## 计算个体总分
d$score <- as.integer(pts[["年龄（每 5 岁）"]]) * pmax(floor((d$age - 45)/5), 0) +
  as.integer(pts[["教育：小学（vs 未上学）"]]) * as.integer(d$edu_group == "Primary school") +
  as.integer(pts[["教育：初中（vs 未上学）"]]) * as.integer(d$edu_group == "Middle school") +
  as.integer(pts[["教育：高中及以上（vs 未上学）"]]) * as.integer(d$edu_group == "High school or above") +
  as.integer(pts[["抑郁"]]) * d$dep +
  as.integer(pts[["ADL 受限"]]) * d$adl
cat("\n总分范围:", range(d$score, na.rm = TRUE), "｜分布:\n"); print(table(d$score, useNA = "ifany"))

## ---------- 3) 风险分层（三分位 + 概率门槛） ----------
q <- quantile(d$score, c(1/3, 2/3), na.rm = TRUE)
d$stratum <- cut(d$score, c(-Inf, q[1], q[2], Inf), labels = c("低", "中", "高"), right = TRUE)
strata <- d %>% group_by(stratum) %>% summarise(n = n(), 事件 = sum(endpoint_2020),
  事件率 = sprintf("%.1f%%", 100*mean(endpoint_2020)), 平均分 = round(mean(score), 1), .groups = "drop")
print(strata, row.names = FALSE)

d$p_hat <- predict(m_u, newdata = d, type = "response")
dd <- d[!is.na(d$p_hat), ]
thr <- do.call(rbind, lapply(c(0.10, 0.20, 0.30), function(t) {
  pos <- dd$p_hat >= t
  data.frame(门槛 = sprintf("预测概率 ≥%d%%", t*100), 阳性数 = sum(pos), 阳性占比 = sprintf("%.1f%%", 100*mean(pos)),
             阳性者事件率 = sprintf("%.1f%%", 100*mean(dd$endpoint_2020[pos])),
             灵敏度 = sprintf("%.1f%%", 100*sum(dd$endpoint_2020[pos])/sum(dd$endpoint_2020)),
             特异度 = sprintf("%.1f%%", 100*mean(dd$endpoint_2020[!pos] == 0)))
}))
print(thr, row.names = FALSE)
write.csv(strata, file.path(DIR, "screening_index_v3_risk_strata_2026-09-26.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(thr, file.path(DIR, "screening_index_v3_thresholds_2026-09-26.csv"), row.names = FALSE, fileEncoding = "UTF-8")

## ---------- 4) 初步性能（正式版见 A08/A09） ----------
ok <- complete.cases(d[, c("endpoint_2020","p_hat","score")])
r1 <- roc(d$endpoint_2020[ok], d$p_hat[ok], quiet = TRUE)
r2 <- roc(d$endpoint_2020[ok], d$score[ok], quiet = TRUE)
cal <- glm(endpoint_2020 ~ p_hat, data = d[ok, ], family = binomial())
brier <- mean((d$p_hat[ok] - d$endpoint_2020[ok])^2)
perf <- data.frame(指标 = c("AUC（连续预测概率）","AUC（整数总分）","校准斜率（logistic 截距固定 0? 含截距）","校准斜率 β1","Brier"),
                   值 = c(sprintf("%.3f (%.3f-%.3f)", as.numeric(auc(r1)), ci.auc(r1)[1], ci.auc(r1)[3]),
                          sprintf("%.3f (%.3f-%.3f)", as.numeric(auc(r2)), ci.auc(r2)[1], ci.auc(r2)[3]),
                          sprintf("%.3f (%.3f-%.3f)", coef(cal)[1], coef(cal)[1]-1.96*summary(cal)$coef[1,2], coef(cal)[1]+1.96*summary(cal)$coef[1,2]),
                          sprintf("%.2f (%.2f-%.2f)", coef(cal)[2], coef(cal)[2]-1.96*summary(cal)$coef[2,2], coef(cal)[2]+1.96*summary(cal)$coef[2,2]),
                          sprintf("%.4f", brier)))
print(perf, row.names = FALSE)
write.csv(perf, file.path(DIR, "screening_index_v3_performance_2026-09-26.csv"), row.names = FALSE, fileEncoding = "UTF-8")

saveRDS(list(model_u = m_u, model_w = m_w, coef_table = cf, score_table = score_tab,
             strata = strata, thresholds = thr, performance = perf, data = d),
        file.path(DIR, "screening_index_v3_2026-09-26.rds"))
cat("\n已保存至", DIR, "\n")
