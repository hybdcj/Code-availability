# ============================================================================
# 19_screening_index_memory_v3.R —— A07（选项 2）：记忆分结局版筛查指数
# 结局：记忆分（即时+延迟回忆，0–20）；固定切点 = 2011 全队列 mean−1SD（跨波沿用）
# 设计：纵向 CHARLS 2011 → 2015（两波均为 1 次学习试次、相同词表，可比）
# 次要：2011 同时性（横断面）同切点
# 变量：年龄（每5岁）、教育（4类）、抑郁（CES-D-10 ≥10）、ADL（任一受限），均 2011 基线
# 主分析不加权；加权（2011 权重 + PSU=communityID）作敏感性
# 生成：2026-09-26
# ============================================================================
invisible(try(Sys.setlocale("LC_ALL", "zh_CN.UTF-8"), silent = TRUE))
.libPaths(c("D:/R/library", .libPaths()))
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(survival); library(survey); library(pROC) })
P <- "D:/BaiduSyncdisk/科研/徐中豪"; OUT <- file.path(P, "output", "04_纵向队列")
DIR <- file.path(P, "output", "10_筛查指数"); if (!dir.exists(DIR)) dir.create(DIR, recursive = TRUE)

coh <- readRDS(file.path(OUT, "cohort_change_2026-09-13.rds"))$cohort
w <- coh %>% select(KEY, year, imrc, dlrc) %>%
  pivot_wider(names_from = year, values_from = c(imrc, dlrc)) %>%
  mutate(mem2011 = imrc_2011 + dlrc_2011, mem2015 = imrc_2015 + dlrc_2015)

b <- coh %>% filter(year == 2011) %>%
  transmute(KEY, age, female, edu_group, dep = as.integer(cesd10 >= 10), adl, iadl2, ms10_2011 = ms10,
            weight, communityID) %>%
  left_join(w, by = "KEY")
b$edu_group <- factor(b$edu_group, levels = c("Less than primary","Primary school","Middle school","High school or above"))

mu <- mean(b$mem2011, na.rm = TRUE); sdd <- sd(b$mem2011, na.rm = TRUE); cut <- mu - sdd
cat(sprintf("记忆分 2011：mean=%.2f SD=%.2f → 固定切点 mean−1SD=%.2f（取 <= %d）\n", mu, sdd, cut, floor(cut)))
b$impaired_2011 <- as.integer(b$mem2011 <= floor(cut))
# 事件定义：仅对「基线未受损」者评定新发（基线受损者不进风险人群）
b$incident_2015 <- ifelse(is.na(b$mem2011) | is.na(b$mem2015) | b$mem2011 <= floor(cut),
                          NA_integer_,
                          as.integer(b$mem2015 <= floor(cut)))
cat(sprintf("同时性受损率 %.1f%%｜纵向新发率 %.1f%%（随访 4 年）\n\n",
    100*mean(b$impaired_2011, na.rm=TRUE), 100*mean(b$incident_2015, na.rm=TRUE)))

## 分析样本（主：纵向）
d <- b %>% filter(!is.na(incident_2015)) %>% filter(complete.cases(age, edu_group, dep, adl))
cat("纵向分析样本 n =", nrow(d), "｜事件 =", sum(d$incident_2015), "\n")
d_conc <- b %>% filter(complete.cases(impaired_2011, age, edu_group, dep, adl))

## ---------- 模型 ----------
f <- incident_2015 ~ age + edu_group + dep + adl
m_u <- glm(f, data = d, family = binomial())
m_w <- svyglm(f, design = svydesign(ids = ~communityID, weights = ~weight, data = d), family = quasibinomial())
cf <- data.frame(term = names(coef(m_u)),
  OR_unweighted = exp(coef(m_u)), ci_u = sprintf("%.2f (%.2f-%.2f)", exp(coef(m_u)),
    exp(coef(m_u) - 1.96*sqrt(diag(vcov(m_u)))), exp(coef(m_u) + 1.96*sqrt(diag(vcov(m_u))))),
  OR_weighted = exp(coef(m_w)[names(coef(m_u))]), ci_w = sprintf("%.2f (%.2f-%.2f)",
    exp(coef(m_w)[names(coef(m_u))]), exp(coef(m_w)[names(coef(m_u))] - 1.96*sqrt(diag(vcov(m_w)))[names(coef(m_u))]),
    exp(coef(m_w)[names(coef(m_u))] + 1.96*sqrt(diag(vcov(m_w)))[names(coef(m_u))])))
cat("\n== 主模型（未加权 / 加权）==\n"); print(cf, row.names = FALSE)
write.csv(cf, file.path(DIR, "screening_index_memory_v3_models_2026-09-26.csv"), row.names = FALSE, fileEncoding = "UTF-8")

## ---------- 整数化（SAP §7.3）----------
bb <- coef(m_u); b_age5 <- 5*bb[["age"]]
terms_beta <- c(age5 = b_age5, edu_primary = bb[["edu_groupPrimary school"]],
                edu_middle = bb[["edu_groupMiddle school"]], edu_high = bb[["edu_groupHigh school or above"]],
                dep = bb[["dep"]], adl = bb[["adl"]])
# 修订规则（SAP v3.1）：以折算后 |β| 最大者为 10 分基准，其余按比例四舍五入
# （原「最小 |β| 为 1 分基准」在某个变量无效时会导致分值爆炸，见 output/10_筛查指数/ 说明）
bmax <- max(abs(terms_beta)); pts <- round(10 * terms_beta / bmax)
score_tab <- data.frame(项目 = c("年龄（每 5 岁）","教育：小学","教育：初中","教育：高中及以上","抑郁","ADL 受限"),
  beta_5岁折算 = round(terms_beta, 3), 每单位分值 = as.integer(pts),
  编码 = c("45–49=0；每增 5 岁 +2 分","vs 未上学","vs 未上学","vs 未上学","CES-D-10 ≥10","任一受限（β≈0 → 0 分，见 SAP 条款）"))
print(score_tab, row.names = FALSE)
write.csv(score_tab, file.path(DIR, "screening_index_memory_v3_score_table_2026-09-26.csv"), row.names = FALSE, fileEncoding = "UTF-8")

d$score <- as.integer(pts[["age5"]]) * pmax(floor((d$age - 45)/5), 0) +
  as.integer(pts[["edu_primary"]]) * as.integer(d$edu_group == "Primary school") +
  as.integer(pts[["edu_middle"]]) * as.integer(d$edu_group == "Middle school") +
  as.integer(pts[["edu_high"]]) * as.integer(d$edu_group == "High school or above") +
  as.integer(pts[["dep"]]) * d$dep + as.integer(pts[["adl"]]) * d$adl
cat("\n总分范围:", range(d$score), "\n")

## ---------- 分层 ----------
q <- quantile(d$score, c(1/3, 2/3))
d$stratum <- cut(d$score, c(-Inf, q[1], q[2], Inf), labels = c("低","中","高"), right = TRUE)
strata <- d %>% group_by(stratum) %>% summarise(n = n(), 事件 = sum(incident_2015),
  事件率 = sprintf("%.1f%%", 100*mean(incident_2015)), 平均分 = round(mean(score),1), .groups="drop")
print(strata, row.names = FALSE)
write.csv(strata, file.path(DIR, "screening_index_memory_v3_risk_strata_2026-09-26.csv"), row.names = FALSE, fileEncoding = "UTF-8")

## ---------- 性能 ----------
d$p_hat <- predict(m_u, type = "response")
r1 <- roc(d$incident_2015, d$p_hat, quiet = TRUE); r2 <- roc(d$incident_2015, d$score, quiet = TRUE)
lp_hat <- qlogis(pmin(pmax(d$p_hat, 1e-6), 1 - 1e-6))      # 校准斜率须在 logit 尺度计算
cal <- glm(incident_2015 ~ lp_hat, data = d, family = binomial())
perf <- data.frame(指标 = c("AUC（连续预测概率）","AUC（整数总分）","校准截距","校准斜率","Brier"),
  值 = c(sprintf("%.3f (%.3f-%.3f)", as.numeric(auc(r1)), ci.auc(r1)[1], ci.auc(r1)[3]),
         sprintf("%.3f (%.3f-%.3f)", as.numeric(auc(r2)), ci.auc(r2)[1], ci.auc(r2)[3]),
         sprintf("%.3f (%.3f-%.3f)", coef(cal)[1], coef(cal)[1]-1.96*summary(cal)$coef[1,2], coef(cal)[1]+1.96*summary(cal)$coef[1,2]),
         sprintf("%.2f (%.2f-%.2f)", coef(cal)[2], coef(cal)[2]-1.96*summary(cal)$coef[2,2], coef(cal)[2]+1.96*summary(cal)$coef[2,2]),
         sprintf("%.4f", mean((d$p_hat - d$incident_2015)^2))))
print(perf, row.names = FALSE)
write.csv(perf, file.path(DIR, "screening_index_memory_v3_performance_2026-09-26.csv"), row.names = FALSE, fileEncoding = "UTF-8")

## 阈值
thr <- do.call(rbind, lapply(c(0.10, 0.15, 0.20), function(t) {
  pos <- d$p_hat >= t
  data.frame(门槛 = sprintf("预测概率 ≥%d%%", t*100), 阳性数 = sum(pos), 阳性占比 = sprintf("%.1f%%", 100*mean(pos)),
    阳性者事件率 = sprintf("%.1f%%", 100*mean(d$incident_2015[pos])),
    灵敏度 = sprintf("%.1f%%", 100*sum(d$incident_2015[pos])/sum(d$incident_2015)),
    特异度 = sprintf("%.1f%%", 100*mean(d$incident_2015[!pos] == 0)))
}))
print(thr, row.names = FALSE)
write.csv(thr, file.path(DIR, "screening_index_memory_v3_thresholds_2026-09-26.csv"), row.names = FALSE, fileEncoding = "UTF-8")

## 次要：同时性（2011 横断面）
m_c <- glm(impaired_2011 ~ age + edu_group + dep + adl, data = d_conc, family = binomial())
cat("\n== 次要：2011 同时性 AUC =", sprintf("%.3f", as.numeric(auc(roc(d_conc$impaired_2011,
    predict(m_c, newdata = d_conc, type = "response"), quiet = TRUE)))), "（n =", nrow(d_conc), "）\n")

saveRDS(list(model_u = m_u, model_w = m_w, coef = cf, score_table = score_tab, strata = strata,
             perf = perf, thresholds = thr, cut = floor(cut), data = d, model_concurrent = m_c),
        file.path(DIR, "screening_index_memory_v3_2026-09-26.rds"))
cat("\n已保存至", DIR, "\n")
