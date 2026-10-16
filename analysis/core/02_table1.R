# =============================================================================
# 02_table1.R  --  CHARLS 2020 基线特征表（Table 1）+ 患病率（v4.4，2026-09-12）
# 口径：主分析 = 认知 9 分项完整（教育已补齐，零缺失）；切点 = 教育分层（<本层均值−1SD）
# 版本：v4.4（2026-09-12）——① TableS1 放回 edu_f（选择偏倚教育梯度，回填后两组教育零缺失）；
#                      ② S2/S4/S5 的 CSV 末行脚注改为合法带引号字段（原含逗号会被解析成多列）；
#                      ③ 教育分层（B/A 两口径）统一为加权估计 + 95%CI，并另列未加权点估计；
#                      ④ TableS2/S3/S4 输出新增 n 列（各模型实际样本量）
# 输入：output/01_2020单期基线/analysis_2020.rds
# 输出：output/01_2020单期基线/Table1_baseline.docx / .html、Table1_and_prevalence.xlsx、
#       TableS1_complete_vs_incomplete.csv、TableS2_pain_sleep_varform.csv（加权）、
#       TableS3_comorbidity_adj.csv（慢病/功能/行为变量加权校正关联）、
#       TableS4_sensitivityA_adj.csv（敏感性A口径下的加权校正关联，与 S3 对称）
# =============================================================================
Sys.setlocale("LC_ALL", "zh_CN.UTF-8")
.libPaths(c("D:/R/library", .libPaths()))
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(gtsummary); library(survey)
  library(openxlsx); library(flextable); library(srvyr)
})
proj <- "D:/BaiduSyncdisk/科研/徐中豪"; out <- file.path(proj, "output", "01_2020单期基线")
dat_all <- readRDS(file.path(out, "analysis_2020.rds"))

## ---- 1. 样本切分 ----
dat    <- dat_all %>% filter(primary == 1, !is.na(cog_impair))   # 主分析样本
# 注：noedu（认知完整但教育缺失）在 2026-09-12 教育补齐后恒为空集；
#     该变量仅为兼容旧流程保留，下方 S1 分组构造已自动跳过空组。
noedu  <- dat_all %>% filter(primary == 1,  is.na(cog_impair))
incomp <- dat_all %>% filter(primary == 0)                       # 未测全（敏感性）
cat(sprintf("主分析 n=%d | 认知完整但教育缺失 n=%d | 未测全 n=%d\n",
            nrow(dat), nrow(noedu), nrow(incomp)))

## ---- 2. 变量整理与标签 ----
prep <- function(d) d %>% mutate(
  age_cat   = cut(age, breaks = c(44, 54, 64, 74, Inf),
                  labels = c("45-54", "55-64", "65-74", "≥75")),
  female_f  = factor(female, c(0,1), c("Male","Female")),
  edu_f     = factor(edu_group,
                     levels = c("Less than primary","Primary school","Middle school","High school or above")),
  married_f = factor(ifelse(marital %in% c(1,2), "Married/partnered", "Not partnered"),
                     levels = c("Married/partnered","Not partnered")),
  urban_f   = factor(urban3, c(0,1), c("Rural","Urban")),
  hukou_f   = factor(case_when(hukou == 1 ~ "Agricultural",
                               hukou %in% c(2,3) ~ "Non-agricultural/Unified"),
                     levels = c("Agricultural","Non-agricultural/Unified")),
  pain_f    = factor(pain_any, c(0,1), c("No pain","Any pain")),
  pain3_f   = factor(pain_3cat, c(0,1,2), c("No pain","Mild","Severe")),
  sleep_f   = factor(sleep_cat, levels = c("<6h","6-8h",">8h")),
  sleep_ref = factor(sleep_cat, levels = c("6-8h","<6h",">8h")),   # 建模时以 6-8h 为参照
  dep_f     = factor(depression, c(0,1), c("No","Yes")),
  adl_f     = factor(adl_any,  c(0,1), c("No","Yes")),
  iadl_f    = factor(iadl_any, c(0,1), c("No","Yes")),
  multi_f   = factor(multimorbidity, c(0,1), c("No","Yes")),
  cog_f     = factor(cog_impair, levels = c(0,1), labels = c("Normal","Impairment")),
  hypertension_f   = factor(hypertension,   c(0,1), c("No","Yes")),
  dyslipidemia_f   = factor(dyslipidemia,   c(0,1), c("No","Yes")),
  diabetes_f       = factor(diabetes,       c(0,1), c("No","Yes")),
  heart_disease_f  = factor(heart_disease,  c(0,1), c("No","Yes")),
  stroke_f         = factor(stroke,         c(0,1), c("No","Yes")),
  lung_disease_f   = factor(lung_disease,   c(0,1), c("No","Yes")),
  arthritis_f      = factor(arthritis,      c(0,1), c("No","Yes")),
  kidney_disease_f = factor(kidney_disease, c(0,1), c("No","Yes")),
  smoke_f   = factor(smoke, c(0,1,2), c("Never","Former","Current")),
  drink_f   = factor(drink_freq, c(3,2,1),
                     c("None in past year","Less than monthly","At least monthly")),
  pa_any_f  = factor(pa_any, c(0,1), c("No","Yes")),
  sleep_restless_f = factor(sleep_restless, c(0,1), c("No","Yes"))
)
dat <- prep(dat); noedu <- prep(noedu); incomp <- prep(incomp)

## ---- 3. Table 1（主分析样本，按认知状态分组） ----
t1_vars <- c("age", "age_cat", "female_f", "edu_f", "married_f", "urban_f", "hukou_f",
             "hypertension_f","dyslipidemia_f","diabetes_f","heart_disease_f","stroke_f",
             "lung_disease_f","arthritis_f","kidney_disease_f","multi_f","n_chronic",
             "pain3_f","sleep_h","sleep_f","sleep_restless_f",
             "smoke_f","drink_f","pa_any_f",
             "dep_f","cesd10","adl_f","iadl_f",
             "imrc","dlrc","memory","orient","ser7","draw","mental","global")
labels1 <- list(
  age ~ "Age, years", age_cat ~ "Age group, n (%)", female_f ~ "Sex", edu_f ~ "Education",
  married_f ~ "Marital status", urban_f ~ "Residence", hukou_f ~ "Hukou type",
  hypertension_f ~ "Hypertension", dyslipidemia_f ~ "Dyslipidemia", diabetes_f ~ "Diabetes",
  heart_disease_f ~ "Heart disease", stroke_f ~ "Stroke", lung_disease_f ~ "Chronic lung disease",
  arthritis_f ~ "Arthritis/rheumatism", kidney_disease_f ~ "Kidney disease",
  multi_f ~ "Multimorbidity (≥2)", n_chronic ~ "Number of chronic diseases",
  pain_f ~ "Body pain", pain3_f ~ "Body pain severity", sleep_h ~ "Night sleep duration, h", sleep_f ~ "Sleep duration category",
  sleep_restless_f ~ "Restless sleep (≥3-4 days/week)",
  smoke_f ~ "Smoking status", drink_f ~ "Alcohol use, past year",
  pa_any_f ~ "Any regular physical activity (≥10 min)",
  dep_f ~ "Depression (CES-D-10 ≥ 10)", cesd10 ~ "CES-D-10 score",
  adl_f ~ "ADL limitation", iadl_f ~ "IADL limitation",
  imrc ~ "Immediate word recall (0-10)", dlrc ~ "Delayed word recall (0-10)",
  memory ~ "Episodic memory (0-20)", orient ~ "Orientation (0-5)", ser7 ~ "Serial 7s (0-5)",
  draw ~ "Figure drawing (0-1)", mental ~ "Mental status (0-11)", global ~ "Global cognition (0-31)")
# 按实际入表变量动态筛选标签（避免 label 引用不存在的列）
labels_for <- function(vars) labels1[vapply(labels1, function(f) as.character(f[[2]]), character(1)) %in% vars]

tab1 <- dat %>% select(all_of(c("cog_f", t1_vars))) %>%
  tbl_summary(by = cog_f,
    statistic = list(all_continuous() ~ "{mean} ({sd})", all_categorical() ~ "{n} ({p}%)"),
    digits = list(all_continuous() ~ 1, all_categorical() ~ c(0, 1)),
    missing = "ifany", label = labels_for(t1_vars)) %>%
  add_p(test = list(all_continuous() ~ "wilcox.test", all_categorical() ~ "chisq.test")) %>%
  add_overall(last = TRUE) %>% bold_labels() %>% modify_header(label ~ "**Characteristic**") %>%
  modify_table_styling(columns = label,
    footnote = "Cognitive impairment was defined with education-stratified cutoffs; between-group comparisons for education and sex are by design and not interpretable as tests of association.")

## ---- 4. Table S1：主分析 vs 未测全（缺失偏倚评估） ----
# 2026-09-12（v4.2）：edu_f 已放回 s1_vars。
#   旧理由（"教育缺失本就是第 2 组的定义"）对应回填前的分组：当时第 2 组 = 认知完整但教育缺失。
#   回填后第 2 组已是"认知未测全"，且两组教育均零缺失，教育梯度（小学以下 71.7% vs 36.8%）
#   是选择偏倚的关键指标，必须随表提交，否则正文数字在附表里查不到。
s1_vars <- c("age", "age_cat", "female_f", "edu_f", "urban_f", "married_f",
             "hypertension_f","diabetes_f","stroke_f","arthritis_f",
             "pain3_f","dep_f","adl_f","iadl_f","sleep_h")
# 2026-09-12：教育经 zredu+ba010 合并后已无缺失，空组自动剔除
grp_list <- list(
  "Analytic sample (complete cognition)"       = dat,
  "Complete cognition, education missing"      = noedu,
  "Incomplete cognition (sensitivity only)"    = incomp
)
grp_list <- grp_list[vapply(grp_list, nrow, integer(1)) > 0]
grp <- bind_rows(lapply(names(grp_list), function(g)
  grp_list[[g]] %>% select(all_of(s1_vars)) %>% mutate(grp = g))) %>%
  mutate(grp = factor(grp, levels = names(grp_list)))
tabS1 <- grp %>% tbl_summary(by = grp,
    statistic = list(all_continuous() ~ "{mean} ({sd})", all_categorical() ~ "{n} ({p}%)"),
    digits = list(all_categorical() ~ c(0,1)), missing = "no",
    label = labels_for(s1_vars)) %>%
  add_p(test = list(all_continuous() ~ "kruskal.test", all_categorical() ~ "chisq.test")) %>%
  bold_labels()

## ---- 4b. Table S2：疼痛/睡眠"变量形式"敏感性分析（粗模型 + 校正模型 OR） ----
# 目的：同时呈现二分类与有序分组两种口径，避免"只报显著结果"的质疑
# 口径：survey 加权 logistic（INDV_weight，社区聚类）；校正 年龄+性别+教育(四组)+城乡
des_main <- dat %>% filter(!is.na(INDV_weight), !is.na(communityID)) %>%
  as_survey_design(weights = INDV_weight, ids = communityID)
# 注：只汇报暴露项系数；协变量（含 edu_f）系数不汇报——教育正是本结局分层切点的依据，
#     其系数不可解释，若输出易被误读为"学历越高、认知受损风险越高"。
or_block <- function(expo, label) {
  adj <- " + age + female_f + edu_f + urban_f"
  grab <- function(f, tag) {
    m <- svyglm(f, design = des_main, family = quasibinomial)
    s <- summary(m)$coefficients
    s <- s[grepl(paste0("^", expo), rownames(s)), , drop = FALSE]   # 仅暴露项
    data.frame(exposure = label, term = rownames(s), model = tag,
               OR_95CI = sprintf("%.2f (%.2f-%.2f)", exp(s[,1]),
                                 exp(s[,1] - 1.96*s[,2]), exp(s[,1] + 1.96*s[,2])),
               p = sprintf("%.3g", s[,4]),
               # 该模型实际进入分析的样本量（2026-09-12 补：各模型缺失不同，n 不同）
               n = nrow(model.frame(m)), row.names = NULL, stringsAsFactors = FALSE)
  }
  rbind(grab(as.formula(paste0("cog_impair ~ ", expo)), "Crude"),
        grab(as.formula(paste0("cog_impair ~ ", expo, adj)), "Adjusted"))
}
NOTE_S2 <- paste0("Note: all models are survey-weighted logistic regressions (individual weight INDV_weight, ",
                  "community as cluster), adjusted for age, sex, education (4 groups) and residence. ",
                  "Only exposure coefficients are reported; covariate coefficients (including education) are ",
                  "deliberately omitted because education defines the stratified cut-off used for the outcome ",
                  "and is therefore not interpretable.")
tabS2 <- bind_rows(
  or_block("pain_f",        "Body pain, binary (any vs none)"),
  or_block("pain3_f",       "Body pain severity (3 levels, ref = no pain)"),
  or_block("sleep_h",       "Sleep duration per 1 h (continuous)"),
  or_block("sleep_ref",     "Sleep duration grouped (ref = 6-8 h)"),
  or_block("sleep_restless_f", "Restless sleep")
)

## ---- 4c. Table S3：慢病及功能/行为变量的加权校正关联 ----
# 与 Table S2 同一口径：svyglm 加权，校正 年龄+性别+教育(四组)+城乡
run_adj <- function(v, label){
  f <- as.formula(paste0("cog_impair ~ ", v, " + age + female_f + edu_f + urban_f"))
  m <- svyglm(f, design = des_main, family = quasibinomial)
  s <- summary(m)$coefficients
  i <- grep(paste0("^", v), rownames(s))
  data.frame(exposure = label, term = rownames(s)[i],
             OR_95CI = sprintf("%.2f (%.2f-%.2f)", exp(s[i,1]),
                               exp(s[i,1] - 1.96*s[i,2]), exp(s[i,1] + 1.96*s[i,2])),
             p = sprintf("%.3g", s[i,4]),
             n = nrow(model.frame(m)), row.names = NULL, stringsAsFactors = FALSE)
}
tabS3 <- bind_rows(
  run_adj("hypertension",   "Hypertension"),
  run_adj("stroke",         "Stroke"),
  run_adj("arthritis",      "Arthritis/rheumatism"),
  run_adj("diabetes",       "Diabetes"),
  run_adj("lung_disease",   "Chronic lung disease"),
  run_adj("heart_disease",  "Heart disease"),
  run_adj("kidney_disease", "Kidney disease"),
  run_adj("n_chronic",      "Number of chronic diseases"),
  run_adj("multimorbidity", "Multimorbidity (>=2)"),
  run_adj("depression",     "Depression (CES-D-10 >= 10)"),
  run_adj("adl_any",        "ADL limitation"),
  run_adj("iadl_any",       "IADL limitation"),
  run_adj("pa_any",         "Regular physical activity"),
  run_adj("dyslipidemia",   "Dyslipidemia"),
  run_adj("drink_f",        "Alcohol use, past year (ref = none)")
)

## ---- 4d. Table S4：敏感性 A 口径（整体均值−1SD 切点）下的加权校正关联 ----
# 与 S3 对称：同一批暴露、同一加权与校正集，仅结局换成 cog_low_1sd。
# 目的：卒中/多重共病/饮酒等结论存在切点口径敏感性，须与主定义并列报告。
run_adj_A <- function(v, label){
  f <- as.formula(paste0("cog_low_1sd ~ ", v, " + age + female_f + edu_f + urban_f"))
  m <- svyglm(f, design = des_main, family = quasibinomial)
  s <- summary(m)$coefficients
  i <- grep(paste0("^", v), rownames(s))
  data.frame(exposure = label, term = rownames(s)[i],
             OR_95CI = sprintf("%.2f (%.2f-%.2f)", exp(s[i,1]),
                               exp(s[i,1] - 1.96*s[i,2]), exp(s[i,1] + 1.96*s[i,2])),
             p = sprintf("%.3g", s[i,4]),
             n = nrow(model.frame(m)), row.names = NULL, stringsAsFactors = FALSE)
}
tabS4 <- bind_rows(
  run_adj_A("hypertension",   "Hypertension"),
  run_adj_A("stroke",         "Stroke"),
  run_adj_A("arthritis",      "Arthritis/rheumatism"),
  run_adj_A("diabetes",       "Diabetes"),
  run_adj_A("lung_disease",   "Chronic lung disease"),
  run_adj_A("heart_disease",  "Heart disease"),
  run_adj_A("kidney_disease", "Kidney disease"),
  run_adj_A("n_chronic",      "Number of chronic diseases"),
  run_adj_A("multimorbidity", "Multimorbidity (>=2)"),
  run_adj_A("depression",     "Depression (CES-D-10 >= 10)"),
  run_adj_A("adl_any",        "ADL limitation"),
  run_adj_A("iadl_any",       "IADL limitation"),
  run_adj_A("pa_any",         "Regular physical activity"),
  run_adj_A("dyslipidemia",   "Dyslipidemia"),
  run_adj_A("drink_f",        "Alcohol use, past year (ref = none)")
)

## ---- 4e. Table S5：三种切点口径（B / A / C）下核心暴露的校正关联 ----
# 目的：SAP §6 的 S1（口径 C）与 S2（口径 A）敏感性分析需要核心暴露在两口径下的模型结果；
#       此前仅对慢病/功能/行为变量跑过 A 口径，**睡眠分组与疼痛三档未在 A/C 口径下建模**。
# 口径：survey 加权 logistic；校正 年龄+性别+教育(四组)+城乡；结局分别取 B/A/C 三种切点定义。
core <- list(
  list(v = "depression", lab = "Depression (CES-D-10 >= 10)"),
  list(v = "adl_any",    lab = "ADL limitation"),
  list(v = "iadl_any",   lab = "IADL limitation"),
  list(v = "sleep_ref",  lab = "Sleep duration grouped (ref = 6-8 h)"),
  list(v = "pa_any",     lab = "Regular physical activity"),
  list(v = "pain3_f",    lab = "Body pain severity (ref = no pain)")
)
outc <- c(B = "cog_impair", A = "cog_low_1sd", C = "cog_low_q1")
one_oc <- function(v, oc) {
  f <- as.formula(paste0(oc, " ~ ", v, " + age + female_f + edu_f + urban_f"))
  m <- svyglm(f, design = des_main, family = quasibinomial)
  s <- summary(m)$coefficients
  i <- grep(paste0("^", v), rownames(s))
  data.frame(term = rownames(s)[i],
             OR = sprintf("%.2f (%.2f-%.2f)", exp(s[i,1]),
                          exp(s[i,1] - 1.96*s[i,2]), exp(s[i,1] + 1.96*s[i,2])),
             p  = sprintf("%.3g", s[i,4]), row.names = NULL, stringsAsFactors = FALSE)
}
tabS5 <- do.call(rbind, lapply(core, function(cx) {
  res <- lapply(names(outc), function(k) one_oc(cx$v, outc[[k]]))
  base <- res[[1]][, "term", drop = FALSE]
  for (j in seq_along(outc)) {
    base[[paste0(names(outc)[j], "_OR_95CI")]] <- res[[j]]$OR
    base[[paste0(names(outc)[j], "_p")]]       <- res[[j]]$p
  }
  cbind(exposure = cx$lab, base, row.names = NULL)
}))
NOTE_S5 <- paste0("Note: B = education-stratified cut-off (primary), A = overall mean-1SD, ",
                  "C = lowest quartile. All models are survey-weighted logistic regressions ",
                  "(individual weight INDV_weight, community as cluster), adjusted for age, sex, ",
                  "education (4 groups) and residence. Only exposure coefficients are reported.")

## ---- 5. 患病率（主定义 + 敏感性） ----
# des_main 已在 4b 节创建（主分析样本的 survey 设计）
mk_design <- function(d) d %>% filter(!is.na(INDV_weight), !is.na(communityID)) %>%
  as_survey_design(weights = INDV_weight, ids = communityID)
prev_row <- function(var, label, d, design = NULL) {
  y <- d[[var]]; n <- sum(!is.na(y)); k <- sum(y == 1, na.rm = TRUE)
  ci_u <- prop.test(k, n)$conf.int
  if (!is.null(design)) {
    f <- svymean(as.formula(paste0("~", var)), design, na.rm = TRUE); ci_w <- confint(f)
    w <- sprintf("%.1f%% (%.1f-%.1f)", 100*coef(f)[1], 100*ci_w[1], 100*ci_w[2])
    nw <- nrow(design$variables)
    yw <- design$variables[[var]]; nwv <- sum(!is.na(yw)); kw <- sum(yw == 1, na.rm = TRUE)
    ci_uw <- prop.test(kw, nwv)$conf.int
    uw <- sprintf("%.1f%% (%.1f-%.1f)", 100*kw/nwv, 100*ci_uw[1], 100*ci_uw[2])
  } else { w <- NA_character_; uw <- NA_character_; nw <- NA_integer_ }
  data.frame(definition = label, n = n, unweighted = sprintf("%.1f%% (%.1f-%.1f)", 100*k/n, 100*ci_u[1], 100*ci_u[2]),
             unweighted_wsub = uw, weighted = w, n_weighted = nw, stringsAsFactors = FALSE)
}
des_all  <- mk_design(dat_all %>% filter(!is.na(cog_impair)))
prev_tab <- bind_rows(
  prev_row("cog_impair",   "【主定义】教育分层切点，认知完整者",        dat, des_main),
  prev_row("cog_low_1sd",  "敏感性A：整体均值−1SD，认知完整者",        dat, des_main),
  prev_row("cog_low_q1",   "敏感性C：最低四分位，认知完整者",          dat, des_main),
  prev_row("cog_impair",   "敏感性D：教育分层切点，纳入未测全者",      dat_all %>% filter(!is.na(cog_impair)), des_all)
)

## ---- 6. 分层患病率 ----
strat <- function(d, byvar, label, design, def = "cog_impair", tag = "primary") {
  dd <- d %>% filter(!is.na(.data[[byvar]]))
  ds <- dd %>% filter(!is.na(INDV_weight), !is.na(communityID)) %>%
    as_survey_design(weights = INDV_weight, ids = communityID)
  ds %>% group_by(.data[[byvar]]) %>%
    summarise(p = survey_mean(.data[[def]], na.rm = TRUE, vartype = "ci"), n = unweighted(n())) %>%
    mutate(stratum = paste0(label, " [", tag, "]"), level = as.character(.data[[byvar]]),
           prevalence = sprintf("%.1f%% (%.1f-%.1f)", 100*p, 100*p_low, 100*p_upp)) %>%
    select(stratum, level, prevalence, n) %>% as.data.frame()
}
strat_tab <- bind_rows(
  strat(dat, "age_cat",  "Age group",  des_main, "cog_impair",  "主定义B"),
  strat(dat, "female_f", "Sex",        des_main, "cog_impair",  "主定义B"),
  strat(dat, "urban_f",  "Residence",  des_main, "cog_impair",  "主定义B"),
  strat(dat, "age_cat",  "Age group",  des_main, "cog_low_1sd", "敏感性A"),
  strat(dat, "female_f", "Sex",        des_main, "cog_low_1sd", "敏感性A"),
  strat(dat, "urban_f",  "Residence",  des_main, "cog_low_1sd", "敏感性A")
) %>% mutate(unweighted = NA_character_)
# 教育梯度（v4.3，2026-09-12）：与其余分层统一为【加权 + 95%CI】，另列【未加权】点估计（全样本 n=14,655）。
#   主定义 B 的切点按教育层内定义，B 列教育梯度被按构造压缩（加权约 13-17%，组间 CI 重叠）；
#   敏感性 A（整体切点）下教育梯度即真实梯度图谱（加权 34.8% -> 1.8%）。
edu_unw <- bind_rows(lapply(c("主定义B", "敏感性A"), function(tag) {
  v <- if (tag == "主定义B") "cog_low_edu" else "cog_low_1sd"
  dat %>% filter(!is.na(edu_f)) %>% group_by(edu_f) %>%
    summarise(p = mean(.data[[v]], na.rm = TRUE), .groups = "drop") %>%
    transmute(tag = tag, level = as.character(edu_f), unweighted = sprintf("%.1f%%", 100*p))
}))
strat_edu <- bind_rows(lapply(c("主定义B", "敏感性A"), function(tag) {
  v <- if (tag == "主定义B") "cog_low_edu" else "cog_low_1sd"
  ds <- dat %>% filter(!is.na(edu_f), !is.na(INDV_weight), !is.na(communityID)) %>%
    as_survey_design(weights = INDV_weight, ids = communityID)
  ds %>% group_by(edu_f) %>%
    summarise(p = survey_mean(.data[[v]], na.rm = TRUE, vartype = "ci"), n = unweighted(n())) %>%
    mutate(stratum = paste0("Education [", tag, ", 加权]"), level = as.character(edu_f),
           prevalence = sprintf("%.1f%% (%.1f-%.1f)", 100*p, 100*p_low, 100*p_upp), tag = tag) %>%
    select(tag, stratum, level, prevalence, n) %>% as.data.frame()
})) %>% left_join(edu_unw, by = c("tag", "level")) %>%
  select(stratum, level, prevalence, n, unweighted) %>% as.data.frame()

strat_all <- bind_rows(strat_tab, strat_edu)
NOTE_STRAT <- paste0("Education rows: prevalence = survey-weighted estimate (95% CI) on the weighted-eligible ",
  "subsample; unweighted = point estimate on the full analytic sample (n=14,655). Other strata: weighted ",
  "estimates only; n = weighted-eligible respondents.")

## ---- 7. 输出 ----
dir.create(out, showWarnings = FALSE)
print(tab1)
as_gt(tab1) %>% gt::gtsave(file.path(out, "Table1_baseline.html"))
as_gt(tabS1) %>% gt::gtsave(file.path(out, "TableS1_complete_vs_incomplete.html"))
write.csv(as.data.frame(tabS1$table_body), file.path(out, "TableS1_complete_vs_incomplete.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
write.csv(tabS2, file.path(out, "TableS2_pain_sleep_varform.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
cat(sprintf('"# %s"\n', NOTE_S2), file = file.path(out, "TableS2_pain_sleep_varform.csv"), append = TRUE)
write.csv(tabS3, file.path(out, "TableS3_comorbidity_adj.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
write.csv(tabS4, file.path(out, "TableS4_sensitivityA_adj.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
cat(sprintf('"# %s"\n', NOTE_S2), file = file.path(out, "TableS4_sensitivityA_adj.csv"), append = TRUE)
write.csv(tabS5, file.path(out, "TableS5_three_cutpoints_core.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
cat(sprintf('"# %s"\n', NOTE_S5), file = file.path(out, "TableS5_three_cutpoints_core.csv"), append = TRUE)

w1 <- createWorkbook()
addWorksheet(w1, "Table1_raw");            writeData(w1, "Table1_raw", as_tibble(tab1$table_body))
addWorksheet(w1, "TableS1_raw");           writeData(w1, "TableS1_raw", as_tibble(tabS1$table_body))
addWorksheet(w1, "TableS2_varform");       writeData(w1, "TableS2_varform", tabS2)
writeData(w1, "TableS2_varform", NOTE_S2, startRow = nrow(tabS2) + 3, startCol = 1)
addWorksheet(w1, "TableS3_comorbidity");   writeData(w1, "TableS3_comorbidity", tabS3)
writeData(w1, "TableS3_comorbidity", NOTE_S2, startRow = nrow(tabS3) + 3, startCol = 1)
addWorksheet(w1, "TableS4_sensitivityA");  writeData(w1, "TableS4_sensitivityA", tabS4)
writeData(w1, "TableS4_sensitivityA", NOTE_S2, startRow = nrow(tabS4) + 3, startCol = 1)
addWorksheet(w1, "TableS5_three_cutpoints"); writeData(w1, "TableS5_three_cutpoints", tabS5)
writeData(w1, "TableS5_three_cutpoints", NOTE_S5, startRow = nrow(tabS5) + 3, startCol = 1)
addWorksheet(w1, "Prevalence");            writeData(w1, "Prevalence", prev_tab)
addWorksheet(w1, "Prevalence_stratified"); writeData(w1, "Prevalence_stratified", strat_all)
writeData(w1, "Prevalence_stratified", NOTE_STRAT, startRow = nrow(strat_all) + 3, startCol = 1)
saveWorkbook(w1, file.path(out, "Table1_and_prevalence.xlsx"), overwrite = TRUE)

# docx 最后写：若被 Word 占用，不影响 html/xlsx 产出
save_as_docx(as_flex_table(tab1), path = file.path(out, "Table1_baseline.docx"))

cat("\n=== 认知受损患病率 ===\n"); print(prev_tab)
cat("\n=== 分层 ===\n"); print(strat_all)
cat("\n=== Table S2：变量形式敏感性分析 ===\n"); print(tabS2)
cat("\n=== Table S3：慢病及功能/行为变量加权校正关联 ===\n"); print(tabS3)
cat("\n=== Table S4：敏感性A口径加权校正关联 ===\n"); print(tabS4)
cat("\n=== Table S5：三口径核心暴露 ===\n"); print(tabS5)
cat("\n完成：Table1_baseline.docx / .html、TableS1_*.html/.csv、TableS2/S3/S4/S5_*.csv、Table1_and_prevalence.xlsx\n")
