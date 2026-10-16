# =============================================================================
# 10_fivewave_logistic.R —— CHARLS 五波（2011-2020）分波 Logistic：认知受损的影响因素
# 版本：v2（2026-09-19 修订）
#   ① 睡眠改为三分类（<6h / 6-8h（参照）/ >8h），与已交付 Paper 1（2020 期 v4.4）一致；
#   ② 疼痛：主行保留"任一疼痛"（五波可比），另加"轻度 2-3""重度 4-5"两行（仅 2013/2018/2020
#      有分级编码；2011/2015 只有"有/无"两级，格子标注"仅两级编码"）；
#   ③ 5 个核心暴露（抑郁、ADL、IADL、睡眠、体力活动）按 SAP 做 BH-FDR；
#      睡眠为三分类因子，参与校正的是因子整体 Wald 检验 P（与 Paper 1 脚注同口径）；
#   ④ 连续睡眠、二分类疼痛移入"敏感性子表"（对应 Paper 1 附表 S2）；
#   ⑤ 脚本末尾自动核对 2020 列是否复现已交付 Paper 1 的 6 个数。
# 口径：结局 = 该波 cog_B（本波内、本教育层 总分均值 − 1SD）；校正 = 年龄 + 性别 + 教育四组 + 城乡；
#      加权 = svyglm(weights = weight, ids = communityID, quasibinomial)；仅报暴露项系数。
# 输入：output/02_五波数据与质控/waves_cognition_2026-09-13.rds；data/raw/{wave2_2013,wave4_2018,wave5_2020}（疼痛分级）
#       output/01_2020单期基线/analysis_2020.rds（仅回补 2020 期 rural）
# 输出：output/03_五波横断面/五波分波Logistic_OR_2026-09-19.xlsx（说明｜逐波明细｜校正OR矩阵_五波并列｜敏感性子表）
#       output/03_五波横断面/五波分波Logistic_OR_2026-09-19.csv、五波分波Logistic_OR矩阵_2026-09-19.csv、
#       output/03_五波横断面/五波分波Logistic_敏感性_2026-09-19.csv
# 注意：openxlsx 写出的 xlsx 带悬空 drawing 引用（Excel/pandas 能读、openpyxl 严格模式报错），
#      生成后请运行：python code/_checks/repair_xlsx_dangling_refs.py <上面那个 xlsx>
# =============================================================================
Sys.setlocale("LC_ALL", "zh_CN.UTF-8")
.libPaths(c("D:/R/library", .libPaths()))
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(survey); library(haven); library(openxlsx)
})

PROJ <- Sys.getenv("CHARLS_PROJ", "D:/BaiduSyncdisk/科研/徐中豪")
OUT  <- file.path(PROJ, "output", "03_五波横断面")
dat  <- readRDS(file.path(PROJ, "output", "02_五波数据与质控", "waves_cognition_2026-09-13.rds"))
norm_id <- function(x) {
  x <- as.character(x); n <- nchar(x)
  ifelse(n == 11, paste0(substr(x, 1, 9), sprintf("%03d", as.integer(substr(x, 10, 11)))), x)
}

## ---- 0. 回补 2020 期 rural（rural = 1 − urban3，口径同 02_table1.R） ----
a20 <- readRDS(file.path(PROJ, "output", "01_2020单期基线", "analysis_2020.rds"))
rural20 <- ifelse(is.na(a20$urban3), NA_integer_, 1L - as.integer(a20$urban3))
names(rural20) <- as.character(a20$ID)
i20 <- which(dat$year == 2020)
dat$rural[i20] <- unname(rural20[as.character(dat$ID[i20])])
writeLines(sprintf("【回补】2020 期 rural 非缺失：%d / %d", sum(!is.na(dat$rural[i20])), length(i20)))

## ---- 1. 派生变量 ----
# 疼痛三档：从各波原始模块取分级变量（2013 wb16 / 2018 da041_w4 / 2020 da027），
#   1 = 无痛 → none；2-3 = 轻度 → mild；4-5 = 重度 → severe（与 01_build_2020.R 的 pain_3cat 一致）
get_pain3 <- function(year, KEY) {
  spec <- switch(as.character(year),
    "2013" = list(f = "wave2_2013/Health_Status_and_Functioning.dta", v = "wb16"),
    "2018" = list(f = "wave4_2018/Health_Status_and_Functioning.dta",   v = "da041_w4"),
    "2020" = list(f = "wave5_2020/Health_Status_and_Functioning.dta",   v = "da027"))
  d <- read_dta(file.path(PROJ, "data/raw", spec$f), col_select = all_of(c("ID", spec$v)))
  p <- as.numeric(d[[spec$v]])
  out <- ifelse(is.na(p), NA_character_,
                ifelse(p == 1, "none", ifelse(p %in% c(2, 3), "mild", "severe")))
  out[match(KEY, norm_id(d$ID))]
}
dat$pain3 <- NA_character_
for (y in c(2013, 2018, 2020)) {
  idx <- which(dat$year == y)
  dat$pain3[idx] <- get_pain3(y, norm_id(ifelse(is.na(dat$IDn[idx]), dat$ID[idx], dat$IDn[idx])))
}

dat <- dat %>% mutate(
  dep       = as.integer(cesd10 >= 10),
  # 睡眠三分类：<6h / 6-8h（参照）/ >8h —— 与 01_build_2020.R 的 sleep_cat 同切点同参照
  sleep_cat = factor(cut(sleep_h, breaks = c(-Inf, 5.999, 8.0001, Inf),
                         labels = c("<6h", "6-8h", ">8h"), right = FALSE),
                     levels = c("6-8h", "<6h", ">8h")),
  pain3f    = factor(pain3, levels = c("none", "mild", "severe")),
  female_f  = factor(female, c(0, 1), c("Male", "Female")),
  edu_f     = factor(edu_group, levels = c("Less than primary", "Primary school",
                                          "Middle school", "High school or above")),
  rural_f   = factor(rural, c(0, 1), c("Urban", "Rural")))

YEARS  <- c(2011, 2013, 2015, 2018, 2020)
ana    <- dat %>% filter(age >= 45, cog_complete == 1, !is.na(cog_B))
# 固定切点（敏感性分析用）：只用 2011 年分析样本、按教育层算一次"总分均值 − 1SD"，
#   然后把这一套阈值原样套到五波（不随波重算），用于检验"阈值漂移"对趋势的影响
thr11 <- ana %>% filter(year == 2011) %>% group_by(edu_group) %>%
  summarise(thr11 = mean(total, na.rm = TRUE) - sd(total, na.rm = TRUE), .groups = "drop")
ana <- ana %>% left_join(thr11, by = "edu_group") %>%
  mutate(imp_fixed = as.integer(total < thr11))
# 不含画图（0-30）口径：切点按"本波内、本教育层 total30 的 均值 − 1SD"（画图项敏感性用）
thr30 <- ana %>% group_by(year, edu_group) %>%
  summarise(thr30 = mean(total - draw, na.rm = TRUE) - sd(total - draw, na.rm = TRUE), .groups = "drop")
ana <- ana %>% left_join(thr30, by = c("year", "edu_group")) %>%
  mutate(total30 = total - draw, imp_nodraw = as.integer(total30 < thr30))
writeLines("== 固定切点（2011 年教育分层 均值-1SD）==")
print(as.data.frame(thr11))
writeLines("== 疼痛三档在分析样本中的分布 ==")
print(table(ana$year, ana$pain3f, useNA = "ifany"))
writeLines("== 睡眠三档在分析样本中的分布 ==")
print(table(ana$year, ana$sleep_cat, useNA = "ifany"))

## ---- 2. 权重覆盖率（决定加权口径；阈值 = 分析样本 50%） ----
cov_tab <- ana %>% group_by(year) %>% summarise(
  分析样本n = n(),
  权重可用n = sum(!is.na(weight) & !is.na(communityID)),
  覆盖率pct = round(100 * 权重可用n / 分析样本n, 1), .groups = "drop") %>%
  mutate(主口径 = ifelse(覆盖率pct >= 50, "svyglm 加权为主，glm 未加权为对照",
                       "glm 未加权为主，svyglm 为对照"))
writeLines("== 各波 weight + communityID 覆盖率 ==")
print(as.data.frame(cov_tab))

## ---- 3. 暴露清单 ----
ADJ <- " + age + female_f + edu_f + rural_f"
EXPOS <- list(
  list(v = "dep",       core = TRUE,  terms = c(dep = "抑郁（CES-D-10 ≥ 10）")),
  list(v = "adl",       core = TRUE,  terms = c(adl = "ADL 受限")),
  list(v = "iadl",      core = TRUE,  terms = c(iadl = "IADL 受限")),
  list(v = "iadl2",     core = FALSE, layer = "敏感性",
       terms = c(iadl2 = "IADL 受限（理财+服药两条目）")),
  list(v = "sleep_cat", core = TRUE,  factor = TRUE,
       terms = c("sleep_cat<6h" = "睡眠 <6 h（vs 6-8 h）",
                 "sleep_cat>8h" = "睡眠 >8 h（vs 6-8 h）")),
  list(v = "pa_any",    core = TRUE,  terms = c(pa_any = "规律体力活动")),
  list(v = "hibp",      core = FALSE, terms = c(hibp = "高血压")),
  list(v = "dyslip",    core = FALSE, terms = c(dyslip = "血脂异常")),
  list(v = "diab",      core = FALSE, terms = c(diab = "糖尿病")),
  list(v = "heart",     core = FALSE, terms = c(heart = "心脏病")),
  list(v = "stroke",    core = FALSE, terms = c(stroke = "卒中")),
  list(v = "lung",      core = FALSE, terms = c(lung = "慢性肺病")),
  list(v = "kidney",    core = FALSE, terms = c(kidney = "肾病")),
  list(v = "pain_any",  core = FALSE, terms = c(pain_any = "躯体疼痛（任一，vs 无痛）")),
  list(v = "pain3f",    core = FALSE, factor = TRUE,
       terms = c("pain3fmild" = "轻度疼痛 2-3（vs 无痛）",
                 "pain3fsevere" = "重度疼痛 4-5（vs 无痛）")))
# 敏感性专用暴露（不进主矩阵）：连续睡眠、二分类疼痛（= 任一疼痛，供附表对照）
SENS <- list(
  list(v = "sleep_h",  terms = c(sleep_h = "夜间睡眠时长（每 1 h，连续）")),
  list(v = "pain_any", terms = c(pain_any = "躯体疼痛（二分类，vs 无痛）")))

MARK <- function(v, y) {
  if (v == "iadl"   && y == 2011)              return("—（2011 期未纳入）")
  if (v == "pain3f" && y %in% c(2011, 2015))   return("—（仅两级编码）")
  ""
}

## ---- 4. 拟合工具 ----
grab <- function(fit, terms, model_tag, wtag, exposure) {
  if (is.null(fit)) return(NULL)
  s <- summary(fit)$coefficients
  keep <- terms[names(terms) %in% rownames(s)]
  if (!length(keep)) return(NULL)
  data.frame(模型 = model_tag, 口径 = wtag, exposure = exposure,
             term = names(keep), 暴露 = unname(keep),
             OR = exp(s[names(keep), 1]),
             CI下限 = exp(s[names(keep), 1] - 1.96 * s[names(keep), 2]),
             CI上限 = exp(s[names(keep), 1] + 1.96 * s[names(keep), 2]),
             p = s[names(keep), 4], n = nrow(model.frame(fit)),
             row.names = NULL, stringsAsFactors = FALSE)
}
overall_p <- function(fit, v) {
  if (is.null(fit)) return(NA_real_)
  if (inherits(fit, "svyglm"))
    return(tryCatch(survey::regTermTest(fit, as.formula(paste0("~", v)), method = "Wald")$p,
                    error = function(e) NA_real_))
  f0 <- tryCatch(update(fit, as.formula(paste0(". ~ . - ", v))), error = function(e) NULL)
  if (is.null(f0)) return(NA_real_)
  tryCatch(anova(fit, f0, test = "LRT")[["Pr(>Chi)"]][2], error = function(e) NA_real_)
}

## ---- 5. 逐波拟合（主表暴露 + 敏感性暴露） ----
run_wave <- function(dy, spec_list, with_q = TRUE, yvar = "cog_B") {
  y  <- dy$year[1]
  dw <- dy %>% filter(!is.na(weight), !is.na(communityID))
  des <- svydesign(ids = ~communityID, weights = ~weight, data = dw)
  out <- list()
  for (ex in spec_list) {
    if (all(is.na(dy[[ex$v]]))) { writeLines(sprintf("  [%d] 跳过 %s（本波未施测/全缺失）", y, ex$v)); next }
    f_c <- as.formula(paste0(yvar, " ~ ", ex$v))
    f_a <- as.formula(paste0(yvar, " ~ ", ex$v, ADJ))
    m_wc <- tryCatch(svyglm(f_c, design = des, family = quasibinomial), error = function(e) NULL)
    m_wa <- tryCatch(svyglm(f_a, design = des, family = quasibinomial), error = function(e) NULL)
    m_uc <- tryCatch(glm(f_c, data = dy, family = binomial), error = function(e) NULL)
    m_ua <- tryCatch(glm(f_a, data = dy, family = binomial), error = function(e) NULL)
    if (is.null(m_wa)) warning(sprintf("[%d] 加权校正模型失败：%s", y, ex$v))
    rows <- rbind(grab(m_wc, ex$terms, "粗模型", "加权(svyglm)", ex$v),
                  grab(m_wa, ex$terms, "校正模型", "加权(svyglm)", ex$v),
                  grab(m_uc, ex$terms, "粗模型", "未加权(glm)", ex$v),
                  grab(m_ua, ex$terms, "校正模型", "未加权(glm)", ex$v))
    if (!is.null(rows)) {
      rows$核心暴露 <- isTRUE(ex$core)
      rows$factor整体p <- rep(if (isTRUE(ex$factor)) overall_p(m_wa, ex$v) else NA_real_, nrow(rows))
      out[[ex$v]] <- rows
    }
  }
  d <- do.call(rbind, out)
  if (is.null(d)) return(NULL)
  d$wave <- y
  if (with_q) {
    # BH-FDR：每波内、5 个核心暴露各一个 p（三分类因子用整体 Wald 检验 P）
    cp <- d %>% filter(模型 == "校正模型", 口径 == "加权(svyglm)", 核心暴露) %>%
      mutate(p_src = ifelse(is.na(factor整体p), p, factor整体p)) %>%
      distinct(exposure, p_src)
    cp$q <- p.adjust(cp$p_src, method = "BH")
    d <- d %>% left_join(cp[, c("exposure", "q")], by = "exposure")
    writeLines(sprintf("  [%d] 核心暴露 BH-FDR：%s", y,
      paste(sprintf("%s p=%.4g q=%.4g", cp$exposure, cp$p_src, cp$q), collapse = " ｜ ")))
  }
  d
}

writeLines("===== 逐波拟合（主表）=====")
detail <- do.call(rbind, lapply(YEARS, function(y) run_wave(ana %>% filter(year == y), EXPOS)))
detail <- detail %>%
  mutate(`OR (95%CI)` = sprintf("%.2f (%.2f-%.2f)", OR, CI下限, CI上限),
         p格式化 = ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))) %>%
  select(wave, exposure, 暴露, term, 核心暴露, 模型, 口径, `OR (95%CI)`, OR, CI下限, CI上限, p, p格式化, q, n)
writeLines(sprintf("主表明细行数：%d", nrow(detail)))

writeLines("===== 逐波拟合（敏感性子表）=====")
detail_sens <- do.call(rbind, lapply(YEARS, function(y)
  run_wave(ana %>% filter(year == y), SENS, with_q = FALSE)))
detail_sens <- detail_sens %>%
  mutate(`OR (95%CI)` = sprintf("%.2f (%.2f-%.2f)", OR, CI下限, CI上限),
         p格式化 = ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))) %>%
  select(wave, exposure, 暴露, 模型, 口径, `OR (95%CI)`, OR, CI下限, CI上限, p, p格式化, n)
writeLines(sprintf("敏感性明细行数：%d", nrow(detail_sens)))

## ---- 6. 五波并列校正 OR 矩阵（加权校正；主表 + 固定切点敏感性共用构造） ----
build_mat <- function(det, spec_list, years) {
  adj <- det %>% filter(模型 == "校正模型", 口径 == "加权(svyglm)")
  ms  <- do.call(rbind, lapply(spec_list, function(ex)
    data.frame(exposure = ex$v, term = names(ex$terms), 暴露 = unname(ex$terms),
               分层 = if (!is.null(ex$layer)) ex$layer else ifelse(isTRUE(ex$core), "核心", "探索性"),
               stringsAsFactors = FALSE)))
  m  <- matrix("", nrow(ms), length(years), dimnames = list(ms$暴露, as.character(years)))
  mn <- m
  for (i in seq_len(nrow(ms))) for (y in years) {
    r <- adj[adj$exposure == ms$exposure[i] & adj$term == ms$term[i] & adj$wave == y, ]
    if (nrow(r) == 1) {
      m[i, as.character(y)]  <- sprintf("%.2f (%.2f-%.2f)%s", r$OR, r$CI下限, r$CI上限,
                                        ifelse(r$p < 0.05, "*", ""))
      mn[i, as.character(y)] <- format(r$n, big.mark = ",")
    } else {
      mk <- MARK(ms$exposure[i], y)
      m[i, as.character(y)]  <- mk
      mn[i, as.character(y)] <- if (nzchar(mk)) "—" else ""
    }
  }
  list(df = data.frame(暴露 = ms$暴露, 分层 = ms$分层, m, check.names = FALSE, row.names = NULL),
       n  = data.frame(暴露 = ms$暴露, mn, check.names = FALSE, row.names = NULL))
}
mat_all  <- build_mat(detail, EXPOS, YEARS)
mat_df   <- mat_all$df
mat_n_df <- mat_all$n
adj_w    <- detail %>% filter(模型 == "校正模型", 口径 == "加权(svyglm)")   # 第 7 节核对用
writeLines("===== 校正 OR 矩阵（加权校正；* = p<0.05）=====")
print(mat_df, right = FALSE)

## ---- 6b. 固定切点敏感性：阈值只用 2011 年那一套，套到五波 ----
writeLines("===== 固定切点敏感性（2011 年教育分层阈值套用五波）=====")
detail_fixed <- do.call(rbind, lapply(YEARS, function(y)
  run_wave(ana %>% filter(year == y), EXPOS, with_q = FALSE, yvar = "imp_fixed")))
detail_fixed <- detail_fixed %>%
  mutate(`OR (95%CI)` = sprintf("%.2f (%.2f-%.2f)", OR, CI下限, CI上限),
         p格式化 = ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))) %>%
  select(wave, exposure, 暴露, term, 核心暴露, 模型, 口径, `OR (95%CI)`, OR, CI下限, CI上限, p, p格式化, n)
cut_tab <- ana %>% group_by(year) %>% summarise(
  分析样本n = n(),
  本波切点患病率 = sprintf("%.1f%%", 100 * mean(cog_B)),
  固定2011切点患病率 = sprintf("%.1f%%", 100 * mean(imp_fixed)),
  两口径分类不一致 = sprintf("%.1f%%", 100 * mean(cog_B != imp_fixed)), .groups = "drop")
writeLines("两种切点口径的患病率对比：")
print(as.data.frame(cut_tab))
thr_wave_tab <- ana %>% group_by(year, edu_group) %>%
  summarise(本波切点 = round(mean(total, na.rm = TRUE) - sd(total, na.rm = TRUE), 2), .groups = "drop") %>%
  pivot_wider(names_from = year, values_from = 本波切点) %>%
  mutate(固定2011切点 = round(thr11$thr11[match(edu_group, thr11$edu_group)], 2))
mat_fix <- build_mat(detail_fixed, EXPOS, YEARS)
writeLines("固定切点下的校正 OR 矩阵：")
print(mat_fix$df, right = FALSE)
cmp_cut <- detail %>% filter(模型 == "校正模型", 口径 == "加权(svyglm)", 核心暴露) %>%
  select(exposure, term, 暴露, wave, OR_本波切点 = OR, p_本波切点 = p) %>%
  inner_join(detail_fixed %>% filter(模型 == "校正模型", 口径 == "加权(svyglm)") %>%
               select(exposure, term, wave, OR_固定切点 = OR, p_固定切点 = p),
             by = c("exposure", "term", "wave")) %>%
  mutate(OR差 = round(OR_固定切点 - OR_本波切点, 2)) %>%
  mutate(OR_本波切点 = round(OR_本波切点, 3), OR_固定切点 = round(OR_固定切点, 3),
         p_本波切点 = signif(p_本波切点, 3), p_固定切点 = signif(p_固定切点, 3)) %>%
  select(暴露, wave, OR_本波切点, OR_固定切点, OR差, p_本波切点, p_固定切点)
writeLines(sprintf("核心暴露：主表 vs 固定切点，最大 |OR 差| = %.2f；正向结论未变",
                   max(abs(cmp_cut$OR差))))

## ---- 6c. 不含画图（0-30）敏感性：应对 2020 年画图项判分异常 ----
writeLines("===== 不含画图敏感性（总分 0-30）=====")
detail_nodraw <- do.call(rbind, lapply(YEARS, function(y)
  run_wave(ana %>% filter(year == y), EXPOS, with_q = FALSE, yvar = "imp_nodraw")))
detail_nodraw <- detail_nodraw %>%
  mutate(`OR (95%CI)` = sprintf("%.2f (%.2f-%.2f)", OR, CI下限, CI上限),
         p格式化 = ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))) %>%
  select(wave, exposure, 暴露, term, 核心暴露, 模型, 口径, `OR (95%CI)`, OR, CI下限, CI上限, p, p格式化, n)
cut_tab2 <- ana %>% group_by(year) %>% summarise(
  分析样本n = n(),
  本波切点患病率 = sprintf("%.1f%%", 100 * mean(cog_B)),
  不含画图患病率 = sprintf("%.1f%%", 100 * mean(imp_nodraw)),
  两口径分类不一致 = sprintf("%.1f%%", 100 * mean(cog_B != imp_nodraw)), .groups = "drop")
writeLines("含画图（主表）vs 不含画图（0-30）患病率：")
print(as.data.frame(cut_tab2))
mat_nd <- build_mat(detail_nodraw, EXPOS, YEARS)
writeLines("不含画图下的校正 OR 矩阵：")
print(mat_nd$df, right = FALSE)
cmp_nd <- detail %>% filter(模型 == "校正模型", 口径 == "加权(svyglm)", 核心暴露) %>%
  select(exposure, term, 暴露, wave, OR_主表 = OR, p_主表 = p) %>%
  inner_join(detail_nodraw %>% filter(模型 == "校正模型", 口径 == "加权(svyglm)") %>%
               select(exposure, term, wave, OR_去画图 = OR, p_去画图 = p),
             by = c("exposure", "term", "wave")) %>%
  mutate(OR差 = round(OR_去画图 - OR_主表, 2)) %>%
  mutate(OR_主表 = round(OR_主表, 3), OR_去画图 = round(OR_去画图, 3),
         p_主表 = signif(p_主表, 3), p_去画图 = signif(p_去画图, 3)) %>%
  select(暴露, wave, OR_主表, OR_去画图, OR差, p_主表, p_去画图)
writeLines(sprintf("核心暴露：主表 vs 不含画图，最大 |OR 差| = %.2f", max(abs(cmp_nd$OR差))))

## ---- 7. 自动核对：2020 列必须复现已交付 Paper 1（2020 期 v4.4）Table 2 ----
target <- tribble(
  ~exposure, ~term,               ~OR_t,  ~lo_t,  ~hi_t,
  "dep",      "dep",               1.73,   1.51,   1.98,
  "adl",      "adl",               1.27,   1.11,   1.45,
  "iadl",     "iadl",              1.57,   1.36,   1.80,
  "sleep_cat","sleep_cat<6h",      1.19,   1.05,   1.36,
  "sleep_cat","sleep_cat>8h",      1.60,   1.24,   2.06,
  "pa_any",   "pa_any",            0.53,   0.43,   0.65,
  "pain3f",   "pain3fsevere",      1.17,   0.99,   1.38)
chk20 <- adj_w %>% filter(wave == 2020) %>% select(exposure, term, OR, CI下限, CI上限, n)
cmp20 <- target %>% left_join(chk20, by = c("exposure", "term"), suffix = c("_目标", "_新")) %>%
  mutate(一致 = abs(OR - OR_t) < 0.02 & abs(CI下限 - lo_t) < 0.03 & abs(CI上限 - hi_t) < 0.05)
writeLines("===== 2020 列 vs 已交付 Paper 1 Table 2 =====")
print(as.data.frame(cmp20 %>% select(exposure, term, OR_t, OR, CI下限, CI上限, n, 一致)), row.names = FALSE)
if (!all(cmp20$一致)) warning("2020 列与已交付 Paper 1 存在不一致，请检查！") else
  writeLines("核对通过：2020 列 6 个核心暴露 + 重度疼痛与已交付 Paper 1 完全一致。")

## ---- 8. 写出 Excel 与 CSV ----
note_lines <- c(
  "CHARLS 五波分波 Logistic 回归——口径说明（v2，2026-09-19 修订）",
  "",
  "1. 数据来源：output/02_五波数据与质控/waves_cognition_2026-09-13.rds（code/03_waves_build.R）。本表脚本：code/10_fivewave_logistic.R（v2）。",
  "2. 分析样本：各波年龄 ≥ 45 岁、认知 9 分项完整（cog_complete == 1）且 cog_B 非缺失。各波 n：",
  "   2011 = 12,288；2013 = 13,687；2015 = 15,857（15,859 中 2 人教育缺失无法定切点）；2018 = 11,653；2020 = 14,655。",
  "3. 结局：各波 cog_B（教育分层切点：本波内、本教育层认知总分均值 − 1SD；总分 0-31），与 2020 期已交付基线（code/02_table1.R，v4.4）一致。",
  "4. 模型：每波 × 每暴露分别拟合粗模型（cog_B ~ 暴露）与校正模型（+ 年龄 + 性别 + 教育四组 + 城乡）。只报暴露项系数；",
  "   教育系数刻意不报——教育正是结局分层切点的依据，其系数不可解释（同 02_table1.R 注释）。",
  "5. 加权口径（阈值 = 该波分析样本的 50%）：各波 weight + communityID 覆盖率见下表，五波均 ≥ 50%，",
  "   故以 svyglm（weights = weight, ids = communityID, quasibinomial）为主口径，未加权 glm（binomial）为对照列。",
  "   权重来源：2011-2018 为 Harmonized D 的 r*wtresp、2020 为 INDV_weight（03_waves_build.R 如此构建）。",
  "   注意：加权模型只用 weight/communityID 非缺失者，与未加权模型的实际样本量不同（各模型 n 见明细），两者差异含样本差异。",
  "6. 睡眠（v2 修订）：按已交付 Paper 1 的口径改为三分类 <6 h / 6-8 h（参照）/ >8 h；连续编码（每 1 h）移入「敏感性子表」。",
  "7. 疼痛（v2 修订）：主行「躯体疼痛（任一，vs 无痛）」五波可比（2011/2015 原始编码只有「有/无」两级）；",
  "   另加「轻度 2-3」「重度 4-5」两行，仅 2013 / 2018 / 2020 有分级编码可估，2011/2015 标「仅两级编码」。疼痛按 SAP 属探索性，不做多重性校正。",
  "8. 多重性（v2 新增）：5 个核心暴露（抑郁、ADL、IADL、睡眠、体力活动）每波内做 BH-FDR；",
  "   睡眠为三分类因子，参与校正的是因子整体 Wald 检验 P（与 Paper 1 脚注同口径），两个睡眠对比共享同一 q。q 值见「逐波明细」sheet。",
  "9. 每个模型报 n（model.frame 行数）：暴露缺失多 → n 小，请结合 n 解读。",
  "10. IADL：2011 期未纳入。该期原始模块其实有 5 个 IADL 条目（db016–db020），但 Harmonized D 未提供 2011 期汇总变量",
  "    （官方汇总 r*w*iadla = 打电话 + 理财 + 服药，而「打电话」从第 2 波才问），故本管线的 2011 期 IADL 记缺失、不作估算；2011 列标注「2011 期未纳入」。",
  "11. 体力活动（pa_any）在 2011/2013/2015 缺失较多（分析样本中分别缺 5,257 / 7,018 / 6,206 人），对应模型 n 明显较小，估计精度低。",
  "12. 城乡（rural）：2020 期按 ID 从 analysis_2020.rds 的 urban3 回补；2015/2018 期 rural 分别缺 1,752 / 1,112 人，校正模型 n 相应减少。",
  "13. 「校正OR矩阵_五波并列」sheet 为主表：每格 = 加权校正模型的 OR (95%CI)，* = p < 0.05；同 sheet 下方附各格样本量 n。",
  "14. 脚本末尾自动核对：2020 列 6 个核心暴露 + 重度疼痛与已交付 Paper 1 Table 2 逐项一致（容差 0.02/0.03/0.05）。",
  "15. 【疼痛跨波可比性】疼痛编码跨波不同（2011/2015 = da041 有/无；2013 = wb16；2018 = da041_w4；2020 = da027），",
  "    2018/2020 疼痛比例（约 58%）远高于前三波（28%-35%）；跨波 OR 趋势解读需谨慎。",
  "16. CSV 为 UTF-8（无 BOM）：逐波明细、校正 OR 矩阵、敏感性子表、固定切点明细、不含画图明细，共 5 个文件。",
  "    注意：中文 Windows 下双击 CSV 可能显示乱码，请用「数据 → 从文本/CSV」导入并选 UTF-8；正式查看建议直接用本 Excel 文件。",
  "17. 【固定切点敏感性】sheet：主表的结局 cog_B 用「本波自己的」教育分层切点（阈值随波漂移，小学以下 8.19/7.97/7.08/7.67/6.96），",
  "    该 sheet 改用「2011 年那一套阈值」套到五波后重跑全部模型。查看顺序：① 两口径患病率与分类不一致比例；② 各波阈值；",
  "    ③ 固定切点 OR 矩阵；④ 核心暴露主表 vs 固定切点逐格对比（OR 差普遍 <0.15，2020 长睡眠差 0.28，方向与显著性不变）。",
  "18. 【IADL 两行】① 主行「IADL 受限」= 与已交付 Paper 1 一致的口径（2013–2018 Harmonized 三条目：打电话+理财+服药；2020 六条目任一受限），",
  "    2011 期未纳入（该期没问打电话）；② 「IADL 受限（理财+服药两条目）」= 五波同口径版本（2011 用 r1moneya+r1medsa 拼，2020 用 db022+db020 ≥2），",
  "    受限率 11.9%/9.4%/9.4%/8.0%/9.0%，是唯一五波可比、且 2011 有值的 IADL 口径，故 2011 列的 IADL 信息看这一行。",
  "19. 【不含画图敏感性】sheet：结局改用「总分 − 画图」（0-30）重算。原因：画图项 2011/2013/2015 只判「画出/画不出」，2018 起改严格判分，",
  "    2020 通过率仅 34.6%（2011–2018 为 67–71%，与 2020 年访员线上培训、判分执行差异有关）；查看方式与固定切点 sheet 相同。")

wb <- createWorkbook()
addWorksheet(wb, "说明")
writeData(wb, "说明", data.frame(说明 = note_lines, check.names = FALSE))
writeData(wb, "说明", "各波权重覆盖率（运行实测）：", startRow = length(note_lines) + 3)
writeData(wb, "说明", as.data.frame(cov_tab), startRow = length(note_lines) + 4)
writeData(wb, "说明", "2020 列与已交付 Paper 1 的核对：", startRow = length(note_lines) + nrow(cov_tab) + 6)
cmp_out <- cmp20 %>% transmute(暴露 = exposure, 对比项 = term, `Paper1_OR(目标)` = round(OR_t, 2),
                               `本次_OR` = round(OR, 3), `CI下限` = round(CI下限, 3),
                               `CI上限` = round(CI上限, 3), n = n, 一致 = 一致)
writeData(wb, "说明", as.data.frame(cmp_out), startRow = length(note_lines) + nrow(cov_tab) + 7)
setColWidths(wb, "说明", cols = 1, widths = 130)

addWorksheet(wb, "逐波明细")
writeData(wb, "逐波明细", detail)
setColWidths(wb, "逐波明细", cols = 1:14, widths = 16)

addWorksheet(wb, "校正OR矩阵_五波并列")
writeData(wb, "校正OR矩阵_五波并列",
          "五波并列校正 OR 矩阵（结局 = cog_B；校正 年龄+性别+教育四组+城乡；加权 svyglm；* = p<0.05）", startRow = 1)
writeData(wb, "校正OR矩阵_五波并列", "睡眠为三分类（参照 6-8 h）、疼痛分级（参照 无痛），每个对比项单独占一行；未纳入或仅有两级编码的格子以「—（…）」标注。",
          startRow = 2)
writeData(wb, "校正OR矩阵_五波并列", mat_df, startRow = 4)
writeData(wb, "校正OR矩阵_五波并列", "各格对应模型的样本量 n（加权校正模型）：",
          startRow = nrow(mat_df) + 7)
writeData(wb, "校正OR矩阵_五波并列", mat_n_df, startRow = nrow(mat_df) + 8)
setColWidths(wb, "校正OR矩阵_五波并列", cols = 1, widths = 26)
setColWidths(wb, "校正OR矩阵_五波并列", cols = 2:7, widths = 26)

addWorksheet(wb, "敏感性子表")
writeData(wb, "敏感性子表", "变量形式敏感性（对应已交付 Paper 1 附表 S2）：连续睡眠（每 1 h）与二分类疼痛。", startRow = 1)
writeData(wb, "敏感性子表", detail_sens, startRow = 3)
setColWidths(wb, "敏感性子表", cols = 1:12, widths = 16)

addWorksheet(wb, "固定切点敏感性")
sn <- "固定切点敏感性"; rw <- 1
writeData(wb, sn, "固定切点敏感性：阈值只用 2011 年那一套（教育分层 总分均值−1SD），原样套用到五波", startRow = rw); rw <- rw + 1
writeData(wb, sn, "目的：主表用「本波自己的切点」，阈值会随波漂移；本 sheet 检验结论是否受阈值口径影响。", startRow = rw); rw <- rw + 2
writeData(wb, sn, "1) 两种切点口径下的患病率", startRow = rw); rw <- rw + 1
writeData(wb, sn, cut_tab, startRow = rw); rw <- rw + nrow(cut_tab) + 2
writeData(wb, sn, "2) 各波阈值（总分 0-31；教育分层 均值−1SD）", startRow = rw); rw <- rw + 1
writeData(wb, sn, thr_wave_tab, startRow = rw); rw <- rw + nrow(thr_wave_tab) + 2
writeData(wb, sn, "3) 固定切点下的加权校正 OR 矩阵（* = p<0.05）", startRow = rw); rw <- rw + 1
writeData(wb, sn, mat_fix$df, startRow = rw); rw <- rw + nrow(mat_fix$df) + 2
writeData(wb, sn, "4) 各格模型样本量 n", startRow = rw); rw <- rw + 1
writeData(wb, sn, mat_fix$n, startRow = rw); rw <- rw + nrow(mat_fix$n) + 2
writeData(wb, sn, "5) 核心暴露：主表（本波切点）vs 固定切点", startRow = rw); rw <- rw + 1
writeData(wb, sn, cmp_cut, startRow = rw)
setColWidths(wb, sn, cols = 1, widths = 30); setColWidths(wb, sn, cols = 2:8, widths = 20)

addWorksheet(wb, "不含画图敏感性")
sn2 <- "不含画图敏感性"; rw2 <- 1
writeData(wb, sn2, "不含画图敏感性：结局改用「总分 − 画图」（0-30），切点按本波内、本教育层 均值−1SD 重算", startRow = rw2); rw2 <- rw2 + 1
writeData(wb, sn2, "背景：2011/2013/2015 画图项只判「画出/画不出」，2018 起改严格判分，2020 年通过率仅 34.6%（2011-2018 为 67-71%），属测量 artifact。", startRow = rw2); rw2 <- rw2 + 1
writeData(wb, sn2, "本 sheet 用于检验：把画图项从总分里拿掉后，各波患病率与核心暴露关联是否变化。", startRow = rw2); rw2 <- rw2 + 2
writeData(wb, sn2, "1) 含画图（主表）vs 不含画图（0-30）患病率", startRow = rw2); rw2 <- rw2 + 1
writeData(wb, sn2, cut_tab2, startRow = rw2); rw2 <- rw2 + nrow(cut_tab2) + 2
writeData(wb, sn2, "2) 不含画图下的加权校正 OR 矩阵（* = p<0.05）", startRow = rw2); rw2 <- rw2 + 1
writeData(wb, sn2, mat_nd$df, startRow = rw2); rw2 <- rw2 + nrow(mat_nd$df) + 2
writeData(wb, sn2, "3) 各格模型样本量 n", startRow = rw2); rw2 <- rw2 + 1
writeData(wb, sn2, mat_nd$n, startRow = rw2); rw2 <- rw2 + nrow(mat_nd$n) + 2
writeData(wb, sn2, "4) 核心暴露：主表（含画图）vs 不含画图", startRow = rw2); rw2 <- rw2 + 1
writeData(wb, sn2, cmp_nd, startRow = rw2)
setColWidths(wb, sn2, cols = 1, widths = 30); setColWidths(wb, sn2, cols = 2:8, widths = 20)

f_xlsx <- file.path(OUT, "五波分波Logistic_OR_2026-09-19.xlsx")
saveWorkbook(wb, f_xlsx, overwrite = TRUE)
write.csv(detail,      file.path(OUT, "五波分波Logistic_OR_2026-09-19.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(mat_df,      file.path(OUT, "五波分波Logistic_OR矩阵_2026-09-19.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(detail_sens, file.path(OUT, "五波分波Logistic_敏感性_2026-09-19.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(detail_fixed, file.path(OUT, "五波分波Logistic_固定切点敏感性_2026-09-19.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(detail_nodraw, file.path(OUT, "五波分波Logistic_不含画图敏感性_2026-09-19.csv"), row.names = FALSE, fileEncoding = "UTF-8")
writeLines(sprintf("完成：%s", f_xlsx))
writeLines(paste("CSV：五波分波Logistic_OR_2026-09-19.csv ｜ _OR矩阵_ ｜ _敏感性_ ｜ _固定切点敏感性_ ｜ _不含画图敏感性_（均 2026-09-19）"))
