# =============================================================================
# 01_build_2020.R  --  CHARLS 2020 (Wave 5) 分析数据集构建
# 项目：中国中老年人群认知功能受损的患病率及术前相关多维脆弱性因素
# 输入：data/raw/wave5_2020/*.dta
# 输出：output/01_2020单期基线/analysis_2020.rds, output/01_2020单期基线/analysis_2020.csv, output/01_2020单期基线/build_log.txt
# 说明：raw/ 只读；本脚本不修改任何原始文件
# 版本：v3.2（2026-09-12）——教育改为 zredu + ba010 合并（2020 期内部互补，合并后零缺失）
# =============================================================================

## ---- 0. 环境 ----
Sys.setlocale("LC_ALL", "zh_CN.UTF-8")          # 修复中文路径
.libPaths(c("D:/R/library", .libPaths()))
suppressPackageStartupMessages({
  library(haven); library(dplyr); library(tidyr)
})

proj <- "D:/BaiduSyncdisk/科研/徐中豪"
raw_dir <- file.path(proj, "data/raw/wave5_2020")
out_dir <- file.path(proj, "output", "01_2020单期基线")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

log_con <- file(file.path(out_dir, "build_log.txt"), open = "wt", encoding = "UTF-8")
say <- function(...) { msg <- paste0(...); cat(msg, "\n"); writeLines(msg, log_con) }

say("== CHARLS 2020 数据构建 ", format(Sys.time()), " ==")
say("口径备注（v3.2，2026-09-12）：定向 / 连续减 7 / 画图中『不知道(997)』一律按答错计 0（算已作答）；")
say("         『拒绝(999)/不适用(-1)』记缺失；画图『未评估(3)』记缺失；")
say("         词语回忆中『没有回忆起任何词』计 0 分，『拒绝回忆/无法测试』记缺失；")
say("         主分析限定认知 9 分项完整且教育已知者；切点采用教育分层（<本层均值−1SD）。")

## ---- 1. 读取原始数据（只选所需列） ----
hs_file   <- file.path(raw_dir, "Health_Status_and_Functioning.dta")
demo_file <- file.path(raw_dir, "Demographic_Background.dta")
wt_file   <- file.path(raw_dir, "Weights.dta")

hs_cols <- c(
  "ID",
  paste0("dc00", 1:5),                                    # 定向：年/季/日/星期/月
  paste0("dc007_", 1:5), paste0("dc007_", 1:5, "_1"),     # 连续减7
  "dc009",                                                # 画图
  paste0("dc012_s", 1:13), paste0("dc013_s", 1:12),       # 即时回忆 试次1、2（含 s11/s12/s13 标记位）
  paste0("dc014_s", 1:12),                                # 即时回忆 试次3
  paste0("dc028_s", 1:12),                                # 延迟回忆
  paste0("dc0", 16:25),                                   # CES-D-10
  paste0("da002_", 1:15, "_"),                            # 慢病（上期已患，比较）
  paste0("da003_", 1:15, "_"),                            # 慢病（医生是否告知过）
  paste0("db", sprintf("%03d", c(1,3,5,7,9,11,12,14,16,18,20,22))),  # ADL/IADL
  "da027", "da030", "da031",                              # 疼痛、睡眠
  "da046", "da047", "zsmoke", "da051",                    # 吸烟、饮酒（da051 为过去一年饮酒频率）
  "da032_1_", "da032_2_", "da032_3_"                      # 体力活动（高/中/低强度 ≥10 分钟）
)
hs <- read_dta(hs_file, col_select = all_of(hs_cols))

demo_cols <- c("ID","householdID","communityID","proxy_2","ba008","ba009","ba010","ba011",
               "xrage","xrgender","zrgender","zredu","xrpartner","zrbirthyear")
demo <- read_dta(demo_file, col_select = all_of(demo_cols))
wt   <- read_dta(wt_file, col_select = all_of(c("ID","HH_weight","HH_weight_ad1","INDV_weight","INDV_weight_ad2")))

say("HS 行数: ", nrow(hs), " | Demo 行数: ", nrow(demo), " | 权重行数: ", nrow(wt))
say("注：Demographic 模块 ", nrow(demo), " 人，其中 ", nrow(demo) - nrow(hs),
    " 人无 Health Status 模块记录；分析以 HS 模块为起点。")
say("教育来源：zredu 非缺失 ", sum(!is.na(demo$zredu)), " ｜ ba010 非缺失 ", sum(!is.na(demo$ba010)),
    " ｜ 合并后缺失 ", sum(is.na(dplyr::coalesce(as.numeric(demo$zredu), as.numeric(demo$ba010)))), "")

## ---- 2. 工具函数 ----
# 特殊码 -> NA（997 不知道 / 999 拒绝 / -1 不适用）
to_na <- function(x, codes = c(997, 999, -1)) { x <- as.numeric(x); x[x %in% codes] <- NA; x }

# 词表回忆计分：s1..s10 中取值为该词在词表中的编号(1-10)，未回忆=0
# "没有回忆起任何词"(s11) -> 0 分；"拒绝回忆"(s12) / "无法测试"(s13) -> 缺失
recall_score <- function(df, prefix, ns = 10, refuse_cols = character(0)) {
  m <- as.matrix(df[paste0(prefix, "_s", seq_len(ns))])
  m <- apply(m, 2, as.numeric)
  cnt <- rowSums(m >= 1 & m <= 10, na.rm = TRUE)
  asked <- rowSums(!is.na(m)) > 0
  cnt[!asked] <- NA
  for (fc in refuse_cols) {
    if (fc %in% names(df)) {
      v <- as.numeric(df[[fc]])
      cnt[!is.na(v) & v > 0] <- NA            # 拒绝 / 无法测试 -> 缺失
    }
  }
  as.numeric(cnt)
}

## ---- 3. 认知维度 ----
# 3.1 定向力 dc001-dc005（1=正确；2 错误 与 997 不知道 均计 0；999 拒绝 / -1 不适用 = 缺失）
orient_items <- lapply(paste0("dc00", 1:5), function(v) {
  x <- as.numeric(hs[[v]])
  x[x %in% c(999, -1)] <- NA                          # 仅拒绝/不适用记缺失
  ifelse(is.na(x), NA, ifelse(x == 1, 1, 0))          # 错误与不知道均计 0（HRS 惯例）
})
orient_df <- as.data.frame(orient_items); names(orient_df) <- paste0("orient", 1:5)

# 3.2 连续减7 dc007_1..5（答案在 _1_ 变量；按官方 Harmonized 逻辑相对前一答案判定）
# 【v3 修正】作答判定 = 有数值 或 明确回答"不知道(997)"；997 计为答错（不加分）。
#            999=拒绝 / -1=不适用 / 未测 -> 缺失。
s7_num <- lapply(1:5, function(i) {
  x <- as.numeric(hs[[paste0("dc007_", i, "_1")]])
  x[!is.na(x) & (x < 0 | x > 1000)] <- NA       # 明显录入异常 / 哨兵值
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
  ok <- cbind(ok1, ok2, ok3, ok4, ok5) & do.call(cbind, s7_eng)   # 未作答的步不计分
  score   <- rowSums(ok, na.rm = TRUE)
  any_ans <- Reduce(`|`, s7_eng)
  score[!any_ans] <- NA
  as.numeric(score)
})

# 3.3 画图 dc009（1 正确→1；2 错误、997 不知道→0（答错）；3 未评估 / 999 拒绝 / 未测→缺失）
draw_x <- as.numeric(hs$dc009); draw_x[draw_x %in% c(999, -1)] <- NA
draw <- ifelse(is.na(draw_x), NA,
               ifelse(draw_x == 1, 1, ifelse(draw_x %in% c(2, 997), 0, NA)))

# 3.4 词语回忆（v3：拒绝回忆 / 无法测试 -> 缺失；"一个都没回忆出" -> 0 分）
imrc1 <- recall_score(hs, "dc012", refuse_cols = c("dc012_s12","dc012_s13"))  # 即时回忆 试次1
imrc2 <- recall_score(hs, "dc013", refuse_cols = c("dc013_s12"))
imrc3 <- recall_score(hs, "dc014", refuse_cols = c("dc014_s12"))
imrc_mean <- rowMeans(cbind(imrc1, imrc2, imrc3), na.rm = TRUE)   # 三次学习试次均值
dlrc  <- recall_score(hs, "dc028", refuse_cols = c("dc028_s12"))   # 延迟回忆

# 3.5 合成分数（遵循 Harmonized CHARLS 惯例：rowtotal 忽略缺失项；全部缺失才记 NA）
sum_avail <- function(mat) {
  mat <- as.matrix(mat)
  out <- rowSums(mat, na.rm = TRUE)
  out[rowSums(!is.na(mat)) == 0] <- NA
  as.numeric(out)
}
orient <- sum_avail(orient_df)                                   # 定向 0-5
memory <- sum_avail(cbind(imrc1, dlrc))                          # 情景记忆 0-20
mental <- sum_avail(cbind(orient_df, ser7, draw))                # 精神状态 0-11
global <- ifelse(is.na(memory) | is.na(mental), NA, memory + mental)   # 总体认知 0-31
# 完整病例标记（9 个认知分项：即时回忆、延迟回忆、定向 5 项、连续减 7、画图，全部非缺失）
cog_mat <- cbind(imrc1, dlrc, orient_df, ser7, draw)
cog_items_avail <- rowSums(!is.na(cog_mat))
cog_complete <- as.integer(cog_items_avail == ncol(cog_mat))

## ---- 4. 抑郁 CES-D-10（dc016-dc025；1-4 -> 0-3；dc020、dc023 反向） ----
cesd_items <- paste0("dc0", 16:25)
cesd_raw <- as.data.frame(lapply(cesd_items, function(v) to_na(hs[[v]])))
names(cesd_raw) <- cesd_items
cesd <- as.data.frame(lapply(names(cesd_raw), function(v) {
  x <- cesd_raw[[v]]
  if (v %in% c("dc020","dc023")) { 4 - x } else { x - 1 }   # 反向题
}))
names(cesd) <- cesd_items
cesd10 <- rowSums(cesd)
cesd10[rowSums(!is.na(cesd)) < 10] <- NA                     # 要求 10 题完整
cesd_avail <- rowSums(!is.na(cesd))
dep <- ifelse(is.na(cesd10), NA, as.integer(cesd10 >= 10))    # 常用界值 ≥10

## ---- 5. 慢病（da002 与 da003 联合判定） ----
dis_labels <- c("hypertension","dyslipidemia","diabetes","cancer","lung_disease",
                "liver_disease","heart_disease","stroke","kidney_disease","digestive_disease",
                "psychiatric","memory_disease","parkinson","arthritis","asthma")
dis <- data.frame(ID = hs$ID)
for (i in 1:15) {
  prev <- to_na(hs[[paste0("da002_", i, "_")]])      # 1/2/3 = 上期已患病；99 = 上期未患
  ever <- to_na(hs[[paste0("da003_", i, "_")]])      # 1 = 医生曾告知；2 = 否
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
  x <- to_na(hs[[v]]); ifelse(is.na(x), NA, as.integer(x >= 2))   # ≥2 表示有困难
}))
iadl_x <- as.data.frame(lapply(iadl_items, function(v) {
  x <- to_na(hs[[v]]); ifelse(is.na(x), NA, as.integer(x >= 2))
}))
names(adl_x) <- adl_items; names(iadl_x) <- iadl_items
adl_count  <- rowSums(adl_x,  na.rm = TRUE);  adl_count[rowSums(!is.na(adl_x))  == 0] <- NA
iadl_count <- rowSums(iadl_x, na.rm = TRUE); iadl_count[rowSums(!is.na(iadl_x)) == 0] <- NA
adl_any  <- ifelse(is.na(adl_count),  NA, as.integer(adl_count  >= 1))
iadl_any <- ifelse(is.na(iadl_count), NA, as.integer(iadl_count >= 1))

## ---- 7. 疼痛 / 睡眠 / 生活方式 ----
pain <- to_na(hs$da027)
pain_any <- ifelse(is.na(pain), NA, as.integer(pain >= 2))
# 疼痛三档（2026-09-10 决策：二分化掩盖 J 形剂量反应）：1 无 / 2-3 轻度 / 4-5 重度
pain_3cat <- ifelse(is.na(pain), NA,
                    ifelse(pain == 1, 0L, ifelse(pain %in% c(2, 3), 1L, 2L)))
sleep_h  <- to_na(hs$da030); sleep_h[!is.na(sleep_h) & (sleep_h < 0 | sleep_h > 24)] <- NA
nap_min  <- to_na(hs$da031); nap_min[!is.na(nap_min) & (nap_min < 0 | nap_min > 600)] <- NA
sleep_cat <- cut(sleep_h, breaks = c(-Inf, 5.999, 8.0001, Inf),
                 labels = c("<6h","6-8h",">8h"), right = FALSE)
# 睡眠质量代理：CES-D 条目 dc022"睡眠不安"（1-4），取 3-4（每周 ≥3-4 天）为睡眠不安
# 注意：dc022 同时计入 CES-D 总分，回归模型中二者并列时需在 SAP 中说明
slp_restless_x <- to_na(hs$dc022)
sleep_restless <- ifelse(is.na(slp_restless_x), NA, as.integer(slp_restless_x >= 3))

## ---- 7b. 吸烟 / 饮酒 / 体力活动 ----
# 吸烟：zsmoke==1（上期吸烟）者本期免问 da046，直接答 da047（1 仍在吸 / 2 已戒 / 3 从未）；
#       其余人先答 da046（1 曾吸 / 2 从未），曾吸者再答 da047。以 da047 为最终状态。
da046v <- to_na(hs$da046); da047v <- to_na(hs$da047)
smoke <- rep(NA_integer_, nrow(hs))
smoke[!is.na(da047v) & da047v == 3] <- 0L          # 从未
smoke[!is.na(da046v) & da046v == 2] <- 0L          # 从未
smoke[!is.na(da047v) & da047v == 2] <- 1L          # 已戒烟
smoke[!is.na(da047v) & da047v == 1] <- 2L          # 当前吸烟
# 饮酒（过去一年，da051）：1 ≥每月一次 / 2 少于每月 / 3 未饮酒
drink_freq <- to_na(hs$da051)
drink_any  <- ifelse(is.na(drink_freq), NA, as.integer(drink_freq %in% c(1, 2)))
# 体力活动（da032，每次 ≥10 分钟）：1 是 / 2 否
pa_vigorous <- { x <- to_na(hs$da032_1_); ifelse(is.na(x), NA, as.integer(x == 1)) }
pa_moderate <- { x <- to_na(hs$da032_2_); ifelse(is.na(x), NA, as.integer(x == 1)) }
pa_light    <- { x <- to_na(hs$da032_3_); ifelse(is.na(x), NA, as.integer(x == 1)) }
pa_mat <- cbind(pa_vigorous, pa_moderate, pa_light)
pa_any <- ifelse(rowSums(!is.na(pa_mat)) == 0, NA,
                 as.integer(rowSums(pa_mat == 1, na.rm = TRUE) > 0))

## ---- 7c. BMI：不纳入 ----
# 2020 期核心模块无身高体重测量；Harmonized D 最近一次实测为 2015 年（r3mbmi），
# 跨 5 年混用存在时间错位与非随机缺失问题。经决策（2026-09-10）：BMI 不纳入本研究，
# 局限性中注明数据不可得。

## ---- 8. 人口学 ----
urba <- to_na(demo$ba008)      # 1 城市中心/镇中心；2 城乡结合部；3 村；4 特殊区域
hukou <- to_na(demo$ba009)     # 1 农业；2 非农业；3 居民户口；4 无户口
df <- data.frame(
  ID          = hs$ID,
  # 认知
  imrc        = imrc1, imrc2 = imrc2, imrc3 = imrc3, imrc_mean = imrc_mean,
  dlrc        = dlrc, memory = memory, orient = orient,
  ser7        = ser7, draw = draw, mental = mental, global = global,
  cog_items_avail = cog_items_avail, cog_complete = cog_complete,
  # 心理/功能
  cesd10 = cesd10, cesd_avail = cesd_avail, depression = dep,
  adl_count = adl_count, adl_any = adl_any, iadl_count = iadl_count, iadl_any = iadl_any,
  pain = pain, pain_any = pain_any, pain_3cat = pain_3cat,
  sleep_h = sleep_h, nap_min = nap_min, sleep_cat = sleep_cat,
  sleep_restless = sleep_restless,
  # 生活方式
  smoke = smoke, drink_freq = drink_freq, drink_any = drink_any,
  pa_vigorous = pa_vigorous, pa_moderate = pa_moderate, pa_light = pa_light, pa_any = pa_any,
  stringsAsFactors = FALSE
)
df <- cbind(df, dis[dis_labels], n_chronic = dis$n_chronic, multimorbidity = dis$multimorbidity)

demo_small <- data.frame(
  ID = demo$ID, age = demo$xrage,
  female = ifelse(is.na(demo$xrgender), NA, as.integer(demo$xrgender == 2)),
  # 教育（2026-09-12 修订）：zredu 与 ba010 在 2020 期完全互补——
  #   zredu 非缺失 17,687（既往波次结转）、ba010 非缺失 1,708（本期新问），零重叠、零缺口，合计 19,395。
  #   合并后教育 0 缺失，故不再需要跨期回填；两变量为同一 11 级编码（已按出生队列分布比对确认）。
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

## ---- 9. 分析样本定义（v3：主分析 = 认知 9 分项完整） ----
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
say("样本流程：全部 ", n0, " -> 年龄≥45 且非代理 ", n1,
    " -> 认知 9 分项完整 ", n2, "（其中教育已知、进入主分析 ", n2b,
    "；教育缺失 ", n2 - n2b, "），认知未测全（仅作敏感性） ", n1 - n2)

## ---- 10. 认知结局定义（切点与 z 分数一律以主分析样本为参照人群） ----
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
  # 主定义（方案 B，2026-09-10 定）：按教育程度分层切点 = 低于本教育层均值 − 1SD
  cog_low_edu = ifelse(is.na(global) | is.na(thr_edu), NA, as.integer(global < thr_edu)),
  # 敏感性定义 A：主分析样本整体均值 − 1SD
  cog_low_1sd = ifelse(is.na(global), NA, as.integer(global < thr_all)),
  # 敏感性定义 C：最低四分位
  cog_low_q1  = ifelse(is.na(global), NA, as.integer(global <= q1_all))
)
dat$cog_impair <- dat$cog_low_edu          # 主结局别名
say("切点：认知完整者整体 均值−1SD = ", round(thr_all, 2), "；最低四分位 = ", round(q1_all, 2))
say("教育分层切点：", paste(sprintf("%s %.2f", thr_edu_tab$edu_group, thr_edu_tab$thr_edu), collapse = " | "))

## ---- 11. 保存 ----
saveRDS(dat, file.path(out_dir, "analysis_2020.rds"))
write.csv(dat, file.path(out_dir, "analysis_2020.csv"), row.names = FALSE, fileEncoding = "UTF-8")

## ---- 12. 构建日志：关键变量分布 ----
say("\n---- 关键变量分布（构建后） ----")
say("分项非缺失数：imrc1=", sum(!is.na(dat$imrc)), " dlrc=", sum(!is.na(dat$dlrc)),
    " orient=", sum(!is.na(dat$orient)), " ser7=", sum(!is.na(dat$ser7)),
    " draw=", sum(!is.na(dat$draw)))
d1 <- dat[dat$primary == 1, ]
say("认知完整者 n=", nrow(d1), "（其中教育已知、进入主分析 ", sum(!is.na(d1$edu_group)),
    "；教育缺失 ", sum(is.na(d1$edu_group)), "）。以下分布统计基于认知完整者（n=", nrow(d1), "）")
say("总体认知 global（认知完整者）: mean=", round(mean(d1$global, na.rm=TRUE),2),
    " sd=", round(sd(d1$global, na.rm=TRUE),2))
say("情景记忆 memory: mean=", round(mean(d1$memory, na.rm=TRUE),2),
    " | 即时(试次1) mean=", round(mean(d1$imrc, na.rm=TRUE),2),
    " | 延迟 mean=", round(mean(d1$dlrc, na.rm=TRUE),2))
say("精神状态 mental: mean=", round(mean(d1$mental, na.rm=TRUE),2),
    " | 定向 mean=", round(mean(d1$orient, na.rm=TRUE),2),
    " | 连续减7 mean=", round(mean(d1$ser7, na.rm=TRUE),2),
    " | 画图 mean=", round(mean(d1$draw, na.rm=TRUE),2))
say("CES-D-10: n=", sum(!is.na(d1$cesd10)), " mean=", round(mean(d1$cesd10, na.rm=TRUE),2),
    " | 抑郁(≥10) 患病率=", round(100*mean(d1$depression, na.rm=TRUE),1), "%")
say("疼痛(任何): ", round(100*mean(d1$pain_any, na.rm=TRUE),1), "%",
    " | 三档[无/轻/重]: ", sum(d1$pain_3cat==0, na.rm=TRUE), " / ",
    sum(d1$pain_3cat==1, na.rm=TRUE), " / ", sum(d1$pain_3cat==2, na.rm=TRUE),
    " | 睡眠时长 mean=", round(mean(d1$sleep_h, na.rm=TRUE),2), "h")
say("ADL受限: ", round(100*mean(d1$adl_any, na.rm=TRUE),1), "%",
    " | IADL受限: ", round(100*mean(d1$iadl_any, na.rm=TRUE),1), "%")
say("吸烟: 从未 ", sum(d1$smoke==0, na.rm=TRUE), " / 已戒 ", sum(d1$smoke==1, na.rm=TRUE),
    " / 当前 ", sum(d1$smoke==2, na.rm=TRUE), " / 缺失 ", sum(is.na(d1$smoke)))
say("饮酒(过去一年): ≥每月 ", sum(d1$drink_freq==1, na.rm=TRUE), " / <每月 ", sum(d1$drink_freq==2, na.rm=TRUE),
    " / 未饮 ", sum(d1$drink_freq==3, na.rm=TRUE), " / 缺失 ", sum(is.na(d1$drink_freq)))
say("体力活动(任何,≥10分钟/次): ", round(100*mean(d1$pa_any, na.rm=TRUE),1), "%",
    " | 高强度: ", round(100*mean(d1$pa_vigorous, na.rm=TRUE),1), "%")
say("睡眠不安(≥3-4天/周): ", round(100*mean(d1$sleep_restless, na.rm=TRUE),1), "%")
say("认知受损（认知完整者；主定义按教育已知者计）：主定义[教育分层切点] ", round(100*mean(d1$cog_low_edu, na.rm=TRUE),1), "%",
    " | 敏感性A[整体均值−1SD] ", round(100*mean(d1$cog_low_1sd, na.rm=TRUE),1), "%",
    " | 敏感性C[最低四分位] ", round(100*mean(d1$cog_low_q1, na.rm=TRUE),1), "%")
say("未测全者（n=", nrow(dat) - nrow(d1), "）若并入：主定义患病率 ",
    round(100*mean(dat$cog_low_edu, na.rm=TRUE),1), "%（敏感性分析用）")
say("权重可用(INDV_weight 非缺失): ", sum(!is.na(dat$INDV_weight)), " / ", nrow(dat))

close(log_con)
cat("\n完成。输出：output/01_2020单期基线/analysis_2020.rds\n")
