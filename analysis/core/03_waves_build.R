# ============================================================================
# 03_waves_build.R —— 五波认知数据集构建（路线 B：逐波从原始模块按 2020 口径重建）
# 版本：v1.0（2026-09-13）
# 输入：data/raw/wave1_2011、wave2_2013、wave3_2015、wave4_2018 的原始模块
#       （2020 期直接复用已建好的 output/01_2020单期基线/analysis_2020.rds）
# 输出：output/02_五波数据与质控/waves_cognition_2026-09-13.rds（五波合并的纵向认知数据）
#       output/02_五波数据与质控/waves_cognition_table_2026-09-13.csv（分波表）
# 口径要点：
#   ① 9 分项：定向 5 项、连续减 7（5 步）、画图、即时回忆、延迟回忆
#   ② 记忆分三套：原始分（raw）、比例分（比例化，跨波可比）、波内 z 分
#   ③ 患病率：每波用"本波自己的教育分层切点"（主定义 B）与整体/四分位切点（A/C）
#   ④ 2018 期口径（2026-09-13 修正）：主口径用 DC028_W4（10 词表·即时试次1）+ DC047_W4（延迟）
#      与 DC024_W4（复制重叠五边形）；这三项按问卷 PROCEDURE 是 45 岁以上全样本施测。
#      60 岁以上的 MMSE 子集另有 3 词表（DC013_W4_1 / DC015_W4）与地点定向等，不作纵向主口径。
#      另：第 4 轮词表内容与施测方式（3 次学习试次）和 2011–2015 不同，官方建议纵向用等百分位等分。
# ============================================================================
invisible(try(Sys.setlocale("LC_ALL", "zh_CN.UTF-8"), silent = TRUE))
if (dir.exists("D:/R/library")) .libPaths(c("D:/R/library", .libPaths()))
suppressPackageStartupMessages({ library(haven); library(dplyr) })

PROJ <- Sys.getenv("CHARLS_PROJ", "D:/BaiduSyncdisk/科研/徐中豪")
RAW  <- file.path(PROJ, "data/raw")
OUT  <- file.path(PROJ, "output", "02_五波数据与质控")
if (!dir.exists(OUT)) dir.create(OUT, recursive = TRUE)

## ---------- 通用小工具 ----------
num <- function(x) as.numeric(x)
# 词表回忆：前 ns 个词位中取值 1..ns 视为"回忆出该词"；标记位另行处理
#   none_var   : "一个词也没回忆出" 标记（11 / 4），用于判定"确实问过"
#   refuse_var : "拒绝回忆" 标记（12）→ 记为缺失，与 2020 期 01_build_2020.R 口径一致
recall_cnt <- function(d, vars, ns, none_var = NULL, refuse_var = NULL) {
  m <- sapply(vars, function(v) num(d[[v]]))
  m <- matrix(m, ncol = length(vars))
  cnt <- rowSums(m >= 1 & m <= ns, na.rm = TRUE)
  any_asked <- rowSums(!is.na(m)) > 0
  if (!is.null(none_var) && none_var %in% names(d)) {
    none <- num(d[[none_var]]); any_asked <- any_asked | (!is.na(none) & none != 0)
  }
  if (!is.null(refuse_var) && refuse_var %in% names(d)) {
    rv <- num(d[[refuse_var]]); cnt[!is.na(rv) & rv == 12] <- NA
  }
  cnt[!any_asked] <- NA
  as.numeric(cnt)
}
# 连续减 7：答案按"前一步 − 7"判对；第一步以 93（100−7）为正确
ser7_score <- function(d, ans_vars) {
  a <- sapply(ans_vars, function(v) { x <- num(d[[v]]); x[!is.na(x) & (x < 0 | x > 100)] <- NA; x })
  a <- matrix(a, ncol = length(ans_vars))
  n <- ncol(a); ok <- matrix(FALSE, nrow = nrow(a), ncol = n)
  ok[, 1] <- !is.na(a[, 1]) & a[, 1] == 93
  for (i in 2:n) ok[, i] <- !is.na(a[, i-1]) & !is.na(a[, i]) & a[, i] == a[, i-1] - 7
  ans <- rowSums(!is.na(a)) > 0
  s <- rowSums(ok, na.rm = TRUE); s[!ans] <- NA
  as.numeric(s)
}
# 教育四组（与 2020 口径一致：<小学 / 小学 / 初中 / 高中及以上）
edu4_raw <- function(x) {
  x <- num(x)
  dplyr::case_when(x %in% 1:3 ~ "Less than primary", x == 4 ~ "Primary school",
                   x == 5 ~ "Middle school", x >= 6 ~ "High school or above", TRUE ~ NA_character_)
}

## ---------- 各波配置 ----------
## 人口学统一取自 Harmonized D（命名统一、编码已核对）；认知一律取原始模块
cfg_year <- c(2011, 2013, 2015, 2018)     # 与 cfg 顺序一致，用于定位 Harmonized 的 r1–r4 字段
norm_id <- function(x) {
  x <- as.character(x); n <- nchar(x)
  ifelse(n == 11, paste0(substr(x, 1, 9), sprintf("%03d", as.integer(substr(x, 10, 11)))), x)
}
hd <- read_dta(file.path(RAW, "harmonized_d/H_CHARLS_D_Data.dta"),
               col_select = any_of(c("ID", "communityID",
                                     paste0("r", 1:4, "agey"), "ragender", "raeduc_c", "raeducl",
                                     paste0("r", 1:4, "hukou"), paste0("r", 1:4, "rural2"),
                                     paste0("r", 1:4, "adlab_c"), paste0("r", 2:4, "iadla"),
                                     paste0("r", 1:4, "cesd10"),
                                     paste0("r", 1:4, "hibpe"), paste0("r", 1:4, "diabe"),
                                     paste0("r", 1:4, "hearte"), paste0("r", 1:4, "stroke"),
                                     paste0("r", 1:4, "lunge"), paste0("r", 1:4, "kidneye"),
                                     paste0("r", 1:4, "dyslipe"),
                                     paste0("r", 1:4, "moneya"), paste0("r", 1:4, "medsa"),
                                     paste0("r", 1:4, "wtresp"))))
hd$IDn <- norm_id(hd$ID)

cfg <- list(
  list(year = 2011, dir = "wave1_2011", hsf = "health_status_and_functioning.dta",
       demo = "demographic_background.dta", wt = "weight.dta", wt_var = "ind_weight",
       orient = c("dc001s1","dc001s2","dc001s3","dc002","dc003"),
       orient_ok = c(1, 2, 3, 1, 1),
       ser7 = paste0("dc0", 19:23), draw = "dc025",
       pain = "da041", pain_ge2 = FALSE, pa_prefix = "da051_",
       imm = paste0("dc006s", 1:10), imm_none = "dc006s11",
       dlr = paste0("dc027s", 1:10), dlr_none = "dc027s11", draw_ok = 1),
  list(year = 2013, dir = "wave2_2013", hsf = "Health_Status_and_Functioning.dta",
       pain = "wb16", pain_ge2 = TRUE, pa_prefix = "da051_",
       demo = "Demographic_Background.dta", wt = "Weights.dta", wt_var = "INDV_weight",
       orient = c("dc001s1","dc001s2","dc001s3","dc002","dc003"),
       orient_ok = c(1, 2, 3, 1, 1),
       ser7 = paste0("dc0", 19:23), draw = "dc025",
       imm = paste0("dc006_1_s", 1:10), imm_none = "dc006_1_s11",
       dlr = paste0("dc027s", 1:10), dlr_none = "dc027s11", draw_ok = 1),
   list(year = 2015, dir = "wave3_2015", hsf = "Health_Status_and_Functioning.dta",
        pain = "da041", pain_ge2 = FALSE, pa_prefix = "da051_",
        demo = "Demographic_Background.dta", wt = "Weights.dta", wt_var = "INDV_weight",
        orient = c("dc001s1","dc001s2","dc001s3","dc002","dc003"),
        orient_ok = c(1, 2, 3, 1, 1),
        ser7 = paste0("dc0", 19:23), draw = "dc025",
        imm = paste0("dc006s", 1:10), imm_none = "dc006s11", imm_refuse = "dc006s12",
        dlr = paste0("dc027s", 1:10), dlr_none = "dc027s11", dlr_refuse = "dc027s12", draw_ok = 1),
  list(year = 2018, dir = "wave4_2018", hsf = "Cognition.dta",
       extra = "Health_Status_and_Functioning.dta",     # 2018 的睡眠/慢病在健康模块里
       demo = "Demographic_Background.dta", wt = "Weights.dta", wt_var = "INDV_weight",
       orient = c("dc001_w4","dc006_w4","dc005_w4","dc003_w4","dc002_w4"),
       orient_ok = c(1, 1, 1, 1, 1),                       # 1=对 / 5=错 / 97=不知道
       ser7 = paste0("dc014_w4_", 1:5, "_1"),
       draw = "dc024_w4", draw_err = 5,                     # 复制重叠五边形：1=对 / 5=错 / 97=未评估
       pain = "da041_w4", pain_ge2 = TRUE, pa_prefix = "da051_",
       imm = paste0("dc028_w4_s", 1:10), imm_none = "dc028_w4_s11", imm_refuse = "dc028_w4_s12",
       dlr = paste0("dc047_w4_s", 1:10), dlr_none = "dc047_w4_s11", dlr_refuse = "dc047_w4_s12",
       draw_ok = 1)
)

build_wave <- function(cf) {
  hs  <- read_dta(file.path(RAW, cf$dir, cf$hsf))
  dm  <- read_dta(file.path(RAW, cf$dir, cf$demo))
  wt  <- read_dta(file.path(RAW, cf$dir, cf$wt))
  cat("\n== ", cf$year, " ==\n  HSF 行数 ", nrow(hs), " ｜ 人口学 ", nrow(dm), " ｜ 权重 ", nrow(wt), "\n", sep = "")

  # 定向：取值等于"该题的正确码"记 1，其余记 0（NA 视为未答对 = 0；与 2020 把"不知道"记 0 一致）
  orient_m <- sapply(seq_along(cf$orient), function(i) {
    x <- num(hs[[cf$orient[i]]]); as.integer(!is.na(x) & x == cf$orient_ok[i])
  })
  orient_m <- matrix(orient_m, ncol = length(cf$orient))
  orient <- rowSums(orient_m)

  ser7 <- ser7_score(hs, cf$ser7)

  draw <- if (is.na(cf$draw)) rep(NA_integer_, nrow(hs)) else {
    x <- num(hs[[cf$draw]])
    err_code <- if (!is.null(cf$draw_err)) cf$draw_err else 2      # 各波"画错"的码不同：2011/2013/2015=2，2018=5
    x[!is.na(x) & !(x %in% c(1, err_code))] <- NA                 # 97/3/997/999 等"未评估/不知道"记缺失
    ifelse(is.na(x), NA_integer_, as.integer(x == cf$draw_ok))
  }
  imm <- recall_cnt(hs, cf$imm, length(cf$imm), cf$imm_none, cf$imm_refuse)
  dlr <- recall_cnt(hs, cf$dlr, length(cf$dlr), cf$dlr_none, cf$dlr_refuse)

  df <- data.frame(ID = as.character(hs$ID), year = cf$year,
                   orient = orient, ser7 = ser7, draw = draw, imrc = imm, dlrc = dlr)
  df$list_len_imm <- length(cf$imm); df$list_len_dlr <- length(cf$dlr)
  df$memory_raw <- df$imrc + df$dlrc
  df$memory_prop <- (df$imrc + df$dlrc) / (df$list_len_imm + df$list_len_dlr)   # A：比例化
  # 精神状态：有画图项时 3 项合成（0–11），无画图项时 2 项（0–10）
  df$mental <- if (all(is.na(draw))) rowSums(cbind(df$orient, df$ser7), na.rm = FALSE) else
    rowSums(cbind(df$orient, df$ser7, df$draw), na.rm = FALSE)
  df$draw_available <- !all(is.na(draw))
  df$cog_complete <- as.integer(!is.na(df$orient) & !is.na(df$ser7) &
                                  (!df$draw_available | !is.na(df$draw)) &
                                  !is.na(df$imrc) & !is.na(df$dlrc))

  # 人口学：自 Harmonized D 对应波次字段取（w = 1..4）
  w <- which(sapply(cfg_year, function(y) y == cf$year))[1]
  dm2 <- data.frame(
    IDn = hd$IDn,
    age = num(hd[[paste0("r", w, "agey")]]),
    female = as.integer(num(hd$ragender) == 2),
    edu = num(hd$raeduc_c),
    rural = num(hd[[paste0("r", w, "rural2")]]),
    communityID = hd$communityID,
    weight = num(hd[[paste0("r", w, "wtresp")]]),
    adl    = as.integer(num(hd[[paste0("r", w, "adlab_c")]]) >= 1),
    iadl   = if (paste0("r", w, "iadla") %in% names(hd)) as.integer(num(hd[[paste0("r", w, "iadla")]]) >= 1) else NA_integer_,
    # IADL 两条目版（理财 + 服药）：官方汇总 r*iadla 含"打电话"，而 2011 期没问打电话，
    # 故 r1iadla 不存在；这两条目 2011 期也有（r1moneya / r1medsa），可拼出五波可比的两条目口径
    iadl2  = {
      mo <- num(hd[[paste0("r", w, "moneya")]]); me <- num(hd[[paste0("r", w, "medsa")]])
      ifelse(is.na(mo) & is.na(me), NA_integer_,
             as.integer((!is.na(mo) & mo == 1) | (!is.na(me) & me == 1)))
    },
    cesd10 = num(hd[[paste0("r", w, "cesd10")]]),
    hibp   = num(hd[[paste0("r", w, "hibpe")]]),
    diab   = num(hd[[paste0("r", w, "diabe")]]),
    heart  = num(hd[[paste0("r", w, "hearte")]]),
    stroke = num(hd[[paste0("r", w, "stroke")]]),
    lung   = num(hd[[paste0("r", w, "lunge")]]),
    kidney = num(hd[[paste0("r", w, "kidneye")]]),
    dyslip = num(hd[[paste0("r", w, "dyslipe")]]))
  dm2$edu_group <- edu4_raw(dm2$edu)
  df$IDn <- norm_id(df$ID)
  slp_var <- if ("da049" %in% names(hs)) "da049" else if ("da030" %in% names(hs)) "da030" else NA
  if (is.na(slp_var) && !is.null(cf$extra)) {         # 2018：睡眠在另一个模块
    ex <- read_dta(file.path(RAW, cf$dir, cf$extra))
    if ("da049" %in% names(ex)) {
      x <- num(ex$da049); x[!is.na(x) & (x < 0 | x > 24)] <- NA
      df$sleep_h <- x[match(df$ID, as.character(ex$ID))]
      slp_var <- "da049"
    }
  }
  if (!"sleep_h" %in% names(df)) df$sleep_h <- if (is.na(slp_var)) NA_real_ else {
    x <- num(hs[[slp_var]]); x[!is.na(x) & (x < 0 | x > 24)] <- NA; x }
  # 疼痛（各波编码不同：da041 为 1=有痛/2=无痛；da041_w4、wb16、da027 为 1=无痛、2-5=疼痛等级）
  px <- hs                                   # 疼痛/体力活动所在模块（2018 在 extra 模块里）
  if (!is.null(cf$extra)) px <- read_dta(file.path(RAW, cf$dir, cf$extra))
  idx <- match(df$ID, as.character(px$ID))
  if (!is.null(cf$pain) && cf$pain %in% names(px)) {
    x <- num(px[[cf$pain]]); x[!is.na(x) & x %in% c(997, 999, -1)] <- NA
    v <- if (isTRUE(cf$pain_ge2)) as.integer(x >= 2) else as.integer(x == 1)
    df$pain_any <- v[idx]
  } else df$pain_any <- NA_integer_
  pa_vars <- intersect(paste0(cf$pa_prefix, 1:3, "_"), names(px))
  if (length(pa_vars)) {
    m <- sapply(pa_vars, function(v) num(px[[v]])); m <- matrix(m, ncol = length(pa_vars))
    v <- as.integer(rowSums(m == 1, na.rm = TRUE) > 0); v[rowSums(!is.na(m)) == 0] <- NA
    df$pa_any <- v[idx]
  } else df$pa_any <- NA_integer_
  df <- left_join(df, dm2, by = "IDn")

  # 切点（每波自己算）：教育分层 B、整体 −1SD 的 A、最低四分位 C
  prim <- df %>% filter(cog_complete == 1, !is.na(memory_raw), !is.na(mental))
  gm <- mean(prim$memory_raw + prim$mental, na.rm = TRUE); gs <- sd(prim$memory_raw + prim$mental, na.rm = TRUE)
  thrA <- gm - gs; thrC <- quantile(prim$memory_raw + prim$mental, 0.25, na.rm = TRUE)
  edu_thr <- prim %>% filter(!is.na(edu_group)) %>% group_by(edu_group) %>%
    summarise(thr = mean(memory_raw + mental, na.rm = TRUE) - sd(memory_raw + mental, na.rm = TRUE), .groups = "drop")
  df <- df %>% left_join(edu_thr, by = "edu_group") %>%
    mutate(total = memory_raw + mental,
           cog_B = as.integer(total < thr), cog_A = as.integer(total < thrA),
           cog_C = as.integer(total <= thrC))
  # 波内 z（B 方案：用于跨波比较与轨迹模型）
  df$total_z <- (df$total - gm) / gs
  df$memory_z <- (df$memory_raw - mean(prim$memory_raw, na.rm = TRUE)) / sd(prim$memory_raw, na.rm = TRUE)
  df$wave <- as.character(cf$year)
  df
}

waves <- lapply(cfg, build_wave)

## 2020 期：复用已建好的数据集，映射成同一套字段
a20 <- readRDS(file.path(PROJ, "output", "01_2020单期基线", "analysis_2020.rds"))
w2020 <- data.frame(
  ID = as.character(a20$ID), year = 2020, wave = "2020",
  # 【2026-10-11 修复】此前漏写 IDn，bind_rows 会把 2020 行的 IDn 填成 NA，
  # 导致该波无法与其他波做个体纵向链接（详见 code/33_fix_waves_cognition_idn.R）。
  # 2020 的 ID 本身已是 12 位规范 ID，norm_id() 会原样返回。
  IDn = norm_id(as.character(a20$ID)),
  orient = a20$orient, ser7 = a20$ser7, draw = a20$draw, imrc = a20$imrc, dlrc = a20$dlrc,
  list_len_imm = 10, list_len_dlr = 10,
  memory_raw = a20$memory, memory_prop = a20$memory / 20,
  mental = a20$mental, draw_available = TRUE, cog_complete = a20$cog_complete,
  age = a20$age, female = a20$female, edu_group = a20$edu_group,
  communityID = a20$communityID, weight = a20$INDV_weight,
  total = a20$global, cog_B = a20$cog_low_edu, cog_A = a20$cog_low_1sd,
  cog_C = a20$cog_low_q1,
  adl = a20$adl_any, iadl = a20$iadl_any, cesd10 = a20$cesd10, sleep_h = a20$sleep_h,
  hibp = a20$hypertension, diab = a20$diabetes, heart = a20$heart_disease,
  stroke = a20$stroke, lung = a20$lung_disease, kidney = a20$kidney_disease,
  dyslip = a20$dyslipidemia, pain_any = a20$pain_any, pa_any = a20$pa_any)
# 2020 期 IADL 两条目（理财 db022 + 服药 db020）：与 2011-2018 的 Harmonized 单条目口径一致（>=2 记为有困难）
hs20 <- read_dta(file.path(RAW, "wave5_2020/Health_Status_and_Functioning.dta"),
                 col_select = any_of(c("ID", "db020", "db022")))
mo20 <- num(hs20$db022); me20 <- num(hs20$db020)
mo20 <- ifelse(!is.na(mo20), as.integer(mo20 >= 2), NA_integer_)
me20 <- ifelse(!is.na(me20), as.integer(me20 >= 2), NA_integer_)
iadl2_20 <- ifelse(is.na(mo20) & is.na(me20), NA_integer_, as.integer(mo20 == 1 | me20 == 1))
w2020$iadl2 <- iadl2_20[match(as.character(a20$ID), as.character(hs20$ID))]
prim20 <- w2020 %>% filter(cog_complete == 1, !is.na(total))
w2020$total_z <- (w2020$total - mean(prim20$total, na.rm = TRUE)) / sd(prim20$total, na.rm = TRUE)
w2020$memory_z <- (w2020$memory_raw - mean(prim20$memory_raw, na.rm = TRUE)) / sd(prim20$memory_raw, na.rm = TRUE)
waves <- c(waves, list(w2020))

allw <- bind_rows(lapply(waves, function(d) d %>% mutate(across(where(is.numeric), ~ .))))
saveRDS(allw, file.path(OUT, "waves_cognition_2026-09-13.rds"))

## ---------- 分波表 ----------
tab <- do.call(rbind, lapply(waves, function(d) {
  inw <- d[!is.na(d$age) & d$age >= 45, ]
  comp <- inw[inw$cog_complete == 1, ]
  data.frame(
    年份 = d$year[1],
    可用样本 = nrow(inw),
    认知完整 = nrow(comp),
    记忆完整 = sum(!is.na(comp$memory_raw)),
    平均年龄 = round(mean(inw$age, na.rm = TRUE), 1),
    女性 = sprintf("%.1f%%", 100*mean(inw$female, na.rm = TRUE)),
    记忆分原始均值 = round(mean(comp$memory_raw, na.rm = TRUE), 2),
    记忆分比例均值 = round(mean(comp$memory_prop, na.rm = TRUE), 3),
    精神状态均值 = round(mean(comp$mental, na.rm = TRUE), 2),
    总体认知均值 = round(mean(comp$total, na.rm = TRUE), 2),
    患病率B = sprintf("%.1f%%", 100*mean(comp$cog_B, na.rm = TRUE)),
    患病率A = sprintf("%.1f%%", 100*mean(comp$cog_A, na.rm = TRUE)),
    患病率C = sprintf("%.1f%%", 100*mean(comp$cog_C, na.rm = TRUE)),
    ADL受限 = sprintf("%.1f%%", 100*mean(comp$adl, na.rm = TRUE)),
    IADL受限 = if (all(is.na(comp$iadl))) "-" else sprintf("%.1f%%", 100*mean(comp$iadl, na.rm = TRUE)),
    抑郁 = sprintf("%.1f%%", 100*mean(comp$cesd10 >= 10, na.rm = TRUE)),
    高血压 = sprintf("%.1f%%", 100*mean(comp$hibp == 1, na.rm = TRUE)),
    糖尿病 = sprintf("%.1f%%", 100*mean(comp$diab == 1, na.rm = TRUE)),
    睡眠小时 = round(mean(comp$sleep_h, na.rm = TRUE), 2),
    疼痛 = sprintf("%.1f%%", 100*mean(comp$pain_any, na.rm = TRUE)),
    体力活动 = sprintf("%.1f%%", 100*mean(comp$pa_any, na.rm = TRUE)),
    含画图项 = unique(comp$draw_available)[1])
}))
write.csv(tab, file.path(OUT, "waves_cognition_table_2026-09-13.csv"), row.names = FALSE, fileEncoding = "UTF-8")
cat("\n================ 分波表（认知） ================\n"); print(tab, row.names = FALSE)
cat("\n已存：output/02_五波数据与质控/waves_cognition_2026-09-13.rds、output/02_五波数据与质控/waves_cognition_table_2026-09-13.csv\n")
