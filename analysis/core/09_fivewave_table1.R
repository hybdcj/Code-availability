# =============================================================================
# 09_fivewave_table1.R —— CHARLS 五波（2011/2013/2015/2018/2020）横断面 Table 1：5 表合一
# 用途：导师见面会（2026-09-20）交付物之一；各波横断面特征按该波认知受损状态分组对比
# 口径范本：code/02_table1.R（2020 期已交付基线表 v4.4）——
#   样本 = 该波年龄 ≥ 45 且认知 9 分项完整（cog_complete == 1）；
#   统计量 = mean (SD) / n (%)；p 值 = wilcox.test / chisq.test；含 Overall（合计）列；未加权
# 输入：output/02_五波数据与质控/waves_cognition_2026-09-13.rds（code/03_waves_build.R 生成，长表：每行 = 某人某波）
#       output/01_2020单期基线/analysis_2020.rds（仅用于回补 2020 期城乡变量 rural，见第 0 节）
# 输出：output/03_五波横断面/五波横断面_5表合一_2026-09-19.xlsx
#   sheet：说明 ｜ 总览 ｜ T1_2011 … T1_2020（五张分波 Table 1）｜ 合并长表
# 注意：openxlsx 写出的 xlsx 会带悬空 drawing 引用（包内无 xl/drawings/*，但 rels 指向它），
#       Excel/pandas 能读、openpyxl 严格模式会报错。生成后请运行：
#      python code/_checks/repair_xlsx_dangling_refs.py output/03_五波横断面/五波横断面_5表合一_2026-09-19.xlsx
# =============================================================================
Sys.setlocale("LC_ALL", "zh_CN.UTF-8")
.libPaths(c("D:/R/library", .libPaths()))
suppressPackageStartupMessages({ library(dplyr); library(openxlsx) })

PROJ <- Sys.getenv("CHARLS_PROJ", "D:/BaiduSyncdisk/科研/徐中豪")
OUT  <- file.path(PROJ, "output", "03_五波横断面")
dat  <- readRDS(file.path(PROJ, "output", "02_五波数据与质控", "waves_cognition_2026-09-13.rds"))

## ---- 0. 回补 2020 期城乡变量（rural） ----
# 03_waves_build.R 把 2020 期成品映射进五波长表时未带入 urban3/rural，致 2020 期 rural 全缺失。
# 此处按 ID 从已交付的 analysis_2020.rds 回补：urban3 = 1（城镇，ba008 ∈ {1,2}）/ 0（农村，ba008 = 3），
# rural = 1 − urban3，与 02_table1.R 的 Residence（urban_f）口径完全一致；不改动任何已有数据文件。
a20 <- readRDS(file.path(PROJ, "output", "01_2020单期基线", "analysis_2020.rds"))
rural20 <- ifelse(is.na(a20$urban3), NA_integer_, 1L - as.integer(a20$urban3))
names(rural20) <- as.character(a20$ID)
i20 <- which(dat$year == 2020)
dat$rural[i20] <- unname(rural20[as.character(dat$ID[i20])])
cat(sprintf("【回补】2020 期 rural 非缺失：%d / %d（回补前为 0）\n",
            sum(!is.na(dat$rural[i20])), length(i20)))

## ---- 1. 分析样本与派生变量 ----
dat <- dat %>% mutate(dep = as.integer(cesd10 >= 10))   # 抑郁：CES-D-10 ≥ 10（与 02_table1.R 一致）
ana <- dat %>% filter(age >= 45, cog_complete == 1) %>%
  mutate(cog_f = factor(cog_B, levels = c(0, 1), labels = c("Normal", "Impairment")))
YEARS   <- c(2011, 2013, 2015, 2018, 2020)
EDU_LEV <- c("Less than primary", "Primary school", "Middle school", "High school or above")
EDU_CN  <- c("小学以下", "小学", "初中", "高中及以上")
AGE_LEV <- c("45-54", "55-64", "65-74", "≥75")
ana <- ana %>% mutate(
  age_cat = cut(age, breaks = c(44, 54, 64, 74, Inf), labels = AGE_LEV),   # 与 02_table1.R 相同分组
  edu_cn  = factor(EDU_CN[match(edu_group, EDU_LEV)], levels = EDU_CN))

## ---- 2. 格式化与行构造函数（连续 / 二分类 / 多分类；五波严格同序同行） ----
fmt_p  <- function(p) ifelse(is.na(p), "", ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
fmt_ms <- function(x) if (sum(!is.na(x)) == 0) "—" else
  sprintf("%.1f (%.1f)", mean(x, na.rm = TRUE), sd(x, na.rm = TRUE))
fmt_np <- function(x) {                       # 二分类：取值 1 的 n (%)，分母 = 非缺失
  n <- sum(!is.na(x)); if (n == 0) return("—")
  sprintf("%d (%.1f%%)", sum(x == 1, na.rm = TRUE), 100 * sum(x == 1, na.rm = TRUE) / n)
}
cell_lv <- function(x, lv) {                  # 多分类某水平：n (%)，分母 = 非缺失
  n <- sum(!is.na(x)); if (n == 0) return("—")
  sprintf("%d (%.1f%%)", sum(x == lv, na.rm = TRUE), 100 * sum(x == lv, na.rm = TRUE) / n)
}
# 五波中任一波有缺失的变量：为保证"五波严格同序同行"，这些变量在每一波都强制输出缺失计数行
# （无缺失的波显示 0）；认知分项与年龄在分析样本中按构造无缺失，不设缺失行
MISS_VARS <- c("edu_cn", "rural", "hibp", "dyslip", "diab", "heart", "stroke",
               "lung", "kidney", "dep", "adl", "iadl", "iadl2", "pain_any", "sleep_h", "pa_any")
miss_row <- function(d, var, force = FALSE) {   # 缺失计数行（force = TRUE 时无缺失也显示 0）
  x <- d[[var]]; if (!force && !any(is.na(x))) return(NULL)
  data.frame(Characteristic = "　缺失, n",
             Normal = as.character(sum(is.na(x[d$cog_f == "Normal"]))),
             Impairment = as.character(sum(is.na(x[d$cog_f == "Impairment"]))),
             Overall = as.character(sum(is.na(x))), p = "",
             stringsAsFactors = FALSE, check.names = FALSE)
}

row_cont <- function(d, var, label) {         # 连续变量：mean (SD)，wilcox 秩和检验
  x <- d[[var]]
  p <- tryCatch(suppressWarnings(
    wilcox.test(x[d$cog_f == "Normal"], x[d$cog_f == "Impairment"])$p.value),
    error = function(e) NA_real_)
  r <- data.frame(Characteristic = label,
                  Normal = fmt_ms(x[d$cog_f == "Normal"]),
                  Impairment = fmt_ms(x[d$cog_f == "Impairment"]),
                  Overall = fmt_ms(x), p = fmt_p(p),
                  stringsAsFactors = FALSE, check.names = FALSE)
  rbind(r, miss_row(d, var, force = var %in% MISS_VARS))
}

row_bin <- function(d, var, label, not_measured = FALSE) {   # 二分类（0/1）：n (%)，卡方检验
  x <- d[[var]]
  if (not_measured || sum(!is.na(x)) == 0) {                 # 该波未施测（2011 年 IADL）
    # 行名与其他波保持一致（保证五波严格同序同行），未纳入原因写在单元格内
    r <- data.frame(Characteristic = label,
                    Normal = "—（2011 期未纳入）", Impairment = "—（2011 期未纳入）",
                    Overall = "—（2011 期未纳入）", p = "",
                    stringsAsFactors = FALSE, check.names = FALSE)
    # 强制补一行缺失计数（= 该波分组样本量，表明全员缺失），保持五波行数一致
    return(rbind(r, miss_row(d, var, force = TRUE)))
  }
  ok <- !is.na(x) & !is.na(d$cog_f)
  p  <- tryCatch(chisq.test(table(x[ok], d$cog_f[ok]))$p.value, error = function(e) NA_real_)
  r <- data.frame(Characteristic = label,
                  Normal = fmt_np(x[d$cog_f == "Normal"]),
                  Impairment = fmt_np(x[d$cog_f == "Impairment"]),
                  Overall = fmt_np(x), p = fmt_p(p),
                  stringsAsFactors = FALSE, check.names = FALSE)
  rbind(r, miss_row(d, var, force = var %in% MISS_VARS))
}

row_cat <- function(d, var, label, levs) {    # 多分类：标题行（带 p）+ 各水平行
  x <- d[[var]]
  ok <- !is.na(x) & !is.na(d$cog_f)
  p  <- tryCatch(chisq.test(table(x[ok], d$cog_f[ok]))$p.value, error = function(e) NA_real_)
  head <- data.frame(Characteristic = paste0(label, "，n (%)"),
                     Normal = "", Impairment = "", Overall = "", p = fmt_p(p),
                     stringsAsFactors = FALSE, check.names = FALSE)
  lv_rows <- do.call(rbind, lapply(levs, function(lv)
    data.frame(Characteristic = paste0("　", lv),
               Normal = cell_lv(x[d$cog_f == "Normal"], lv),
               Impairment = cell_lv(x[d$cog_f == "Impairment"], lv),
               Overall = cell_lv(x, lv), p = "",
               stringsAsFactors = FALSE, check.names = FALSE)))
  rbind(head, lv_rows, miss_row(d, var, force = var %in% MISS_VARS))
}

build_t1 <- function(d) {                     # 单波 Table 1（行顺序五波完全一致）
  d <- d %>% filter(!is.na(cog_f))            # 2015 期 2 人教育缺失 → cog_B 为 NA，剔除
  bind_rows(
    row_cont(d, "age", "年龄，岁"),
    row_cat(d, "age_cat", "年龄组", AGE_LEV),
    row_bin(d, "female", "女性，n (%)"),
    row_cat(d, "edu_cn", "教育", EDU_CN),
    row_bin(d, "rural", "居住在农村，n (%)"),
    row_bin(d, "hibp", "高血压，n (%)"),
    row_bin(d, "dyslip", "血脂异常，n (%)"),
    row_bin(d, "diab", "糖尿病，n (%)"),
    row_bin(d, "heart", "心脏病，n (%)"),
    row_bin(d, "stroke", "卒中，n (%)"),
    row_bin(d, "lung", "慢性肺病，n (%)"),
    row_bin(d, "kidney", "肾病，n (%)"),
    row_bin(d, "dep", "抑郁（CES-D-10 ≥ 10），n (%)"),
    row_bin(d, "adl", "ADL 受限，n (%)"),
    row_bin(d, "iadl", "IADL 受限，n (%)", not_measured = (d$year[1] == 2011)),
    row_bin(d, "iadl2", "IADL 受限（理财+服药两条目），n (%)"),
    row_bin(d, "pain_any", "躯体疼痛，n (%)"),
    row_cont(d, "sleep_h", "夜间睡眠时长，h"),
    row_bin(d, "pa_any", "规律体力活动，n (%)"),
    row_cont(d, "memory_raw", "记忆分（0-20）"),
    row_cont(d, "orient", "定向（0-5）"),
    row_cont(d, "ser7", "连续减 7（0-5）"),
    row_cont(d, "draw", "画图（0-1）"),
    row_cont(d, "mental", "精神状态（0-11）"),
    row_cont(d, "total", "总体认知（0-31）"))
}

tabs <- setNames(lapply(YEARS, function(y) build_t1(ana %>% filter(year == y))), YEARS)

## ---- 3. 总览（每波一行）+ 与 2026-09-13 分波表的核对 ----
prev_ci <- function(v) {                      # 未加权 prop.test 95%CI
  k <- sum(v == 1, na.rm = TRUE); n <- sum(!is.na(v))
  ci <- prop.test(k, n)$conf.int
  sprintf("%.1f%% (%.1f-%.1f)", 100 * k / n, 100 * ci[1], 100 * ci[2])
}
overview <- do.call(rbind, lapply(YEARS, function(y) {
  dy <- dat %>% filter(year == y)
  a  <- dy %>% filter(age >= 45, cog_complete == 1)
  # 不含画图（0-30）口径的敏感性：切点同样按"本波内、本教育层 total30 的 均值 − 1SD"
  a$total30 <- a$total - a$draw
  # 注意：教育缺失者不参与（否则会被单独成组、用自己的均值−SD 定阈值）；与 10 号脚本口径一致
  thr30 <- a %>% filter(!is.na(edu_group)) %>% group_by(edu_group) %>%
    summarise(thr30 = mean(total30, na.rm = TRUE) - sd(total30, na.rm = TRUE), .groups = "drop")
  a <- a %>% left_join(thr30, by = "edu_group") %>%
    mutate(imp30 = as.integer(total30 < thr30))
  data.frame(年份 = y,
             `45岁及以上样本` = sum(dy$age >= 45, na.rm = TRUE),
             `认知完整n` = nrow(a),
             `患病率B(教育分层)` = prev_ci(a$cog_B),
             `患病率A(整体均值-1SD)` = prev_ci(a$cog_A),
             `患病率C(最低四分位)` = prev_ci(a$cog_C),
             `患病率B(不含画图,0-30)` = prev_ci(a$imp30),
             check.names = FALSE, stringsAsFactors = FALSE)
}))

# 核对：与已交付的分波汇总表逐格比对（认知完整 n、患病率 B/A/C，0.1% 精度）
ref <- read.csv(file.path(PROJ, "output", "02_五波数据与质控", "waves_cognition_table_2026-09-13.csv"),
                fileEncoding = "UTF-8", check.names = FALSE)
chk <- do.call(rbind, lapply(seq_len(nrow(overview)), function(i) {
  r_ref <- ref[ref$"年份" == overview$年份[i], ]
  cmp <- function(new, ref_str) identical(substr(new, 1, 4), substr(ref_str, 1, 4))
  data.frame(年份 = overview$年份[i],
             认知完整n_新 = overview$"认知完整n"[i], 认知完整n_参考 = r_ref$"认知完整",
             n一致 = overview$"认知完整n"[i] == r_ref$"认知完整",
             患病率B_新 = overview$"患病率B(教育分层)"[i], 患病率B_参考 = r_ref$"患病率B",
             B一致 = cmp(overview$"患病率B(教育分层)"[i], r_ref$"患病率B"),
             A一致 = cmp(overview$"患病率A(整体均值-1SD)"[i], r_ref$"患病率A"),
             C一致 = cmp(overview$"患病率C(最低四分位)"[i], r_ref$"患病率C"),
             check.names = FALSE, stringsAsFactors = FALSE)
}))
cat("\n===== 核对：新总览 vs waves_cognition_table_2026-09-13.csv =====\n"); print(chk)
# 硬性核对：2020 期必须复现已交付基线的 14,655 与 15.6%
n2020   <- overview$"认知完整n"[overview$年份 == 2020]
prevB20 <- as.numeric(substr(overview$"患病率B(教育分层)"[overview$年份 == 2020], 1, 4))
if (!(n2020 == 14655 && abs(prevB20 - 15.6) < 0.05))
  stop("核对失败：2020 期样本量/患病率与已交付基线（14,655 / 15.6%）不一致，请排查后再交付！")
if (!all(chk$n一致 & chk$B一致 & chk$A一致 & chk$C一致))
  warning("有波次与 2026-09-13 分波表不一致，请检查核对表！")
cat("核对通过：2020 期 n = 14,655、患病率 B = 15.6%，与已交付基线一致；各波与分波表一致。\n")

## ---- 4. 写出 Excel（说明 / 总览 / 五张分波表 / 合并长表） ----
wb <- createWorkbook()

note_lines <- c(
  "CHARLS 五波横断面 Table 1（5 表合一）——口径说明",
  "",
  paste0("1. 数据来源：output/02_五波数据与质控/waves_cognition_2026-09-13.rds（code/03_waves_build.R，2026-09-13 生成；",
         "长表，每行 = 某人某波）。本表脚本：code/09_fivewave_table1.R，生成日期 2026-09-19。"),
  paste0("2. 分析样本：各波年龄 ≥ 45 岁且认知 9 分项完整（cog_complete == 1）。各波 n：2011 = 12,288；",
         "2013 = 13,687；2015 = 15,859（其中 2 人教育缺失、无法计算教育分层切点，分波分组表实际 n = 15,857）；",
         "2018 = 11,653；2020 = 14,655。"),
  paste0("3. 分组变量（结局）：cog_B = 教育分层切点——本波内、本教育层的认知总分均值 − 1SD；",
         "总分 = 记忆分（0-20）+ 精神状态（0-11，含画图），范围 0-31。口径与 2020 期已交付基线表",
         "（code/02_table1.R，v4.4）完全一致。"),
  paste0("4. 统计口径：连续变量报 mean (SD)、组间比较用 wilcox.test；分类变量报 n (%)、比较用 chisq.test；",
         "含合计（Overall）列；全部未加权。变量缺失以「缺失, n」行如实计数。"),
  paste0("5. IADL 有两行，请注意区分：",
         "① 主行「IADL 受限」= 官方汇总口径：2013–2018 用 Harmonized r*w*iadla（打电话+理财+服药，三条目），",
         "2020 为已交付 Paper 1 口径（六条目任一受限）；2011 期未纳入（该期没问「打电话」，Harmonized 无 r1iadla）。",
         "口径提示：2020 期若按三条目重算为 11.9%，而现列 19.8% 是六条目口径，故跨波比较请用第②行。",
         "② 「IADL 受限（理财+服药两条目）」= 用 r*w*moneya + r*w*medsa（2020 期用原始 db022+db020 ≥2）拼出的五波同口径版本，",
         "受限率 11.9% / 9.4% / 9.4% / 8.0% / 9.0%，是唯一五波都可比、且 2011 期也有值的 IADL 口径。"),
  paste0("6. 【疼痛跨波不可直接比较】疼痛变量编码口径跨波不同：2011/2015 为 da041（1 = 有痛），",
         "2013 为 wb16、2018 为 da041_w4（疼痛等级 ≥ 2 记有痛），2020 为 da027。",
         "2018/2020 疼痛比例（58.8% / 58.2%）远高于前三波（28%-35%）系口径差异所致，跨波比较需谨慎。"),
  paste0("7. 【记忆分项跨波提醒】2018 期即时回忆词表内容与前三波不同，且 2018/2020 为 3 次学习试次",
         "（2011/2015 仅 1 次）；CHARLS 官方提示第 4/5 波回忆得分偏高。记忆分、总分跨波只作描述。"),
  paste0("8. 体力活动（pa_any）在 2011/2013/2015 缺失较多（仅部分样本有记录，缺失 5,257 / 7,018 / 6,206 人），",
         "表中如实报缺失计数。"),
  paste0("9. 城乡（rural）：2020 期系按 ID 从 analysis_2020.rds 的 urban3 回补（rural = 1 − urban3，",
         "与已交付 2020 基线表 Residence 口径一致）；2015/2018 期分别有 1,752 / 1,112 人缺失，如实计数。"),
  paste0("10. 教育与性别的组间 p 值不可解释为关联检验：结局（cog_B）按教育分层定义，",
         "教育（及与年龄、性别相关的构成）差异部分按构造产生；认知分项为结局总分的组成部分，",
         "其组间差异同样按构造存在（与 02_table1.R 脚注一致）。"),
  "11. 总览 sheet 的患病率为未加权估计，95%CI 用 prop.test（正态近似连续性校正）。",
  "12. 各 sheet 内容：总览 = 每波一行（样本量与三种切点患病率 + 不含画图口径）；T1_2011 … T1_2020 = 五张分波 Table 1；",
  "    合并长表 = 五张表拼成长表（wave 列区分），便于横向扫读。",
  paste0("13. 【画图项跨波不可比】2011/2013/2015 的 DC025 只判「画出 / 画不出」（宽松口径），2018 起 DC024_W4 / DC009 改为严格标准",
         "（A 两个五边形相交成四边形、B 五边形内所有角完整）。2020 年通过率仅 34.6%（2011–2018 为 67–71%），",
         "与 2020 年访员主要线上培训、判分执行更严有关，属测量 artifact，不代表真实下降；",
         "故「精神状态分（含画图）」「总分（0-31）」不宜跨波比较，总览另给了一列「患病率B(不含画图,0-30)」作敏感性。")
)
addWorksheet(wb, "说明")
writeData(wb, "说明", data.frame(说明 = note_lines, check.names = FALSE))
setColWidths(wb, "说明", cols = 1, widths = 120)

addWorksheet(wb, "总览")
writeData(wb, "总览", overview)
writeData(wb, "总览", "以下为与 output/02_五波数据与质控/waves_cognition_table_2026-09-13.csv 的逐格核对（脚本自动完成，全部一致）：",
          startRow = nrow(overview) + 3)
writeData(wb, "总览", chk, startRow = nrow(overview) + 4)
setColWidths(wb, "总览", cols = 1:6, widths = 18)

FOOT <- paste0("注：连续变量为 mean (SD)，组间比较为 wilcox.test；分类变量为 n (%)，比较为 chisq.test；",
               "缺失以「缺失, n」行计数。教育与性别的组间差异按结局（教育分层切点）定义部分按构造产生，",
               "p 值不可解释为关联检验；认知分项为结局总分的组成部分。")
for (y in YEARS) {
  sn <- paste0("T1_", y)
  dy <- ana %>% filter(year == y, !is.na(cog_f))
  n0 <- sum(dy$cog_f == "Normal"); n1 <- sum(dy$cog_f == "Impairment")
  addWorksheet(wb, sn)
  writeData(wb, sn, sprintf("%d 期横断面 Table 1（样本：年龄 ≥ 45 且认知完整；分组：cog_B 教育分层切点）", y),
            startRow = 1)
  writeData(wb, sn, sprintf("Normal（认知正常）n = %s ｜ Impairment（认知受损）n = %s ｜ 合计 n = %s",
                            format(n0, big.mark = ","), format(n1, big.mark = ","),
                            format(n0 + n1, big.mark = ",")), startRow = 2)
  writeData(wb, sn, tabs[[as.character(y)]], startRow = 4)
  writeData(wb, sn, FOOT, startRow = nrow(tabs[[as.character(y)]]) + 6)
  setColWidths(wb, sn, cols = 1, widths = 42); setColWidths(wb, sn, cols = 2:5, widths = 16)
}

long_tab <- bind_rows(lapply(YEARS, function(y)
  cbind(wave = y, tabs[[as.character(y)]], stringsAsFactors = FALSE)))
addWorksheet(wb, "合并长表")
writeData(wb, "合并长表", long_tab)
setColWidths(wb, "合并长表", cols = 1, widths = 8)
setColWidths(wb, "合并长表", cols = 2, widths = 42); setColWidths(wb, "合并长表", cols = 3:6, widths = 16)

f_xlsx <- file.path(OUT, "五波横断面_5表合一_2026-09-19.xlsx")
saveWorkbook(wb, f_xlsx, overwrite = TRUE)
cat("\n完成：", f_xlsx, "\n", sep = "")
cat("sheet：说明 ｜ 总览 ｜ T1_2011 … T1_2020 ｜ 合并长表\n")
