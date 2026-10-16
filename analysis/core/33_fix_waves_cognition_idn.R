# 修复 waves_cognition_2026-09-13.rds 中 2020 波 IDn 全缺的问题
# ---------------------------------------------------------------------------
# 问题：03_waves_build.R 里 2020 波是单独拼的 w2020 数据框，列清单漏了 IDn，
#       bind_rows 时把 2020 行的 IDn 全部填成 NA（17,730 行），
#       导致该波无法与其他波做个体纵向链接（"五波齐全者"退化为 0 人）。
# 修复：2020 的 ID 列本身已是 12 位规范 ID，直接用 norm_id(ID) 补 IDn。
#       本脚本幂等——已修复的数据再跑一次不会有任何变化。
# 说明：R 在 Windows 下读写中文路径不稳，脚本只处理英文临时路径，
#       由外部（PowerShell）负责备份与回填。
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({ library(dplyr) })

TMP     <- "D:/codex_tmp/ckm_stage"
IN_RDS  <- file.path(TMP, "waves_cognition_orig.rds")
OUT_RDS <- file.path(TMP, "waves_cognition_fixed.rds")
OUT_LOG <- file.path(TMP, "out_log_idn_fix.txt")

norm_id <- function(x) {
  x <- as.character(x); n <- nchar(x)
  ifelse(n == 11,
         paste0(substr(x, 1, 9), sprintf("%03d", as.integer(substr(x, 10, 11)))),
         x)
}

L <- c()
p <- function(...) L <<- c(L, paste0(...))

dat <- readRDS(IN_RDS)
p("输入：", IN_RDS, "（", nrow(dat), " 行 × ", ncol(dat), " 列）")

## ---------- 1) 修复前诊断 ----------
p("\n== 修复前：各波 IDn 缺失情况 ==")
before <- dat %>% group_by(year) %>%
  summarise(rows = n(), idn_na = sum(is.na(IDn)),
            id_ok = sum(!is.na(ID)), .groups = "drop")
p(capture.output(print(as.data.frame(before))))

nw_before <- dat %>% group_by(IDn) %>% summarise(nw = n_distinct(year), .groups = "drop")
p("修复前个体数（按 IDn）：", nrow(nw_before))
p("修复前五波齐全者：", sum(nw_before$nw == 5))
p("修复前最长随访波数：", max(nw_before$nw))

## ---------- 2) 修复 ----------
stopifnot("ID" %in% names(dat))
need <- is.na(dat$IDn) & !is.na(dat$ID)
dat$IDn[need] <- norm_id(dat$ID[need])
p("\n== 修复动作 ==")
p("按 ID 补齐 IDn 的行数：", sum(need))

## ---------- 3) 修复后校验 ----------
p("\n== 修复后校验 ==")
p("IDn 仍缺失：", sum(is.na(dat$IDn)))
p("IDn 长度分布：", paste(names(table(nchar(dat$IDn))), table(nchar(dat$IDn)),
                          sep = ":", collapse = " "))
after <- dat %>% group_by(year) %>%
  summarise(rows = n(), idn_na = sum(is.na(IDn)), .groups = "drop")
p(capture.output(print(as.data.frame(after))))

nw_after <- dat %>% group_by(IDn) %>% summarise(nw = n_distinct(year), .groups = "drop")
p("修复后个体数（按 IDn）：", nrow(nw_after))
p("修复后五波齐全者：", sum(nw_after$nw == 5))
p("修复后最长随访波数：", max(nw_after$nw))
p("各随访波数人数：", paste(names(table(nw_after$nw)), table(nw_after$nw),
                          sep = ":", collapse = " "))

## 关键一致性检查：2020 与 2011 的 IDn 必须落在同一套编码上
ids11 <- dat$IDn[dat$year == 2011]; ids20 <- dat$IDn[dat$year == 2020]
p("2011 与 2020 交集人数：", length(intersect(ids11, ids20)))
p("2020 中不在 2011 的人数：", length(setdiff(ids20, ids11)))

## 2020 内部一致性：IDn 应等于 ID（12 位原样）
p("2020 行中 IDn == ID 的比例：",
  sprintf("%.4f", mean(dat$IDn[dat$year == 2020] == dat$ID[dat$year == 2020])))

saveRDS(dat, OUT_RDS)
p("\n输出：", OUT_RDS)
writeLines(L, OUT_LOG)
cat(paste(L, collapse = "\n"), "\n")
