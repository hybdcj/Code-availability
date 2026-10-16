# ============================================================================
# 17_hrs_merge_inventory.R —— HRS 两文件按 ID 合并 + 各波可分析样本盘点
# 目的：决定外部验证用哪些波次（SAP §9.0 方案 A/B）并为 A11 提供样本量依据
# 生成：2026-09-26
# 输出：external_data/HRS/logs/hrs_wave_inventory_2026-09-26.csv / .md
# ============================================================================
invisible(try(Sys.setlocale("LC_ALL", "zh_CN.UTF-8"), silent = TRUE))
.libPaths(c("D:/R/library", .libPaths()))
suppressPackageStartupMessages({ library(haven); library(dplyr) })
HRS <- "D:/BaiduSyncdisk/科研/徐中豪/external_data/HRS"
OUT <- file.path(HRS, "logs"); if (!dir.exists(OUT)) dir.create(OUT, recursive = TRUE)
RD <- file.path(HRS, "harmonized/randhrs1992_2022v1_STATA.zip")
GW <- file.path(HRS, "harmonized/H_HRS_d_stata.zip")
rd_tmp <- "D:/codex_tmp/hrs/randhrs1992_2022v1.dta"
gw_tmp <- "D:/codex_tmp/hrs/H_HRS_d.dta"

waves <- 2:13
rd_cols <- c("hhidpn", "raedyrs", paste0("r", waves, "agey_b"), paste0("r", waves, "cesd"),
             paste0("r", waves, "adl5a"), paste0("r", waves, "ser7"))
gw_cols <- c("hhidpn", "raeducl", paste0("r", waves, "orient"), paste0("r", waves, "adlfive"))
rd <- read_dta(rd_tmp, col_select = any_of(rd_cols))
gw <- read_dta(gw_tmp, col_select = any_of(gw_cols))
rd <- as.data.frame(lapply(rd, function(x) { x <- as.numeric(x); x[x < 0] <- NA; x }))
gw <- as.data.frame(lapply(gw, function(x) { x <- as.numeric(x); x[x < 0] <- NA; x }))
cat("RAND 行数:", nrow(rd), " Gateway 行数:", nrow(gw), "\n")
m <- inner_join(rd, gw, by = "hhidpn", suffix = c(".rd", ".gw"))
cat("按 hhidpn 合并后:", nrow(m), " 行（RAND 匹配率 ", sprintf("%.1f%%", 100*nrow(m)/nrow(rd)), "）\n\n")

inv <- do.call(rbind, lapply(waves, function(w) {
  g <- function(p) paste0("r", w, p)
  o <- m[[g("orient")]]; s <- m[[g("ser7")]]
  # 定义：4 个预测变量（年龄 agey_b、教育 raedyrs、抑郁 cesd、ADL adl5a）+ 结局两成分（orient + ser7）
  ok <- complete.cases(m[, c(g("agey_b"), g("cesd"), g("adl5a"), g("ser7"), g("orient"))]) & !is.na(m$raedyrs)
  data.frame(wave = w,
    pct_orient = round(100*mean(!is.na(o)), 1), pct_ser7 = round(100*mean(!is.na(s)), 1),
    pct_both = round(100*mean(!is.na(o) & !is.na(s)), 1),
    n_analysis = sum(ok), pct_analysis = round(100*mean(ok), 1))
}))
print(inv, row.names = FALSE)
write.csv(inv, file.path(OUT, "hrs_wave_inventory_2026-09-26.csv"), row.names = FALSE, fileEncoding = "UTF-8")

best <- inv[which.max(inv$n_analysis), ]
cat("\n可用样本最大的波次: w", best$wave, "（n=", best$n_analysis, "）\n", sep = "")
# 两波纵向
for (pair in list(c(8, 13), c(9, 13), c(10, 13), c(best$wave, 13))) {
  a <- pair[1]; b <- pair[2]
  ok <- function(w) { g <- function(p) paste0("r", w, p)
    complete.cases(m[, c(g("agey_b"), g("cesd"), g("adl5a"), g("ser7"), g("orient"))]) & !is.na(m$raedyrs) }
  cat(sprintf("  纵向 w%d → w%d：两波完整 n = %d\n", a, b, sum(ok(a) & ok(b))))
}
md <- c("# HRS 两文件合并与波次样本盘点（2026-09-26）", "",
  sprintf("- RAND HRS 2022：%d 行；Gateway Harmonized HRS D：%d 行；按 `hhidpn` 合并后 **%d 行**。", nrow(rd), nrow(gw), nrow(m)),
  "- 负值缺失码已统一转为 NA。", "",
  "## 各波可分析样本（合并后）", "", paste(capture.output(print(inv, row.names = FALSE)), collapse = "\n"), "",
  sprintf("## 结论"), "", sprintf("- 单波最大可用样本：**w%d（n=%d）**。", best$wave, best$n_analysis),
  "- 纵向（两波都完整）样本明显更小，见上方打印；**A11 需据此在 SAP §9.0 的方案 A（纵向）与方案 B（同时性）之间定夺**。", "",
  "*脚本：code/17_hrs_merge_inventory.R*")
writeLines(md, file.path(OUT, "hrs_wave_inventory_2026-09-26.md"), useBytes = TRUE)
