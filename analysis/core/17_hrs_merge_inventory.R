# ============================================================================
# 17_hrs_merge_inventory.R —— HRS two files joined by ID combine + analyzable samples by wave
# purpose: select waves for external validation(SAP §9.0 scheme A/B)and for A11 provide sample size evidence
# generate: 2026-09-26
# output: external_data/HRS/logs/hrs_wave_inventory_2026-09-26.csv / .md
# ============================================================================
invisible(try(Sys.setlocale("LC_ALL", "zh_CN.UTF-8"), silent = TRUE))
.libPaths(c("D:/R/library", .libPaths()))
suppressPackageStartupMessages({ library(haven); library(dplyr) })
HRS <- "D:/BaiduSyncdisk/research/Xu_Zhonghao/external_data/HRS"
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
cat("RAND row count:", nrow(rd), " Gateway row count:", nrow(gw), "\n")
m <- inner_join(rd, gw, by = "hhidpn", suffix = c(".rd", ".gw"))
cat("by hhidpn after combining:", nrow(m), " rows(RAND match rate ", sprintf("%.1f%%", 100*nrow(m)/nrow(rd)), ")\n\n")

inv <- do.call(rbind, lapply(waves, function(w) {
  g <- function(p) paste0("r", w, p)
  o <- m[[g("orient")]]; s <- m[[g("ser7")]]
  # definition: 4 predictor variables(age agey_b, education raedyrs, depression cesd, ADL adl5a)+ two outcome components(orient + ser7)
  ok <- complete.cases(m[, c(g("agey_b"), g("cesd"), g("adl5a"), g("ser7"), g("orient"))]) & !is.na(m$raedyrs)
  data.frame(wave = w,
    pct_orient = round(100*mean(!is.na(o)), 1), pct_ser7 = round(100*mean(!is.na(s)), 1),
    pct_both = round(100*mean(!is.na(o) & !is.na(s)), 1),
    n_analysis = sum(ok), pct_analysis = round(100*mean(ok), 1))
}))
print(inv, row.names = FALSE)
write.csv(inv, file.path(OUT, "hrs_wave_inventory_2026-09-26.csv"), row.names = FALSE, fileEncoding = "UTF-8")

best <- inv[which.max(inv$n_analysis), ]
cat("\n wave with largest available sample: w", best$wave, "(n=", best$n_analysis, ")\n", sep = "")
# two wave longitudinal
for (pair in list(c(8, 13), c(9, 13), c(10, 13), c(best$wave, 13))) {
  a <- pair[1]; b <- pair[2]
  ok <- function(w) { g <- function(p) paste0("r", w, p)
    complete.cases(m[, c(g("agey_b"), g("cesd"), g("adl5a"), g("ser7"), g("orient"))]) & !is.na(m$raedyrs) }
  cat(sprintf("  longitudinal w%d → w%d: complete at both waves n = %d\n", a, b, sum(ok(a) & ok(b))))
}
md <- c("# HRS merge two files and inventory wave samples(2026-09-26)", "",
  sprintf("- RAND HRS 2022: %d rows; Gateway Harmonized HRS D: %d rows; by `hhidpn` after combining **%d rows**.", nrow(rd), nrow(gw), nrow(m)),
  "- negative missing codes converted consistently to NA.", "",
  "## analyzable sample by wave(after combining)", "", paste(capture.output(print(inv, row.names = FALSE)), collapse = "\n"), "",
  sprintf("## conclusion"), "", sprintf("- largest single wave sample: **w%d(n=%d)**.", best$wave, best$n_analysis),
  "- longitudinal(complete at both waves)substantially smaller sample, see printed output above; **A11 use this to choose between SAP §9.0 scheme of A(longitudinal)and scheme B(concurrent)make selection**.", "",
  "*scripts: code/17_hrs_merge_inventory.R*")
writeLines(md, file.path(OUT, "hrs_wave_inventory_2026-09-26.md"), useBytes = TRUE)
