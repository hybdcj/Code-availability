invisible(try(Sys.setlocale("LC_ALL", "zh_CN.UTF-8"), silent = TRUE))
suppressPackageStartupMessages({ library(haven); library(dplyr) })
if (dir.exists("D:/R/library")) .libPaths(c("D:/R/library", .libPaths()))
OUT <- "C:/y03tmp/out"
con <- file(file.path(OUT, "F16_diag.txt"), open = "wt", encoding = "UTF-8")
say <- function(...) { s <- paste0(..., collapse = ""); cat(s, "\n"); cat(s, "\n", file = con) }
neg2na <- function(x) { x <- as.numeric(x); x[x < 0] <- NA; x }
step <- function(d, nm, conds) {
  say("  ", nm, " = ", nrow(d))
  for (c in names(conds)) say("    +", c, " -> ", sum(conds[[c]]))
}

say("======== A. HRS ========")
rd <- as.data.frame(read_dta("D:/codex_tmp/hrs/randhrs1992_2022v1.dta",
  col_select = any_of(c("hhidpn","raedyrs","r10agey_b","r10cesd","r10adl5a",
                        "r10imrc","r10dlrc","r12imrc","r12dlrc"))))
say("原始行数 ", nrow(rd))
rd <- as.data.frame(lapply(rd, neg2na))
say("age 非缺失 ", sum(!is.na(rd$r10agey_b)),
    "| age>=45 ", sum(!is.na(rd$r10agey_b) & rd$r10agey_b >= 45))
say("raedyrs 非缺失 ", sum(!is.na(rd$raedyrs)),
    "| r10cesd 非缺失 ", sum(!is.na(rd$r10cesd)),
    "| r10imrc 非缺失 ", sum(!is.na(rd$r10imrc)),
    "| r10dlrc 非缺失 ", sum(!is.na(rd$r10dlrc)),
    "| r12imrc 非缺失 ", sum(!is.na(rd$r12imrc)),
    "| r12dlrc 非缺失 ", sum(!is.na(rd$r12dlrc)))
k <- !is.na(rd$r10agey_b) & rd$r10agey_b >= 45
say("A11 序列（先 age>=45）：", sum(k),
    " -> +edu ", sum(k & !is.na(rd$raedyrs)),
    " -> +dep ", sum(k & !is.na(rd$raedyrs) & !is.na(rd$r10cesd)),
    " -> +mem_base ", sum(k & !is.na(rd$raedyrs) & !is.na(rd$r10cesd) &
                          !is.na(rd$r10imrc) & !is.na(rd$r10dlrc)))
say("A11 记录：age>=45 = 21697；baseline eligible = 20232；fu memory = 16345")

say("\n======== B. ELSA ========")
el <- as.data.frame(read_dta("D:/codex_tmp/elsa/gh_elsa_h.dta",
  col_select = any_of(c("idauniq","raeducl","r4agey","r4cesd","r4adlfive",
                        "r4imrc","r4dlrc","r6imrc","r6dlrc"))))
say("原始行数 ", nrow(el))
el <- as.data.frame(lapply(el, neg2na))
k <- !is.na(el$r4agey) & el$r4agey >= 45
say("A11 序列：", sum(k),
    " -> +edu ", sum(k & !is.na(el$raeducl)),
    " -> +dep ", sum(k & !is.na(el$raeducl) & !is.na(el$r4cesd)),
    " -> +mem_base ", sum(k & !is.na(el$raeducl) & !is.na(el$r4cesd) &
                          !is.na(el$r4imrc) & !is.na(el$r4dlrc)))
say("A11 记录：age>=45 = 10956；baseline eligible = 9620；fu memory = 7523")

say("\n======== C. SHARE ========")
sh <- as.data.frame(read_dta("D:/codex_tmp/share/easySHARE_rel9-0-0.dta"))
say("原始行数 ", nrow(sh), "｜含 wave 变量：", "wave" %in% names(sh))
say("wave 取值频数：")
print(table(as.numeric(sh$wave), useNA = "ifany"))
say("mergeid 在 wave4 的重复数：", sum(duplicated(sh$mergeid[as.numeric(sh$wave) == 4])))
say("mergeid 在 wave6 的重复数：", sum(duplicated(sh$mergeid[as.numeric(sh$wave) == 6])))
sh$mergeid <- as.character(sh$mergeid)
for (v in setdiff(names(sh), "mergeid")) sh[[v]] <- neg2na(sh[[v]])
sh$isced1997_r[!is.na(sh$isced1997_r) & sh$isced1997_r > 6] <- NA
a <- sh[sh$wave == 4, ]; b <- sh[sh$wave == 6, ]
say("wave4 行 = ", nrow(a), "（A11 记录 57982）｜wave6 行 = ", nrow(b))
say("wave4 中 age>=45 = ", sum(!is.na(a$age) & a$age >= 45), "（A11 记录 57645）")
mm <- merge(a, b, by = "mergeid")
say("merge 后行 = ", nrow(mm))
k <- !is.na(mm$age) & mm$age >= 45
say("merge 后 age>=45 = ", sum(k))
say("  -> +edu ", sum(k & !is.na(mm$isced1997_r)),
    " -> +dep ", sum(k & !is.na(mm$isced1997_r) & !is.na(mm$eurod)),
    " -> +mem_base ", sum(k & !is.na(mm$isced1997_r) & !is.na(mm$eurod) &
                          !is.na(mm$recall_1) & !is.na(mm$recall_2)))
say("A11 记录：baseline eligible = 54034；fu memory = 31452；risk = 26879")
k2 <- !is.na(a$age) & a$age >= 45
say("（仅 wave4 口径）age>=45 ", sum(k2),
    " -> +edu ", sum(k2 & !is.na(a$isced1997_r)),
    " -> +dep ", sum(k2 & !is.na(a$isced1997_r) & !is.na(a$eurod)),
    " -> +mem_base ", sum(k2 & !is.na(a$isced1997_r) & !is.na(a$eurod) &
                          !is.na(a$recall_1) & !is.na(a$recall_2)))
close(con)
cat("done\n")
