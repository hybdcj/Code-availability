# repair waves_cognition_2026-09-13.rds in 2020 wave IDn all missing issue
# ---------------------------------------------------------------------------
# issue: 03_waves_build.R within 2020 wave assembled separately w2020 data frame, column list omitted IDn, 
#       bind_rows at that point converted 2020 rows with IDn all filled as NA(17,730 rows), 
#       prevents longitudinal linkage to other waves("participants complete at five waves"degenerates to 0 participants).
# repair: 2020 of ID column already is 12 digit standardization ID, use directly norm_id(ID) supplement IDn.
#       script is idempotent——rerunning repaired data changes nothing.
# notes: R at Windows Chinese path I/O unstable under, script handles English temporary paths only, 
#       external process(PowerShell)handles backup and backfill.
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
p("input: ", IN_RDS, "(", nrow(dat), " rows × ", ncol(dat), " columns)")

## ---------- 1) pre repair diagnostics ----------
p("\n== before repair: each wave IDn missingness ==")
before <- dat %>% group_by(year) %>%
  summarise(rows = n(), idn_na = sum(is.na(IDn)),
            id_ok = sum(!is.na(ID)), .groups = "drop")
p(capture.output(print(as.data.frame(before))))

nw_before <- dat %>% group_by(IDn) %>% summarise(nw = n_distinct(year), .groups = "drop")
p("individual count before repair(by IDn): ", nrow(nw_before))
p("five wave complete count before repair: ", sum(nw_before$nw == 5))
p("maximum follow up waves before repair: ", max(nw_before$nw))

## ---------- 2) repair ----------
stopifnot("ID" %in% names(dat))
need <- is.na(dat$IDn) & !is.na(dat$ID)
dat$IDn[need] <- norm_id(dat$ID[need])
p("\n== repair action ==")
p("by ID complete IDn row count of: ", sum(need))

## ---------- 3) post repair verification ----------
p("\n== post repair verification ==")
p("IDn still missing: ", sum(is.na(dat$IDn)))
p("IDn length distribution: ", paste(names(table(nchar(dat$IDn))), table(nchar(dat$IDn)),
                          sep = ":", collapse = " "))
after <- dat %>% group_by(year) %>%
  summarise(rows = n(), idn_na = sum(is.na(IDn)), .groups = "drop")
p(capture.output(print(as.data.frame(after))))

nw_after <- dat %>% group_by(IDn) %>% summarise(nw = n_distinct(year), .groups = "drop")
p("individual count after repair(by IDn): ", nrow(nw_after))
p("five wave complete count after repair: ", sum(nw_after$nw == 5))
p("maximum follow up waves after repair: ", max(nw_after$nw))
p("participants by follow up wave count: ", paste(names(table(nw_after$nw)), table(nw_after$nw),
                          sep = ":", collapse = " "))

## key consistency checks: 2020 and 2011 of IDn must use identical coding
ids11 <- dat$IDn[dat$year == 2011]; ids20 <- dat$IDn[dat$year == 2020]
p("2011 and 2020 intersection count: ", length(intersect(ids11, ids20)))
p("2020 not present in 2011 participant count: ", length(setdiff(ids20, ids11)))

## 2020 internal consistency: IDn should equal ID(12 digits unchanged)
p("2020 within rows IDn == ID proportion of: ",
  sprintf("%.4f", mean(dat$IDn[dat$year == 2020] == dat$ID[dat$year == 2020])))

saveRDS(dat, OUT_RDS)
p("\n output: ", OUT_RDS)
writeLines(L, OUT_LOG)
cat(paste(L, collapse = "\n"), "\n")
