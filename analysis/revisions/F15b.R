invisible(try(Sys.setlocale("LC_ALL", "zh_CN.UTF-8"), silent = TRUE))
suppressPackageStartupMessages({ library(data.table) })
if (dir.exists("D:/R/library")) .libPaths(c("D:/R/library", .libPaths()))
TMP <- "C:/y03tmp/tmp"; OUT <- "C:/y03tmp/out"
con <- file(file.path(OUT, "F15b_probe.txt"), open = "wt", encoding = "UTF-8")
say <- function(...) { s <- paste0(..., collapse = ""); cat(s, "\n"); cat(s, "\n", file = con) }
w <- as.data.table(readRDS(file.path(TMP, "waves_cognition.rds")))
say("列名：", paste(names(w), collapse = ", "))
for (y in c(2011, 2015, 2018)) {
  d <- w[year == y]
  say("\n-- 年份 ", y, "：n = ", nrow(d), " --")
  a <- as.numeric(d$age)
  say("年龄缺失 ", sum(is.na(a)), " | min ", min(a, na.rm = TRUE), " max ", max(a, na.rm = TRUE))
  say("分位：", paste(names(quantile(a, c(0, .01, .05, .25, .5, .75, .95, .99, 1), na.rm = TRUE)),
                      quantile(a, c(0, .01, .05, .25, .5, .75, .95, .99, 1), na.rm = TRUE),
                      sep = "=", collapse = "  "))
  say("<45 岁人数：", sum(a < 45, na.rm = TRUE), " | 45-49：", sum(a >= 45 & a < 50, na.rm = TRUE),
      " | >=50：", sum(a >= 50, na.rm = TRUE))
  low <- d[a < 45]
  if (nrow(low)) {
    say("  <45 岁者的记忆分可用数：", sum(!is.na(as.numeric(low$imrc)) & !is.na(as.numeric(low$dlrc))),
        " | 教育缺失：", sum(is.na(low$edu_group)))
    say("  <45 岁者 id 样例：", paste(head(as.character(low$IDn), 5), collapse = ", "))
  }
  hh <- d[a >= 45 & a < 50]
  if (nrow(hh)) say("  45-49 岁者教育分布：",
                    paste(names(table(hh$edu_group)), as.integer(table(hh$edu_group)), sep = "=", collapse = "  "))
}
close(con)
cat("done\n")
