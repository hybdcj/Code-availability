invisible(Sys.setlocale('LC_CTYPE', 'Chinese_China.utf8'))
.libPaths(c('D:/R/library', .libPaths()))
suppressPackageStartupMessages(library(dplyr))
options(stringsAsFactors = FALSE)
P <- 'D:/BaiduSyncdisk/科研/徐中豪'
ROOT <- file.path(P, '16.论文修订第三轮')
OUT <- file.path(ROOT, '_核查/输出')
BACK <- file.path(ROOT, '_核查/历史归档/F29_修正前_2026-10-15')
source_path <- file.path(P, 'code/29_A15_sensitivity_analyses.R')
rds_rel <- 'output/13_A15_敏感性分析/A15_results_2026-09-26.rds'
csv_rel <- unlist(lapply(c('output/13_A15_敏感性分析', 'output/14_A16_图表/补充表_A15敏感性'),
 function(d) file.path(d, c('A15_结局定义_性能比较.csv', 'A15_内部敏感性汇总.csv'))))
fig_rel <- unlist(lapply(c('output/09_图表', 'output/14_A16_图表/补充图'),
 function(d) file.path(d, paste0('FigureS8_A15_sensitivity_2026-09-26.', c('png','tif')))))
rels <- c(rds_rel, csv_rel, fig_rel)
for (rel in rels) {
  dest <- file.path(BACK, rel)
  dir.create(dirname(dest), recursive=TRUE, showWarnings=FALSE)
  if (!file.exists(dest)) stopifnot(file.copy(file.path(P, rel), dest))
}
write.csv(data.frame(path=c('code/29_A15_sensitivity_analyses.R', rels),
 md5=unname(tools::md5sum(file.path(BACK, c('code/29_A15_sensitivity_analyses.R', rels))))),
 file.path(OUT, 'F29_修正前文件校验.csv'), row.names=FALSE, fileEncoding='UTF-8')
# Load only reusable functions, never execute the full original pipeline.
expr <- parse(source_path, encoding='UTF-8')
wanted <- c('auc_ci_unweighted', 'calibration_glm', 'model_metrics', 'fixed_lp_manual')
loaded <- character()
for (e in expr) if (is.call(e) && identical(e[[1]], as.name('<-')) &&
 is.symbol(e[[2]]) && as.character(e[[2]]) %in% wanted) {
  eval(e); loaded <- c(loaded, as.character(e[[2]]))
}
stopifnot(setequal(loaded, wanted))
old <- readRDS(file.path(BACK, rds_rel))
z <- old
a07 <- readRDS('C:/y03tmp/tmp/a07_memory_v3.rds')
d <- z$main
stopifnot(!anyDuplicated(d$KEY), !anyDuplicated(a07$data$KEY),
 setequal(d$KEY, a07$data$KEY))
idx <- match(d$KEY, a07$data$KEY)
bad <- as.numeric(predict(z$model_main, type='response'))
good <- as.numeric(predict(z$model_main, newdata=d, type='response'))
stopifnot(max(abs(good - bad[idx])) < 1e-12,
 max(abs(good - plogis(fixed_lp_manual(d, coef(z$model_main))))) < 1e-12,
 all(d$event == a07$data$incident_2015[idx]))
scenario <- '主结局：4年跨阈值 ≤4'
corrected <- model_metrics('结局定义', scenario, 'unweighted', d$event, good,
 score=d$score, n_note=nrow(d), B=0)
reproduced <- model_metrics('结局定义', scenario, 'unweighted', d$event, bad,
 score=d$score, n_note=nrow(d), B=0)
oldrow <- old$met_outcome_alt[old$met_outcome_alt$scenario == scenario, ]
stopifnot(isTRUE(all.equal(reproduced, oldrow, check.attributes=FALSE, tolerance=1e-10)))
# Regression check: reverse row order, predict on that data, and compare all metrics.
rev_d <- d[rev(seq_len(nrow(d))), ]
rev_p <- as.numeric(predict(z$model_main, newdata=rev_d, type='response'))
rev_m <- model_metrics('结局定义', scenario, 'unweighted', rev_d$event, rev_p,
 score=rev_d$score, n_note=nrow(rev_d), B=0)
same_metrics <- function(a,b,tol=1e-10) {
  numeric_cols <- names(a)[vapply(a,is.numeric,logical(1))]
  other_cols <- setdiff(names(a),numeric_cols)
  isTRUE(all.equal(a[other_cols],b[other_cols],check.attributes=FALSE)) &&
    max(abs(as.matrix(a[numeric_cols])-as.matrix(b[numeric_cols]))) < tol
}
stopifnot(same_metrics(corrected, rev_m))
# Adjacent rows: reconstruct from saved coefficients; no model refitting.
for (sc in setdiff(z$met_outcome_alt$scenario, scenario)) {
  co <- z$outcome_alt[z$outcome_alt$scenario == sc, ]
  f <- if (grepl('mem2011', sc)) ~ mem2011 + age + edu_group + dep + adl else ~ age + edu_group + dep + adl
  mm <- model.matrix(f, d)
  b <- setNames(co$beta, co$term)
  stopifnot(setequal(colnames(mm), names(b)))
  p <- plogis(drop(mm %*% b[colnames(mm)]))
  row <- z$met_outcome_alt[z$met_outcome_alt$scenario == sc, ]
  y <- if (row$model_type == 'alt_binary') d$decline2 else d$event
  m <- model_metrics('结局定义', sc, row$model_type, y, p, score=d$score, n_note=nrow(d), B=0)
  stopifnot(same_metrics(m, row, tol=1e-9))
}
replace_row <- function(tab) {
  ii <- which(tab$scenario == scenario)
  stopifnot(length(ii)==1)
  cols <- intersect(names(tab), names(corrected))
  tab[ii, cols] <- corrected[1, cols]
  tab
}
z$met_outcome_alt <- replace_row(z$met_outcome_alt)
z$summary_internal <- replace_row(z$summary_internal)
unchanged <- setdiff(names(z), c('met_outcome_alt','summary_internal'))
stopifnot(all(vapply(unchanged, function(n) identical(z[[n]], old[[n]]), logical(1))))
for (n in c('met_outcome_alt','summary_internal')) {
  keep <- old[[n]]$scenario != scenario
  stopifnot(identical(z[[n]][keep, ], old[[n]][keep, ]))
}
for (rel in csv_rel) {
  tab <- read.csv(file.path(BACK, rel), fileEncoding='UTF-8', check.names=FALSE)
  repaired <- replace_row(tab)
  stopifnot(identical(tab[tab$scenario!=scenario, ], repaired[repaired$scenario!=scenario, ]))
  write.csv(repaired, file.path(P, rel), row.names=FALSE, fileEncoding='UTF-8')
}
saveRDS(z, file.path(P, rds_rel))
write.csv(z$met_outcome_alt, file.path(ROOT, 'F29_S29结局定义性能_修正版.csv'),
 row.names=FALSE, fileEncoding='UTF-8')
fields <- names(corrected)[vapply(corrected, is.numeric, logical(1))]
write.csv(data.frame(metric=fields, before=as.numeric(oldrow[1,fields]),
 after=as.numeric(corrected[1,fields])), file.path(OUT, 'F29_主模型性能修正前后.csv'),
 row.names=FALSE, fileEncoding='UTF-8')
# Regenerate only the affected plot section with corrected data and wider AUC axis.
summary_internal <- z$summary_internal
key_or <- z$key_or
FIG <- file.path(P, 'output/09_图表')
src <- readLines(source_path, encoding='UTF-8', warn=FALSE)
start <- grep('^plot_tab <- summary_internal', src)
end <- grep('^# 11[.] 保存 RDS', src) - 2L
stopifnot(length(start)==1, length(end)==1, end>start)
eval(parse(text=src[start:end]))
for (ext in c('png','tif')) {
  nm <- paste0('FigureS8_A15_sensitivity_2026-09-26.', ext)
  stopifnot(file.copy(file.path(FIG,nm), file.path(P,'output/14_A16_图表/补充图',nm), overwrite=TRUE))
  stopifnot(file.copy(file.path(FIG,nm), file.path(ROOT,paste0('F29_内部敏感性修正版.',ext)), overwrite=TRUE))
}
checks <- c(paste0('n=', nrow(d), '; events=',sum(d$event)),
 paste0('different row positions=',sum(idx != seq_along(idx))),
 paste0('outcome disagreements before matching=',sum(d$event != a07$data$incident_2015)),
 paste0('max prediction error after KEY matching=', max(abs(good-bad[idx]))),
 'Reversed-row invariance: PASS', 'Adjacent alternative outcomes from saved coefficients: PASS',
 'All other saved RDS objects unchanged: PASS', 'Other rows in 4 CSV files unchanged: PASS',
 paste0('AUC plot limits: ',paste(auc_xlim,collapse=', ')), capture.output(print(corrected, digits=12)),
 capture.output(sessionInfo()))
writeLines(checks, file.path(OUT,'F29_核验日志.txt'), useBytes=TRUE)
cat(paste(checks[seq_len(9)],collapse='\n'), '\n')
print(corrected, digits=10)
