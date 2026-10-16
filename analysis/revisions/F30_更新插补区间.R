invisible(Sys.setlocale('LC_CTYPE', 'Chinese_China.utf8'))
.libPaths(c('D:/R/library',.libPaths()))
P <- 'D:/BaiduSyncdisk/科研/徐中豪'
ROOT <- file.path(P,'16.论文修订第三轮')
OUT <- file.path(ROOT,'_核查/输出')
BACK <- file.path(ROOT,'_核查/历史归档/F30_区间修正前_2026-10-15')
src <- file.path(P,'code/29_A15_sensitivity_analyses.R')
for(e in parse(src,encoding='UTF-8')) if(is.call(e) && identical(e[[1]],as.name('<-')) &&
 is.symbol(e[[2]]) && as.character(e[[2]])=='summarise_mi_performance') eval(e)
perf <- readRDS(file.path(OUT,'F30_逐次性能.rds'))
fixed <- lapply(perf,summarise_mi_performance)
metrics <- c('auc_model','auc_score','cal_slope','cal_intercept')
bounds <- unlist(lapply(c('auc','auc_score','cal_slope','cal_intercept'),
 function(x)paste0(x,c('_low','_high'))))
for(kind in names(perf)) {
  # Quantiles independently calculated outside the helper; row-order invariance.
  for(metric in metrics) {
    prefix <- if(metric=='auc_model') 'auc' else metric
    stopifnot(identical(fixed[[kind]][[paste0(prefix,'_low')]],
      unname(stats::quantile(perf[[kind]][[metric]],.025,type=7))),
      identical(fixed[[kind]][[paste0(prefix,'_high')]],
      unname(stats::quantile(perf[[kind]][[metric]],.975,type=7))))
  }
  shuffled <- summarise_mi_performance(perf[[kind]][rev(seq_len(20)),])
  num <- vapply(shuffled,is.numeric,logical(1))
  stopifnot(max(abs(as.matrix(shuffled[num])-as.matrix(fixed[[kind]][num])))<1e-12)
}
stopifnot(fixed$cov$auc_high>fixed$cov$auc_low,
 fixed$out$auc_high>fixed$out$auc_low,
 fixed$cov$auc_score_low==fixed$cov$auc_score_high)
dirs <- c('output/13_A15_敏感性分析','output/14_A16_图表/补充表_A15敏感性')
csv_names <- c('A15_缺失_MICE_协变量_性能汇总.csv','A15_缺失_MICE_结局_性能汇总.csv','A15_内部敏感性汇总.csv')
csv_rel <- unlist(lapply(dirs,function(d)file.path(d,csv_names)))
md_rel <- file.path(dirs,'A15_敏感性分析报告_2026-09-26.md')
rds_rel <- 'output/13_A15_敏感性分析/A15_results_2026-09-26.rds'
rels <- c(csv_rel,md_rel,rds_rel)
for(rel in rels) {
  dest <- file.path(BACK,rel)
  dir.create(dirname(dest),recursive=TRUE,showWarnings=FALSE)
  if(!file.exists(dest)) stopifnot(file.copy(file.path(P,rel),dest))
}
old <- readRDS(file.path(BACK,rds_rel))
z <- old
kind_types <- c(cov='mice_cov',out='mice_outcome')
repair <- function(tab) {
  before <- tab
  for(kind in names(kind_types)) {
    ii <- which(tab$model_type==kind_types[[kind]])
    if(!length(ii)) next
    stopifnot(length(ii)==1)
    fields <- intersect(bounds,names(tab))
    tab[ii,fields] <- fixed[[kind]][1,fields]
  }
  stopifnot(identical(tab[setdiff(names(tab),bounds)],before[setdiff(names(tab),bounds)]))
  keep <- !before$model_type %in% kind_types
  stopifnot(identical(tab[keep,],before[keep,]))
  tab
}
for(nm in c('perf_cov_sum','perf_out_sum','summary_internal')) z[[nm]] <- repair(z[[nm]])
other <- setdiff(names(old),c('perf_cov_sum','perf_out_sum','summary_internal'))
stopifnot(all(vapply(other,function(n)identical(old[[n]],z[[n]]),logical(1))))
for(rel in csv_rel) {
  tab <- read.csv(file.path(BACK,rel),fileEncoding='UTF-8',check.names=FALSE)
  write.csv(repair(tab),file.path(P,rel),row.names=FALSE,fileEncoding='UTF-8')
}
saveRDS(z,file.path(P,rds_rel))
note <- 'MICE 连续模型在每个插补数据集中重新拟合并在同一数据上评价，属于表观性能；校准斜率约1、截距约0不代表外部校准。整数评分固定应用原评分规则。结局缺失分析的事件数为20次插补的平均值，不是实际观察到的整数事件数。'
for(rel in md_rel) {
  lines <- readLines(file.path(BACK,rel),encoding='UTF-8',warn=FALSE)
  for(kind in c('cov','out')) {
    tag <- if(kind=='cov') '协变量缺失；结局观察者' else '结局缺失，MAR；含固定队列辅助变量'
    ii <- which(grepl(tag,lines,fixed=TRUE) & grepl('（插补间 ',lines,fixed=TRUE))
    stopifnot(length(ii)==1)
    values <- fixed[[kind]]
    replacement <- sprintf('%.6f（插补间 %.6f–%.6f）',values$auc_model,values$auc_low,values$auc_high)
    lines[ii] <- sub('[0-9.]+（插补间 [0-9.]+–[0-9.]+）',replacement,lines[ii])
    if(kind=='out') lines[ii] <- sub('|2,451 |','|2,450.9（插补均值） |',lines[ii],fixed=TRUE)
  }
  anchor <- which(grepl('MICE 行括号内为 20 次插补',lines,fixed=TRUE))
  stopifnot(length(anchor)==1)
  lines <- append(lines,c('',note),after=anchor)
  writeLines(lines,file.path(P,rel),useBytes=TRUE)
}
delivery <- rbind(as.data.frame(z$perf_cov_sum),as.data.frame(z$perf_out_sum))
delivery$performance_scope <- '连续模型逐次重新拟合的表观性能；整数评分规则固定'
delivery$imputations <- 20L
delivery$events_min <- c(min(perf$cov$events_model),min(perf$out$events_model))
delivery$events_max <- c(max(perf$cov$events_model),max(perf$out$events_model))
delivery$events_type <- c('实际观察事件数','20次插补事件数均值')
write.csv(delivery,file.path(ROOT,'F30_多重插补性能_修正版.csv'),row.names=FALSE,fileEncoding='UTF-8')
# Read back all updated files and verify that F29's corrected model row survives.
disk <- readRDS(file.path(P,rds_rel))
stopifnot(isTRUE(all.equal(disk,z,check.environment=FALSE)),
 identical(disk$summary_internal[disk$summary_internal$model_type=='unweighted',],
 old$summary_internal[old$summary_internal$model_type=='unweighted',]))
for(rel in csv_rel) {
  disk_csv <- read.csv(file.path(P,rel),fileEncoding='UTF-8',check.names=FALSE)
  expected <- repair(read.csv(file.path(BACK,rel),fileEncoding='UTF-8',check.names=FALSE))
  stopifnot(isTRUE(all.equal(disk_csv,expected,tolerance=1e-12,check.attributes=FALSE)))
}
all_rel <- c('code/29_A15_sensitivity_analyses.R',rels)
write.csv(data.frame(path=all_rel,before_md5=unname(tools::md5sum(file.path(BACK,all_rel))),
 after_md5=unname(tools::md5sum(file.path(P,all_rel)))),
 file.path(OUT,'F30_更新文件校验.csv'),row.names=FALSE,fileEncoding='UTF-8')
checks <- c('Original imputation draws preserved: PASS; no MICE rerun',
 'Corrected quantiles versus independent stats::quantile: PASS',
 'Imputation-order invariance: PASS',
 'All point estimates and non-MICE rows unchanged: PASS',
 'All other RDS objects unchanged, including frozen model, coefficient tables and F29 results: PASS',
 'Six CSVs and saved RDS read-back: PASS',
 'Two report tables and interpretation notes updated; figures unchanged',
 capture.output(print(delivery[c('model_type','auc_model','auc_low','auc_high','auc_score','auc_score_low','auc_score_high')],digits=12)),
 capture.output(sessionInfo()))
writeLines(checks,file.path(OUT,'F30_更新核验日志.txt'),useBytes=TRUE)
cat(paste(head(checks,7),collapse='\n'),'\n')
print(delivery[c('model_type','auc_model','auc_low','auc_high','auc_score','auc_score_low','auc_score_high')],digits=10)
