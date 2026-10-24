invisible(Sys.setlocale('LC_CTYPE', 'Chinese_China.utf8'))
.libPaths(c('D:/R/library', .libPaths()))
suppressPackageStartupMessages({library(dplyr); library(mice)})
P <- 'D:/BaiduSyncdisk/research/Xu_Zhonghao'
ROOT <- file.path(P,'16.manuscript_revision_round_three')
OUT <- file.path(ROOT,'_audit/output')
src <- file.path(P,'code/29_A15_sensitivity_analyses.R')
archive <- file.path(ROOT,'_audit/historical_archive/F30_before_interval_correction_2026_10_15/output/13_A15_sensitivity_analysis/A15_results_2026-09-26.rds')
z <- readRDS(if(file.exists(archive)) archive else
 file.path(P,'output/13_A15_sensitivity_analysis/A15_results_2026-09-26.rds'))
needed <- c('auc_ci_unweighted','calibration_glm','model_metrics','fixed_score')
for(e in parse(src,encoding='UTF-8')) if(is.call(e) && identical(e[[1]],as.name('<-')) &&
 is.symbol(e[[2]]) && as.character(e[[2]]) %in% needed) eval(e)
all_perf <- list(); comparisons <- list(); diagnoses <- list(); coef_checks <- list()
for(kind in c('cov','out')) {
  imp <- z[[paste0('imp_',kind)]]
  stopifnot(imp$m==20)
  fits <- vector('list',imp$m); perf <- vector('list',imp$m)
  for(i in seq_len(imp$m)) {
    d <- mice::complete(imp,i)
    if(kind=='out') d$event <- as.integer(d$mem2015<=4)
    stopifnot(!anyNA(d[c('event','age','edu_group','dep','adl')]))
    mm <- glm(event ~ age + edu_group + dep + adl, data=d, family=binomial())
    stopifnot(mm$converged)
    fits[[i]] <- mm
    perf[[i]] <- model_metrics('missing',paste0(kind,'_imputation_',i),
      if(kind=='cov') 'mice_cov' else 'mice_outcome',d$event,
      predict(mm,newdata=d,type='response'),score=fixed_score(d),n_note=nrow(d),B=0)
    perf[[i]]$imputation <- i
  }
  perf <- bind_rows(perf)
  all_perf[[kind]] <- perf
  old <- z[[paste0('perf_',kind,'_sum')]]
  for(v in c('auc_model','auc_score','cal_slope','cal_intercept','brier')) {
    stopifnot(abs(mean(perf[[v]])-old[[v]]) < 1e-10)
    prefix <- if(v=='auc_model') 'auc' else v
    if(v!='brier') stopifnot(abs(old[[paste0(prefix,'_low')]]-old[[v]])<1e-12,
      abs(old[[paste0(prefix,'_high')]]-old[[v]])<1e-12)
    comparisons[[paste(kind,v)]] <- data.frame(analysis=kind,metric=v,
      mean=mean(perf[[v]]),old_low=if(v=='brier') NA_real_ else old[[paste0(prefix,'_low')]],
      old_high=if(v=='brier') NA_real_ else old[[paste0(prefix,'_high')]],
      p025=unname(quantile(perf[[v]],.025,type=7)),p975=unname(quantile(perf[[v]],.975,type=7)),
      min=min(perf[[v]]),max=max(perf[[v]]),sd=sd(perf[[v]]),distinct=n_distinct(perf[[v]]))
  }
  # Verify that existing coefficient intervals are pooled, and unaffected by the performance bug.
  pooled <- summary(mice::pool(fits),conf.int=TRUE)
  saved <- z[[if(kind=='cov') 'coef_cov_mi' else 'coef_out_mi']]
  jj <- match(saved$term,pooled$term)
  stopifnot(!anyNA(jj))
  diffs <- c(beta=max(abs(saved$beta-pooled$estimate[jj])),
    se=max(abs(saved$se-pooled$std.error[jj])),
    OR_low=max(abs(saved$OR_low-exp(pooled$`2.5 %`[jj]))),
    OR_high=max(abs(saved$OR_high-exp(pooled$`97.5 %`[jj]))))
  stopifnot(max(diffs)<1e-9)
  coef_checks[[kind]] <- data.frame(analysis=kind,metric=names(diffs),max_abs_error=unname(diffs))
  diagnoses[[kind]] <- data.frame(analysis=kind,variable=names(imp$data),
    missing=vapply(imp$data,function(x)sum(is.na(x)),integer(1)),
    method=unname(imp$method),stringsAsFactors=FALSE)
  cat(kind,': n=',nrow(imp$data),'; m=',imp$m,'; iterations=',imp$iteration,
    '; events range=',range(perf$events_model),'; distinct AUC=',n_distinct(perf$auc_model),'\n')
}
perf_all <- bind_rows(all_perf,.id='analysis')
comparison <- bind_rows(comparisons)
write.csv(perf_all,file.path(OUT,'F30_40_per_imputation_performance.csv'),row.names=FALSE,fileEncoding='UTF-8')
write.csv(comparison,file.path(OUT,'F30_audit_performance_means_and_quantiles.csv'),row.names=FALSE,fileEncoding='UTF-8')
write.csv(bind_rows(diagnoses),file.path(OUT,'F30_imputation_object_missingness_and_methods.csv'),row.names=FALSE,fileEncoding='UTF-8')
write.csv(bind_rows(coef_checks),file.path(OUT,'F30_audit_original_coefficient_intervals.csv'),row.names=FALSE,fileEncoding='UTF-8')
saveRDS(all_perf,file.path(OUT,'F30_per_imputation_performance.rds'))
print(comparison,digits=11,row.names=FALSE)
# Minimal reproduction of scalar overwrite, independent of study results.
bug_demo <- data.frame(auc_model=c(.6,.7,.8)) %>% summarise(
 auc_model=mean(auc_model),auc_low=quantile(auc_model,.025),auc_high=quantile(auc_model,.975))
stopifnot(all(unlist(bug_demo)==.7))
writeLines(c('Original performance means reproduced: PASS',
 'Original pooled coefficient estimates and intervals reproduced: PASS',
 'Original interval endpoints equal means: PASS; scalar-overwrite bug reproduced',
 capture.output(print(comparison,digits=15,row.names=FALSE)),capture.output(sessionInfo())),
 file.path(OUT,'F30_interval_audit_log.txt'),useBytes=TRUE)
