# Recompute screening metrics only; retain the original cohorts, outcomes and scores.
# Usage: Rscript F24_recompute_screening.R [project_directory]
invisible(Sys.setlocale('LC_CTYPE','Chinese_China.utf8'))
.libPaths(c('D:/R/library', .libPaths()))
suppressPackageStartupMessages({library(data.table); library(pROC)})
args <- commandArgs(trailingOnly=TRUE)
root <- if(length(args)) args[1] else stop('Provide project directory')
src <- file.path(root, '_audit', 'scripts'); out <- file.path(root, '_audit', 'output')
archive <- file.path(root, '_audit', 'historical_archive', 'F24_before_screening_metric_correction_2026_10_15', 'output')
dir.create(archive, recursive=TRUE, showWarnings=FALSE)

# Reuse only the already-reviewed preparation expressions; suppress report/file I/O.
prefix_env <- function(filename, marker) {
  lines <- readLines(file.path(src, filename), encoding='UTF-8', warn=FALSE)
  idx <- grep(marker, lines, fixed=TRUE)[1]
  stopifnot(!is.na(idx))
  expr <- parse(text=lines[seq_len(idx-1)])
  e <- new.env(parent=globalenv()); e$say <- function(...) invisible(NULL); e$pt <- e$say
  for(ex in expr) {
    head <- if(is.call(ex)) as.character(ex[[1]]) else ''
    if(head %in% c('fwrite','dir.create','say','pt','invisible')) next
    if(head %in% c('<-','=')) {
      lhs <- paste(deparse(ex[[2]]),collapse='')
      if(lhs %in% c('con','say','pt')) next
    }
    eval(ex,e)
  }
  e
}
ratio <- function(a,b) if(b==0) NA_real_ else a/b
confusion <- function(y,pos) {
  stopifnot(length(y)==length(pos), !anyNA(y), !anyNA(pos), all(y %in% 0:1))
  tp<-sum(pos & y==1); fp<-sum(pos & y==0); tn<-sum(!pos & y==0); fn<-sum(!pos & y==1)
  data.table(n=length(y), events=sum(y), cases_detected=sum(pos),
    actual_screening_fraction=100*mean(pos), TP=tp, FP=fp, TN=tn, FN=fn,
    sensitivity=100*ratio(tp,tp+fn), specificity=100*ratio(tn,tn+fp),
    positive_predictive_value=100*ratio(tp,tp+fp), negative_predictive_value=100*ratio(tn,tn+fn))
}
# Known 2x2 example plus empty-positive/negative edge cases.
test <- confusion(c(1,1,1,0,0,0,0),c(TRUE,TRUE,FALSE,TRUE,FALSE,FALSE,FALSE))
stopifnot(test$TP==2,test$FP==1,test$TN==3,test$FN==1,
          abs(test$specificity-75)<1e-10,abs(test$negative_predictive_value-75)<1e-10,
          is.na(confusion(c(0,1),c(FALSE,FALSE))$positive_predictive_value),
          is.na(confusion(c(0,1),c(TRUE,TRUE))$negative_predictive_value))
test2<-confusion(c(1,1,1,1,0,0,0,0),c(TRUE,TRUE,FALSE,FALSE,TRUE,FALSE,FALSE,FALSE))
stopifnot(test2$specificity==75,test2$negative_predictive_value==60)

rows<-list()
add <- function(source, cohort, outcome, model, y, s, fr=NA_real_, low=NULL, fixed=NULL) {
  stopifnot(length(y)==length(s), !anyNA(s), all(is.finite(s)))
  if(!is.null(fixed)) { th<-fixed; pos<-s>=th; label<-'fixed cutoff'; strategy<-'fixed score threshold, include all ties' }
  else if(is.na(fr)) {
    best<-as.data.frame(coords(roc(y,s,quiet=TRUE,direction='<'),'best',best.method='youden',ret=c('threshold','sensitivity','specificity')))
    best<-best[which.max(best$specificity),,drop=FALSE]
    th<-best$threshold; pos<-s>=th; label<-'optimal Youden cutoff'; strategy<-'apparent optimal Youden cutoff in same sample; choose highest specificity if tied'
  } else {
    th<-sort(s,decreasing=TRUE)[max(1L,floor(fr*length(s)))]; pos<-s>=th
    label<-paste0(fr*100,'%'); strategy<-'score cutoff for target fraction; include all ties'
  }
  z<-cbind(data.table(source=source,cohort_or_window=cohort,outcome_definition=outcome,model=model,
    target_screening_fraction=label,threshold=th,positive_rule='score>=threshold',tie_handling=strategy),confusion(y,pos))
  if(is.na(fr) && is.null(fixed)) stopifnot(abs(z$sensitivity/100-best$sensitivity)<1e-10,abs(z$specificity/100-best$specificity)<1e-10)
  if(!is.null(low)) {
    known<-pos & !is.na(low); baseline_low<-sum(pos & low==1,na.rm=TRUE)
    incident<-sum(pos & y==1 & low==0,na.rm=TRUE)
    z[, `:=`(screen_positive_participants_with_known_baseline_status=sum(known),screen_positive_participants_already_low_scoring_at_baseline=baseline_low,
      fraction_already_low_scoring_at_baseline_among_detected_cases=100*ratio(baseline_low,sum(known)),
      true_incident_events=incident,per_detected_1_incident_case_number_needed_to_screen=ratio(length(y),incident))]
  }
  rows[[length(rows)+1L]]<<-z
}

# F14 / F22: use saved predictions, without modifying BDSI or Tu implementation.
d<-as.data.table(readRDS('C:/y03tmp/tmp/F22_d.rds'))
d1<-copy(d[mem15>=cutoff]); d1[,event:=as.integer(mem18<cutoff)]
d2<-copy(d[mem15>4]); d2[,event:=as.integer(mem18<=4)]
stopifnot(nrow(d1)==6436,sum(d1$event)==1116,nrow(d2)==5872,sum(d2$event)==1157)
arms<-c(score_us='our three item scorecard',score_bdsi='BDSI(existing mapping)',lp_tu='Tu2026(existing reconstruction)')
for(v in names(arms)) for(fr in c(.1,.2,.3,NA_real_))
  add('F14','CHARLS 2015→2018','education strata',arms[[v]],d1$event,d1[[v]],fr)
for(lab in c('T1 education strata','T2 fixed≤4')) {
  x<-if(lab=='T1 education strata') d1 else d2
  for(v in names(arms)) add('F22','CHARLS 2015→2018',lab,arms[[v]],x$event,x[[v]],.3)
}

# F15 preserves the original common-complete-case requirement, including ADL.
e15<-prefix_env('F15.R','## ---------------- 4. sample flow')
for(nm in c('w1','w2')) {
  x<-get(nm,envir=e15); x<-x[e15$core(x)]
  stopifnot(nrow(x)==if(nm=='w1')11753 else 12460)
  for(fr in c(.05,.1,.2,.3,NA_real_)) add('F15',if(nm=='w1')'CHARLS 2011→2015' else 'CHARLS 2015→2018',
    'education strata','our three item scorecard',x$event,x$score,fr,x$low0)
}

for(task in c('F16','F19')) {
  ee<-if(task=='F16') prefix_env('F16_analysis.R','say("cutoff verification') else prefix_env('F19.R','say("threshold: ')
  for(co in ee$CO) {
    x<-ee$D[[co]][!is.na(mem_fu)]
    stopifnot(nrow(x)==c(HRS=16345,ELSA=7529,SHARE=31650)[[co]])
    for(fr in if(task=='F16') c(.1,.2,.3,NA_real_) else c(.1,.2,.3))
      add(task,co,if(task=='F16')paste0('cohort specific≤',unique(x$cut_int)) else 'fixed≤4',
        'our three item scorecard',x$event,x$score,fr,x$low0)
  }
}

# F03: reproduce the two existing integer score definitions; no bootstrap rerun.
e3<-prefix_env('F03_headtohead.R','## ---------------- A: temporal validation of frozen model')
for(nm in names(e3$defs)) {
  def<-e3$defs[[nm]]; x<-copy(e3$m[mem18>def$k18]); x[,event:=as.integer(mem20<=def$k20)]
  b<-coef(e3$fit_new(x)); terms<-c(age5=5*b[['age18']],eduP=b[['edu_groupPrimary school']],
    eduM=b[['edu_groupMiddle school']],eduH=b[['edu_groupHigh school or above']],dep=b[['dep']],adl=b[['adl18']])
  pts<-round(10*terms/max(abs(terms))); edu<-as.integer(x$edu_group)
  old<-2*pmax(floor((x$age18-45)/5),0)+c(0,-2,-7,-10)[edu]+3*x$dep
  new<-pts[['age5']]*pmax(floor((x$age18-45)/5),0)+c(0,pts[['eduP']],pts[['eduM']],pts[['eduH']])[edu]+pts[['dep']]*x$dep
  add('F03','CHARLS 2018→2020',def$lab,'old integer weights',x$event,old,fixed=5)
  add('F03','CHARLS 2018→2020',def$lab,'original analysis reestimated integer weights',x$event,new,fixed=5)
}
all<-rbindlist(rows,fill=TRUE)
stopifnot(nrow(all)==53,all$TP+all$FP==all$cases_detected,all$TP+all$FN==all$events,
          all$TP+all$FP+all$TN+all$FN==all$n)
# Verify the untouched operating quantities against the old exported results.
files<-c(F14='F14_operating.csv',F22='F22_operating.csv',F15='F15_screening.csv',F16='F16_screening.csv',F19='F19_screening.csv')
metrics<-c('actual_screening_fraction','sensitivity','specificity','positive_predictive_value','negative_predictive_value','fraction_already_low_scoring_at_baseline_among_detected_cases','per_detected_1_incident_case_number_needed_to_screen')
for(task in names(files)) {
  path<-file.path(out,files[[task]]); backup<-file.path(archive,files[[task]])
  if(!file.exists(backup)) stopifnot(file.copy(path,backup))
  prior<-fread(backup,encoding='UTF-8'); now<-copy(all[source==task])
  stopifnot(nrow(prior)==nrow(now),all(abs(prior$sensitivity-round(now$sensitivity,1))<1e-8))
  ppv<-if('event_rate_among_screen_positives'%in%names(prior)) prior$event_rate_among_screen_positives else prior$positive_predictive_value
  stopifnot(all(abs(ppv-round(now$positive_predictive_value,1))<1e-8))
  if('cases_detected'%in%names(prior)) stopifnot(identical(as.integer(prior$cases_detected),as.integer(now$cases_detected)))
  now[,old_table_specificity:=prior$specificity]
  now[,specificity_correction_difference:=round(round(specificity,1)-old_table_specificity,1)]
  for(col in intersect(metrics,names(now))) set(now,j=col,value=round(now[[col]],if(col=='actual_screening_fraction')2 else 1))
  fwrite(now,path,bom=TRUE)
}
fwrite(all,file.path(out,'F24_complete_screening_metric_audit.csv'),bom=TRUE)
fwrite(all[source=='F03'],file.path(out,'F03_screening_corrected.csv'),bom=TRUE)
cat('PASS: 53 operating points; confusion-matrix identities; original sensitivity/PPV/selected counts reproduced.\n')
print(all[source=='F22',.(outcome_definition,model,n,cases_detected,actual_screening_fraction,sensitivity,specificity,positive_predictive_value,negative_predictive_value)])
print(all[source=='F19' & target_screening_fraction=='20%',.(cohort_or_window,n,cases_detected,actual_screening_fraction,sensitivity,specificity,positive_predictive_value,negative_predictive_value)])
# All modified source scripts must remain syntactically valid.
for(f in c('F03_headtohead.R','F14_bdsi_vs_ours_vs_tu.R','F15.R','F16_analysis.R','F19.R','F22b.R')) parse(file.path(src,f),encoding='UTF-8')
# Execute only patched screening blocks, with file writers disabled, then compare
# against independently constructed 2x2 metrics. This avoids rerunning AUC/bootstrap.
check_block <- function(file,start,end,e,object,task) {
  lines<-readLines(file.path(src,file),encoding='UTF-8',warn=FALSE)
  lo<-grep(start,lines,fixed=TRUE)[1]; hi<-grep(end,lines,fixed=TRUE)[1]
  stopifnot(!is.na(lo),!is.na(hi),hi>lo)
  e$fwrite<-function(...) invisible(NULL); e$say<-function(...) invisible(NULL); e$pt<-e$say
  eval(parse(text=lines[lo:(hi-1)]),e)
  z<-get(object,envir=e); ref<-all[source==task]
  if(any(abs(z$negative_predictive_value-round(ref$negative_predictive_value,1))>1e-8)) print(data.table(task=task,row=seq_len(nrow(z)),script=z$negative_predictive_value,independent=round(ref$negative_predictive_value,1)))
  stopifnot(nrow(z)==nrow(ref),all(abs(z$specificity-round(ref$specificity,1))<1e-8),
    all(abs(z$negative_predictive_value-round(ref$negative_predictive_value,1))<1e-8),
    all(abs(z$actual_screening_fraction-round(ref$actual_screening_fraction,2))<1e-8))
}
e14<-new.env(parent=globalenv()); e14$d<-d1
check_block('F14_bdsi_vs_ours_vs_tu.R','## ---------------- 7.','## ---------------- 8.',e14,'op','F14')
e22<-new.env(parent=globalenv()); e22$d1<-d1; e22$d2<-d2; e22$ARMS<-arms
check_block('F22b.R','say("\\n===== 4.','close(con)',e22,'op','F22')
e15$W1L<-'2011→2015'; e15$W2L<-'2015→2018'
e15$getc<-eval(quote(function(nm) { x<-get(nm,envir=environment(getc)); x[core(x)] }),e15)
check_block('F15.R','## ---------------- 8.','## ---------------- 9.',e15,'scrT','F15')
e16<-prefix_env('F16_analysis.R','say("cutoff verification')
check_block('F16_analysis.R','scr <- rbindlist','strat_auc <-',e16,'scr','F16')
e19<-prefix_env('F19.R','say("threshold: ')
check_block('F19.R','scr <- rbindlist','say("\\n========== 5.',e19,'scr','F19')
cat('PASS: six patched R scripts parsed; five patched screening blocks match the independent 2x2 calculations.\n')
cat('No AUC/bootstrap rerun; F03 original refit reproduced only to recover its integer scores.\n')
