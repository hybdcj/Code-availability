invisible(Sys.setlocale('LC_CTYPE','Chinese_China.utf8'))
.libPaths(c('D:/R/library',.libPaths()))
suppressPackageStartupMessages({library(data.table);library(pROC)})
P<-'D:/BaiduSyncdisk/research/Xu_Zhonghao';ROOT<-file.path(P,'16.manuscript_revision_round_three')
OUT<-file.path(ROOT,'_audit/output')
src<-file.path(ROOT,'_audit/scripts/F15.R')
# Read original preparation and frozen formulas only; suppress legacy writes.
e<-new.env(parent=globalenv());e$TMP<-'C:/y03tmp/tmp';e$OUT<-OUT
e$say<-function(...)invisible(NULL);e$pt<-function(...)invisible(NULL);e$fwrite<-function(...)invisible(NULL)
for(x in parse(src,encoding='UTF-8')) if(is.call(x)&&identical(x[[1]],as.name('<-'))&&
 is.symbol(x[[2]])&&as.character(x[[2]])%in%c('B','LV','G3')) eval(x,e)
lines<-readLines(src,encoding='UTF-8',warn=FALSE)
start<-grep('^w <- as.data.table',lines);end<-grep('^# same inclusion criteria as',lines)-1L
stopifnot(length(start)==1,length(end)==1)
eval(parse(text=lines[start:end]),e)
a07<-readRDS('C:/y03tmp/tmp/a07_memory_v3.rds');ids<-as.character(a07$data$KEY)
tables<-list();flow<-list();dec<-list();screen<-list();ages<-list();sets<-list();audit<-list()
auc_met<-function(y,p) {
 stopifnot(!anyNA(y),!anyNA(p),all(y%in%c(0,1)))
 r<-roc(y,p,quiet=TRUE,direction='<',levels=c(0,1));ci<-as.numeric(ci.auc(r))
 c(auc=ci[2],auc_low=ci[1],auc_high=ci[3])
}
armset<-function(d) list(U=d,N=d[!is.na(low0)&low0==0],I=d[!is.na(low0)&low0==1],
 M=d[is.na(low0)],B=d[!is.na(low0)])
add_perf<-function(d,win,def,arm,model,pop) {
 v<-if(model=='integer score') d$score else d$phat
 a<-auc_met(d$event,v)
 ca<-c(slope=NA_real_,offset_intercept=NA_real_,joint_intercept=NA_real_,brier=NA_real_)
 if(model=='continuous model') {
   joint<-glm(event~lp,data=d,family=binomial());off<-glm(event~offset(lp),data=d,family=binomial())
   ca<-c(slope=unname(coef(joint)[2]),offset_intercept=unname(coef(off)[1]),
    joint_intercept=unname(coef(joint)[1]),brier=mean((d$phat-d$event)^2))
 }
 tables[[length(tables)+1L]]<<-cbind(data.frame(window=win,definition=def,arm=arm,
  model=model,population=pop,n=nrow(d),events=sum(d$event),event_rate=sum(d$event)/nrow(d),
  development_id_overlap=sum(as.character(d$ID)%in%ids),
  baseline_rule=switch(arm,U='no baseline memory restriction',N='baseline memory above low score boundary',
    I='baseline meets low score definition',M='missing baseline memory',B='baseline memory measurable'),
  outcome_rule=if(def=='fixed≤4') 'follow up memory≤4; N baseline>4, I baseline≤4' else
    'follow up memory<2011 education strata mean-SD; N baseline≥stratification boundaries, I baseline<stratification boundaries',
  role='CHARLS post hoc broader population extension; includes original development participants',ci_method='DeLong'),
  as.data.frame(as.list(c(a,ca))))
}
add_screen<-function(d,win,def,target,threshold) {
 pos<-d$score>=threshold;y<-d$event;low<-d$low0
 tp<-sum(pos&y==1);fp<-sum(pos&y==0);tn<-sum(!pos&y==0);fn<-sum(!pos&y==1)
 newtp<-sum(pos&y==1&!is.na(low)&low==0,na.rm=TRUE)
 oldtp<-sum(pos&y==1&!is.na(low)&low==1,na.rm=TRUE)
 unknowntp<-sum(pos&y==1&is.na(low))
 knownpos<-sum(pos&!is.na(low));lowpos<-sum(pos&!is.na(low)&low==1,na.rm=TRUE)
 stopifnot(tp+fp+tn+fn==nrow(d),newtp+oldtp+unknowntp==tp)
 screen[[length(screen)+1L]]<<-data.frame(window=win,definition=def,population='full score available sample U',
  target_fraction=target,threshold=threshold,tie_rule='score>=threshold; include all ties',
  n=nrow(d),events=sum(y),selected=sum(pos),actual_fraction=mean(pos),TP=tp,FP=fp,TN=tn,FN=fn,
  sensitivity=tp/(tp+fn),specificity=tn/(tn+fp),PPV=tp/(tp+fp),NPV=tn/(tn+fn),
  selected_baseline_known=knownpos,selected_baseline_low=lowpos,
  selected_baseline_low_fraction_all=lowpos/sum(pos),selected_baseline_low_fraction_known=lowpos/knownpos,
  TP_baseline_above=newtp,TP_baseline_low=oldtp,TP_baseline_unknown=unknowntp,
  total_scored_per_detected_crossing=if(newtp>0)nrow(d)/newtp else NA_real_,
  selected_per_detected_crossing=if(newtp>0)sum(pos)/newtp else NA_real_,
  threshold_status=if(is.na(target))'Youden selection in same sample; not a prespecified action threshold' else 'score cutoff corresponding to target fraction; actual fraction may exceed target')
}
for(win in c('2011-2015','2015-2018')) {
 raw<-copy(if(win=='2011-2015')e$w1 else e$w2)
 stopifnot(!anyDuplicated(raw$ID))
 nd<-data.frame(age=raw$age0,edu_group=factor(raw$edu0,levels=e$LV),
   dep=as.integer(raw$cesd0>=10),adl=raw$adl0)
 pred<-as.numeric(predict(a07$model_u,newdata=nd,type='response'))
 stopifnot(all(is.na(pred)==is.na(raw$phat)),max(abs(pred-raw$phat),na.rm=TRUE)<1e-12)
 for(def in c('fixed≤4','education strata')) {
  d<-copy(raw)
  if(def=='fixed≤4') d[,`:=`(event=as.integer(mem1<=4),low0=as.integer(mem0<=4))]
  score_ok<-!is.na(d$score)&!is.na(d$event)&!is.na(d$age0)&d$age0>=45
  model_ok<-score_ok&!is.na(d$lp)
  s<-d[score_ok];m<-d[model_ok];extra<-d[score_ok&!model_ok]
  stopifnot(all(is.na(extra$adl0)),setequal(s$ID,c(m$ID,extra$ID)),!anyNA(s$score),!anyNA(m$lp))
  key<-paste(win,def,sep=' | ');sets[[key]]<-list(score=s,model=m)
  flow[[length(flow)+1L]]<-data.frame(window=win,definition=def,matched_ids=nrow(d),
   followup_memory_observed=sum(!is.na(d$mem1)),score_sample=nrow(s),model_sample=nrow(m),
   restored_ADL_missing=nrow(extra),restored_events=sum(extra$event),
   score_development_overlap=sum(as.character(s$ID)%in%ids),model_development_overlap=sum(as.character(m$ID)%in%ids))
  ss<-armset(s);ms<-armset(m)
  for(grp in names(ss)) {
    add_perf(ss[[grp]],win,def,grp,'integer score','score specific complete sample')
    add_perf(ms[[grp]],win,def,grp,'integer score','common complete sample for continuous model')
    add_perf(ms[[grp]],win,def,grp,'continuous model','common complete sample for continuous model')
  }
  for(pop in c('score specific complete sample','common complete sample for continuous model')) {
    aa<-if(pop=='score specific complete sample')ss else ms
    narm<-vapply(aa,nrow,integer(1));earm<-vapply(aa,function(x)sum(x$event),numeric(1))
    stopifnot(narm['U']==sum(narm[c('N','I','M')]),earm['U']==sum(earm[c('N','I','M')]),
     narm['B']==sum(narm[c('N','I')]),earm['B']==sum(earm[c('N','I')]))
    dec[[length(dec)+1L]]<-data.frame(window=win,definition=def,population=pop,
      n_U=narm['U'],n_N=narm['N'],n_I=narm['I'],n_M=narm['M'],
      events_U=earm['U'],events_N=earm['N'],events_I=earm['I'],events_M=earm['M'],
      baseline_low_fraction_of_ALL_events=earm['I']/earm['U'],
      baseline_low_fraction_of_KNOWN_baseline_events=earm['I']/earm['B'],
      known_baseline_events=earm['B'],row.names=NULL)
    x<-aa$U
    for(band in c('45-59','60-69','70-79','80-89','90+')) {
      bounds<-switch(band,'45-59'=c(45,60),'60-69'=c(60,70),'70-79'=c(70,80),'80-89'=c(80,90),'90+'=c(90,Inf))
      ab<-x[age0>=bounds[1]&age0<bounds[2]]
      ages[[length(ages)+1L]]<-data.frame(window=win,definition=def,population=pop,band=band,
        n=nrow(ab),events=sum(ab$event),event_rate=if(nrow(ab)>0)mean(ab$event) else NA_real_)
    }
  }
  for(fr in c(.05,.10,.20,.30)) {
    th<-sort(s$score,decreasing=TRUE)[floor(fr*nrow(s))]
    add_screen(s,win,def,fr,th)
  }
  co<-coords(roc(s$event,s$score,quiet=TRUE,direction='<'),'best',best.method='youden',
     ret=c('threshold','specificity'))
  add_screen(s,win,def,NA_real_,co$threshold[which.max(co$specificity)])
 }
}
perf<-rbindlist(tables);ft<-rbindlist(flow);dc<-rbindlist(dec);sc<-rbindlist(screen);ag<-rbindlist(ages)
# Reproduce every old F15 AUC on its actual common-complete population.
old<-fread(file.path(OUT,'F15_auc.csv'),encoding='UTF-8')
old_cut<-fread(file.path(OUT,'F15_cutoffs.csv'),encoding='UTF-8')
stopifnot(max(abs(e$cut3$cutoff-old_cut$cutoff[match(e$cut3$g3,old_cut$g3)]))<1e-12)
for(i in seq_len(nrow(old))) {
  oo<-old[i];win<-if(grepl('2011',oo$window)) '2011-2015' else '2015-2018'
  arm_key<-substr(oo$definition,1,1);mdl<-if(oo$predictors=='three item scorecard')'integer score' else 'continuous model'
  rr<-perf[window==win&definition=='education strata'&arm==arm_key&model==mdl&population=='common complete sample for continuous model']
  stopifnot(nrow(rr)==1,rr$n==oo$n,rr$events==oo$events,round(rr$auc,4)==oo$AUC)
}
oldc<-fread(file.path(OUT,'F15c_auc_globalcut.csv'),encoding='UTF-8')
for(i in seq_len(nrow(oldc))) {
 oo<-oldc[i];arm_key<-substr(oo$definition,1,1)
 for(mdl in c('integer score','continuous model')) {
  rr<-perf[window=='2011-2015'&definition=='fixed≤4'&arm==arm_key&model==mdl&population=='common complete sample for continuous model']
  expected<-if(mdl=='integer score')oo$AUC_scorecard else oo$AUC_continuous
  stopifnot(nrow(rr)==1,rr$n==oo$n,rr$events==oo$events,round(rr$auc,4)==expected)
 }
}
stopifnot(ft[window=='2011-2015',all(score_sample==11869&model_sample==11753&restored_ADL_missing==116)])
fwrite(perf,file.path(ROOT,'F32_broader_population_performance_and_denominators_corrected_version.csv'),bom=TRUE)
fwrite(sc,file.path(ROOT,'F32_broader_population_screening_metrics_corrected_version.csv'),bom=TRUE)
fwrite(dc,file.path(OUT,'F32_denominator_for_event_composition.csv'),bom=TRUE)
fwrite(ft,file.path(OUT,'F32_sample_flow.csv'),bom=TRUE)
fwrite(ag,file.path(OUT,'F32_age_composition.csv'),bom=TRUE)
fwrite(e$cut3,file.path(OUT,'F32_retain_education_specific_thresholds.csv'),bom=TRUE)
saveRDS(sets,file.path(OUT,'F32_verification_sample.rds'))
for(pair in list(list(perf,file.path(ROOT,'F32_broader_population_performance_and_denominators_corrected_version.csv')),
 list(sc,file.path(ROOT,'F32_broader_population_screening_metrics_corrected_version.csv')))) {
 disk<-fread(pair[[2]],encoding='UTF-8')
 stopifnot(nrow(disk)==nrow(pair[[1]]),identical(names(disk),names(pair[[1]])))
 numeric_cols<-names(pair[[1]])[vapply(pair[[1]],is.numeric,logical(1))]
 for(nm in numeric_cols) {
  stopifnot(identical(is.na(disk[[nm]]),is.na(pair[[1]][[nm]])))
  delta<-abs(disk[[nm]]-pair[[1]][[nm]])
  if(any(!is.na(delta))) stopifnot(max(delta,na.rm=TRUE)<1e-10)
 }
}
checks<-c('Original F15 and F15c common-complete AUC/n/events reproduced: PASS',
 'Frozen predictions verified against saved A07 model: PASS',
 'Score population = common-complete population + ADL-missing restored participants: PASS',
 'U = N + I + M for participants and events; B = N + I: PASS',
 'TP+FP+TN+FN = n; TP split by baseline status sums to total: PASS',
 capture.output(print(ft)),capture.output(print(dc,digits=10)),capture.output(sessionInfo()))
writeLines(checks,file.path(OUT,'F32_verification_log.txt'),useBytes=TRUE)
cat(paste(head(checks,5),collapse='\n'),'\n');print(ft);print(dc,digits=8)
print(perf[arm=='U',.(window,definition,model,population,n,events,auc)],digits=8)
