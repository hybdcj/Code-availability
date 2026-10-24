# F27: participant-disjoint evaluation of frozen scores after F26 mapping.
# Rscript F27_nonoverlap_validation.R <third-round directory>
invisible(Sys.setlocale('LC_CTYPE','Chinese_China.utf8'))
.libPaths(c('D:/R/library',.libPaths()))
suppressPackageStartupMessages({library(data.table); library(pROC)})
root<-commandArgs(trailingOnly=TRUE)[1]
stopifnot(!is.na(root),dir.exists(root))
out<-file.path(root,'_audit','output'); tmp<-'C:/y03tmp/tmp'
model<-readRDS(file.path(tmp,'a07_memory_v3.rds'))
training<-as.data.table(model$data)
parent<-as.data.table(readRDS(file.path(tmp,'cohort_change.rds'))$cohort)[year==2011]
d<-as.data.table(readRDS(file.path(tmp,'F26_bdsi_audit.rds')))
norm_id<-function(x) {
  x<-as.character(x)
  stopifnot(!anyNA(x),all(grepl('^[0-9]{11,12}$',x)))
  j<-nchar(x)==11
  x[j]<-paste0(substr(x[j],1,9),sprintf('%03d',as.integer(substr(x[j],10,11))))
  x
}
dev_ids<-norm_id(training$KEY); parent_ids<-norm_id(parent$KEY)
d[,ID:=norm_id(ID)]
stopifnot(length(dev_ids)==4441,length(parent_ids)==5672,
          !anyDuplicated(dev_ids),!anyDuplicated(parent_ids),!anyDuplicated(d$ID),
          all(dev_ids%in%parent_ids),nobs(model$model_u)==nrow(training))
# Verify the saved development IDs actually index the fitted model's rows.
mf<-model.frame(model$model_u)
for(v in names(mf)) stopifnot(isTRUE(all.equal(as.character(mf[[v]]),as.character(training[[v]]))))
stopifnot(max(abs(as.numeric(fitted(model$model_u))-training$p_hat))<1e-12)
stopifnot(identical(as.integer(model$score_table[['per_score_unit']]),c(2L,-2L,-7L,-10L,3L,0L)))
expected_score<-2*pmax(floor((d$age15-45)/5),0)-2*(d$edu15=='Primary school')-
  7*(d$edu15=='Middle school')-10*(d$edu15=='High school or above')+3*(d$cesd15>=10)
stopifnot(all(expected_score==d$score_us))
d[,`:=`(overlap_train=ID%in%dev_ids,overlap_parent=ID%in%parent_ids)]
# Existing F26 results are reused, not recomputed. Check identity of the populations.
f26_cmp<-fread(file.path(out,'F26_BDSI_paired_comparison.csv'))[
  scope=='65 to 79 years old' & stage=='paired sample with scorable items' & comparator=='BDSI item correction mapping']

ratio<-function(a,b) if(b==0) NA_real_ else a/b
cm<-function(y,w) {
  tp<-sum(w*y); fp<-sum(w*(1-y)); tn<-sum((1-w)*(1-y)); fn<-sum((1-w)*y)
  stopifnot(abs(tp+fp+tn+fn-length(y))<1e-9)
  data.table(n=length(y),events=sum(y),selected=sum(w),screen_pct=100*mean(w),TP=tp,FP=fp,TN=tn,FN=fn,
    sensitivity=100*ratio(tp,tp+fn),specificity=100*ratio(tn,tn+fp),PPV=100*ratio(tp,tp+fp),NPV=100*ratio(tn,tn+fn))
}
budget_weights<-function(s) {
  k<-floor(.3*length(s)); th<-sort(s,decreasing=TRUE)[k]
  p<-(k-sum(s>th))/sum(s==th)
  w<-as.numeric(s>th)+p*as.numeric(s==th)
  stopifnot(abs(sum(w)-k)<1e-8,p>=0,p<=1)
  list(w=w,threshold=th,tie_probability=p)
}
budget_sens<-function(y,s) sum(y*budget_weights(s)$w)/sum(y)
rankauc<-function(y,s) { n1<-sum(y); n0<-length(y)-n1; (sum(rank(s)[y==1])-n1*(n1+1)/2)/(n1*n0) }
flows<-list(); aucs<-list(); comparisons<-list(); screening<-list(); budgets<-list(); characteristics<-list(); qc<-list()
reps_all<-list(); saved_samples<-list()
for(def in c('T1','T2')) {
  z<-if(def=='T1') copy(d[mem15>=cutoff]) else copy(d[mem15>4])
  z[,event:=if(def=='T1')as.integer(mem18<cutoff) else as.integer(mem18<=4)]
  stopifnot(nrow(z)==if(def=='T1')6436 else 5872,
            sum(z$overlap_train)==if(def=='T1')3001 else 2901)
  mapped<-z[!is.na(score_map)]
  eligible<-mapped[age15>=65 & age15<=79]
  stopifnot(nrow(eligible)==f26_cmp[outcome==def,n],sum(eligible$event)==f26_cmp[outcome==def,events])
  for(set_name in c('F22_original_sample','F26_item_complete_sample','F26_complete_items_and_65_to_79_years_old')) {
    s<-switch(set_name,'F22_original_sample'=z,'F26_item_complete_sample'=mapped,'F26_complete_items_and_65_to_79_years_old'=eligible)
    flows[[length(flows)+1]]<-data.table(outcome=def,population=set_name,n=nrow(s),events=sum(s$event),
      overlap_train=sum(s$overlap_train),overlap_train_events=sum(s$event[s$overlap_train]),
      after_train_exclusion=sum(!s$overlap_train),after_train_exclusion_events=sum(s$event[!s$overlap_train]),
      overlap_parent=sum(s$overlap_parent),overlap_parent_events=sum(s$event[s$overlap_parent]),
      after_parent_exclusion=sum(!s$overlap_parent),after_parent_exclusion_events=sum(s$event[!s$overlap_parent]))
  }
  for(rule in c('exclude 4441 model fitting participants','exclude 5672 development source participants')) {
    keep<-if(rule=='exclude 4441 model fitting participants') !eligible$overlap_train else !eligible$overlap_parent
    for(group in c('retain','exclude')) {
      g<-eligible[if(group=='retain')keep else !keep]
      characteristics[[length(characteristics)+1]]<-data.table(outcome=def,exclusion=rule,group=group,n=nrow(g),events=sum(g$event),
        event_pct=100*mean(g$event),mean_age=mean(g$age15),sd_age=sd(g$age15),female_pct=100*mean(g$female15),
        less_primary_pct=100*mean(g$edu15=='Less than primary'),cesd10_ge10_pct=100*mean(g$cesd15>=10),
        baseline_memory_mean=mean(g$mem15),baseline_memory_sd=sd(g$mem15),score_us_mean=mean(g$score_us),score_bdsi_mean=mean(g$score_map))
    }
    a<-eligible[keep]
    stopifnot(nrow(a)>0,length(unique(a$event))==2,!any(a$ID%in%dev_ids))
    if(rule=='exclude 5672 development source participants') stopifnot(!any(a$ID%in%parent_ids))
    rs<-list()
    for(v in c('score_us','score_map')) {
      label<-if(v=='score_us')'our three item score' else 'BDSI item correction mapping'
      rr<-roc(a$event,a[[v]],direction='<',levels=c(0,1),quiet=TRUE); rs[[v]]<-rr
      ci<-as.numeric(ci.auc(rr,method='delong'))
      stopifnot(abs(as.numeric(auc(rr))-rankauc(a$event,a[[v]]))<1e-12)
      aucs[[length(aucs)+1]]<-data.table(outcome=def,exclusion=rule,scope='65 to 79 years old',model=label,n=nrow(a),events=sum(a$event),AUC=as.numeric(auc(rr)),CI_low=ci[1],CI_high=ci[3])
      bw<-budget_weights(a[[v]])
      rules<-c('target 30%_all ties','fixed quota 30%_expected value under random tie breaking','retain F26 cutoff')
      if(v=='score_map') rules<-c(rules,'original publication 22 points_exploratory cross outcome application')
      for(sr in rules) {
        th<-bw$threshold; p<-1
        if(sr=='retain F26 cutoff') th<-if(v=='score_map')19 else if(def=='T1')10 else 9
        if(sr=='original publication 22 points_exploratory cross outcome application') th<-22
        w<-as.numeric(a[[v]]>=th)
        if(sr=='fixed quota 30%_expected value under random tie breaking'){w<-bw$w;p<-bw$tie_probability}
        screening[[length(screening)+1]]<-cbind(data.table(outcome=def,exclusion=rule,model=label,rule=sr,threshold=th,tie_probability=p),cm(a$event,w))
      }
    }
    test<-roc.test(rs$score_us,rs$score_map,paired=TRUE,method='delong',conf.int=TRUE)
    comparisons[[length(comparisons)+1]]<-data.table(outcome=def,exclusion=rule,n=nrow(a),events=sum(a$event),
      AUC_us=as.numeric(auc(rs$score_us)),AUC_BDSI=as.numeric(auc(rs$score_map)),delta=as.numeric(auc(rs$score_us)-auc(rs$score_map)),
      CI_low=as.numeric(test$conf.int[1]),CI_high=as.numeric(test$conf.int[2]),P=test$p.value)
    seed<-202610160L+if(def=='T1')1L else 2L
    if(rule=='exclude 5672 development source participants') seed<-seed+10L
    set.seed(seed)
    reps<-replicate(2000,{ii<-sample.int(nrow(a),replace=TRUE); budget_sens(a$event[ii],a$score_us[ii])-budget_sens(a$event[ii],a$score_map[ii])})
    stopifnot(all(is.finite(reps)))
    ci<-quantile(reps,c(.025,.975)); sens_u<-budget_sens(a$event,a$score_us); sens_b<-budget_sens(a$event,a$score_map)
    budgets[[length(budgets)+1]]<-data.table(outcome=def,exclusion=rule,n=nrow(a),events=sum(a$event),budget=floor(.3*nrow(a)),
      sensitivity_us=100*sens_u,sensitivity_BDSI=100*sens_b,delta_pp=100*(sens_u-sens_b),CI_low_pp=100*ci[1],CI_high_pp=100*ci[2],bootstrap_reps=2000,seed=seed)
    key<-paste(def,rule,sep='_'); reps_all[[key]]<-reps
    saved_samples[[key]]<-a[,.(ID,event,score_us,score_map)]
    qc[[key]]<-data.table(outcome=def,exclusion=rule,n=nrow(a),training_overlap=sum(a$ID%in%dev_ids),parent_overlap=sum(a$ID%in%parent_ids),unique_IDs=uniqueN(a$ID),valid_bootstrap=sum(is.finite(reps)))
  }
}
results<-list('sample flow'=rbindlist(flows),'AUC'=rbindlist(aucs),'AUC paired comparison'=rbindlist(comparisons),
  'screening metrics'=rbindlist(screening),'sensitivity at equal quota'=rbindlist(budgets),'sample characteristics'=rbindlist(characteristics),'verify'=rbindlist(qc))
for(nm in names(results)) fwrite(results[[nm]],file.path(out,paste0('F27_',nm,'.csv')),bom=TRUE)
saveRDS(list(samples=saved_samples,sensitivity_bootstrap=reps_all),file.path(tmp,'F27_nonoverlap_audit.rds'))
writeLines(c(capture.output(sessionInfo()),'All ID, score, AUC rank-sum, confusion matrix, budget, and bootstrap checks passed.'),file.path(out,'F27_sessionInfo.txt'))
print(results[['sample flow']]);print(results[['AUC paired comparison']]);print(results[['sensitivity at equal quota']]);print(results[['verify']])
