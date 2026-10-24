# F26: audit BDSI on frozen F22 participants; never refit the study score.
# Rscript F26_bdsi_audit.R <third-round directory>
invisible(Sys.setlocale('LC_CTYPE','Chinese_China.utf8'))
.libPaths(c('D:/R/library',.libPaths()))
suppressPackageStartupMessages({library(data.table); library(haven); library(pROC)})
root <- commandArgs(trailingOnly=TRUE)[1]
stopifnot(!is.na(root), dir.exists(root))
out <- file.path(root,'_audit','output')
src <- file.path(root,'_audit','scripts')
project <- dirname(root)
tmp <- 'C:/y03tmp/tmp'
norm_id <- function(x) { x<-as.character(x); ifelse(nchar(x)==11,paste0(substr(x,1,9),sprintf('%03d',as.integer(substr(x,10,11)))),x) }
clean <- function(x, valid) { x<-as.numeric(x); x[!x %in% valid]<-NA_real_; x }
# Three-valued OR: one positive establishes positivity; two known negatives establish negativity.
or_known <- function(a,b) ifelse((!is.na(a)&a==1)|(!is.na(b)&b==1),1L,ifelse(!is.na(a)&!is.na(b),0L,NA_integer_))
stopifnot(identical(or_known(c(1,0,0,NA),c(NA,0,NA,NA)),c(1L,0L,NA_integer_,NA_integer_)))
d <- as.data.table(readRDS(file.path(tmp,'F22_d.rds')))
stopifnot(!anyDuplicated(d$ID))
vars <- c('ID','db019','db020','dc012',paste0('da020s',1:4),'da007_11_','da007_w2_2_11_')
r <- as.data.table(read_dta(file.path(tmp,'w3_hs.dta'),col_select=all_of(vars)))
labels <- rbindlist(lapply(setdiff(vars,'ID'),function(v) data.table(variable=v,label=attr(r[[v]],'label'),values=paste(names(attr(r[[v]],'labels')),collapse=' | '))))
fwrite(labels,file.path(out,'F26_raw_field_labels.csv'),bom=TRUE)
r[,ID:=as.character(ID)]
h <- as.data.table(read_dta(file.path(project,'data/raw/harmonized_d/H_CHARLS_D_Data.dta'),col_select=all_of(c('ID','r3psyche','raeduc_c'))))
h[,ID:=norm_id(ID)]
stopifnot(!anyDuplicated(r$ID),!anyDuplicated(h$ID),all(d$ID %in% r$ID))
d <- merge(d,r,by='ID',all.x=TRUE,sort=FALSE)
d <- merge(d,h,by='ID',all.x=TRUE,sort=FALSE)
stopifnot(nrow(d)==7839)
for(v in c('db019','db020','dc012')) set(d,j=v,value=clean(d[[v]],1:4))
d[, func_new := or_known(as.integer(db019>=3),as.integer(db020>=3))]
d[, effort := as.integer(dc012>=3)]
# Reconstruct diagnosis only to identify questionnaire structural skips.
# Harmonized r3psyche has missing .x for some new entrants: use their direct answer as fallback.
d[, psych := clean(r3psyche,0:1)]
direct <- clean(d$da007_11_,1:2); update <- clean(d$da007_w2_2_11_,1:2)
idx_direct <- which(is.na(d$psych)&!is.na(direct))
d[idx_direct,psych:=as.integer(direct[idx_direct]==1)]
# Explicit current-wave affirmative answer takes precedence over a recorded negative.
d[(!is.na(direct)&direct==1)|(!is.na(update)&update==1),psych:=1L]
d[is.na(psych)&!is.na(update)&update==2,psych:=0L]
tx <- sapply(paste0('da020s',1:4),function(v) clean(d[[v]],1:4))
any_tx <- rowSums(!is.na(tx))>0
anti_yes <- rowSums(tx==2,na.rm=TRUE)>0
d[, anti := ifelse(anti_yes,1L,ifelse(any_tx,0L,ifelse(!is.na(psych)&psych==0,0L,NA_integer_)))]
d[, anti_status := fifelse(anti_yes,'explicit antidepressant selection',fifelse(any_tx,'responded DA020 but did not select antidepressants',fifelse(!is.na(psych)&psych==0,'no psychiatric problem coded by skip logic 0','unknown medication status')))]
d[, dep_new := or_known(effort,anti)]
d[, age_new := pmax(0,pmin(age15,79)-65)] # outside 65-79: audit extrapolation ONLY
d[, score_old := score_bdsi]
d[, score_age := score_old-b_age+age_new]
d[, score_map := age_new+b_edu+b_stroke+b_diab+b_bmi+10L*func_new+6L*dep_new]
d[, score_effort := age_new+b_edu+b_stroke+b_diab+b_bmi+10L*func_new+6L*effort]
stopifnot(all(d$score_old==d$b_age+d$b_edu+d$b_stroke+d$b_diab+d$b_bmi+d$b_func+d$b_dep),
          all(d[age15>=65 & age15<=79,score_age==score_old-1]),
          identical(pmax(0,pmin(c(65,66,69,79),79)-65),c(0,1,4,14)))
fwrite(d[,.(n=.N),by=.(anti_status,psych,anti)],file.path(out,'F26_medication_skip_logic_audit.csv'),bom=TRUE)
# No subject-level data are written to the report directory.
saveRDS(d,file.path(tmp,'F26_bdsi_audit.rds'))

aucs<-list(); comparisons<-list(); screens<-list(); flows<-list(); changes<-list()
ratio<-function(a,b) if(b==0) NA_real_ else a/b
metrics<-function(y,w) {
  tp<-sum(w*y); fp<-sum(w*(1-y)); tn<-sum((1-w)*(1-y)); fn<-sum((1-w)*y)
  stopifnot(abs(tp+fp+tn+fn-length(y))<1e-8)
  data.table(n=length(y),events=sum(y),selected=sum(w),screen_pct=100*mean(w),TP=tp,FP=fp,TN=tn,FN=fn,
    sensitivity=100*ratio(tp,tp+fn),specificity=100*ratio(tn,tn+fp),PPV=100*ratio(tp,tp+fp),NPV=100*ratio(tn,tn+fn))
}
labels_model<-c(score_us='our three item score',score_old='BDSI old mapping',score_age='BDSI age correction only',score_map='BDSI item correction mapping',score_effort='BDSI effort item alone as proxy')
for(def in c('T1','T2')) {
  dd<-if(def=='T1') copy(d[mem15>=cutoff]) else copy(d[mem15>4])
  dd[,event:=if(def=='T1') as.integer(mem18<cutoff) else as.integer(mem18<=4)]
  stopifnot(nrow(dd)==if(def=='T1')6436L else 5872L)
  for(scope in c('full sample extrapolation audit','65 to 79 years old')) {
    z<-if(scope=='65 to 79 years old') dd[age15>=65 & age15<=79] else copy(dd)
    flows[[length(flows)+1]]<-data.table(outcome=def,scope=scope,n_original=nrow(z),events_original=sum(z$event),
      age_under65=sum(z$age15<65),age_over79=sum(z$age15>79),function_unknown=sum(is.na(z$func_new)),
      anti_unknown=sum(is.na(z$anti)),depression_unknown=sum(is.na(z$dep_new)),
      mapping_complete=sum(!is.na(z$score_map)),events_complete=sum(z$event[!is.na(z$score_map)]))
    zz<-z[!is.na(score_map)]
    changes[[length(changes)+1]]<-data.table(outcome=def,scope=scope,n=nrow(zz),
      function_old_positive=sum(zz$b_func==10),function_new_positive=sum(zz$func_new==1),
      function_changed=sum(zz$b_func!=10*zz$func_new),depression_old_positive=sum(zz$b_dep==6),
      depression_new_positive=sum(zz$dep_new==1),depression_changed=sum(zz$b_dep!=6*zz$dep_new),
      anti_positive=sum(zz$anti==1,na.rm=TRUE),anti_additional=sum(zz$dep_new==1 & zz$effort==0,na.rm=TRUE))
    for(stage in c('historical complete sample_isolate age correction','paired sample with scorable items')) {
      a<-if(stage=='paired sample with scorable items') zz else z
      models<-if(stage=='paired sample with scorable items') c('score_us','score_old','score_age','score_map') else c('score_us','score_old','score_age')
      rs<-list()
      for(v in models) {
        rs[[v]]<-roc(a$event,a[[v]],direction='<',levels=c(0,1),quiet=TRUE)
        ci<-as.numeric(ci.auc(rs[[v]],method='delong'))
        aucs[[length(aucs)+1]]<-data.table(outcome=def,scope=scope,stage=stage,model=labels_model[[v]],n=nrow(a),events=sum(a$event),AUC=as.numeric(auc(rs[[v]])),CI_low=ci[1],CI_high=ci[3])
        for(rule in c('target 30%_all ties','fixed quota 30%_expected value under random tie breaking')) {
          k<-floor(.30*nrow(a)); th<-sort(a[[v]],decreasing=TRUE)[k]
          w<-as.numeric(a[[v]]>=th); tie_prob<-1
          if(rule=='fixed quota 30%_expected value under random tie breaking') {
            tie_prob<-(k-sum(a[[v]]>th))/sum(a[[v]]==th)
            w<-as.numeric(a[[v]]>th)+as.numeric(a[[v]]==th)*tie_prob
            stopifnot(abs(sum(w)-k)<1e-8)
          }
          screens[[length(screens)+1]]<-cbind(data.table(outcome=def,scope=scope,stage=stage,model=labels_model[[v]],rule=rule,threshold=th,tie_probability=tie_prob),metrics(a$event,w))
        }
        if(v=='score_map' & scope=='65 to 79 years old') screens[[length(screens)+1]]<-cbind(data.table(outcome=def,scope=scope,stage=stage,model=labels_model[[v]],rule='original publication 22 points_exploratory cross outcome application',threshold=22,tie_probability=1),metrics(a$event,as.numeric(a[[v]]>=22)))
      }
      for(v in setdiff(models,'score_us')) {
        tst<-roc.test(rs$score_us,rs[[v]],paired=TRUE,method='delong',conf.int=TRUE,conf.level=.95)
        comparisons[[length(comparisons)+1]]<-data.table(outcome=def,scope=scope,stage=stage,comparator=labels_model[[v]],n=nrow(a),events=sum(a$event),
          AUC_us=as.numeric(auc(rs$score_us)),AUC_comparator=as.numeric(auc(rs[[v]])),delta=as.numeric(auc(rs$score_us)-auc(rs[[v]])),
          delta_low=tst$conf.int[1],delta_high=tst$conf.int[2],P=tst$p.value)
      }
      if(scope=='65 to 79 years old') stopifnot(abs(as.numeric(auc(rs$score_old)-auc(rs$score_age)))<1e-12)
    }
  }
}
aa<-rbindlist(aucs); cc<-rbindlist(comparisons); ss<-rbindlist(screens)
fwrite(aa,file.path(out,'F26_BDSI_AUC.csv'),bom=TRUE)
fwrite(cc,file.path(out,'F26_BDSI_paired_comparison.csv'),bom=TRUE)
fwrite(ss,file.path(out,'F26_BDSI_screening_metrics.csv'),bom=TRUE)
fwrite(rbindlist(flows),file.path(out,'F26_BDSI_sample_flow.csv'),bom=TRUE)
fwrite(rbindlist(changes),file.path(out,'F26_BDSI_item_changes.csv'),bom=TRUE)
# Independent AUC verification by Mann-Whitney ranks, including ties.
rankauc<-function(y,s) { n1<-sum(y); n0<-length(y)-n1; (sum(rank(s,ties.method='average')[y==1])-n1*(n1+1)/2)/(n1*n0) }
for(def in c('T1','T2')) {
  z<-if(def=='T1') d[mem15>=cutoff & !is.na(score_map)] else d[mem15>4 & !is.na(score_map)]
  yy<-if(def=='T1') as.integer(z$mem18<z$cutoff) else as.integer(z$mem18<=4)
  check<-aa[outcome==def & scope=='full sample extrapolation audit' & stage=='paired sample with scorable items' & model==labels_model[['score_map']],AUC]
  stopifnot(abs(rankauc(yy,z$score_map)-check)<1e-12)
}
cat('AUC rows:',nrow(aa),'Paired comparisons:',nrow(cc),'Screening rows:',nrow(ss),'\n')
print(rbindlist(flows))
print(cc[stage=='paired sample with scorable items' & comparator=='BDSI item correction mapping'])
print(ss[stage=='paired sample with scorable items' & scope=='65 to 79 years old' & model %in% c('our three item score','BDSI item correction mapping')])
cat('All arithmetic checks passed. Frozen F22 input and study coefficients unchanged.\n')
writeLines(c(capture.output(sessionInfo()),'All arithmetic checks passed.'),file.path(out,'F26_sessionInfo.txt'))

# New matched-budget question: paired participant bootstrap, thresholds reselected
# without outcomes. This estimates sampling uncertainty in expected random-tie
# sensitivity, not the variance of executing one particular random tie allocation.
budget_sens <- function(y,s) {
  k<-floor(.30*length(y)); th<-sort(s,decreasing=TRUE)[k]
  p<-(k-sum(s>th))/sum(s==th)
  sum(y*(as.numeric(s>th)+p*as.numeric(s==th)))/sum(y)
}
budget_results<-list()
set.seed(20261015)
for(def in c('T1','T2')) {
  z<-if(def=='T1') d[mem15>=cutoff & age15>=65 & age15<=79 & !is.na(score_map)] else d[mem15>4 & age15>=65 & age15<=79 & !is.na(score_map)]
  y<-if(def=='T1') as.integer(z$mem18<z$cutoff) else as.integer(z$mem18<=4)
  est<-budget_sens(y,z$score_us)-budget_sens(y,z$score_map)
  reps<-replicate(2000,{ii<-sample.int(nrow(z),replace=TRUE); budget_sens(y[ii],z$score_us[ii])-budget_sens(y[ii],z$score_map[ii])})
  ci<-quantile(reps,c(.025,.975))
  budget_results[[def]]<-data.table(outcome=def,scope='65 to 79 years old',n=nrow(z),events=sum(y),budget=floor(.3*nrow(z)),
    sensitivity_us=100*budget_sens(y,z$score_us),sensitivity_BDSI=100*budget_sens(y,z$score_map),delta_pp=100*est,
    CI_low_pp=100*ci[1],CI_high_pp=100*ci[2],bootstrap_reps=length(reps),seed=20261015)
}
bb<-rbindlist(budget_results)
fwrite(bb,file.path(out,'F26_BDSI_sensitivity_comparison_at_equal_quota.csv'),bom=TRUE)
print(bb)
# Verify cached functional proxy against raw difficulty items where both observed.
raw_difficulty<-or_known(as.integer(d$db019>=2),as.integer(d$db020>=2))
qc<-data.table(check=c('ID duplicate','original functional proxy differs from raw difficulty item_both known','antidepressants and none of above both selected','subtract within age range 1 point identity failed'),
  count=c(anyDuplicated(d$ID),sum(raw_difficulty!=d$iadl2_15,na.rm=TRUE),sum(anti_yes & tx[,4]==4,na.rm=TRUE),sum(d[age15>=65 & age15<=79,score_age!=score_old-1])))
fwrite(qc,file.path(out,'F26_mapping_quality_audit.csv'),bom=TRUE)
print(qc)
