# F28: lock outcome hierarchy and create a source-traceable numeric registry.
# Rscript F28_outcome_registry.R <third-round directory>
invisible(Sys.setlocale('LC_CTYPE','Chinese_China.utf8'))
.libPaths(c('D:/R/library',.libPaths()))
suppressPackageStartupMessages({library(data.table); library(pROC)})
root<-commandArgs(trailingOnly=TRUE)[1]; project<-dirname(root)
out<-file.path(root,'_audit','output'); tmp<-'C:/y03tmp/tmp'
source_a12<-file.path(project,'output/11_external_validation/A12_primary_results_three_cohorts_S0_S1_S2.csv')
perf<-fread(source_a12)
frozen<-readRDS(file.path(tmp,'a07_memory_v3.rds'))
models<-frozen$model_u
windows<-c(HRS='2010 to 2014',ELSA='2008 to 2012',SHARE='2011 to 2015')
defs<-c(S0='baseline_all_mean_sd',S1='complete_fu_mean_sd',S2='fixed_charls_le4')
rows<-list(); checks<-list(); calibration<-list()
predict_frozen<-function(d) {
  nd<-data.frame(age=d$age,edu_group=factor(as.character(d$edu_group),levels=models$xlevels$edu_group),dep=d$dep,adl=d$adl)
  as.numeric(predict(models,newdata=nd,type='response'))
}
calc_cal<-function(y,p) {
  lp<-qlogis(p)
  joint<-glm(y~lp,family=binomial())
  offset_fit<-glm(y~offset(lp),family=binomial())
  c(offset_intercept=unname(coef(offset_fit)[1]),joint_intercept=unname(coef(joint)[1]),slope=unname(coef(joint)[2]),brier=mean((p-y)^2))
}
auc_val<-function(y,s) as.numeric(auc(roc(y,s,direction='<',levels=c(0,1),quiet=TRUE)))
addrow<-function(co,win,scenario,priority,threshold,rawcut,baseline,kind,population,n,events,a,lo,hi,ci_source,source,cal=NULL) {
  rows[[length(rows)+1]]<<-data.table(cohort=co,window=win,scenario=scenario,priority=priority,
    outcome='follow up memory at or below threshold',cut_int=threshold,cut_raw=rawcut,baseline_rule=baseline,
    model=kind,population=population,n=n,events=events,event_rate=events/n,AUC=a,AUC_low=lo,AUC_high=hi,
    calibration_offset_intercept=if(is.null(cal))NA_real_ else cal[['offset_intercept']],
    calibration_slope=if(is.null(cal))NA_real_ else cal[['slope']],
    Brier=if(is.null(cal))NA_real_ else cal[['brier']],CI_source=ci_source,source=source)
}
for(co in names(windows)) {
  risks<-readRDS(file.path(tmp,paste0('ext_',tolower(co),'_sens.rds')))
  for(sc in names(defs)) {
    cd<-defs[[sc]]; d<-as.data.table(risks[[paste(co,cd,sep='__')]])
    d[,p:=predict_frozen(d)]
    row<-perf[cohort==co & cut_definition==cd]
    stopifnot(nrow(row)==1,!anyDuplicated(d$ID),all(d$mem_base>d$cut_int),
              all(d$event==as.integer(d$mem_fu<=d$cut_int)),nrow(d)==row$n_total,sum(d$event)==row$events_total)
    m<-d[!is.na(p)]
    stopifnot(nrow(m)==row$n_model,sum(m$event)==row$events_model)
    ac<-c(auc_val(m$event,m$p),auc_val(d$event,d$score),auc_val(m$event,m$score))
    expected<-as.numeric(row[,.(auc_model,auc_score_full,auc_score_same)])
    stopifnot(max(abs(ac-expected))<1e-10)
    ca<-calc_cal(m$event,m$p)
    stopifnot(abs(ca[['offset_intercept']]-row$cal_intercept)<1e-8,abs(ca[['slope']]-row$cal_slope)<1e-8,
              abs(ca[['brier']]-row$brier)<1e-10)
    priority<-if(sc=='S0')'V2 established external primary analysis' else 'V2 established threshold sensitivity'
    rawcut<-unique(d$cut); stopifnot(length(rawcut)==1)
    for(j in 1:3) {
      kind<-if(j==1)'continuous model' else 'integer score'
      pop<-if(j==2)'all at risk participants with calculable scores' else 'continuous model complete risk population'
      nn<-if(j==2)nrow(d) else nrow(m); ee<-if(j==2)sum(d$event) else sum(m$event)
      prefix<-c('auc_model','auc_score_full','auc_score_same')[j]
      addrow(co,windows[[co]],sc,priority,row$cut_int,rawcut,'baseline memory score>threshold',kind,pop,nn,ee,
        row[[prefix]],row[[paste0(prefix,'_low')]],row[[paste0(prefix,'_high')]],'A12 existing DeLong interval',source_a12,if(j==1)ca else NULL)
    }
    calibration[[length(calibration)+1]]<-data.table(cohort=co,scenario=sc,population='baseline>threshold',n=nrow(m),
      offset_intercept=ca[['offset_intercept']],joint_intercept=ca[['joint_intercept']],slope=ca[['slope']],Brier=ca[['brier']],source='verify frozen predictions A12')
    checks[[length(checks)+1]]<-data.table(cohort=co,scenario=sc,n_score=nrow(d),n_model=nrow(m),max_auc_difference=max(abs(ac-expected)),
      offset_intercept_difference=ca[['offset_intercept']]-row$cal_intercept,slope_difference=ca[['slope']]-row$cal_slope)
  }
  # Keep the full-population post-hoc analyses separate from incident risk sets.
  base<-as.data.table(readRDS(file.path(tmp,paste0('F16_base_',co,'.rds'))))[!is.na(age)&age>=45]
  base[,edu_group:=factor(fcase(edu_pts==0,'Less than primary',edu_pts== -2,'Primary school',edu_pts== -7,'Middle school',edu_pts== -10,'High school or above'),levels=models$xlevels$edu_group)]
  base[,score:=2*pmax(floor((age-45)/5),0)+edu_pts+3*dep]
  base[,p:=predict_frozen(base)]
  for(sc in c('S0','S2')) {
    cutpoint<-perf[cohort==co & cut_definition==defs[[sc]],cut_int]
    old<-as.data.table(risks[[paste(co,defs[[sc]],sep='__')]])
    same<-base[!is.na(mem_base)&mem_base>cutpoint&!is.na(mem_fu)]
    stopifnot(setequal(as.character(same$ID),as.character(old$ID)))
    a<-base[!is.na(mem_fu)];a[,event:=as.integer(mem_fu<=cutpoint)]
    source_csv<-file.path(out,if(sc=='S0')'F16_auc.csv' else 'F19_auc.csv')
    oldtable<-fread(source_csv)
    for(kind in c('integer score','continuous model')) {
      predlabel<-if(kind=='integer score')'three item scorecard' else 'continuous model'
      rr<-oldtable[cohort==co & definition=='U_full_population' & predictors==predlabel]
      b<-if(kind=='continuous model')a[!is.na(p)] else a[!is.na(score)]
      av<-auc_val(b$event,if(kind=='continuous model')b$p else b$score)
      stopifnot(nrow(rr)==1,nrow(b)==rr$n,sum(b$event)==rr$events,abs(av-rr$AUC)<0.000051)
      ca<-if(kind=='continuous model')calc_cal(b$event,b$p) else NULL
      addrow(co,windows[[co]],paste0('U_',sc),'post hoc broader population follow up low score analysis',cutpoint,NA_real_,
        'baseline memory availability or above threshold status not required',kind,if(kind=='continuous model')'continuous predictors complete and follow up measurable' else 'score variables complete and follow up measurable',
        nrow(b),sum(b$event),rr$AUC,rr$AUC_low,rr$AUC_high,'F16/F19_existing_DeLong_interval_four_decimal_places',source_csv,ca)
      if(!is.null(ca)) calibration[[length(calibration)+1]]<-data.table(cohort=co,scenario=paste0('U_',sc),population='no baseline memory restriction',n=nrow(b),
        offset_intercept=ca[['offset_intercept']],joint_intercept=ca[['joint_intercept']],slope=ca[['slope']],Brier=ca[['brier']],source='complete and harmonize frozen predictions offset intercept definition')
    }
  }
}
dev<-as.data.table(frozen$data)
for(kind in c('continuous model','integer score')) {
  pred<-if(kind=='continuous model')dev$p_hat else dev$score
  r<-roc(dev$incident_2015,pred,direction='<',quiet=TRUE);ci<-as.numeric(ci.auc(r,method='delong'))
  ca<-if(kind=='continuous model')calc_cal(dev$incident_2015,pred) else NULL
  addrow('CHARLS','2011 to 2015','DEV','apparent development sample performance',4,NA_real_,'complete five wave common items and baseline memory>4',kind,'frozen development sample',
    nrow(dev),sum(dev$incident_2015),as.numeric(auc(r)),ci[1],ci[3],'saved predictions DeLong audit; not bootstrap corrected interval','C:/y03tmp/tmp/a07_memory_v3.rds',ca)
}
registry<-rbindlist(rows);registry[,record_id:=sprintf('F28-%02d',.I)];setcolorder(registry,c('record_id',setdiff(names(registry),'record_id')))
fwrite(registry,file.path(out,'F28_outcome_and_numerical_source_table.csv'),bom=TRUE)
fwrite(rbindlist(checks),file.path(out,'F28_A12_independent_verification.csv'),bom=TRUE)
fwrite(rbindlist(calibration),file.path(out,'F28_calibration_intercept_definition_audit.csv'),bom=TRUE)
stopifnot(sum(registry[scenario=='S0' & model=='integer score' & population=='all at risk participants with calculable scores',n])==47655,
          sum(registry[scenario=='S0' & model=='continuous model',n])==47571)
writeLines(c(capture.output(sessionInfo()),'All numerical, eligibility, ID-set and calibration-definition checks passed.'),file.path(out,'F28_sessionInfo.txt'))
cat('Registry records:',nrow(registry),'\n')
print(registry[scenario%in%c('S0','S1','S2') & population!='continuous model complete risk population'|scenario%in%c('S0','S1','S2') & model=='continuous model',.(cohort,scenario,model,n,events,AUC)])
print(rbindlist(calibration)[scenario%in%c('S0','S2')])
