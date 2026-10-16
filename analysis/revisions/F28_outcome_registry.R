# F28: lock outcome hierarchy and create a source-traceable numeric registry.
# Rscript F28_outcome_registry.R <third-round directory>
invisible(Sys.setlocale('LC_CTYPE','Chinese_China.utf8'))
.libPaths(c('D:/R/library',.libPaths()))
suppressPackageStartupMessages({library(data.table); library(pROC)})
root<-commandArgs(trailingOnly=TRUE)[1]; project<-dirname(root)
out<-file.path(root,'_核查','输出'); tmp<-'C:/y03tmp/tmp'
source_a12<-file.path(project,'output/11_外部验证/A12_主结果_三队列_S0_S1_S2.csv')
perf<-fread(source_a12)
frozen<-readRDS(file.path(tmp,'a07_memory_v3.rds'))
models<-frozen$model_u
windows<-c(HRS='2010至2014',ELSA='2008至2012',SHARE='2011至2015')
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
    outcome='随访记忆分达到或低于阈值',cut_int=threshold,cut_raw=rawcut,baseline_rule=baseline,
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
    priority<-if(sc=='S0')'V2既定外部主分析' else 'V2既定阈值敏感性分析'
    rawcut<-unique(d$cut); stopifnot(length(rawcut)==1)
    for(j in 1:3) {
      kind<-if(j==1)'连续模型' else '整数评分'
      pop<-if(j==2)'全部评分可计算风险人群' else '连续模型完整风险人群'
      nn<-if(j==2)nrow(d) else nrow(m); ee<-if(j==2)sum(d$event) else sum(m$event)
      prefix<-c('auc_model','auc_score_full','auc_score_same')[j]
      addrow(co,windows[[co]],sc,priority,row$cut_int,rawcut,'基线记忆分>阈值',kind,pop,nn,ee,
        row[[prefix]],row[[paste0(prefix,'_low')]],row[[paste0(prefix,'_high')]],'A12既有DeLong区间',source_a12,if(j==1)ca else NULL)
    }
    calibration[[length(calibration)+1]]<-data.table(cohort=co,scenario=sc,population='基线>阈值',n=nrow(m),
      offset_intercept=ca[['offset_intercept']],joint_intercept=ca[['joint_intercept']],slope=ca[['slope']],Brier=ca[['brier']],source='冻结模型预测值核对A12')
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
    for(kind in c('整数评分','连续模型')) {
      predlabel<-if(kind=='整数评分')'三项评分卡' else '连续模型'
      rr<-oldtable[队列==co & 口径=='U_全人口' & 预测子==predlabel]
      b<-if(kind=='连续模型')a[!is.na(p)] else a[!is.na(score)]
      av<-auc_val(b$event,if(kind=='连续模型')b$p else b$score)
      stopifnot(nrow(rr)==1,nrow(b)==rr$n,sum(b$event)==rr$事件,abs(av-rr$AUC)<0.000051)
      ca<-if(kind=='连续模型')calc_cal(b$event,b$p) else NULL
      addrow(co,windows[[co]],paste0('U_',sc),'事后全人群随访低分状态分析',cutpoint,NA_real_,
        '不要求基线记忆可测或高于阈值',kind,if(kind=='连续模型')'连续模型变量完整且随访可测' else '评分变量完整且随访可测',
        nrow(b),sum(b$event),rr$AUC,rr$AUC低,rr$AUC高,'F16/F19既有DeLong区间（四位小数）',source_csv,ca)
      if(!is.null(ca)) calibration[[length(calibration)+1]]<-data.table(cohort=co,scenario=paste0('U_',sc),population='不限制基线记忆状态',n=nrow(b),
        offset_intercept=ca[['offset_intercept']],joint_intercept=ca[['joint_intercept']],slope=ca[['slope']],Brier=ca[['brier']],source='冻结预测值补齐统一offset截距定义')
    }
  }
}
dev<-as.data.table(frozen$data)
for(kind in c('连续模型','整数评分')) {
  pred<-if(kind=='连续模型')dev$p_hat else dev$score
  r<-roc(dev$incident_2015,pred,direction='<',quiet=TRUE);ci<-as.numeric(ci.auc(r,method='delong'))
  ca<-if(kind=='连续模型')calc_cal(dev$incident_2015,pred) else NULL
  addrow('CHARLS','2011至2015','DEV','开发样本表观性能',4,NA_real_,'五波共同题集完整且基线记忆分>4',kind,'冻结开发样本',
    nrow(dev),sum(dev$incident_2015),as.numeric(auc(r)),ci[1],ci[3],'保存预测值DeLong复核；非bootstrap校正区间','C:/y03tmp/tmp/a07_memory_v3.rds',ca)
}
registry<-rbindlist(rows);registry[,record_id:=sprintf('F28-%02d',.I)];setcolorder(registry,c('record_id',setdiff(names(registry),'record_id')))
fwrite(registry,file.path(out,'F28_结局与数字底表.csv'),bom=TRUE)
fwrite(rbindlist(checks),file.path(out,'F28_A12独立核对.csv'),bom=TRUE)
fwrite(rbindlist(calibration),file.path(out,'F28_校准截距定义核对.csv'),bom=TRUE)
stopifnot(sum(registry[scenario=='S0' & model=='整数评分' & population=='全部评分可计算风险人群',n])==47655,
          sum(registry[scenario=='S0' & model=='连续模型',n])==47571)
writeLines(c(capture.output(sessionInfo()),'All numerical, eligibility, ID-set and calibration-definition checks passed.'),file.path(out,'F28_sessionInfo.txt'))
cat('Registry records:',nrow(registry),'\n')
print(registry[scenario%in%c('S0','S1','S2') & population!='连续模型完整风险人群'|scenario%in%c('S0','S1','S2') & model=='连续模型',.(cohort,scenario,model,n,events,AUC)])
print(rbindlist(calibration)[scenario%in%c('S0','S2')])
