invisible(Sys.setlocale('LC_CTYPE','Chinese_China.utf8'))
.libPaths(c('D:/R/library',.libPaths()))
suppressPackageStartupMessages({library(data.table);library(pROC);library(splines)})
P <- 'D:/BaiduSyncdisk/research/Xu_Zhonghao'
ROOT <- file.path(P,'16.manuscript_revision_round_three')
OUT <- file.path(ROOT,'_audit/output')
src <- file.path(ROOT,'_audit/scripts/F04_incremental.R')
# Evaluate only the existing data preparation and formula definitions.
# No bootstrap, plot, report-writing, or original output code is executed.
e <- new.env(parent=globalenv())
e$TMP <- 'C:/y03tmp/tmp';e$say <- function(...) invisible(NULL);e$pt <- function(...) invisible(NULL)
parsed <- parse(src,encoding='UTF-8')
for(x in parsed) if(is.call(x) && identical(x[[1]],as.name('<-')) &&
 is.symbol(x[[2]]) && as.character(x[[2]]) %in% c('BETA_OLD','LV')) eval(x,e)
lines <- readLines(src,encoding='UTF-8',warn=FALSE)
start <- grep('^w <- as.data.table',lines)
end <- grep('^## ---------------- apparent results',lines)-1L
stopifnot(length(start)==1,length(end)==1)
eval(parse(text=lines[start:end]),e)
main <- e$m;risk <- e$d_all;inc <- e$d
w <- e$w
paired <- merge(w[year==2018,.(ID=IDn,mem18=as.numeric(imrc)+as.numeric(dlrc))],
 w[year==2020,.(ID=IDn,mem20=as.numeric(imrc)+as.numeric(dlrc))],by='ID')
paired <- paired[!is.na(mem18)&!is.na(mem20)]
stopifnot(nrow(paired)==12571,nrow(main)==12331,nrow(risk)==10163,
 nrow(inc)==9178,sum(risk$event)==423,sum(inc$event)==397,
 e$cut18==3,e$cut20==6,!anyDuplicated(main$ID),!anyDuplicated(paired$ID))
a07 <- readRDS('C:/y03tmp/tmp/a07_memory_v3.rds')
dev_ids <- as.character(a07$data$KEY)
oldb <- coef(a07$model_u)
saved_coef <- oldb[c('(Intercept)','age','edu_groupPrimary school',
 'edu_groupMiddle school','edu_groupHigh school or above','dep','adl')]
stopifnot(max(abs(unname(e$BETA_OLD)-unname(saved_coef)))<1e-12)
prediction <- function(d) {
 b<-e$BETA_OLD;lv<-e$LV
 plogis(b[1]+b['age']*d$age18+b['dep']*d$dep+b['adl']*d$adl18+
 b['eduP']*(d$edu_group==lv[2])+b['eduM']*(d$edu_group==lv[3])+b['eduH']*(d$edu_group==lv[4]))
}
met <- function(y,p) {
 r<-roc(y,p,direction='<',quiet=TRUE);ci<-as.numeric(ci.auc(r));lp<-qlogis(p)
 joint<-glm(y~lp,family=binomial());offsetfit<-glm(y~offset(lp),family=binomial())
 c(n=length(y),events=sum(y),auc=ci[2],auc_low=ci[1],auc_high=ci[3],
 slope=unname(coef(joint)[2]),joint_intercept=unname(coef(joint)[1]),
 offset_intercept=unname(coef(offsetfit)[1]),brier=mean((p-y)^2))
}
out <- list();checks <- character()
add <- function(task,analysis,sample,outcome,model,y,p,d) {
 v<-met(y,p)
 out[[length(out)+1L]] <<- cbind(data.frame(task,analysis,sample,outcome,model,
   baseline_year=2018,outcome_year=2020,development_id_overlap=sum(as.character(d$ID)%in%dev_ids),
   role='post hoc later window exploration in same cohort; not independent external validation'),as.data.frame(as.list(v)))
 v
}
f03old<-fread(file.path(OUT,'F03_A_temporal_validation.csv'),encoding='UTF-8')
for(k in c(3,6)) {
 d<-copy(risk);d$event<-as.integer(d$mem20<=k);p<-prediction(d);lp<-qlogis(p)
 sc<-if(k==3) 'D1' else 'D2'
 v<-add('F03',sc,'2018 memory>3 and age≥45 and complete four variable data',paste0('2020 memory≤',k),
  'frozen continuous model',d$event,p,d)
 row<-if(k==3) 1L else 2L
 stopifnot(round(v['auc'],3)==f03old$AUC[row],round(v['joint_intercept'],3)==f03old$calibration_intercept[row])
 f<-glm(event~age18+edu_group+dep+adl18,data=d,family=binomial())
 add('F03',sc,'2018 memory>3 and age≥45 and complete four variable data',paste0('2020 memory≤',k),'four variable refit(apparent)',d$event,predict(f,type='response'),d)
 fi<-glm(event~offset(lp),data=d,family=binomial())
 add('F03',sc,'2018 memory>3 and age≥45 and complete four variable data',paste0('2020 memory≤',k),'intercept only recalibration(same sample evaluation)',d$event,
  plogis(lp+coef(fi)[1]),d)
}
apparent<-fread(file.path(OUT,'F04_apparent.csv'),encoding='UTF-8')
for(nm in names(e$FORMS)) {
 f<-glm(e$FORMS[[nm]],data=inc,family=binomial());p<-predict(f,type='response')
 v<-add('F04',nm,'complete common candidate variable sample; 2018 memory>3','2020 memory≤3',
  paste0('refit(apparent): ',e$LBL[[nm]]),inc$event,p,inc)
 row<-match(nm,apparent$model)
 stopifnot(round(v['auc'],4)==apparent$apparent_AUC[row],round(v['brier'],4)==apparent$Brier[row])
}
# F05's standardization uses all 12,571 paired observations, not the model subset.
mu18<-mean(paired$mem18);mu20<-mean(paired$mem20);sd18<-sd(paired$mem18);sd20<-sd(paired$mem20)
d<-copy(main)
d[,dz:=(mem20-mu20)/sd20-(mem18-mu18)/sd18]
d[,draw:=mem20-mem18]
f05old<-fread(file.path(OUT,'F05_outcome_definitions.csv'),encoding='UTF-8')
targets<-list(list('O1',TRUE,function(x)x$mem20<=3),list('O2',TRUE,function(x)x$mem20<=6),
 list('O6',TRUE,function(x)x$mem20<=3 & x$draw<=-5),
 list('O3',FALSE,function(x)x$dz<=-1),list('O4',FALSE,function(x)x$dz<=-.5),
 list('O3',TRUE,function(x)x$dz<=-1),list('O4',TRUE,function(x)x$dz<=-.5))
for(i in seq_along(targets)) {
 t<-targets[[i]];dd<-if(t[[2]]) copy(d[mem18>3]) else copy(d)
 dd$event<-as.integer(t[[3]](dd));p<-prediction(dd)
 sc<-paste(t[[1]],if(t[[2]]) 'risk set' else 'all complete cases')
 v<-add('F05',sc,if(t[[2]]) '2018 memory>3' else '≥45 years old and complete four variable data',f05old$definition[i],
   'frozen continuous model',dd$event,p,dd)
 f<-glm(event~age18+edu_group+dep+adl18,data=dd,family=binomial())
 v2<-add('F05',sc,if(t[[2]]) '2018 memory>3' else '≥45 years old and complete four variable data',f05old$definition[i],
   'four variable refit(apparent)',dd$event,predict(f,type='response'),dd)
 stopifnot(nrow(dd)==f05old$n[i],sum(dd$event)==f05old$events[i],
  round(v['auc'],4)==f05old$frozen_AUC[i],round(v2['auc'],4)==f05old$reestimate_AUC[i])
}
registry<-rbindlist(out,fill=TRUE)
fwrite(registry,file.path(OUT,'F31_available_performance_source_table.csv'),bom=TRUE)
flow<-data.frame(sample=c('memory measurable at both waves(age and predictor completeness restrictions not yet applied)',
 'age≥45 and four variables complete','baseline memory>3 risk set','common complete candidate risk set'),
 n=c(nrow(paired),nrow(main),nrow(risk),nrow(inc)),
 development_id_overlap=c(sum(as.character(paired$ID)%in%dev_ids),sum(as.character(main$ID)%in%dev_ids),
 sum(as.character(risk$ID)%in%dev_ids),sum(as.character(inc$ID)%in%dev_ids)))
write.csv(flow,file.path(OUT,'F31_sample_flow_and_overlap.csv'),row.names=FALSE,fileEncoding='UTF-8')
paired[,change:=mem20-mem18]
desc<-data.frame(n=nrow(paired),mean2018=mu18,sd2018=sd18,mean2020=mu20,sd2020=sd20,
 mean_change=mean(paired$change),correlation=cor(paired$mem18,paired$mem20),
 baseline_low=sum(paired$mem18<=3),both_low=sum(paired$mem18<=3&paired$mem20<=3),
 incident_events=sum(paired$mem18>3&paired$mem20<=3),
 cut_anchor_population='12571 paired sample(F05); 12331 complete cases(F03/F04), all exported≤3 and≤6')
write.csv(desc,file.path(OUT,'F31_paired_distribution_description.csv'),row.names=FALSE,fileEncoding='UTF-8')
# Audit the bootstrap estimand from existing code, without repeating bootstrap.
stopifnot(any(grepl('optimism <- app[nms] - colMeans(optAUC',lines,fixed=TRUE)),
 any(grepl('dA <- sweep(bootAUC',lines,fixed=TRUE)))
checks<-c('F03 sample sizes, thresholds, frozen AUC and joint calibration intercept: PASS',
 'F04 all 12 apparent AUC/Brier rows reproduced: PASS',
 'F05 all 7 definitions, sample sizes, event counts and frozen/refit AUC: PASS',
 'Frozen coefficients equal A07 saved coefficients: PASS',
 'No bootstrap rerun; F04 correction formula is simple-bootstrap evaluation, not Harrell optimism subtraction',
 capture.output(print(flow)),capture.output(print(registry[task=='F03'],digits=10)),
 capture.output(print(desc,digits=10)),capture.output(sessionInfo()))
writeLines(checks,file.path(OUT,'F31_verification_log.txt'),useBytes=TRUE)
cat(paste(head(checks,5),collapse='\n'),'\n');print(flow);print(registry[task=='F03'],digits=8)
