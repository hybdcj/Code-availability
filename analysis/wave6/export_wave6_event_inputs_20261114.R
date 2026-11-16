invisible(try(Sys.setlocale('LC_ALL','English_United States.utf8'), silent=TRUE))
out <- 'output/wave6_events_20261114'
dir.create(out,showWarnings=FALSE)
f <- list.files('output',pattern='^A14_final_model_2026-09-26.rds$',recursive=TRUE,full.names=TRUE)
stopifnot(length(f)==1)
m <- readRDS(f)$model
dir.create('output/12_final_model',recursive=TRUE,showWarnings=FALSE)
write.csv(data.frame(term=names(coef(m)),beta=as.numeric(coef(m))), 'output/12_final_model/A14_final_model_coefficients.csv',row.names=FALSE)
d <- m$data
stopifnot(nrow(d)==4441,sum(d$incident_2015)==658)
d$p_frozen <- as.numeric(predict(m,newdata=d,type='response'))
write.csv(d[,c('KEY','communityID','age','edu_group','dep','adl','incident_2015','p_frozen')],file.path(out,'development_internal.csv'),row.names=FALSE,fileEncoding='UTF-8')
r<-read.csv('output/wave6_qc_20261113/riskset2018_audit_internal.csv',colClasses=c(ID='character',communityID='character'))
r$dep<-as.integer(r$cesd10>=10)
write.csv(data.frame(ID=r$ID,p_R=as.numeric(predict(m,newdata=r,type='response'))),file.path(out,'prediction_reference_internal.csv'),row.names=FALSE)
