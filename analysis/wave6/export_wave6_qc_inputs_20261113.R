find_one <- function(pattern) {
 f <- list.files('output',pattern=pattern,recursive=TRUE,full.names=TRUE)
 if(length(f)!=1L) stop(paste('Expected exactly one input:',pattern))
 f
}
invisible(try(Sys.setlocale('LC_ALL','English_United States.utf8'), silent=TRUE))
w<-readRDS(find_one('^waves_cognition_2026-09-13.rds$'))
dir.create('output/wave6_qc_20261113',showWarnings=FALSE)
write.csv(w[w$year==2018,c('IDn','year','age','female','rural','communityID','edu_group','cesd10','adl','imrc','dlrc')],'output/wave6_qc_20261113/baseline2018_internal.csv',row.names=FALSE,fileEncoding='UTF-8')
m<-readRDS(find_one('^A14_final_model_2026-09-26.rds$'))
stopifnot(nrow(m$model$data)==4441)
write.csv(data.frame(KEY=m$model$data$KEY),'output/wave6_qc_20261113/development_ids_internal.csv',row.names=FALSE)
co<-readRDS(find_one('^cohort_change_2026-09-13.rds$'))$cohort
stopifnot(length(unique(co$KEY))==5672)
write.csv(data.frame(KEY=unique(co$KEY)),'output/wave6_qc_20261113/source_ids_internal.csv',row.names=FALSE)
