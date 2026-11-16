args <- commandArgs(trailingOnly=TRUE)
if(length(args)) setwd(args[1])
o <- 'output/wave6_intercept_validation_20261116'
d <- read.csv('output/wave6_validation_20261114/fouryear_predictions_internal.csv',colClasses=c(ID='character',communityID='character'))
g <- read.csv(file.path(o,'cluster_index_internal.csv'),colClasses=c(communityID='character'))
cnt <- as.matrix(read.csv(file.path(o,'first50_counts_internal.csv')))
ix <- match(d$communityID,g$communityID)
b <- read.csv(file.path(o,'bootstrap_results.csv'))
s <- read.csv(file.path(o,'performance.csv'))
stopifnot(nrow(d)==1188,sum(d$event)==269,length(unique(d$communityID))==212,!anyNA(ix))
fit <- function(w) glm(event~1+offset(lp),data=d,weights=w,family=binomial(),control=glm.control(epsilon=1e-12,maxit=100))
full <- fit(rep(1,nrow(d)))
a <- unname(coef(full)[1]); app <- mean((plogis(d$lp+a)-d$event)^2)
frozen <- mean((plogis(d$lp)-d$event)^2)
stopifnot(abs(a-s$alpha)<1e-9,abs(app-s$brier_updated_apparent)<1e-12)
checks <- list()
for(i in seq_len(nrow(cnt))) {
 w <- cnt[i,ix]; f <- fit(w); ab <- unname(coef(f)[1])
 loss <- (plogis(d$lp+ab)-d$event)^2
 tr <- weighted.mean(loss,w); te <- mean(loss)
 checks[[i]] <- data.frame(replicate=i,alpha_difference=ab-b$alpha[i],train_difference=tr-b$brier_train[i],test_difference=te-b$brier_test[i])
}
ch <- do.call(rbind,checks)
stopifnot(max(abs(as.matrix(ch[,-1])))<1e-9)
correction <- mean(b$brier_test-b$brier_train)
corrected <- app+correction
stopifnot(abs(corrected-s$brier_updated_corrected)<1e-12,abs(corrected-frozen-s$delta_corrected)<1e-12)
auc <- function(v) {r<-rank(v); n1<-sum(d$event); n0<-nrow(d)-n1;(sum(r[d$event==1])-n1*(n1+1)/2)/(n1*n0)}
stopifnot(abs(auc(d$lp)-auc(d$lp+a))<1e-12)
write.csv(ch,file.path(o,'independent_R_checks.csv'),row.names=FALSE)
cat('PASS: full fit, first 50 cluster resamples, correction sign/mean, corrected delta, invariant AUC.\n')
cat(sprintf('Corrected Brier %.10f; corrected delta %.10f\n',corrected,corrected-frozen))
