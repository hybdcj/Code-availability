# Redraw labels from saved estimates/data. No model fitting or resampling.
library(grid);library(ragg)
figdir <- 'C:/y03tmp/F34_delivery/TIFF'
box <- function(x,y,w,h,txt,fill='#F3F6F8',size=12) {
  grid.rect(x,y,w,h,gp=gpar(fill=fill,col='#7F909D',lwd=1))
  grid.text(txt,x,y,gp=gpar(fontfamily='SimSun',fontsize=size,lineheight=1.35))
}
arr <- function(x1,y1,x2,y2)grid.lines(c(x1,x2),c(y1,y2),arrow=arrow(length=unit(.09,'in'),type='closed'),gp=gpar(col='#667885'))
agg_tiff(file.path(figdir,'Figure1.tif'),width=10,height=7,units='in',res=600,compression='lzw')
grid.newpage()
grid.text('CHARLS development and external validation in three international cohorts',.5,.96,gp=gpar(fontfamily='SimSun',fontsize=17,fontface='bold'))
box(.25,.82,.44,.13,'CHARLS 2011—2020 year five waves\n fixed cohort complete on common items',size=13)
box(.25,.62,.44,.17,'2011 age in year≥45 years old, memory score>4\n2015 year measurable memory\n age, education, depression, ADL complete',size=13)
box(.25,.41,.44,.13,'development sample: 4,441 participants, 658 events\n2011—2015 year outcome: memory score≤4',fill='#DCEAF4',size=13)
box(.25,.20,.44,.15,'four variable continuous model and three variable integer score\n1,000 times Bootstrap internal validation\n frozen formula and risk cutoffs',size=13)
arr(.25,.755,.25,.705);arr(.25,.535,.25,.475);arr(.25,.345,.25,.275)
box(.76,.82,.43,.13,'HRS 2010—2014, S0≤6\n score 14,051 participants; continuous model 13,972 participants',fill='#FFF1DA',size=12)
box(.76,.62,.43,.13,'ELSA 2008—2012, S0≤6\n score and continuous model both 6,725 participants',fill='#FFF1DA',size=12)
box(.76,.42,.43,.13,'SHARE 2011—2015, S0≤5\n score 26,879 participants; continuous model 26,874 participants',fill='#FFF1DA',size=12)
box(.76,.20,.43,.15,'cohort risk sets: baseline memory score>corresponding threshold\n follow up outcome available, required predictors calculable\nS1/S2 threshold sensitivity listed separately',fill='#E7F2E9',size=12)
arr(.47,.20,.545,.20)
grid.text('common items define development source cohort; primary prediction outcome is four year memory threshold crossing.post hoc broader population extension in supplement.',.5,.055,gp=gpar(fontfamily='SimSun',fontsize=10))
dev.off()

agg_tiff(file.path(figdir,'Figure3.tif'),width=8,height=5.2,units='in',res=600,compression='lzw')
par(family='SimSun',mar=c(4.5,4.5,3,1),cex=1)
obs<-c(89/1581,273/1775,296/1085); pred<-c(.059,.148,.278)
xx<-barplot(obs,names.arg=c('low risk(≤−3 points)','intermediate risk(−2 to+4 points)','high risk(≥+5 points)'),
  col=c('#CEE0ED','#83B5D6','#2D729F'),ylim=c(0,.35),ylab='four year memory threshold crossing event rate',border=NA,
  main='CHARLS fixed score risk strata in development sample')
points(xx,pred,pch=18,col='#B13A38',cex=1.5)
text(xx,obs+.025,labels=c('5.6%\n(n=1,581)','15.4%\n(n=1,775)','27.3%\n(n=1,085)'),cex=1)
legend('topleft',c('observed_event_proportion','mean predicted risk'),fill=c('#83B5D6',NA),pch=c(NA,18),col=c(NA,'#B13A38'),bty='n',cex=.9)
dev.off()

proj<-'D:/BaiduSyncdisk/research/Xu_Zhonghao'
cal<-read.csv(file.path(proj,'output/10_screening_index/A09_calibration_deciles_2026-09-26.csv'),check.names=FALSE)
dca<-read.csv(file.path(proj,'output/10_screening_index/A09_dca_2026-09-26.csv'),check.names=FALSE)
a<-readRDS('C:/y03tmp/tmp/a07_memory_v3.rds')
dd<-a$data
lp<-predict(a$model_u,newdata=dd,type='response')
yy<-model.response(model.frame(a$model_u))
# The saved data were already verified to match model frame in F27; use its observed outcome.
stopifnot(length(yy)==nrow(dd),sum(yy)==658)
rocobj<-pROC::roc(yy,lp,direction='<',quiet=TRUE)
stopifnot(abs(as.numeric(pROC::auc(rocobj))-.701519234585697)<1e-10)
agg_tiff(file.path(figdir,'Figure4.tif'),width=12,height=5.5,units='in',res=600,compression='lzw')
par(mfrow=c(1,3),family='SimSun',mar=c(4.5,4.5,3,1),cex=.9)
plot(1-rocobj$specificities,rocobj$sensitivities,type='l',col='#2C7FB8',lwd=2,xlim=c(0,1),ylim=c(0,1),xlab='1 − specificity',ylab='sensitivity',main='(a) apparent ROC',bty='l');abline(0,1,lty=3,col='grey60')
text(.52,.14,'apparent AUC 0.702\n95%CI 0.681–0.722\n optimism corrected point estimate 0.698',cex=.8)
plot(cal[['mean_prediction']],cal[['observed_rate']],xlim=c(0,.4),ylim=c(0,.4),pch=16,col='#279D63',xlab='predicted probability',ylab='observed_event_proportion',main='(b) apparent calibration',bty='l')
abline(0,1,lty=3,col='grey60');segments(cal[['mean_prediction']],cal[['observed_CI_low']],cal[['mean_prediction']],cal[['observed_CI_high']],col='#279D63');lines(cal[['mean_prediction']],cal[['observed_rate']],col='#7971AE')
text(.09,.34,'slope β = 1.000\noffset intercept α ≈ 0',adj=0,cex=.8)
plot(dca[['threshold']],dca[['net_benefit_model']],type='l',col='#2C7FB8',lwd=2,xlab='threshold probability',ylab='net_benefit',ylim=c(-.02,.16),main='(c) apparent DCA',bty='l')
lines(dca[['threshold']],dca[['net_benefit_intervene_in_all']],lty=2,col='#C45A39',lwd=2);abline(h=0,lty=3,col='grey60')
legend('topright',c('continuous model','intervene in all','no intervention'),col=c('#2C7FB8','#C45A39','grey60'),lty=c(1,2,3),bty='n',cex=.8)
dev.off()
