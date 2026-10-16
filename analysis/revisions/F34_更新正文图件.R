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
grid.text('CHARLS开发及三个国际队列外部验证',.5,.96,gp=gpar(fontfamily='SimSun',fontsize=17,fontface='bold'))
box(.25,.82,.44,.13,'CHARLS 2011—2020年五波\n共同题集完整的固定队列',size=13)
box(.25,.62,.44,.17,'2011年年龄≥45岁、记忆分>4\n2015年记忆可测\n年龄、教育、抑郁、ADL完整',size=13)
box(.25,.41,.44,.13,'开发样本：4,441人，658事件\n2011—2015年结局：记忆分≤4',fill='#DCEAF4',size=13)
box(.25,.20,.44,.15,'四变量连续模型与三变量整数评分\n1,000次Bootstrap内部验证\n冻结公式及风险层切点',size=13)
arr(.25,.755,.25,.705);arr(.25,.535,.25,.475);arr(.25,.345,.25,.275)
box(.76,.82,.43,.13,'HRS 2010—2014，S0≤6\n评分14,051人；连续模型13,972人',fill='#FFF1DA',size=12)
box(.76,.62,.43,.13,'ELSA 2008—2012，S0≤6\n评分及连续模型均6,725人',fill='#FFF1DA',size=12)
box(.76,.42,.43,.13,'SHARE 2011—2015，S0≤5\n评分26,879人；连续模型26,874人',fill='#FFF1DA',size=12)
box(.76,.20,.43,.15,'各队列风险集：基线记忆分>相应阈值\n随访结局可得，所需预测变量可计算\nS1/S2阈值敏感性分析另列',fill='#E7F2E9',size=12)
arr(.47,.20,.545,.20)
grid.text('共同题集用于定义开发来源队列；主预测结局为四年记忆分跨阈值。全人群事后扩展另见补充材料。',.5,.055,gp=gpar(fontfamily='SimSun',fontsize=10))
dev.off()

agg_tiff(file.path(figdir,'Figure3.tif'),width=8,height=5.2,units='in',res=600,compression='lzw')
par(family='SimSun',mar=c(4.5,4.5,3,1),cex=1)
obs<-c(89/1581,273/1775,296/1085); pred<-c(.059,.148,.278)
xx<-barplot(obs,names.arg=c('低风险（≤−3分）','中风险（−2至+4分）','高风险（≥+5分）'),
  col=c('#CEE0ED','#83B5D6','#2D729F'),ylim=c(0,.35),ylab='四年记忆分跨阈值事件率',border=NA,
  main='CHARLS开发样本的固定评分风险层')
points(xx,pred,pch=18,col='#B13A38',cex=1.5)
text(xx,obs+.025,labels=c('5.6%\n(n=1,581)','15.4%\n(n=1,775)','27.3%\n(n=1,085)'),cex=1)
legend('topleft',c('观察事件率','平均预测风险'),fill=c('#83B5D6',NA),pch=c(NA,18),col=c(NA,'#B13A38'),bty='n',cex=.9)
dev.off()

proj<-'D:/BaiduSyncdisk/科研/徐中豪'
cal<-read.csv(file.path(proj,'output/10_筛查指数/A09_calibration_deciles_2026-09-26.csv'),check.names=FALSE)
dca<-read.csv(file.path(proj,'output/10_筛查指数/A09_dca_2026-09-26.csv'),check.names=FALSE)
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
plot(1-rocobj$specificities,rocobj$sensitivities,type='l',col='#2C7FB8',lwd=2,xlim=c(0,1),ylim=c(0,1),xlab='1 − 特异度',ylab='灵敏度',main='(a) 表观ROC',bty='l');abline(0,1,lty=3,col='grey60')
text(.52,.14,'表观AUC 0.702\n95%CI 0.681–0.722\n乐观度校正点估计 0.698',cex=.8)
plot(cal[['预测均值']],cal[['观测率']],xlim=c(0,.4),ylim=c(0,.4),pch=16,col='#279D63',xlab='预测概率',ylab='观察事件率',main='(b) 表观校准',bty='l')
abline(0,1,lty=3,col='grey60');segments(cal[['预测均值']],cal[['观测CI低']],cal[['预测均值']],cal[['观测CI高']],col='#279D63');lines(cal[['预测均值']],cal[['观测率']],col='#7971AE')
text(.09,.34,'斜率β = 1.000\noffset截距α ≈ 0',adj=0,cex=.8)
plot(dca[['阈值']],dca[['净获益_模型']],type='l',col='#2C7FB8',lwd=2,xlab='阈值概率',ylab='净获益',ylim=c(-.02,.16),main='(c) 表观DCA',bty='l')
lines(dca[['阈值']],dca[['净获益_全部干预']],lty=2,col='#C45A39',lwd=2);abline(h=0,lty=3,col='grey60')
legend('topright',c('连续模型','全部干预','不干预'),col=c('#2C7FB8','#C45A39','grey60'),lty=c(1,2,3),bty='n',cex=.8)
dev.off()
