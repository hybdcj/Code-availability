invisible(Sys.setlocale('LC_CTYPE','Chinese_China.utf8'))
.libPaths(c('D:/R/library',.libPaths()))
suppressPackageStartupMessages({library(dplyr);library(pROC);library(ragg)})
P<-'D:/BaiduSyncdisk/科研/徐中豪'
ROOT<-file.path(P,'16.论文修订第三轮')
FIG<-'C:/y03tmp/F36_delivery/TIFF'
dir.create(FIG,recursive=TRUE,showWarnings=FALSE)
oldfig<-file.path(ROOT,'F35_图表文内排版修订版_2026-10-16/TIFF图片')
file.copy(list.files(oldfig,full.names=TRUE),FIG,overwrite=TRUE)
labels_used<-character()
load_named<-function(path,names) {
  for(ex in parse(path,encoding='UTF-8')) {
    if(is.call(ex)&&as.character(ex[[1]])%in%c('<-','=')&&is.symbol(ex[[2]])&&as.character(ex[[2]])%in%names)
      eval(ex,envir=.GlobalEnv)
  }
}
vertical<-function(fun,n=3,left=5) {
  txt<-paste(deparse(fun),collapse='\n')
  txt<-gsub('mfrow = c(1, 3)','mfrow = c(3, 1)',txt,fixed=TRUE)
  txt<-gsub('mfrow = c(1, 2)','mfrow = c(2, 1)',txt,fixed=TRUE)
  txt<-gsub('family = "sans"','family = "SimSun"',txt,fixed=TRUE)
  f<-eval(parse(text=txt));environment(f)<-.GlobalEnv;f
}
savefig<-function(id,fun,w=6.65,h=8.2) {
  agg_tiff(file.path(FIG,paste0(id,'.tif')),width=w,height=h,units='in',res=600,compression='lzw',background='white',pointsize=10)
  par(family='SimSun',cex=.85,cex.main=1.05)
  fun();dev.off();cat('Saved ',id,'\n',sep='')
}

# Reuse drawing functions and saved numerical results only. Do not run analyses.
s32<-file.path(P,'code/32_A16_figures.R')
load_named(s32,c('parse_leading_num','parse_est','orm','risk','iadl2','or_stats','risk_est','risk_hr','iadl_est','plot_candidate_stability'))
plot_candidate_stability<-vertical(plot_candidate_stability,2)
savefig('Figure2',function()plot_candidate_stability('ZH'))

# Internal-performance plot from the already checked frozen model and saved curves.
src<-readLines(file.path(ROOT,'_核查/脚本/F34_更新正文图件.R'),encoding='UTF-8',warn=FALSE)
src<-src[grep("^proj<-",src):length(src)]
txt<-paste(src,collapse='\n')
txt<-gsub('width=12,height=5.5','width=6.65,height=8.2',txt,fixed=TRUE)
txt<-gsub('mfrow=c(1,3)','mfrow=c(3,1)',txt,fixed=TRUE)
txt<-gsub("res=600,compression='lzw'","res=600,compression='lzw',pointsize=10",txt,fixed=TRUE)
figdir<-FIG;eval(parse(text=txt))

res<-readRDS(file.path(P,'output/11_外部验证/A12_results_2026-09-26.rds'))
perf<-res$perf;meta<-res$meta;deciles<-res$deciles;dca<-res$dca
a12<-res;a07<-readRDS(file.path(P,'output/10_筛查指数/screening_index_memory_v3_2026-09-26.rds'))
model_u<-a07$model_u;model_levels<-model_u$xlevels$edu_group
cohorts<-c('HRS','ELSA','SHARE');cohort_cols<-c(HRS='#1F77B4',ELSA='#2CA02C',SHARE='#D62728')
cut_defs<-c('baseline_all_mean_sd','complete_fu_mean_sd','fixed_charls_le4')
scenario_lab<-setNames(c('S0：基线符合条件人群阈值','S1：随访可测人群阈值','S2：统一固定阈值≤4'),cut_defs)
risk_sets<-list()
for(co in cohorts) {
 x<-readRDS(file.path(P,'external_data',co,'analysis',paste0(co,'_sensitivity_risk_sets_2026-09-26.rds')))[[paste0(co,'__baseline_all_mean_sd')]]
 x$edu_group<-factor(as.character(x$edu_group),levels=model_levels)
 x$p_hat<-as.numeric(predict(model_u,newdata=data.frame(age=x$age,edu_group=x$edu_group,dep=x$dep,adl=x$adl),type='response'))
 risk_sets[[co]]<-x
}
load_named(s32,'draw_external');draw_external<-vertical(draw_external)
body(draw_external)<-parse(text=gsub('mar = c(4.8, 4.8, 3, 1.5)','mar = c(4.8, 8, 3, 1.5)',paste(deparse(body(draw_external)),collapse='\n'),fixed=TRUE))[[1]]
savefig('Figure5',function()draw_external('ZH'))
s25<-file.path(P,'code/25_a12_figures_report.R')
load_named(s25,c('draw_sens','draw_cal_sens','draw_dca'))
draw_sens<-vertical(draw_sens);draw_cal_sens<-vertical(draw_cal_sens);draw_dca<-vertical(draw_dca)
savefig('FigureS3',draw_sens);savefig('FigureS5',draw_cal_sens);savefig('FigureS6',draw_dca)
s26<-file.path(P,'code/26_a13_cross_cohort_heterogeneity.R')
load_named(s26,c('draw_score_forest','draw_cal_het'))
cal_het<-read.csv(file.path(P,'output/11_外部验证/A13_校准异质性.csv'),check.names=FALSE,fileEncoding='UTF-8-BOM')
cal_het$scenario<-c('S0','S1','S2')[match(cal_het$cut_definition,cut_defs)]
draw_score_forest<-vertical(draw_score_forest);draw_cal_het<-vertical(draw_cal_het,2)
# The original margin is too narrow for country/scenario labels on portrait pages.
body(draw_cal_het)<-parse(text=gsub('mar = c(5, 4, 3, 1)','mar = c(5, 8, 3, 1)',paste(deparse(body(draw_cal_het)),collapse='\n'),fixed=TRUE))[[1]]
savefig('FigureS7',draw_score_forest);savefig('FigureS8',draw_cal_het)

OUT<-file.path(P,'output/04_纵向队列')
readout<-function(n)read.csv(file.path(OUT,n),check.names=FALSE,fileEncoding='UTF-8-BOM')
retention<-readout('attrition_retention_2026-09-20.csv')
balance_ipw<-readout('attrition_balance_ipw_2026-09-20.csv')
balance_ipw$label<-gsub('ms10','共同认知题集分',balance_ipw$label,fixed=TRUE)
ipw_models<-readout('ipw_longitudinal_models_2026-09-20.csv')
load_named(file.path(P,'code/14_attrition_selection.R'),'draw_attrition_figure')
draw_attrition_figure<-vertical(draw_attrition_figure)
body(draw_attrition_figure)<-parse(text=gsub('mar = c(5, 7, 3, 1)','mar = c(5, 12, 3, 1)',paste(deparse(body(draw_attrition_figure)),collapse='\n'),fixed=TRUE))[[1]]
body(draw_attrition_figure)<-parse(text=gsub('ms10 available (%)','共同认知题集可计分比例（%）',paste(deparse(body(draw_attrition_figure)),collapse='\n'),fixed=TRUE))[[1]]
body(draw_attrition_figure)<-parse(text=gsub('"topright"','"bottomright"',paste(deparse(body(draw_attrition_figure)),collapse='\n'),fixed=TRUE))[[1]]
savefig('FigureS1',draw_attrition_figure)
models<-readout('weighted_longitudinal_models_2026-09-20.csv')
load_named(file.path(P,'code/13_weighted_longitudinal.R'),'draw_weighted_figure')
draw_weighted_figure<-vertical(draw_weighted_figure)
body(draw_weighted_figure)<-parse(text=gsub('mar = c(5, 7, 3, 1)','mar = c(5, 12, 3, 1)',paste(deparse(body(draw_weighted_figure)),collapse='\n'),fixed=TRUE))[[1]]
savefig('FigureS2',draw_weighted_figure)

# Dense sensitivity forest: four labelled portrait panels preserve all plotted rows.
z<-readRDS(file.path(P,'output/13_A15_敏感性分析/A15_results_2026-09-26.rds'))
summary_internal<-z$summary_internal;key_or<-z$key_or
cleanlab<-function(s) {
 s<-gsub('ms10','共同认知题集',s,fixed=TRUE)
 s<-gsub('mem2011','基线记忆分',s,fixed=TRUE)
 s<-gsub('mem2015','随访记忆分',s,fixed=TRUE)
 s<-gsub('A14 ','',s,fixed=TRUE);s<-gsub('原A15','内部',s,fixed=TRUE)
 s<-gsub('基线校正后跨阈值（含 基线记忆分）','基线记忆校正跨阈值',s,fixed=TRUE)
 s<-gsub('宽样本 MICE m=20（协变量缺失；结局观察者）','宽样本：协变量插补（m=20）',s,fixed=TRUE)
 s<-gsub('宽样本 MICE m=20（结局缺失，MAR；含固定队列辅助变量）','宽样本：结局插补（m=20）',s,fixed=TRUE)
 s<-gsub('宽样本（不要求五波 共同认知题集 完整）','宽样本（放宽五波完整性限制）',s,fixed=TRUE)
 s
}
plot_tab<-summary_internal%>%filter(!grepl('插补|MICE',scenario),!is.na(auc_model))%>%distinct(domain,scenario,.keep_all=TRUE)
plot_tab<-plot_tab[order(plot_tab$auc_model),]
plot_tab$label<-cleanlab(paste0(plot_tab$domain,' | ',plot_tab$scenario))
auc_xlim<-range(c(plot_tab$auc_low,plot_tab$auc_high),na.rm=TRUE)+c(-.02,.02)
savefig('FigureS10a',function(){
 par(mar=c(4,19,3,1),cex=.85)
 yy<-seq_len(nrow(plot_tab))
 plot(plot_tab$auc_model,yy,xlim=auc_xlim,ylim=c(.5,nrow(plot_tab)+.5),pch=16,col='#1F4E78',xlab='AUC（连续模型）',ylab='',main='(a) 内部敏感性分析：区分度',yaxt='n')
 segments(plot_tab$auc_low,yy,plot_tab$auc_high,yy,col='#9DC3E6',lwd=2)
 abline(v=summary_internal$auc_model[summary_internal$scenario=='主模型（未加权，冻结）'][1],lty=2,col='#C00000')
 axis(2,at=yy,labels=plot_tab$label,las=1,cex.axis=.75);grid(col='grey90')
},h=8.4)
plot_or<-key_or%>%filter(term!='adl',!is.na(OR),!grepl('插补',scenario))
plot_or$term_plot<-ifelse(plot_or$term=='age','年龄（每5岁）',ifelse(plot_or$term=='dep','抑郁','高中及以上教育'))
plot_or$OR_plot<-ifelse(plot_or$term=='age',plot_or$OR^5,plot_or$OR)
plot_or$low_plot<-ifelse(plot_or$term=='age',plot_or$OR_low^5,plot_or$OR_low)
plot_or$high_plot<-ifelse(plot_or$term=='age',plot_or$OR_high^5,plot_or$OR_high)
plot_or$label<-cleanlab(plot_or$scenario)
terms<-c('年龄（每5岁）','抑郁','高中及以上教育')
cols<-c('#1F4E78','#C00000','#2E7D32')
for(i in seq_along(terms)) {
 dd<-plot_or[plot_or$term_plot==terms[i],];dd<-dd[order(dd$OR_plot),]
 savefig(paste0('FigureS10',letters[i+1]),function(){
  par(mar=c(4,17,3,1),cex=.85);y<-seq_len(nrow(dd))
  plot(dd$OR_plot,y,log='x',xlim=range(c(dd$low_plot,dd$high_plot,1),na.rm=TRUE),ylim=c(.5,nrow(dd)+.5),pch=16,col=cols[i],xlab='OR（对数刻度）',ylab='',main=paste0('(',letters[i+1],') 关联方向：',terms[i]),yaxt='n')
  segments(dd$low_plot,y,dd$high_plot,y,col=cols[i],lwd=2);abline(v=1,lty=2,col='grey40')
  axis(2,at=y,labels=dd$label,las=1,cex.axis=.8);grid(col='grey90')
 },h=7.5)
}
writeLines(c(plot_tab$label,plot_or$label,scenario_lab,cal_het$scenario,balance_ipw$label,ipw_models$exposure_label,models$exposure_label),'C:/y03tmp/F36_delivery/图内标签核对.txt',useBytes=TRUE)
saveRDS(list(plot_tab=plot_tab,plot_or=plot_or),'C:/y03tmp/F36_delivery/敏感性图绘图底表.rds')
cat('PASS: saved estimates reused; no models refitted, no new bootstrap.\n')
