a15dir<-file.path(project,'output/13_A15_敏感性分析')
outdir<-file.path(root,'_核查/输出')
spec <- list()
addspec<-function(id,title,data,note,source,mice=FALSE,raw=FALSE) spec[[length(spec)+1]] <<- list(id=id,title=title,data=data,note=note,source=source,mice=mice,raw=raw)
base_notes <- '注：范围通常为原分析95%CI；校准截距α使用预测logit作为offset，斜率固定为1。模型角色与分析情境分别列示；重拟合指标为同样本表观表现。—表示未报告或不适用。'
mi_note <- '注：MICE性能为20次插补的均值；括号内为插补间2.5%–97.5%经验分位数范围，不是95%抽样CI。连续模型逐次重拟合并在同一插补样本评价；结局事件数为插补平均。回归系数表的区间则按Rubin规则合并。'
map <- c('样本流向_2026-09-26','内部敏感性汇总','外部敏感性汇总','结论稳定性判定','关键OR方向稳定性',
 '权重_设计诊断','权重_模型比较','权重_性能比较','IPW_诊断','IPW_选择模型','IPW_平衡','IPW_模型比较','IPW_性能比较',
 '缺失_固定队列与宽样本比较','缺失_宽样本模型比较','缺失_宽样本性能比较','缺失_外部ADL敏感性',
 '缺失_MICE_协变量_系数','缺失_MICE_协变量_性能汇总','缺失_MICE_结局_系数','缺失_MICE_结局_性能汇总',
 '切点_定义汇总','切点_内部模型比较','切点_内部性能比较','切点_外部汇总','切点_外部补充敏感性',
 '结局定义_连续变化模型','结局定义_模型比较','结局定义_性能比较','2018测量差异_模型比较','2018测量差异_性能比较')
titles<-c('样本流向','内部敏感性汇总','外部敏感性汇总','敏感性结果及解释边界','关键OR方向',
 '调查权重：设计诊断','调查权重：模型系数','调查权重：性能','IPW：诊断','IPW：选择模型','IPW：平衡','IPW：模型系数','IPW：性能',
 '缺失数据：固定队列与宽样本比较','缺失数据：宽样本模型系数','缺失数据：宽样本性能','外部ADL敏感性',
 '协变量MICE：合并系数','协变量MICE：平均表观性能及插补间范围','结局MICE：合并系数','结局MICE：平均表观性能及插补间范围',
 '结局阈值：定义','结局阈值：内部模型系数','结局阈值：内部性能','结局阈值：外部S0/S1/S2','结局阈值：外部补充敏感性',
 '结局定义：连续变化模型','结局定义：模型系数','结局定义：性能','2018测量条件限制：模型系数','2018测量条件限制：性能')
for(i in seq_along(map)) {
  src<-file.path(a15dir,paste0('A15_',map[i],'.csv')); d<-readc(src)
  note<-base_notes
  if(i==2)note<-paste(base_notes,mi_note,'主结局性能已按参与者行序修正；MICE均值和分位数按保存插补对象核验。')
  if(i %in% c(18,20))note<-'注：20次插补后按Rubin规则合并系数、标准误及OR置信区间。与插补性能的经验分位数范围不同。'
  if(i %in% c(19,21))note<-mi_note
  if(i==4) {
    d$conclusion[d$domain=='总体']<-'已报告敏感性分析中的主要关联方向大体一致；AUC、事件率和校准随样本及阈值变化。新增分析包含事后探索，不表示全部预先指定。'
    d$status[d$domain=='总体']<-'限于所报告情境'
    d$conclusion[d$domain=='缺失/选择']<-'宽样本8,296人、事件率22.8%；重拟合AUC 0.729。MICE平均AUC 0.728/0.734为逐次重拟合表观值，插补间范围不是95%CI。'
    note<-'注：该表概括已报告分析，不构成等效性检验；方向相近不能证明选择偏倚消除或绝对风险校准稳定。'
  }
  if(i %in% c(30,31))note<-paste(base_notes,'此处为原A15按2018年测量条件限制的分析，不是2018—2020后续窗口验证；后者见表S36。')
  if(i==25)note<-paste(base_notes,'S0为外部主分析；S1为随访可测样本定义阈值，S2为固定≤4。总n为评分风险集，模型n为ADL完整风险集。')
  if(i==29) {
    d<-readc(file.path(root,'F29_S29结局定义性能_修正版.csv'))
    note<-paste(base_notes,'主模型采用冻结预测，已修正预测与结局行序；下降≥2分和加入基线记忆为另行拟合的表观结果。')
  }
  addspec(paste0('TableS',i),paste0('补充表 S',i,'  ',titles[i]),d,note,src,mice=i %in% c(2,19,21))
}

# S32: every verified F32 performance row, with independent/common sample distinctions.
src<-file.path(root,'F32_全人群性能与分母_修正版.csv');d<-readc(src)
d<-d[,intersect(c('window','definition','arm','model','population','n','events','auc','auc_low','auc_high','slope','offset_intercept','development_id_overlap'),names(d))]
names(d)[names(d)=='development_id_overlap']<-'开发者重叠'
addspec('TableS32a','补充表 S32a  CHARLS放宽基线记忆限制后的性能',d,
 '注：U=N+I+M，B=N+I；N为基线高于低分界值、I为基线已低分、M为基线记忆缺失。评分可用样本不要求ADL完整，共同样本要求ADL完整；两工具直接比较须在共同样本内进行。区间为DeLong 95%CI。教育分层定义属事后探索，使用严格小于界值判低分。',src)
src<-file.path(outdir,'F28_结局与数字底表.csv');d<-readc(src)
u<-d[grepl('^U_',d$scenario),]
stopifnot(nrow(u)>0)
u<-u[,c('cohort','window','scenario','model','population','n','events','AUC','AUC_low','AUC_high','calibration_offset_intercept','calibration_slope')]
names(u)[names(u)=='calibration_offset_intercept']<-'offset_intercept';names(u)[names(u)=='calibration_slope']<-'slope'
addspec('TableS32b','补充表 S32b  外部全人群S0及S2性能',u,
 '注：取消基线记忆状态限制后，预测目标为随访低分状态。S0采用HRS/ELSA≤6、SHARE≤5；S2均≤4。模型与评分分别使用其可计算样本，评分不套用连续模型校准。全人群AUC区间沿用已保存DeLong结果。',src)

src<-file.path(outdir,'F32_样本流程.csv');d<-readc(src)
names(d)<-c('窗口','结局定义','波次匹配人数','随访记忆可测','评分可用n','共同样本n','恢复ADL缺失者','恢复者事件','评分样本开发者重叠','共同样本开发者重叠')
addspec('TableS33a','补充表 S33a  CHARLS全人群入组流程',d,
 '注：恢复者均因ADL缺失被旧共同样本条件排除；不以0填补其ADL。评分仍要求年龄≥45岁、三项变量可计算和随访结局可得。',src,raw=TRUE)
src<-file.path(outdir,'F32_事件构成分母.csv');d<-readc(src)
d$baseline_low_fraction_of_ALL_events<-paste0(fmt(100*d$baseline_low_fraction_of_ALL_events,2),'%')
d$baseline_low_fraction_of_KNOWN_baseline_events<-paste0(fmt(100*d$baseline_low_fraction_of_KNOWN_baseline_events,2),'%')
names(d)<-c('窗口','结局定义','样本','U人数','N人数','I人数','M人数','U事件','N事件','I事件','M事件','I/全部事件','I/基线已知事件','基线已知事件数')
addspec('TableS33b','补充表 S33b  全人群事件构成与两种分母',d,
 '注：I组事件指两次测量均低分，不能推断中间始终低分。M组事件无法判定是否新发生跨阈值；I/全部事件与I/基线状态已知事件的分母不同。',src,raw=TRUE)

screen_specs<-list(c('TableS34a','CHARLS全人群筛查',file.path(root,'F32_全人群筛查指标_修正版.csv')),
 c('TableS34b','外部全人群S0筛查',file.path(outdir,'F16_screening.csv')),
 c('TableS34c','外部全人群S2筛查',file.path(outdir,'F19_screening.csv')))
for(ss in screen_specs) {
  dd<-readc(ss[3])
  if(ss[1]!='TableS34a') {
    dd<-dd[,c('队列或窗口','结局定义','目标筛查比例','阈值','n','事件','检出人数','实际筛查比例','TP','FP','TN','FN','灵敏度','特异度','阳性预测值','阴性预测值','阳性且基线状态可知人数','阳性中基线已低分人数','新发真事件数','每检出1例新发需筛人数')]
    names(dd)[names(dd)=='新发真事件数']<-'可确认跨阈值TP'
    names(dd)[names(dd)=='每检出1例新发需筛人数']<-'评分人数/跨阈值TP'
    addspec(ss[1],paste0('补充表 ',sub('Table','',ss[1]),'  ',ss[2]),dd,
      '注：实际筛查比例、灵敏度、特异度及预测值的单位均为%。目标比例按分数切点纳入全部同分者，实际人数可超过目标。事件为随访低分状态；可确认跨阈值TP要求基线高于相应界值。',ss[3])
    next
  }
  # Explicit whitelist; obsolete specificity fields and historical metadata must not be displayed.
  keep<-intersect(c('window','definition','cohort','target','target_fraction','strategy','threshold','cut','n','events','selected','n_pos','actual_fraction','selected_fraction',
    'TP','FP','TN','FN','tp','fp','tn','fn','sensitivity','specificity','ppv','npv','PPV','NPV','TP_baseline_above','TP_baseline_low','TP_baseline_unknown',
    'detected_known_new','detected_baseline_low','detected_baseline_unknown','tp_N','tp_I','tp_M',
    'total_scored_per_detected_crossing','selected_per_detected_crossing'),names(dd))
  if(!all(c('sensitivity','specificity')%in%keep))stop('Screen fields need inspection')
  dd<-dd[,keep,drop=FALSE]
  dd<-data.frame(dd[,c('window','definition')],选点方式=ifelse(is.na(dd$target_fraction),'约登点','目标比例'),dd[,setdiff(names(dd),c('window','definition'))],check.names=FALSE)
  names(dd)[names(dd)=='sensitivity']<-'灵敏度'
  addspec(ss[1],paste0('补充表 ',sub('Table','',ss[1]),'  ',ss[2]),dd,
   '注：特异度=TN/(TN+FP)，NPV=TN/(TN+FN)。比例按0–1列示；目标比例切点纳入全部同分者，须以实际阳性人数/比例解释。约登点由同样本选取，仅供探索。事件为随访低分状态，不是当前临床认知疾病诊断。',ss[3])
}

src<-file.path(outdir,'F26_BDSI_AUC.csv');d<-readc(src)
d<-d[d$scope=='65至79岁' & d$stage=='条目可计分配对样本' & d$model%in%c('本研究三项评分','BDSI条目修正映射'),]
if(!nrow(d)) {print(unique(readc(src)[,c('scope','stage')]));stop('BDSI filter')}
d$exclusion<-'未排除开发者';d<-d[,c('outcome','exclusion','model','n','events','AUC','CI_low','CI_high')]
p27<-readc(file.path(outdir,'F27_AUC.csv'))
# F27 uses the same column names for model-level AUCs.
d<-rbind(d,p27[,names(d)])
addspec('TableS35a','补充表 S35a  65～79岁映射BDSI比较：AUC',d,
 '注：CHARLS 2015—2018年事后比较。T1为教育分层、T2为固定≤4；仅排除本研究开发参与者不构成独立外部验证。原BDSI预测六年痴呆，本表评价三年记忆低分。年龄65岁计0分；功能项要求理财或服药需帮助；抑郁项为“做事费劲”或抗抑郁药；教育仍为学历代理。完整样本沿用原三臂附加入选条件。',paste(src,'F27_AUC.csv'))
src<-file.path(outdir,'F27_AUC配对比较.csv');d<-readc(src)
addspec('TableS35b','补充表 S35b  排除开发参与者后的配对AUC差异',d,
 '注：配对DeLong 95%CI；ΔAUC=本评分减映射BDSI。P为双侧且未作多重比较校正；两种排除范围并非独立验证样本。',src)
src<-file.path(outdir,'F27_同名额灵敏度.csv');d<-readc(src)
d$bootstrap_reps<-NULL;d$seed<-NULL
addspec('TableS35c','补充表 S35c  排除开发参与者后的同名额期望灵敏度',d,
 '注：名额=floor(0.30×n)；切点同分者按统一概率分配剩余名额，表内为期望灵敏度。差值以百分点表示，区间来自2,000次参与者配对Bootstrap。随机分配未实际实施。T1差值区间跨0，T2不跨0。',src)
src<-file.path(outdir,'F26_BDSI_同名额灵敏度比较.csv');d<-readc(src)
d$bootstrap_reps<-NULL;d$seed<-NULL
addspec('TableS35d','补充表 S35d  排除前65～79岁同名额比较',d,
 '注：修正条目映射后的年龄适用范围内分析；尚含开发参与者。与表S35c并列解释，不能选取较高AUC或较小P值作为独立验证。',src)
src<-file.path(outdir,'F27_筛查指标.csv');d<-readc(src)
drop<-grep('^seed$|^bootstrap_reps$',names(d),value=TRUE);d[drop]<-NULL
names(d)[names(d)=='sensitivity']<-'灵敏度(%)'
names(d)[names(d)=='specificity']<-'特异度(%)'
names(d)[names(d)=='PPV']<-'PPV(%)';names(d)[names(d)=='NPV']<-'NPV(%)'
addspec('TableS35e','补充表 S35e  排除后的筛查策略与实际负担',d,
 '注：区分固定名额期望、全部同分者纳入及沿用F26切点。固定阈值的实际阳性比例不等于30%。比例单位沿表头；期望四格表允许非整数，不能解释为已执行随机分配后的观察人数。',src)

src<-file.path(outdir,'F31_可用性能数字底表.csv');d<-readc(src)
for(task in c('F03','F05')) {
  dd<-d[d$task==task,];dd<-dd[,c('analysis','sample','outcome','model','n','events','development_id_overlap','auc','auc_low','auc_high','slope','offset_intercept','brier')]
  names(dd)[names(dd)=='development_id_overlap']<-'开发者重叠'
  addspec(if(task=='F03')'TableS36a' else 'TableS36b',paste0('补充表 ',if(task=='F03')'S36a  2018—2020冻结、重校准与重拟合' else 'S36b  2018—2020变化结局探索'),dd,
   '注：同队列后续窗口事后分析。冻结预测、同样本截距更新和重拟合分别标示；重校准后的表观改善不是独立验证。波内标准化不保证跨波测量等值；不据此推算MDC、真实错分率或模型性能上限。',src)
}
src<-file.path(outdir,'F31_可用性能数字底表.csv');d<-readc(src)
# Use the independently verified apparent values, never the legacy corrected columns.
d<-d[d$task=='F04',c('analysis','model','n','events','auc')]
names(d)<-c('模型','说明','n','事件','表观AUC')
d$表观AUC<-fmt(d$表观AUC,4)
addspec('TableS36c','补充表 S36c  2018—2020候选变量模型的表观AUC',d,
 '注：共同完整样本9,178人、397事件。全部模型在同一样本拟合和评价，共12个模型含参照模型。仅报告已核验的表观AUC，不报告旧校正增益、差值区间或“达标”判断；不能据此认定新增变量无价值。',src,raw=TRUE)
