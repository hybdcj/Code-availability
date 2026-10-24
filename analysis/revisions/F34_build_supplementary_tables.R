a15dir<-file.path(project,'output/13_A15_sensitivity_analysis')
outdir<-file.path(root,'_audit/output')
spec <- list()
addspec<-function(id,title,data,note,source,mice=FALSE,raw=FALSE) spec[[length(spec)+1]] <<- list(id=id,title=title,data=data,note=note,source=source,mice=mice,raw=raw)
base_notes <- 'note: range generally from original analysis 95%CI; calibration intercept α uses predictions logit as offset, slope fixed at 1.model roles and analysis contexts listed separately; refit metrics are same sample apparent performance.—denotes not reported or not applicable.'
mi_note <- 'note: MICE performance is 20 mean across imputations; parentheses show between imputation 2.5%–97.5%empirical quantile range, is not 95%sampling CI.continuous model refitted and evaluated within each imputed sample; event count is imputation mean.coefficient table intervals follow Rubin pooling rules.'
map <- c('sample flow_2026-09-26','internal sensitivity summary','external sensitivity summary','assessment of conclusion stability','key OR directional stability',
 'weights_design diagnostics','weights_model comparison','weights_performance comparison','IPW_diagnostics','IPW_selection model','IPW_balance','IPW_model comparison','IPW_performance comparison',
 'missing_fixed cohort versus broader sample','missing_broader sample model comparison','missing_broader sample performance comparison','missing_external ADL sensitivity',
 'missing_MICE_covariates_coefficients','missing_MICE_covariates_performance summary','missing_MICE_outcome_coefficients','missing_MICE_outcome_performance summary',
 'cutoff_definition summary','cutoff_internal model comparison','cutoff_internal performance comparison','cutoff_external summary','cutoff_additional external sensitivity',
 'outcome definition_continuous change model','outcome definition_model comparison','outcome definition_performance comparison','2018 measurement differences_model comparison','2018 measurement differences_performance comparison')
titles<-c('sample flow','internal sensitivity summary','external sensitivity summary','sensitivity results and interpretation boundaries','key OR direction',
 'survey weights: design diagnostics','survey weights: model coefficients','survey weights: performance','IPW: diagnostics','IPW: selection model','IPW: balance','IPW: model coefficients','IPW: performance',
 'missing data: fixed cohort versus broader sample','missing data: broader sample model coefficients','missing data: broader sample performance','external ADL sensitivity',
 'covariates MICE: pooled coefficients','covariates MICE: mean apparent performance and between imputation range','outcome MICE: pooled coefficients','outcome MICE: mean apparent performance and between imputation range',
 'outcome threshold: definition','outcome threshold: internal model coefficients','outcome threshold: internal performance','outcome threshold: external S0/S1/S2','outcome threshold: additional external sensitivity',
 'outcome definition: continuous change model','outcome definition: model coefficients','outcome definition: performance','2018 measurement limitations: model coefficients','2018 measurement limitations: performance')
for(i in seq_along(map)) {
  src<-file.path(a15dir,paste0('A15_',map[i],'.csv')); d<-readc(src)
  note<-base_notes
  if(i==2)note<-paste(base_notes,mi_note,'primary outcome performance corrected for participant row alignment; MICE means and quantiles verified against saved imputation objects.')
  if(i %in% c(18,20))note<-'note: 20 after imputations apply Rubin coefficient pooling rules, standard errors and OR confidence interval.distinct from empirical performance quantile range across imputations.'
  if(i %in% c(19,21))note<-mi_note
  if(i==4) {
    d$conclusion[d$domain=='overall']<-'main association directions broadly agree across reported sensitivities; AUC, event rates and calibration vary by sample and threshold.additional analyses include post hoc exploration, not all analyses prespecified.'
    d$status[d$domain=='overall']<-'limited to reported contexts'
    d$conclusion[d$domain=='missing/selection']<-'broader sample 8,296 participants, event rate 22.8%; refit AUC 0.729.MICE mean AUC 0.728/0.734 per imputation refitted apparent values, between imputation range is not 95%CI.'
    note<-'note: table summarizes reported analyses, not an equivalence test; similar directions do not establish elimination of selection bias or stable absolute calibration.'
  }
  if(i %in% c(30,31))note<-paste(base_notes,'here original A15 by 2018 year measurement restriction analysis, is not 2018—2020 later window validation; latter shown in table S36.')
  if(i==25)note<-paste(base_notes,'S0 primary external analysis; S1 threshold defined in follow up measurable sample, S2 is fixed≤4.total n score risk set, model n is ADL complete risk set.')
  if(i==29) {
    d<-readc(file.path(root,'F29_S29_performance_by_outcome_definition_corrected_version.csv'))
    note<-paste(base_notes,'primary model uses frozen predictions, prediction outcome row alignment corrected; decline≥2 points and baseline memory addition are separately fitted apparent results.')
  }
  addspec(paste0('TableS',i),paste0('supplementary table S',i,'  ',titles[i]),d,note,src,mice=i %in% c(2,19,21))
}

# S32: every verified F32 performance row, with independent/common sample distinctions.
src<-file.path(root,'F32_broader_population_performance_and_denominators_corrected_version.csv');d<-readc(src)
d<-d[,intersect(c('window','definition','arm','model','population','n','events','auc','auc_low','auc_high','slope','offset_intercept','development_id_overlap'),names(d))]
names(d)[names(d)=='development_id_overlap']<-'overlap with development participants'
addspec('TableS32a','supplementary table S32a  CHARLS performance after relaxing baseline memory restrictions',d,
 'note: U=N+I+M, B=N+I; N baseline above low score boundary, I baseline low score, M missing baseline memory.score available sample does not require ADL complete, common sample requires ADL complete; direct comparison of two tools requires common sample.interval is DeLong 95%CI.education specific definition is post hoc exploration, low score defined strictly below boundary.',src)
src<-file.path(outdir,'F28_outcome_and_numerical_source_table.csv');d<-readc(src)
u<-d[grepl('^U_',d$scenario),]
stopifnot(nrow(u)>0)
u<-u[,c('cohort','window','scenario','model','population','n','events','AUC','AUC_low','AUC_high','calibration_offset_intercept','calibration_slope')]
names(u)[names(u)=='calibration_offset_intercept']<-'offset_intercept';names(u)[names(u)=='calibration_slope']<-'slope'
addspec('TableS32b','supplementary table S32b  external broader population S0 and S2 performance',u,
 'note: after removing baseline memory restriction, prediction target is follow up low score status.S0 use HRS/ELSA≤6, SHARE≤5; S2 all≤4.model and score use their respective calculable samples, continuous model calibration not applied to score.broader population AUC intervals retain saved DeLong results.',src)

src<-file.path(outdir,'F32_sample_flow.csv');d<-readc(src)
names(d)<-c('window','outcome_definition','matched wave participants','follow up memory measurable','score available n','common sample n','restored ADL participants with missing values','events among restored participants','development overlap in score sample','development overlap in common sample')
addspec('TableS33a','supplementary table S33a  CHARLS broader population enrollment flow',d,
 'note: all restored participants were excluded due to ADL missingness under old common sample restriction; do not use 0 to impute their ADL.score still requires age≥45 years old, three calculable items and available follow up outcome.',src,raw=TRUE)
src<-file.path(outdir,'F32_denominator_for_event_composition.csv');d<-readc(src)
d$baseline_low_fraction_of_ALL_events<-paste0(fmt(100*d$baseline_low_fraction_of_ALL_events,2),'%')
d$baseline_low_fraction_of_KNOWN_baseline_events<-paste0(fmt(100*d$baseline_low_fraction_of_KNOWN_baseline_events,2),'%')
names(d)<-c('window','outcome_definition','sample','U participant count','N participant count','I participant count','M participant count','U events','N events','I events','M events','I/all events','I/events with known baseline status','number of events with known baseline status')
addspec('TableS33b','supplementary table S33b  broader population event composition and two denominators',d,
 'note: I group events mean low scores at both measurements, does not establish persistently low scores between visits.M group events cannot be classified as incident crossing; I/all events and I/known baseline events have different denominators.',src,raw=TRUE)

screen_specs<-list(c('TableS34a','CHARLS broader population screening',file.path(root,'F32_broader_population_screening_metrics_corrected_version.csv')),
 c('TableS34b','external broader population S0 screening',file.path(outdir,'F16_screening.csv')),
 c('TableS34c','external broader population S2 screening',file.path(outdir,'F19_screening.csv')))
for(ss in screen_specs) {
  dd<-readc(ss[3])
  if(ss[1]!='TableS34a') {
    dd<-dd[,c('cohort_or_window','outcome_definition','target_screening_fraction','threshold','n','events','cases_detected','actual_screening_fraction','TP','FP','TN','FN','sensitivity','specificity','positive_predictive_value','negative_predictive_value','screen_positive_participants_with_known_baseline_status','screen_positive_participants_already_low_scoring_at_baseline','true_incident_events','per_detected_1_incident_case_number_needed_to_screen')]
    names(dd)[names(dd)=='true_incident_events']<-'confirmed threshold crossing TP'
    names(dd)[names(dd)=='per_detected_1_incident_case_number_needed_to_screen']<-'number scored/threshold crossing TP'
    addspec(ss[1],paste0('supplementary table ',sub('Table','',ss[1]),'  ',ss[2]),dd,
      'note: actual screening fraction, sensitivity, specificity and predictive values expressed in%.target fraction cutoff includes all ties, actual count may exceed target.event is low score at follow up; confirmed threshold crossing TP baseline must exceed corresponding boundary.',ss[3])
    next
  }
  # Explicit whitelist; obsolete specificity fields and historical metadata must not be displayed.
  keep<-intersect(c('window','definition','cohort','target','target_fraction','strategy','threshold','cut','n','events','selected','n_pos','actual_fraction','selected_fraction',
    'TP','FP','TN','FN','tp','fp','tn','fn','sensitivity','specificity','ppv','npv','PPV','NPV','TP_baseline_above','TP_baseline_low','TP_baseline_unknown',
    'detected_known_new','detected_baseline_low','detected_baseline_unknown','tp_N','tp_I','tp_M',
    'total_scored_per_detected_crossing','selected_per_detected_crossing'),names(dd))
  if(!all(c('sensitivity','specificity')%in%keep))stop('Screen fields need inspection')
  dd<-dd[,keep,drop=FALSE]
  dd<-data.frame(dd[,c('window','definition')],cutoff_selection_method=ifelse(is.na(dd$target_fraction),'Youden point','target fraction'),dd[,setdiff(names(dd),c('window','definition'))],check.names=FALSE)
  names(dd)[names(dd)=='sensitivity']<-'sensitivity'
  addspec(ss[1],paste0('supplementary table ',sub('Table','',ss[1]),'  ',ss[2]),dd,
   'note: specificity=TN/(TN+FP), NPV=TN/(TN+FN).proportions according to 0–1 listed; target fraction cutoff includes all ties, interpret using actual positive count/and fraction.Youden point selected in same sample, exploratory only.event is low score at follow up, not a current clinical cognitive disease diagnosis.',ss[3])
}

src<-file.path(outdir,'F26_BDSI_AUC.csv');d<-readc(src)
d<-d[d$scope=='65 to 79 years old' & d$stage=='paired sample with scorable items' & d$model%in%c('our three item score','BDSI item correction mapping'),]
if(!nrow(d)) {print(unique(readc(src)[,c('scope','stage')]));stop('BDSI filter')}
d$exclusion<-'development participants not excluded';d<-d[,c('outcome','exclusion','model','n','events','AUC','CI_low','CI_high')]
p27<-readc(file.path(outdir,'F27_AUC.csv'))
# F27 uses the same column names for model-level AUCs.
d<-rbind(d,p27[,names(d)])
addspec('TableS35a','supplementary table S35a  65～79 age mapping BDSI comparison: AUC',d,
 'note: CHARLS 2015—2018 year post hoc comparison.T1 education stratified, T2 is fixed≤4; excluding development participants alone is not independent external validation.original BDSI predicts six year dementia, table evaluates three year low memory score.age 65 age scores 0 points; functional item requires help with finances or medication; depression item is“everything was an effort”or antidepressant use; education remains a qualification proxy.complete sample retains original three arm additional eligibility criteria.',paste(src,'F27_AUC.csv'))
src<-file.path(outdir,'F27_AUC_paired_comparison.csv');d<-readc(src)
addspec('TableS35b','supplementary table S35b  paired analysis after excluding development participants AUC difference',d,
 'note: paired DeLong 95%CI; ΔAUC=our score minus mapped BDSI.P two sided without multiplicity adjustment; two exclusion sets are not independent validation samples.',src)
src<-file.path(outdir,'F27_sensitivity_at_equal_quota.csv');d<-readc(src)
d$bootstrap_reps<-NULL;d$seed<-NULL
addspec('TableS35c','supplementary table S35c  expected sensitivity at equal quota after development exclusions',d,
 'note: quota=floor(0.30×n); remaining quota allocated uniformly among cutoff ties, table reports expected sensitivity.differences in percentage points, interval from 2,000 paired participant resamples Bootstrap.random allocation not actually performed.T1 difference interval crosses 0, T2 does not cross 0.',src)
src<-file.path(outdir,'F26_BDSI_sensitivity_comparison_at_equal_quota.csv');d<-readc(src)
d$bootstrap_reps<-NULL;d$seed<-NULL
addspec('TableS35d','supplementary table S35d  before exclusion 65～79 age restricted equal quota comparison',d,
 'note: analysis within age range after corrected item mapping; still includes development participants.and table S35c interpret jointly, do not select higher AUC or smaller P value as independent validation.',src)
src<-file.path(outdir,'F27_screening_metrics.csv');d<-readc(src)
drop<-grep('^seed$|^bootstrap_reps$',names(d),value=TRUE);d[drop]<-NULL
names(d)[names(d)=='sensitivity']<-'sensitivity (%)'
names(d)[names(d)=='specificity']<-'specificity (%)'
names(d)[names(d)=='PPV']<-'PPV(%)';names(d)[names(d)=='NPV']<-'NPV(%)'
addspec('TableS35e','supplementary table S35e  screening strategies and actual burden after exclusion',d,
 'note: distinguish expected fixed quota, all ties included and retained F26 cutoff.actual positive fraction at fixed threshold does not equal 30%.proportion units follow column headings; expected contingency table may contain noninteger counts, not observed counts from implemented random allocation.',src)

src<-file.path(outdir,'F31_available_performance_source_table.csv');d<-readc(src)
for(task in c('F03','F05')) {
  dd<-d[d$task==task,];dd<-dd[,c('analysis','sample','outcome','model','n','events','development_id_overlap','auc','auc_low','auc_high','slope','offset_intercept','brier')]
  names(dd)[names(dd)=='development_id_overlap']<-'overlap with development participants'
  addspec(if(task=='F03')'TableS36a' else 'TableS36b',paste0('supplementary table ',if(task=='F03')'S36a  2018—2020 frozen, recalibration and refitting' else 'S36b  2018—2020 exploratory change outcomes'),dd,
   'note: post hoc later window analysis in same cohort.frozen predictions, same sample intercept updates and refits labeled separately; apparent recalibration improvement is not independent validation.within wave standardization does not guarantee measurement equivalence; do not infer from this MDC, true misclassification rate or performance ceiling.',src)
}
src<-file.path(outdir,'F31_available_performance_source_table.csv');d<-readc(src)
# Use the independently verified apparent values, never the legacy corrected columns.
d<-d[d$task=='F04',c('analysis','model','n','events','auc')]
names(d)<-c('model','notes','n','events','apparent_AUC')
d$apparent_AUC<-fmt(d$apparent_AUC,4)
addspec('TableS36c','supplementary table S36c  2018—2020 candidate model apparent performance AUC',d,
 'note: common complete sample 9,178 participants, 397 events.all models fitted and evaluated in same sample, total 12 models including reference.report verified apparent values only AUC, do not report historical corrected gains, difference intervals or“criterion met”determine; cannot conclude added variables have no value.',src,raw=TRUE)
