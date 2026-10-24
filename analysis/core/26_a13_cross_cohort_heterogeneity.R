# ============================================================================
# 26_a13_cross_cohort_heterogeneity.R —— A13: cross cohort synthesis and heterogeneity analysis
# input: A12_results_2026-09-26.rds, A11 sensitivity risk set, CHARLS model
# output: output/11_external_validation/A13_*, FigureS5–S7
# generate: 2026-09-26
# ============================================================================
invisible(try(Sys.setlocale("LC_ALL", "zh_CN.UTF-8"), silent = TRUE))
.libPaths(c("D:/R/library", .libPaths()))
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(pROC); library(knitr) })

P   <- "D:/BaiduSyncdisk/research/Xu_Zhonghao"
OUT <- file.path(P, "output", "11_external_validation")
FIG <- file.path(P, "output", "09_figures_and_tables")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)

res <- readRDS(file.path(OUT, "A12_results_2026-09-26.rds"))
perf <- res$perf; strata <- res$strata; meta <- res$meta; share_country_s0 <- res$share_country

a07 <- readRDS(file.path(P, "output/10_screening_index/screening_index_memory_v3_2026-09-26.rds"))
model_u <- a07$model_u; lev <- model_u$xlevels$edu_group
cohorts <- c("HRS","ELSA","SHARE")
cut_defs <- c("baseline_all_mean_sd","complete_fu_mean_sd","fixed_charls_le4")
scenario_lab <- c(baseline_all_mean_sd="S0_baseline_all", complete_fu_mean_sd="S1_complete_fu", fixed_charls_le4="S2_fixed_le4")
cohort_cols <- c(HRS="#1f77b4", ELSA="#2ca02c", SHARE="#d62728")

auc_ci <- function(y,x) {
  ok <- !is.na(y) & !is.na(x); y <- y[ok]; x <- x[ok]
  if (length(y) < 10 || length(unique(y)) < 2) return(c(auc=NA,low=NA,high=NA))
  r <- pROC::roc(y, x, direction = "<", quiet = TRUE); ci <- as.numeric(pROC::ci.auc(r))
  c(auc=as.numeric(pROC::auc(r)), low=ci[1], high=ci[3])
}
add_pred <- function(d, model) {
  d <- as.data.frame(d); d$edu_group <- factor(as.character(d$edu_group), levels = lev)
  d$p_hat <- as.numeric(predict(model, newdata = data.frame(age=d$age, edu_group=d$edu_group,
                                                            dep=d$dep, adl=d$adl), type = "response"))
  d
}
meta_dl <- function(auc, low, high, labels) {
  yi <- qlogis(auc); sei <- (qlogis(high) - qlogis(low)) / (2*1.96)
  wi <- 1/sei^2; fixed <- sum(wi*yi)/sum(wi)
  Q <- sum(wi*(yi-fixed)^2); df <- length(yi)-1
  C <- sum(wi) - sum(wi^2)/sum(wi); tau2 <- max(0,(Q-df)/C)
  ws <- 1/(sei^2+tau2); pooled <- sum(ws*yi)/sum(ws); se <- sqrt(1/sum(ws))
  I2 <- if (Q > 0) max(0,(Q-df)/Q)*100 else 0
  data.frame(k=length(yi), pooled_auc=plogis(pooled), pooled_low=plogis(pooled-1.96*se),
             pooled_high=plogis(pooled+1.96*se), tau2=tau2, Q=Q,
             p_Q=pchisq(Q, df, lower.tail=FALSE), I2=I2, cohorts=paste(labels, collapse=", "))
}

# ---------- cross cohort summary ----------
cross <- meta %>% mutate(scenario = scenario_lab[cut_definition]) %>%
  select(cut_definition, scenario, metric, k, pooled_auc, pooled_low, pooled_high,
         tau2, Q, p_Q, I2, cohorts) %>% arrange(factor(cut_definition, levels=cut_defs), metric)

# ---------- calibration heterogeneity ----------
cal_het <- perf %>% select(cohort, cut_definition, scenario, n_model, events_model,
                           cal_slope, cal_slope_low, cal_slope_high,
                           cal_intercept, cal_intercept_low, cal_intercept_high, brier) %>%
  arrange(factor(cut_definition, levels=cut_defs), match(cohort, cohorts))

# ---------- risk gradient ----------
risk_grad <- strata %>% left_join(perf %>% select(cohort, cut_definition, high_low_ratio),
                                  by = c("cohort","cut_definition")) %>%
  select(cohort, cut_definition, stratum, n, events, rate,
         n_model, events_model, pred_mean, oe_model, high_low_ratio) %>%
  mutate(scenario = scenario_lab[cut_definition]) %>% arrange(factor(cut_definition, levels=cut_defs), match(cohort, cohorts), stratum)

# ---------- SHARE country level: S0/S1/S2 + between countries meta ----------
share_rds <- readRDS(file.path(P,"external_data/SHARE/analysis/SHARE_sensitivity_risk_sets_2026-09-26.rds"))
share_country_all <- list(); share_country_meta <- list()
for (cd in cut_defs) {
  d <- add_pred(share_rds[[paste("SHARE", cd, sep="__")]], model_u)
  dm <- d[!is.na(d$p_hat), , drop = FALSE]
  # country level: model AUC use model complete cases; integer score AUC use full risk set
  ctry_model <- bind_rows(lapply(sort(unique(d$country_mod)), function(cc) {
    x <- dm[dm$country_mod == cc, , drop = FALSE]
    a <- auc_ci(x$event, x$p_hat)
    data.frame(cut_definition=cd, country_mod=cc, country_label=unique(x$country_label),
               n_model=nrow(x), events_model=sum(x$event), event_rate_model=mean(x$event),
               auc_model=a["auc"], auc_model_low=a["low"], auc_model_high=a["high"])
  }))
  ctry_score <- bind_rows(lapply(sort(unique(d$country_mod)), function(cc) {
    x <- d[d$country_mod == cc, , drop = FALSE]
    a <- auc_ci(x$event, x$score)
    data.frame(cut_definition=cd, country_mod=cc, country_label=unique(x$country_label),
               n_score=nrow(x), events_score=sum(x$event), event_rate_score=mean(x$event),
               auc_score=a["auc"], auc_score_low=a["low"], auc_score_high=a["high"])
  }))
  ctry <- full_join(ctry_model, ctry_score, by=c("cut_definition","country_mod","country_label"))
  share_country_all[[cd]] <- ctry
  m_model <- meta_dl(ctry$auc_model, ctry$auc_model_low, ctry$auc_model_high, ctry$country_label); m_model$cut_definition <- cd; m_model$metric <- "auc_model"
  m_score <- meta_dl(ctry$auc_score, ctry$auc_score_low, ctry$auc_score_high, ctry$country_label); m_score$cut_definition <- cd; m_score$metric <- "auc_score"
  share_country_meta[[cd]] <- bind_rows(m_model, m_score)
}
share_country_all <- bind_rows(share_country_all)
share_country_meta <- bind_rows(share_country_meta) %>% mutate(scenario = scenario_lab[cut_definition]) %>%
  select(cut_definition, scenario, metric, k, pooled_auc, pooled_low, pooled_high, tau2, Q, p_Q, I2, cohorts)

# ---------- CSV ----------
write.csv(cross, file.path(OUT,"A13_cross_cohort_summary.csv"), row.names=FALSE, fileEncoding="UTF-8")
write.csv(cal_het, file.path(OUT,"A13_calibration_heterogeneity.csv"), row.names=FALSE, fileEncoding="UTF-8")
write.csv(risk_grad, file.path(OUT,"A13_risk_gradient.csv"), row.names=FALSE, fileEncoding="UTF-8")
write.csv(share_country_all, file.path(OUT,"A13_SHARE_country_AUC.csv"), row.names=FALSE, fileEncoding="UTF-8")
write.csv(share_country_meta, file.path(OUT,"A13_SHARE_country_meta.csv"), row.names=FALSE, fileEncoding="UTF-8")
saveRDS(list(cross=cross, cal_het=cal_het, risk_grad=risk_grad,
             share_country_all=share_country_all, share_country_meta=share_country_meta),
        file.path(OUT,"A13_results_2026-09-26.rds"))

writeLines(capture.output(sessionInfo()), file.path(OUT,"A13_sessionInfo_2026-09-26.txt"), useBytes=TRUE)

# ---------- figure ----------
save_both <- function(base,width,height,draw){ png(paste0(base,".png"),width=width,height=height,units="in",res=300); draw(); dev.off()
  tiff(paste0(base,".tif"),width=width,height=height,units="in",res=600,compression="lzw"); draw(); dev.off() }
draw_score_forest <- function(){
  par(mfrow=c(1,3),mar=c(4,3.5,3,1))
  for(cd in cut_defs){
    a <- perf %>% filter(cut_definition==cd) %>% arrange(match(cohort,cohorts))
    m <- meta %>% filter(cut_definition==cd, metric=="score_full")
    est <- c(a$auc_score_full,m$pooled_auc); lo <- c(a$auc_score_full_low,m$pooled_low); hi <- c(a$auc_score_full_high,m$pooled_high)
    labels <- c(a$cohort,"Pooled")
    y <- rev(seq_along(labels))
    plot(NA,xlim=c(0.62,0.88),ylim=c(0.5,length(labels)+0.5),yaxt="n",xlab="AUC",ylab="",main=scenario_lab[[cd]],bty="n")
    abline(v=0.5,lty=3,col="grey70"); segments(lo,y,hi,y,col=c(cohort_cols[a$cohort],"black"),lwd=2)
    points(est,y,pch=19,col=c(cohort_cols[a$cohort],"black"),cex=c(1,1,1,1.3)); axis(2,at=y,labels=labels,las=1,cex.axis=0.8)
  }
}
draw_cal_het <- function(){
  par(mfrow=c(1,2),mar=c(5,4,3,1))
  d <- cal_het
  y <- seq_len(nrow(d))
  plot(NA,xlim=c(0.6,2.4),ylim=c(0.5,nrow(d)+0.5),yaxt="n",xlab="Calibration slope",ylab="",main="Calibration slopes",bty="n")
  abline(v=1,lty=3,col="grey70"); segments(d$cal_slope_low,y,d$cal_slope_high,y,col=cohort_cols[d$cohort],lwd=2)
  points(d$cal_slope,y,pch=19,col=cohort_cols[d$cohort],cex=0.9)
  axis(2,at=y,labels=paste(d$cohort,d$scenario,sep=" / "),las=1,cex.axis=0.7)
  plot(NA,xlim=c(-1.3,0.8),ylim=c(0.5,nrow(d)+0.5),yaxt="n",xlab="Calibration intercept (calibration-in-the-large)",ylab="",main="Calibration intercepts",bty="n")
  abline(v=0,lty=3,col="grey70"); segments(d$cal_intercept_low,y,d$cal_intercept_high,y,col=cohort_cols[d$cohort],lwd=2)
  points(d$cal_intercept,y,pch=19,col=cohort_cols[d$cohort],cex=0.9)
  axis(2,at=y,labels=paste(d$cohort,d$scenario,sep=" / "),las=1,cex.axis=0.7)
}
draw_share_meta <- function(){
  d <- share_country_all %>% filter(cut_definition=="baseline_all_mean_sd") %>% arrange(auc_model)
  y <- seq_len(nrow(d))
  par(mar=c(4,7,3,1))
  plot(NA,xlim=c(0.6,0.9),ylim=c(0.5,nrow(d)+1.5),yaxt="n",xlab="AUC (95% CI)",ylab="",main="SHARE country-level AUC (S0 model, country as unit)",bty="n")
  abline(v=0.5,lty=3,col="grey70")
  segments(d$auc_model_low,y,d$auc_model_high,y,col="#d62728",lwd=2); points(d$auc_model,y,pch=19,col="#d62728",cex=0.9)
  m <- share_country_meta %>% filter(cut_definition=="baseline_all_mean_sd",metric=="auc_model")
  yy <- nrow(d)+1
  segments(m$pooled_low,yy,m$pooled_high,yy,col="black",lwd=3); points(m$pooled_auc,yy,pch=18,col="black",cex=1.5)
  axis(2,at=c(y,yy),labels=c(d$country_label,"Pooled RE"),las=1,cex.axis=0.75)
  mtext(sprintf("Pooled RE AUC %.3f (%.3f-%.3f), I2=%.1f%%",m$pooled_auc,m$pooled_low,m$pooled_high,m$I2),side=3,cex=0.75)
}
save_both(file.path(FIG,"FigureS5_external_score_forest_2026-09-26"),12,4.2,draw_score_forest)
save_both(file.path(FIG,"FigureS6_external_calibration_heterogeneity_2026-09-26"),12,5.5,draw_cal_het)
save_both(file.path(FIG,"FigureS7_SHARE_country_meta_2026-09-26"),8,5.5,draw_share_meta)

# ---------- report ----------
tab_cross <- cross %>% mutate(scenario=scenario_lab[cut_definition], metric=ifelse(metric=="model","continuous model","integer score"),
                              pooled=sprintf("%.3f (%.3f–%.3f)",pooled_auc,pooled_low,pooled_high),
                              I2=sprintf("%.1f%%",I2)) %>% select(cutoff=scenario,metrics=metric,combine_AUC=pooled,I2,k,tau2,Q,p_Q)
tab_cal <- cal_het %>% mutate(scenario=scenario_lab[cut_definition],
                              slope=sprintf("%.2f (%.2f to %.2f)",cal_slope,cal_slope_low,cal_slope_high),
                              intercept=sprintf("%.2f (%.2f to %.2f)",cal_intercept,cal_intercept_low,cal_intercept_high),
                              Brier=sprintf("%.3f",brier)) %>% select(cohort=cohort,cutoff=scenario,model_n=n_model,events=events_model,calibration_slope=slope,calibration_intercept=intercept,Brier)
tab_share_meta <- share_country_meta %>% mutate(scenario=scenario_lab[cut_definition], metric=ifelse(metric=="auc_model","continuous model","integer score"),
                                                pooled=sprintf("%.3f (%.3f–%.3f)",pooled_auc,pooled_low,pooled_high),I2=sprintf("%.1f%%",I2)) %>%
  select(cutoff=scenario,metrics=metric,combine_AUC=pooled,I2,k,tau2,Q,p_Q)
tab_risk <- risk_grad %>% filter(cut_definition=="baseline_all_mean_sd") %>% mutate(stratification=as.character(stratum),observed_event_proportion=sprintf("%.1f%%",100*rate),predicted_risk=sprintf("%.1f%%",100*pred_mean),OE=sprintf("%.2f",oe_model)) %>%
  select(cohort=cohort,stratification,n,events=events,observed_event_proportion,predicted_risk,OE,high_to_low_risk_ratio=high_low_ratio)

md <- c(
"# A13: cross cohort synthesis and heterogeneity analysis(2026-09-26)",
"",
"> task incorporated into A12 formal external validation workflow; report focuses on cross cohort pooling, SHARE country and calibration heterogeneity.",
"> pooling method: DerSimonian–Laird random effects model, logit(AUC) scale; SE from each cohort/country 95%CI back calculate.",
"> sample definition: continuous model AUC uses model complete cases; integer score AUC uses full risk set; SHARE country level analyses also separated(n_model/n_score).",
"",
"## one, cross cohort pooling AUC",
"",
kable(tab_cross,format="pipe",row.names=FALSE,digits=4),
"",
"> S0 continuous model pooling AUC 0.737(0.690–0.779), I² 94.7%; S1 continuous model 0.725(0.689–0.758), I² 91.3%; S2 continuous model 0.777(0.744–0.808), I² 85.4%.",
"> very high cross cohort heterogeneity, pooled values are descriptive only; emphasize cohort results and sources of heterogeneity.additionally, cross cohort meta only k=3, I² itself insufficiently stable, pooled values cannot replace cohort results.",
"",
"## two, SHARE between country heterogeneity",
"",
kable(tab_share_meta,format="pipe",row.names=FALSE,digits=4),
"",
"> S0/S1 country level AUC range approximately 0.691–0.788, event rate approximately 5.4%–22.1%; S2 country level AUC range approximately 0.655–0.830, event rate approximately 3.5%–16.9%.between country values in table meta is**random effects pooling with country as analysis unit**(not participant level pooling; S0 continuous model 0.740, 95%CI 0.725–0.754), and participant level SHARE AUC 0.751(cluster bootstrap 0.734–0.771)different.between countries I²: S0/S1 approximately 50%(moderate), S2 approximately 77%(relatively high), notes SHARE not homogeneous populations; this meta reflects country level differences, not heterogeneity between independent studies.A12 already by country cluster bootstrap report CI, FigureS2/S7 provide country forest plot.",
"",
"## three, calibration heterogeneity",
"",
kable(tab_cal,format="pipe",row.names=FALSE,digits=4),
"",
"> HRS systematic risk underestimation(calibration intercept +0.33), SHARE systematic risk overestimation(intercept −0.31), ELSA calibration slope 1.82.calibration differences cannot be ignored, must accompany AUC report together.",
"",
"## four, risk gradient(S0)",
"",
kable(tab_risk,format="pipe",row.names=FALSE,digits=4),
"",
"> clear risk gradient reproduced in three cohorts: HRS high risk/low risk 4.9×, ELSA 10.5×, SHARE 7.3×.",
"",
"## five, figures tables and outputs",
"",
"- `A13_cross_cohort_summary.csv`",
"- `A13_calibration_heterogeneity.csv`",
"- `A13_risk_gradient.csv`",
"- `A13_SHARE_country_AUC.csv`, `A13_SHARE_country_meta.csv`",
"- `A13_results_2026-09-26.rds`",
"- `A13_independent_verification_2026_09_26.txt`",
"- `A13_sessionInfo_2026-09-26.txt`",
"- `FigureS5_external_score_forest_2026-09-26.png/.tif`(integer score forest plot)",
"- `FigureS6_external_calibration_heterogeneity_2026-09-26.png/.tif`(calibration heterogeneity)",
"- `FigureS7_SHARE_country_meta_2026-09-26.png/.tif`(SHARE country meta forest plot)",
"",
"## six, conclusions and writing notes",
"",
"1. cross cohort AUC similar point estimates, but random effects I² very high(85%–95%), and k=3, heterogeneity partly reflects threshold definitions/measurement/population differences, do not report only pooled values or attribute all differences to true effects; ",
"2. SHARE moderate to high country heterogeneity(S0/S1 I²≈50%, S2 I²≈77%), report country level results or cluster-robust CI; ",
"3. substantial between cohort calibration differences, HRS underestimation, SHARE overestimation, ELSA slope above one; ",
"4. S0/S1 and S2 threshold meanings differ, three complementary result sets; ",
"5. external cohorts are observational aging cohorts, not validated in surgical patients, not a postoperative outcome prediction tool."
)
writeLines(md,file.path(OUT,"A13_cross_cohort_synthesis_and_heterogeneity_2026_09_26.md"),useBytes=TRUE)
cat("A13 outputs written.\n")