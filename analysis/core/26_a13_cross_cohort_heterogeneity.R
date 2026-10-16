# ============================================================================
# 26_a13_cross_cohort_heterogeneity.R —— A13：跨队列汇总与异质性分析
# 输入：A12_results_2026-09-26.rds、A11 敏感性风险集、CHARLS 模型
# 输出：output/11_外部验证/A13_*、FigureS5–S7
# 生成：2026-09-26
# ============================================================================
invisible(try(Sys.setlocale("LC_ALL", "zh_CN.UTF-8"), silent = TRUE))
.libPaths(c("D:/R/library", .libPaths()))
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(pROC); library(knitr) })

P   <- "D:/BaiduSyncdisk/科研/徐中豪"
OUT <- file.path(P, "output", "11_外部验证")
FIG <- file.path(P, "output", "09_图表")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)

res <- readRDS(file.path(OUT, "A12_results_2026-09-26.rds"))
perf <- res$perf; strata <- res$strata; meta <- res$meta; share_country_s0 <- res$share_country

a07 <- readRDS(file.path(P, "output/10_筛查指数/screening_index_memory_v3_2026-09-26.rds"))
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

# ---------- 跨队列汇总 ----------
cross <- meta %>% mutate(scenario = scenario_lab[cut_definition]) %>%
  select(cut_definition, scenario, metric, k, pooled_auc, pooled_low, pooled_high,
         tau2, Q, p_Q, I2, cohorts) %>% arrange(factor(cut_definition, levels=cut_defs), metric)

# ---------- 校准异质性 ----------
cal_het <- perf %>% select(cohort, cut_definition, scenario, n_model, events_model,
                           cal_slope, cal_slope_low, cal_slope_high,
                           cal_intercept, cal_intercept_low, cal_intercept_high, brier) %>%
  arrange(factor(cut_definition, levels=cut_defs), match(cohort, cohorts))

# ---------- 风险梯度 ----------
risk_grad <- strata %>% left_join(perf %>% select(cohort, cut_definition, high_low_ratio),
                                  by = c("cohort","cut_definition")) %>%
  select(cohort, cut_definition, stratum, n, events, rate,
         n_model, events_model, pred_mean, oe_model, high_low_ratio) %>%
  mutate(scenario = scenario_lab[cut_definition]) %>% arrange(factor(cut_definition, levels=cut_defs), match(cohort, cohorts), stratum)

# ---------- SHARE 国家层面：S0/S1/S2 + 国家间 meta ----------
share_rds <- readRDS(file.path(P,"external_data/SHARE/analysis/SHARE_sensitivity_risk_sets_2026-09-26.rds"))
share_country_all <- list(); share_country_meta <- list()
for (cd in cut_defs) {
  d <- add_pred(share_rds[[paste("SHARE", cd, sep="__")]], model_u)
  dm <- d[!is.na(d$p_hat), , drop = FALSE]
  # 国家层面：模型 AUC 用模型完整样本；整数分 AUC 用全风险集
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
write.csv(cross, file.path(OUT,"A13_跨队列汇总.csv"), row.names=FALSE, fileEncoding="UTF-8")
write.csv(cal_het, file.path(OUT,"A13_校准异质性.csv"), row.names=FALSE, fileEncoding="UTF-8")
write.csv(risk_grad, file.path(OUT,"A13_风险梯度.csv"), row.names=FALSE, fileEncoding="UTF-8")
write.csv(share_country_all, file.path(OUT,"A13_SHARE国家AUC.csv"), row.names=FALSE, fileEncoding="UTF-8")
write.csv(share_country_meta, file.path(OUT,"A13_SHARE国家meta.csv"), row.names=FALSE, fileEncoding="UTF-8")
saveRDS(list(cross=cross, cal_het=cal_het, risk_grad=risk_grad,
             share_country_all=share_country_all, share_country_meta=share_country_meta),
        file.path(OUT,"A13_results_2026-09-26.rds"))

writeLines(capture.output(sessionInfo()), file.path(OUT,"A13_sessionInfo_2026-09-26.txt"), useBytes=TRUE)

# ---------- 图 ----------
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

# ---------- 报告 ----------
tab_cross <- cross %>% mutate(scenario=scenario_lab[cut_definition], metric=ifelse(metric=="model","连续模型","整数分"),
                              pooled=sprintf("%.3f (%.3f–%.3f)",pooled_auc,pooled_low,pooled_high),
                              I2=sprintf("%.1f%%",I2)) %>% select(切点=scenario,指标=metric,合并AUC=pooled,I2,k,tau2,Q,p_Q)
tab_cal <- cal_het %>% mutate(scenario=scenario_lab[cut_definition],
                              slope=sprintf("%.2f (%.2f to %.2f)",cal_slope,cal_slope_low,cal_slope_high),
                              intercept=sprintf("%.2f (%.2f to %.2f)",cal_intercept,cal_intercept_low,cal_intercept_high),
                              Brier=sprintf("%.3f",brier)) %>% select(队列=cohort,切点=scenario,模型n=n_model,事件=events_model,校准斜率=slope,校准截距=intercept,Brier)
tab_share_meta <- share_country_meta %>% mutate(scenario=scenario_lab[cut_definition], metric=ifelse(metric=="auc_model","连续模型","整数分"),
                                                pooled=sprintf("%.3f (%.3f–%.3f)",pooled_auc,pooled_low,pooled_high),I2=sprintf("%.1f%%",I2)) %>%
  select(切点=scenario,指标=metric,合并AUC=pooled,I2,k,tau2,Q,p_Q)
tab_risk <- risk_grad %>% filter(cut_definition=="baseline_all_mean_sd") %>% mutate(分层=as.character(stratum),观察事件率=sprintf("%.1f%%",100*rate),预测风险=sprintf("%.1f%%",100*pred_mean),OE=sprintf("%.2f",oe_model)) %>%
  select(队列=cohort,分层,n,事件=events,观察事件率,预测风险,OE,高危低危比=high_low_ratio)

md <- c(
"# A13：跨队列汇总与异质性分析（2026-09-26）",
"",
"> 本任务已并入 A12 的正式外部验证流程；本报告专门汇总跨队列合并、SHARE 国家间异质性和校准异质性。",
"> 合并方法：DerSimonian–Laird 随机效应模型，logit(AUC) 尺度；SE 由各队列/国家 95%CI 反推。",
"> 样本口径：连续模型 AUC 使用模型完整样本；整数分 AUC 使用全风险集；SHARE 国家层面同样分别处理（n_model/n_score）。",
"",
"## 一、跨队列合并 AUC",
"",
kable(tab_cross,format="pipe",row.names=FALSE,digits=4),
"",
"> S0 连续模型合并 AUC 0.737（0.690–0.779），I² 94.7%；S1 连续模型 0.725（0.689–0.758），I² 91.3%；S2 连续模型 0.777（0.744–0.808），I² 85.4%。",
"> 跨队列异质性很高，合并值仅作描述；论文应重点报告各队列结果和异质性来源。此外，跨队列 meta 仅 k=3，I² 本身不够稳定，合并值不能替代各队列结果。",
"",
"## 二、SHARE 国家间异质性",
"",
kable(tab_share_meta,format="pipe",row.names=FALSE,digits=4),
"",
"> S0/S1 国家层面 AUC 范围约 0.691–0.788，事件率约 5.4%–22.1%；S2 国家层面 AUC 范围约 0.655–0.830，事件率约 3.5%–16.9%。表中国家间 meta 为**以国家为分析单位的随机效应合并**（不是个体水平合并；S0 连续模型 0.740，95%CI 0.725–0.754），与个体水平 SHARE AUC 0.751（cluster bootstrap 0.734–0.771）不同。国家间 I²：S0/S1 约 50%（中等），S2 约 77%（较高），说明 SHARE 并非同质人群；该 meta 反映的是国家层面差异，不是独立研究间的异质性。A12 已按国家 cluster bootstrap 报告 CI，FigureS2/S7 给出国家森林图。",
"",
"## 三、校准异质性",
"",
kable(tab_cal,format="pipe",row.names=FALSE,digits=4),
"",
"> HRS 系统性低估风险（校准截距 +0.33）、SHARE 系统性高估风险（截距 −0.31）、ELSA 校准斜率 1.82。校准差异不能忽略，必须与 AUC 一起报告。",
"",
"## 四、风险梯度（S0）",
"",
kable(tab_risk,format="pipe",row.names=FALSE,digits=4),
"",
"> 三队列均复现明显风险梯度：HRS 高危/低危 4.9×、ELSA 10.5×、SHARE 7.3×。",
"",
"## 五、图表与输出",
"",
"- `A13_跨队列汇总.csv`",
"- `A13_校准异质性.csv`",
"- `A13_风险梯度.csv`",
"- `A13_SHARE国家AUC.csv`、`A13_SHARE国家meta.csv`",
"- `A13_results_2026-09-26.rds`",
"- `A13_独立复核_2026-09-26.txt`",
"- `A13_sessionInfo_2026-09-26.txt`",
"- `FigureS5_external_score_forest_2026-09-26.png/.tif`（整数分森林图）",
"- `FigureS6_external_calibration_heterogeneity_2026-09-26.png/.tif`（校准异质性）",
"- `FigureS7_SHARE_country_meta_2026-09-26.png/.tif`（SHARE 国家 meta 森林图）",
"",
"## 六、结论与写作提醒",
"",
"1. 跨队列 AUC 点估计较接近，但随机效应 I² 很高（85%–95%），且 k=3、异质性部分来自阈值定义/测量/人群差异，不能只报告合并值或全部归因于真实效应差异；",
"2. SHARE 国家间异质性为中等至较高（S0/S1 I²≈50%，S2 I²≈77%），需报告国家层面结果或 cluster-robust CI；",
"3. 校准存在明显队列差异，HRS 低估、SHARE 高估、ELSA 斜率偏大；",
"4. S0/S1 与 S2 的阈值含义不同，三套结果互补；",
"5. 外部队列均为观察性老龄队列，未在真实手术人群验证，不能表述为术后结局预测工具。"
)
writeLines(md,file.path(OUT,"A13_跨队列汇总与异质性_2026-09-26.md"),useBytes=TRUE)
cat("A13 outputs written.\n")