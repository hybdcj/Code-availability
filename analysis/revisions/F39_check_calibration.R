.libPaths(c("D:/R/library", .libPaths())); suppressPackageStartupMessages(library(dplyr))
options(digits=15)
b <- c(-4.19524508778947,0.0492039250378111,-0.290105771000328,-0.994711800202235,-1.43588388035287,0.412166484954086,-0.0110119344848014)
a07 <- readRDS("D:/BaiduSyncdisk/research/Xu_Zhonghao/output/10_screening_index/screening_index_memory_v3_2026-09-26.rds"); d <- a07$data
e <- as.character(d$edu_group)
lp <- b[1]+b[2]*d$age+b[3]*(e=="Primary school")+b[4]*(e=="Middle school")+b[5]*(e=="High school or above")+b[6]*d$dep+b[7]*d$adl
p <- plogis(lp); y <- d$incident_2015
cat("max|p - a07 p_hat| =", max(abs(p-d$p_hat)), "\n")
res <- do.call(rbind, lapply(c("all",levels(factor(d$edu_group)),"male","female"), function(lv){
  idx <- if(lv=="all") rep(TRUE,length(y)) else if(lv=="male") d$female==0 else if(lv=="female") d$female==1 else e==lv
  fi <- glm(y[idx]~offset(lp[idx]),family=binomial()); fs <- glm(y[idx]~lp[idx],family=binomial())
  data.frame(level=lv,n=sum(idx),observed=sum(y[idx]),expected=sum(p[idx]),O_minus_E=sum(y[idx])-sum(p[idx]),
    intercept=unname(coef(fi)[1]),converged=fi$converged,iter=fi$iter,slope=unname(coef(fs)[2]))}))
print(res)
write.csv(res,"D:/BaiduSyncdisk/research/Xu_Zhonghao/16.manuscript_revision_round_three/F39_TRIPOD_supplement_2026-10-17/F39_check_calibration_O_vs_E.csv",row.names=FALSE,fileEncoding="UTF-8")
note_fun <- function(f){x <- read.csv(f,fileEncoding="UTF-8",check.names=FALSE)
  x$cal_note <- ifelse(x$cohort=="CHARLS_dev" & x$subgroup %in% c("overall","education"),
    "by construction: apparent performance; education is a model term fitted in this sample (ML score equations force O=E), so intercept=0 (and overall slope=1)","")
  write.csv(x,f,row.names=FALSE,fileEncoding="UTF-8")}
for(f in c("D:/BaiduSyncdisk/research/Xu_Zhonghao/16.manuscript_revision_round_three/F39_TRIPOD_supplement_2026-10-17/F39_23a_CHARLS_dev_subgroups.csv","D:/BaiduSyncdisk/research/Xu_Zhonghao/16.manuscript_revision_round_three/F39_TRIPOD_supplement_2026-10-17/F39_23a_all_subgroups.csv")) note_fun(f)
