# ============================================================================
# 21_internal_validation_figures.R —— A09: internal validation figures and tables(ROC/calibration/DCA/stratification)
# basis: SAP_v3 §8(calibration plot, DCA threshold 5–40%)and §11(Figure 3/4)
# model: incident_2015 ~ age + edu_group + dep + adl(options 2)
# generate: 2026-09-26 ｜ output: output/09_figures_and_tables/Figure3_*, Figure4_*(png 300dpi + tif 600dpi)
#        supplementary table: output/10_screening_index/A09_*(calibration quantile table, DCA table)
# ============================================================================
invisible(try(Sys.setlocale("LC_ALL", "zh_CN.UTF-8"), silent = TRUE))
.libPaths(c("D:/R/library", .libPaths()))
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(pROC) })
P <- "D:/BaiduSyncdisk/research/Xu_Zhonghao"; OUT <- file.path(P, "output", "04_longitudinal_cohort")
DIR <- file.path(P, "output", "10_screening_index"); FIG <- file.path(P, "output", "09_figures_and_tables")
if (!dir.exists(FIG)) dir.create(FIG, recursive = TRUE)

coh <- readRDS(file.path(OUT, "cohort_change_2026-09-13.rds"))$cohort
w <- coh %>% select(KEY, year, imrc, dlrc) %>% pivot_wider(names_from = year, values_from = c(imrc, dlrc)) %>%
  mutate(mem11 = imrc_2011 + dlrc_2011, mem15 = imrc_2015 + dlrc_2015)
b <- coh %>% filter(year == 2011) %>% transmute(KEY, age, edu_group, dep = as.integer(cesd10 >= 10), adl) %>% left_join(w, by = "KEY")
b$edu_group <- factor(b$edu_group, levels = c("Less than primary","Primary school","Middle school","High school or above"))
d <- b[!is.na(b$mem11) & !is.na(b$mem15) & b$mem11 > 4 & complete.cases(b[, c("age","edu_group","dep","adl")]), ]
d$inc <- as.integer(d$mem15 <= 4)
m <- glm(inc ~ age + edu_group + dep + adl, data = d, family = binomial())
d$p <- predict(m, type = "response")
lp <- qlogis(pmin(pmax(d$p, 1e-6), 1 - 1e-6))
cal <- glm(inc ~ lp, data = d, family = binomial())
slope <- coef(cal)[2]; intercept <- coef(cal)[1]
auc_v <- as.numeric(auc(roc(d$inc, d$p, quiet = TRUE)))

## ---------- calibration quantile table(10 quantile)----------
d$dec <- cut(d$p, quantile(d$p, seq(0, 1, 0.1)), include.lowest = TRUE, labels = FALSE)
cal_tab <- d %>% group_by(dec) %>% summarise(n = n(), mean_prediction = mean(p), events = sum(inc),
  observed_rate = mean(inc), .groups = "drop") %>%
  mutate(observed_CI_low = mapply(function(k, n) binom.test(k, n)$conf.int[1], events, n),
         observed_CI_high = mapply(function(k, n) binom.test(k, n)$conf.int[2], events, n))
write.csv(cal_tab, file.path(DIR, "A09_calibration_deciles_2026-09-26.csv"), row.names = FALSE, fileEncoding = "UTF-8")

## ---------- DCA(net benefit, threshold 5–40%)----------
prev <- mean(d$inc); n <- nrow(d)
dca_tab <- do.call(rbind, lapply(seq(0.05, 0.40, 0.01), function(pt) {
  pos <- d$p >= pt
  TP <- sum(pos & d$inc == 1); FP <- sum(pos & d$inc == 0)
  nb_model <- TP/n - FP/n * (pt/(1-pt))
  nb_all   <- prev - (1-prev) * (pt/(1-pt))
  data.frame(threshold = pt, positive_rate = mean(pos), net_benefit_model = nb_model, net_benefit_intervene_in_all = nb_all, net_benefit_no_intervention = 0)
}))
write.csv(dca_tab, file.path(DIR, "A09_dca_2026-09-26.csv"), row.names = FALSE, fileEncoding = "UTF-8")

## ---------- risk stratification(**integer total score**, consistent with tool: SAP v3.1 rules)----------
bb <- coef(m); b_age5 <- 5*bb[["age"]]
tb <- c(age_per_5_years_label = b_age5, primary_school_education = bb[["edu_groupPrimary school"]], middle_school_education = bb[["edu_groupMiddle school"]],
        high_school_education = bb[["edu_groupHigh school or above"]], depression = bb[["dep"]], ADL = bb[["adl"]])
pts <- round(10*tb/max(abs(tb)))
d$score <- as.integer(pts[["age_per_5_years_label"]])*pmax(floor((d$age-45)/5),0) +
  as.integer(pts[["primary_school_education"]])*as.integer(d$edu_group=="Primary school") +
  as.integer(pts[["middle_school_education"]])*as.integer(d$edu_group=="Middle school") +
  as.integer(pts[["high_school_education"]])*as.integer(d$edu_group=="High school or above") +
  as.integer(pts[["depression"]])*d$dep + as.integer(pts[["ADL"]])*d$adl
q <- quantile(d$score, c(1/3, 2/3))   # integer total score tertiles(and A07 consistent)
d$stratum <- cut(d$score, c(-Inf, q[1], q[2], Inf), labels = c("low","in","high"), right = TRUE)
strata <- d %>% group_by(stratum) %>% summarise(n = n(), events = sum(inc), observed_rate = mean(inc),
  predicted_rate = mean(p), mean_score = round(mean(score),1), .groups = "drop")
cat("integer total score tertile cutoffs: ", q, "\n")
write.csv(strata, file.path(DIR, "A09_strata_observed_vs_predicted_2026-09-26.csv"), row.names = FALSE, fileEncoding = "UTF-8")

## ---------- Figure 4: ROC / calibration / DCA ----------
draw4 <- function() {
  par(mfrow = c(1, 3), mar = c(5, 5, 3, 2), cex.lab = 1.2, cex.axis = 1.0)
  # (a) ROC
  ro <- roc(d$inc, d$p, quiet = TRUE)
  plot(1 - ro$specificities, ro$sensitivities, type = "l", lwd = 2.5, col = "#2C7FB8",
       xlab = "1 − specificity", ylab = "sensitivity", main = "(a) ROC curve", xlim = c(0,1), ylim = c(0,1))
  abline(0, 1, lty = 3, col = "grey50")
  text(0.62, 0.18, sprintf("AUC = %.3f\n(95%%CI %.3f–%.3f)", auc_v,
       as.numeric(ci.auc(ro)[1]), as.numeric(ci.auc(ro)[3])), cex = 1.0)
  # (b) calibration plot
  plot(cal_tab$mean_prediction, cal_tab$observed_rate, pch = 19, cex = 1.2, col = "#31A354",
       xlim = c(0, max(cal_tab$mean_prediction, cal_tab$observed_CI_high)*1.05),
       ylim = c(0, max(cal_tab$observed_CI_high)*1.05),
       xlab = "predicted probability", ylab = "observed_event_rate", main = "(b) calibration plot(deciles)")
  segments(cal_tab$mean_prediction, cal_tab$observed_CI_low, cal_tab$mean_prediction, cal_tab$observed_CI_high, col = "#31A354")
  abline(0, 1, lty = 3, col = "grey50")
  lines(lowess(cal_tab$mean_prediction, cal_tab$observed_rate), col = "#756BB1", lwd = 1.5)
  text(max(cal_tab$mean_prediction)*0.45, max(cal_tab$observed_CI_high)*0.95,
       sprintf("calibration slope = %.2f\n intercept = %.2f", slope, intercept), cex = 1.0, adj = 0)
  # (c) DCA
  plot(dca_tab$threshold, dca_tab$net_benefit_model, type = "l", lwd = 2.5, col = "#2C7FB8",
       ylim = c(min(dca_tab$net_benefit_model, -0.02), max(dca_tab$net_benefit_model, prev)*1.05),
       xlab = "threshold probability", ylab = "net_benefit", main = "(c) decision curve analysis")
  lines(dca_tab$threshold, dca_tab$net_benefit_intervene_in_all, lty = 2, col = "#E6550D", lwd = 2)
  abline(h = 0, lty = 3, col = "grey50")
  legend("topright", c("screening index", "intervene in all", "no intervention"), lty = c(1,2,3), lwd = c(2.5,2,1),
         col = c("#2C7FB8","#E6550D","grey50"), bty = "n", cex = 0.95)
}
png(file.path(FIG, "Figure4_internal_validation_2026-09-26.png"), width = 18, height = 6, units = "in", res = 300)
draw4(); dev.off()
tiff(file.path(FIG, "Figure4_internal_validation_2026-09-26.tif"), width = 18, height = 6, units = "in", res = 600, compression = "lzw")
draw4(); dev.off()

## ---------- Figure 3: risk stratification ----------
draw3 <- function() {
  par(mar = c(5, 5, 3, 2), cex.lab = 1.2)
  bp <- barplot(strata$observed_rate, names.arg = c("low risk","intermediate risk","high risk"), col = c("#C6DBEF","#6BAED6","#2171B5"),
                ylim = c(0, max(strata$observed_rate)*1.3), ylab = "4 year incident cognitive impairment rate", xlab = "screening index risk strata",
                main = "screening index risk stratification(development cohort)")
  text(bp, strata$observed_rate, sprintf("%.1f%%\n(n=%d)", 100*strata$observed_rate, strata$n), pos = 3, cex = 1.0)
  points(bp, strata$predicted_rate, pch = 18, cex = 1.6, col = "red")
  legend("topleft", c("observed_event_rate", "model_predicted_rate"), pch = c(15, 18), col = c("#2171B5","red"), bty = "n")
}
png(file.path(FIG, "Figure3_risk_stratification_2026-09-26.png"), width = 8, height = 6, units = "in", res = 300)
draw3(); dev.off()
tiff(file.path(FIG, "Figure3_risk_stratification_2026-09-26.tif"), width = 8, height = 6, units = "in", res = 600, compression = "lzw")
draw3(); dev.off()

cat("== calibration decile table ==\n"); print(as.data.frame(cal_tab), row.names = FALSE)
cat("\n== DCA key thresholds ==\n"); print(dca_tab[seq(1, 36, 5), ], row.names = FALSE)
cat("\n== stratification(observed vs prediction)==\n"); print(as.data.frame(strata), row.names = FALSE)
cat("\n figures saved: ", FIG, "\n")
