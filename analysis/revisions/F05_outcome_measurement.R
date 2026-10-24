# ============================================================================
# F05: outcome measurement issues——extent of misclassification, performance after definition changes, what the model actually predicts
# data: waves_cognition(2018 and 2020 memory score, immediate recall trial 1 + delayed, 0–20)
# output: C:/y03tmp/out/F05_report.txt, F05_*.csv, F05_*.png
# notes: all post hoc analyses; definitions prespecified, no post hoc selection
# ============================================================================
invisible(try(Sys.setlocale("LC_ALL", "zh_CN.UTF-8"), silent = TRUE))
suppressPackageStartupMessages({ library(data.table); library(pROC) })
if (dir.exists("D:/R/library")) .libPaths(c("D:/R/library", .libPaths()))
TMP <- "C:/y03tmp/tmp"; OUT <- "C:/y03tmp/out"
con <- file(file.path(OUT, "F05_report.txt"), open = "wt", encoding = "UTF-8")
say <- function(...) { s <- paste0(..., collapse = ""); cat(s, "\n"); cat(s, "\n", file = con) }
pt  <- function(x) say(paste(capture.output(print(x)), collapse = "\n"))

BETA_OLD <- c("(Intercept)" = -4.19524508778947, age = 0.0492039250378111,
              eduP = -0.290105771000328, eduM = -0.994711800202235, eduH = -1.43588388035287,
              dep = 0.412166484954086, adl = -0.0110119344848014)
LV <- c("Less than primary", "Primary school", "Middle school", "High school or above")

## ---------------- data ----------------
w <- as.data.table(readRDS(file.path(TMP, "waves_cognition.rds")))
a18 <- w[year == 2018, .(ID = IDn, mem18 = as.numeric(imrc) + as.numeric(dlrc),
                         age18 = as.numeric(age), edu18 = edu_group,
                         cesd18 = as.numeric(cesd10), adl18 = as.numeric(adl),
                         wt18 = as.numeric(weight))]
a20 <- w[year == 2020, .(ID = IDn, mem20 = as.numeric(imrc) + as.numeric(dlrc))]
m <- merge(a18, a20, by = "ID")
pt_paired <- m[!is.na(mem18) & !is.na(mem20)]
say("paired sample n = ", nrow(pt_paired))

## within wave standardization(use paired sample within wave mean/SD)
mu18 <- mean(pt_paired$mem18); sd18 <- sd(pt_paired$mem18)
mu20 <- mean(pt_paired$mem20); sd20 <- sd(pt_paired$mem20)
pt_paired[, z18 := (mem18 - mu18)/sd18][, z20 := (mem20 - mu20)/sd20]
pt_paired[, dz := z20 - z18][, draw := mem20 - mem18]
say("2018: mean=", round(mu18, 2), " SD=", round(sd18, 2),
    "｜2020: mean=", round(mu20, 2), " SD=", round(sd20, 2))
say("change score(0–20 points): mean=", round(mean(pt_paired$draw), 2),
    " SD=", round(sd(pt_paired$draw), 2))
say("standardized change dz: mean=", round(mean(pt_paired$dz), 3), " SD=", round(sd(pt_paired$dz), 3))

## ---------------- 1. correlation between waves, reliability and minimal detectable change ----------------
say("\n========== 1. measurement reliability and minimal detectable change ==========")
r <- cor(pt_paired$mem18, pt_paired$mem20)
sem <- sd18 * sqrt(1 - r)
mdc <- 1.96 * sqrt(2) * sem
say("correlation between waves r = ", round(r, 3), "(note: 2018/2020 word list changed, this value will be lower, ")
say("      by ρ≈r thus overestimates measurement error, conservative direction)")
say("estimate SEM = ", round(sem, 2), " points; minimal detectable change MDC95 = ", round(mdc, 2), " points")

## ---------------- 2. transition matrix and movements around cutoff ----------------
say("\n========== 2. transition matrix(paired sample)==========")
cut_pop <- floor(mu18 - sd18)
pt_paired[, g18 := cut(mem18, c(-Inf, cut_pop, 6, 9, Inf),
                       labels = c(paste0("≤", cut_pop), "4-6", "7-9", "≥10"))]
pt_paired[, g20 := cut(mem20, c(-Inf, cut_pop, 6, 9, Inf),
                       labels = c(paste0("≤", cut_pop), "4-6", "7-9", "≥10"))]
tt <- dcast(pt_paired[, .(n = .N), by = .(g18, g20)], g18 ~ g20, value.var = "n", fill = 0)
pt(tt)
say("\n by 2018 score groups, 2020 threshold crossing(≤", cut_pop, ")proportion of: ")
say("(list only 2018 unimpaired participants)")
sub <- pt_paired[mem18 > cut_pop]
sub[, g := cut(mem18, c(cut_pop, 6, 9, Inf), labels = c("4-6(immediately above cutoff)", "7-9", "≥10"))]
pt(sub[, .(n = .N, events = sum(mem20 <= cut_pop),
           crossing_rate = paste0(round(100*mean(mem20 <= cut_pop), 1), "%")), by = g][order(g)])
say("→ near cutoff crossing rate relative to high score group ",
    round(sub[mem18 %in% 4:6, mean(mem20 <= cut_pop)] / sub[mem18 >= 10, mean(mem20 <= cut_pop)], 1),
    " fold")

## decline magnitude among events: how many exceed minimal detectable change
ev <- pt_paired[mem18 > cut_pop & mem20 <= cut_pop]
say("\n threshold crossing events n = ", nrow(ev),
    "; of which decline magnitude ≥ MDC95(", round(mdc, 2), " points)participants ",
    sum(ev$draw <= -mdc), " cases(", round(100*mean(ev$draw <= -mdc), 1), "%)")
say("   median decline = ", round(median(ev$draw), 1), " points; ≥5 scoring participants ",
    sum(ev$draw <= -5), " cases(", round(100*mean(ev$draw <= -5), 1), "%)")
say("→ remaining crossing declines compatible with measurement fluctuation(includes regression to mean)")
say("\n comparison: 2018 among unimpaired participants, median score change in nonevents = ",
    round(median(sub[mem20 > cut_pop]$draw), 1), " points")

## ---------------- 3. model performance under different outcome definitions ----------------
say("\n========== 3. model performance under different outcome definitions(frozen vs reestimate in identical participants)==========")
## m already includes covariates; add directly onto it z score and change, avoid merge produce .x/.y suffix
mm <- m[!is.na(mem18) & !is.na(mem20)]
mm[, `:=`(z18 = (mem18 - mu18)/sd18, z20 = (mem20 - mu20)/sd20)]
mm[, `:=`(dz = z20 - z18, draw = mem20 - mem18)]
mm_full <- copy(mm)
mm <- mm[!is.na(age18) & !is.na(edu18) & !is.na(cesd18) & !is.na(adl18) & age18 >= 45]
mm_full <- mm_full[!is.na(age18) & !is.na(edu18) & !is.na(cesd18) & !is.na(adl18) & age18 >= 45]
mm[, dep := as.integer(cesd18 >= 10)]
mm[, edu_group := factor(edu18, levels = LV)]
mm[, p_frozen := plogis(BETA_OLD[["(Intercept)"]] + BETA_OLD[["age"]]*age18 +
       BETA_OLD[["dep"]]*dep + BETA_OLD[["adl"]]*adl18 +
       BETA_OLD[["eduP"]]*as.integer(edu_group == LV[2]) +
       BETA_OLD[["eduM"]]*as.integer(edu_group == LV[3]) +
       BETA_OLD[["eduH"]]*as.integer(edu_group == LV[4]))]

## equipercentile anchored 2020 threshold
p_anchor <- mean(pt_paired$mem18 <= cut_pop)
cut20p <- as.numeric(quantile(pt_paired$mem20, p_anchor, type = 1))
say("O2 equipercentile anchoring: 2018 ≤", cut_pop, " corresponding percentile ", round(100*p_anchor, 1),
    "% → 2020 cutoff ≤", cut20p)

mm[, o1 := as.integer(mem20 <= cut_pop)]
mm[, o2 := as.integer(mem20 <= cut20p)]
mm[, o3 := as.integer(dz <= -1.0)]
mm[, o4 := as.integer(dz <= -0.5)]
mm[, o6 := as.integer(mem20 <= cut_pop & draw <= -5)]
mm_sub <- mm[mem18 > cut_pop]                 # population A: 2018 unimpaired participants
mm_risk <- copy(mm_sub)                       # population for threshold crossing outcomes
mm_full2 <- mm                               # population B: all paired participants(for decline outcomes, avoid truncation)

aucv <- function(y, p) if (length(unique(y)) < 2) NA_real_ else
  as.numeric(auc(roc(y, p, quiet = TRUE, direction = "<")))
eval_def <- function(col, tag, data = mm_risk) {
  d <- copy(data)[, event := get(col)]
  fit <- glm(event ~ age18 + edu_group + dep + adl18, data = d, family = binomial())
  d[, p_new := predict(fit, newdata = d, type = "response")]
  cal <- glm(event ~ qlogis(pmin(pmax(p_frozen, 1e-9), 1 - 1e-9)), data = d, family = binomial())
  a_old <- aucv(d$event, d$p_frozen); a_new <- aucv(d$event, d$p_new)
  data.table(definition = tag, population = if (nrow(d) == nrow(mm_risk)) "A unimpaired participants" else "B all paired participants",
             n = nrow(d), events = sum(d$event),
             event_rate = round(100*mean(d$event), 2),
             frozen_AUC = round(a_old, 4), reestimate_AUC = round(a_new, 4),
             ΔAUC = round(a_new - a_old, 4),
             frozen_calibration_slope = round(coef(cal)[2], 3), frozen_calibration_intercept = round(coef(cal)[1], 3),
             Brier_frozen = round(mean((d$p_frozen - d$event)^2), 4),
             Brier_reestimate = round(mean((d$p_new - d$event)^2), 4))
}
defs <- rbind(
  eval_def("o1", paste0("O1 absolute cutoff ≤", cut_pop, "(current status)"), mm_risk),
  eval_def("o2", "O2 equipercentile anchoring", mm_risk),
  eval_def("o6", "O6 threshold crossing plus decline ≥5 points", mm_risk),
  eval_def("o3", "O3 within wave standardized decline ≥1 SD", mm_full2),
  eval_def("o4", "O4 within wave standardized decline ≥0.5 SD", mm_full2),
  eval_def("o3", "O3 within wave standardized decline ≥1 SD(restricted to unimpaired participants)", mm_risk),
  eval_def("o4", "O4 within wave standardized decline ≥0.5 SD(restricted to unimpaired participants)", mm_risk))
pt(defs)
fwrite(defs, file.path(OUT, "F05_outcome_definitions.csv"), bom = TRUE)

## ---------------- 4. what the model actually predicts(R²) ----------------
say("\n========== 4. what four variables explain: level vs change ==========")
r2 <- function(form, data) summary(lm(form, data = data))$r.squared
mm_full[, `:=`(dep = as.integer(cesd18 >= 10),
               edu_group = factor(edu18, levels = LV))]
mm_sub <- copy(mm_full)[mem18 > cut_pop]
r2tab <- data.table(
  sample = rep(c("2018 unimpaired participants", "all paired participants"), each = 3),
  outcome = rep(c("2020 memory score(level)", "2 year change(0–20 points)", "2 year standardized change dz"), 2),
  R2 = c(r2(mem20 ~ age18 + edu_group + dep + adl18, mm_sub),
         r2(draw ~ age18 + edu_group + dep + adl18, mm_sub),
         r2(dz ~ age18 + edu_group + dep + adl18, mm_sub),
         r2(mem20 ~ age18 + edu_group + dep + adl18, mm_full),
         r2(draw ~ age18 + edu_group + dep + adl18, mm_full),
         r2(dz ~ age18 + edu_group + dep + adl18, mm_full)))
r2tab[, R2 := round(R2, 4)]
pt(r2tab)
fwrite(r2tab, file.path(OUT, "F05_R2.csv"), bom = TRUE)

say("\n four variables for[2020 level]and[2 year change]linear coefficient of(same scale 0–20 points, all paired participants): ")
lv <- summary(lm(mem20 ~ age18 + edu_group + dep + adl18, data = mm_full2))$coefficients
ch <- summary(lm(draw ~ age18 + edu_group + dep + adl18, data = mm_full2))$coefficients
cm <- data.table(items = rownames(lv), for_2020_level = round(lv[, 1], 3), for_2_year_change = round(ch[, 1], 3))
pt(cm)
fwrite(cm, file.path(OUT, "F05_coefficients_level_vs_change.csv"), bom = TRUE)

## ---------------- figure ----------------
png(file.path(OUT, "F05_fig_change_dist.png"), width = 1150, height = 560, res = 120)
par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 1.5))
hist(pt_paired$draw, breaks = seq(-20.5, 20.5, 1), col = "grey80", border = "white",
     xlab = "Change in memory score 2018 to 2020", ylab = "Count",
     main = "Distribution of 2-year change (n = 12,571)")
abline(v = 0, lty = 3); abline(v = -mdc, col = "red", lwd = 2)
legend("topleft", c("MDC95", "no change"), col = c("red", "black"), lwd = c(2, 1),
       lty = c(1, 3), bty = "n", cex = 0.8)
lvl <- c(paste0("≤", cut_pop), "4-6", "7-9", "≥10")
bandmat <- cbind(w2018 = as.integer(table(factor(pt_paired$g18, levels = lvl))),
                 w2020 = as.integer(table(factor(pt_paired$g20, levels = lvl))))
barplot(t(bandmat), beside = TRUE, col = c("grey60", "grey25"), xlab = "Memory score band",
        ylab = "Count", main = "Band distribution: 2018 vs 2020", las = 1)
legend("topright", c("2018", "2020"), fill = c("grey60", "grey25"), bty = "n", cex = 0.8)
dev.off()

close(con)
cat("\n complete ->", file.path(OUT, "F05_report.txt"), "\n")
