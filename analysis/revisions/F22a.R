# F26 LEGACY: BDSI实现/输入为历史版本；修正分析请运行F26_bdsi_audit.R。
# 原结果数值保留追溯，不可直接作为原版BDSI的有效比较；Tu限制另见F25。
# F22a：构造三臂对比数据（CHARLS 2015->2018，本指数 / BDSI / Tu2026）
invisible(try(Sys.setlocale("LC_ALL", "zh_CN.UTF-8"), silent = TRUE))
suppressPackageStartupMessages({ library(data.table); library(haven) })
if (dir.exists("D:/R/library")) .libPaths(c("D:/R/library", .libPaths()))
TMP <- "C:/y03tmp/tmp"
B_TU <- c(age_lt65 = -0.53, age_ge75 = 0.26, male = 0.38,
          edu_primid = 0.45, edu_none = 0.03, grip = -0.16, memory = -0.33, orient = -0.15)
LV <- c("Less than primary", "Primary school", "Middle school", "High school or above")
G3 <- c("A_none", "B_primid", "C_high")

w <- as.data.table(readRDS(file.path(TMP, "waves_cognition.rds")))
g <- w[year %in% c(2011, 2015, 2018), .(ID = IDn, year = as.numeric(year),
        mem = as.numeric(imrc) + as.numeric(dlrc), orient = as.numeric(orient),
        age = as.numeric(age), female = as.numeric(female), edu_group = edu_group,
        cesd10 = as.numeric(cesd10), adl = as.numeric(adl), iadl2 = as.numeric(iadl2),
        stroke = as.numeric(stroke), diab = as.numeric(diab))]
d15 <- g[year == 2015, .(ID, mem15 = mem, orient15 = orient, age15 = age, female15 = female,
                         edu15 = edu_group, cesd15 = cesd10, adl15 = adl, iadl2_15 = iadl2,
                         stroke15 = stroke, diab15 = diab)]
d18 <- g[year == 2018, .(ID, mem18 = mem)]

b <- as.data.table(read_dta(file.path(TMP, "w3_biomarker.dta"),
      col_select = any_of(c("ID", "qi002", "ql002", "qc003", "qc004", "qc005", "qc006"))))
b[, `:=`(height = fifelse(qi002 >= 120 & qi002 <= 220, qi002, NA_real_),
         weight = fifelse(ql002 >= 25 & ql002 <= 200, ql002, NA_real_))]
for (v in c("qc003", "qc004", "qc005", "qc006"))
  b[[v]] <- fifelse(b[[v]] > 0 & b[[v]] <= 100, b[[v]], NA_real_)
b[, `:=`(bmi = weight / (height / 100)^2, grip_max = pmax(qc003, qc004, qc005, qc006, na.rm = TRUE))]
b <- b[, .(ID = as.character(ID), bmi, grip_max)]

d11 <- g[year == 2011 & !is.na(mem) & !is.na(edu_group) & age >= 50]
d11[, g3 := fcase(edu_group == LV[1], G3[1], edu_group %in% LV[2:3], G3[2], edu_group == LV[4], G3[3])]
cut3 <- d11[!is.na(g3), .(cutoff = mean(mem) - sd(mem)), by = g3]

mm <- merge(merge(d15, d18, by = "ID"), b, by = "ID", all.x = TRUE)
mm[, g3 := fcase(edu15 == LV[1], G3[1], edu15 %in% LV[2:3], G3[2], edu15 == LV[4], G3[3])]
mm <- merge(mm, cut3[, .(g3, cutoff)], by = "g3", all.x = TRUE)
d <- mm[!is.na(mem15) & !is.na(mem18) & age15 >= 50 & !is.na(cutoff) & !is.na(bmi) &
          !is.na(grip_max) & !is.na(orient15) & !is.na(cesd15) & !is.na(iadl2_15) &
          !is.na(stroke15) & !is.na(diab15)]

d[, edu_f := factor(edu15, levels = LV)]
d[, male := as.integer(female15 == 0)]
d[, score_us := 2L * pmax(floor((age15 - 45) / 5), 0) +
      (-2L) * as.integer(edu_f == LV[2]) + (-7L) * as.integer(edu_f == LV[3]) +
      (-10L) * as.integer(edu_f == LV[4]) + 3L * as.integer(cesd15 >= 10)]
d[, `:=`(b_age = pmax(0, pmin(age15, 79) - 64),
         b_edu = 9L * as.integer(edu_f != LV[4]),
         b_stroke = 6L * as.integer(stroke15 == 1),
         b_diab = 3L * as.integer(diab15 == 1),
         b_bmi = 8L * as.integer(bmi < 18.5),
         b_func = 10L * as.integer(iadl2_15 == 1),
         b_dep = 6L * as.integer(cesd15 >= 10))]
d[, score_bdsi := b_age + b_edu + b_stroke + b_diab + b_bmi + b_func + b_dep]
d[, lp_tu := B_TU[["age_lt65"]] * as.integer(age15 < 65) + B_TU[["age_ge75"]] * as.integer(age15 >= 75) +
      B_TU[["male"]] * male + B_TU[["edu_primid"]] * as.integer(g3 == "B_primid") +
      B_TU[["edu_none"]] * as.integer(g3 == "A_none") + B_TU[["grip"]] * grip_max +
      B_TU[["memory"]] * mem15 + B_TU[["orient"]] * orient15]

saveRDS(d, file.path(TMP, "F22_d.rds"))
cat("三臂数据已保存：n =", nrow(d), "\n")
