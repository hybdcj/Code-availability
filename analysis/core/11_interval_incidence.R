# ============================================================================
# 11_interval_incidence.R —— T04：固定队列区间新发认知受损与 95%CI
# 版本：v1.0（2026-09-20）
# 输入：output/04_纵向队列/cohort_change_2026-09-13.rds
# 输出：
#   output/04_纵向队列/interval_incidence_2026-09-20.csv
#   output/04_纵向队列/interval_incidence_summary_2026-09-20.csv
#   output/04_纵向队列/interval_incidence_events_2026-09-20.rds
#   output/04_纵向队列/interval_incidence_2026-09-20.md
#
# 口径：
#   结局：impaired_fixed（2011 年教育分层固定切点，五波均使用同一阈值）
#   随访：2011→2013→2015→2018→2020，对应 2、2、3、2 年，总随访 9 年
#
#   分析 A（primary first episode）：
#       2011 认知正常者中，首次在任一随访波出现 impaired_fixed=1。
#       首次事件后退出风险集（即使之后恢复为正常，也不再次进入）。
#
#   分析 B（confirmed first episode，敏感性）：
#       基线正常者中，首次出现连续两波均受损。事件时间记为连续受损的第一波；
#       2020 年因没有下一波确认，不产生 new confirmed event。
#
#   分析 C（current-state interval incidence）：
#       每个区间开始时当前未受损者，到区间结束时受损；
#       允许既往受损后恢复者再次进入区间风险集，因此不是首次发病率。
#
#   分析 D（endpoint 2020）：
#       2011 正常者到 2020 时仍评估为受损。该指标减少了一过性跨阈值事件的干扰。
#
#   点估计：精确二项 95%CI（Clopper–Pearson）。
#   人年率：事件时间未知，主估计将事件置于区间中点；仅供近似，明确标注。
# ============================================================================

invisible(try(Sys.setlocale("LC_ALL", "zh_CN.UTF-8"), silent = TRUE))
if (dir.exists("D:/R/library")) .libPaths(c("D:/R/library", .libPaths()))
suppressPackageStartupMessages({ library(dplyr); library(tidyr) })

P <- Sys.getenv("CHARLS_PROJ", "D:/BaiduSyncdisk/科研/徐中豪")
OUT <- file.path(P, "output", "04_纵向队列")
if (!dir.exists(OUT)) dir.create(OUT, recursive = TRUE)

obj <- readRDS(file.path(OUT, "cohort_change_2026-09-13.rds"))
coh <- obj$cohort

stopifnot(all(c("KEY", "year", "impaired_fixed") %in% names(coh)))
stopifnot(setequal(sort(unique(coh$year)), c(2011, 2013, 2015, 2018, 2020)))
if (anyNA(coh$impaired_fixed)) stop("impaired_fixed 存在缺失，需先核查")

wide <- coh %>%
  select(KEY, year, impaired_fixed) %>%
  pivot_wider(names_from = year, values_from = impaired_fixed, names_prefix = "im") %>%
  arrange(KEY)

stopifnot(nrow(wide) == length(unique(coh$KEY)))

# ------------------------- 通用函数 -------------------------
binom_ci <- function(events, n) {
  if (is.na(n) || n <= 0) return(c(NA_real_, NA_real_))
  as.numeric(binom.test(events, n)$conf.int)
}

poisson_ci <- function(events, pyears) {
  if (is.na(pyears) || pyears <= 0 || is.na(events)) return(c(NA_real_, NA_real_))
  as.numeric(poisson.test(events, T = pyears)$conf.int)
}

fmt_pct <- function(x) ifelse(is.na(x), "—", sprintf("%.1f%%", 100 * x))
fmt_rate <- function(x) ifelse(is.na(x), "—", sprintf("%.1f", x))
fmt_ci_pct <- function(lo, hi) ifelse(is.na(lo) | is.na(hi), "—", sprintf("%.1f–%.1f", 100 * lo, 100 * hi))
fmt_ci_rate <- function(lo, hi) ifelse(is.na(lo) | is.na(hi), "—", sprintf("%.1f–%.1f", lo, hi))
fmt_int <- function(x) ifelse(is.na(x), "无法确认", format(as.integer(x), big.mark = ",", scientific = FALSE, trim = TRUE))

intervals <- tibble(
  interval = c("2011-2013", "2013-2015", "2015-2018", "2018-2020"),
  start_year = c(2011, 2013, 2015, 2018),
  end_year = c(2013, 2015, 2018, 2020),
  interval_years = c(2, 2, 3, 2),
  midpoint_years = c(1, 3, 5.5, 8)
)

# ------------------------- A. 首次发作 -------------------------
base <- wide %>% filter(im2011 == 0) %>% mutate(
  first_index = case_when(
    im2013 == 1 ~ 1L,
    im2015 == 1 ~ 2L,
    im2018 == 1 ~ 3L,
    im2020 == 1 ~ 4L,
    TRUE ~ NA_integer_
  ),
  first_event = as.integer(!is.na(first_index)),
  endpoint_2020 = im2020
)
baseline_n <- nrow(base)
first_events_total <- sum(base$first_event, na.rm = TRUE)
endpoint_events <- sum(base$endpoint_2020 == 1, na.rm = TRUE)

first_rows <- vector("list", nrow(intervals))
for (i in seq_len(nrow(intervals))) {
  at_risk <- sum(is.na(base$first_index) | base$first_index >= i)
  events <- sum(base$first_index == i, na.rm = TRUE)
  py <- (at_risk - events) * intervals$interval_years[i] + events * intervals$interval_years[i] / 2
  ci <- binom_ci(events, at_risk)
  rci <- poisson_ci(events, py)
  first_rows[[i]] <- tibble(
    analysis = "first_episode",
    analysis_label = "首次发作（首次受损后退出风险集）",
    interval = intervals$interval[i],
    start_year = intervals$start_year[i],
    end_year = intervals$end_year[i],
    interval_years = intervals$interval_years[i],
    at_risk_n = at_risk,
    events = events,
    incidence = events / at_risk,
    ci_lower = ci[1], ci_upper = ci[2],
    cumulative_events = NA_integer_,
    person_years_midpoint = py,
    rate_per_1000py = 1000 * events / py,
    rate_ci_lower = 1000 * rci[1], rate_ci_upper = 1000 * rci[2],
    note = "事件时间置于区间中点；首次事件后退出风险集"
  )
}
first_tab <- bind_rows(first_rows)
first_tab$cumulative_events <- cumsum(first_tab$events)
first_tab$cumulative_incidence <- first_tab$cumulative_events / baseline_n
first_ci <- t(vapply(seq_len(nrow(first_tab)), function(i) binom_ci(first_tab$cumulative_events[i], baseline_n), numeric(2)))
first_tab$cumulative_ci_lower <- first_ci[, 1]
first_tab$cumulative_ci_upper <- first_ci[, 2]

# ------------------------- B. 连续两波确认的首次发作 -------------------------
base <- base %>% mutate(
  confirmed_index = case_when(
    im2013 == 1 & im2015 == 1 ~ 1L,
    im2015 == 1 & im2018 == 1 ~ 2L,
    im2018 == 1 & im2020 == 1 ~ 3L,
    TRUE ~ NA_integer_
  ),
  confirmed_event = as.integer(!is.na(confirmed_index))
)
confirmed_total <- sum(base$confirmed_event, na.rm = TRUE)

confirmed_rows <- vector("list", nrow(intervals))
for (i in seq_len(nrow(intervals))) {
  confirmable <- i < nrow(intervals)   # 2018-2020 需要 2022 波才能确认
  if (confirmable) {
    at_risk <- sum(is.na(base$confirmed_index) | base$confirmed_index >= i)
    events <- sum(base$confirmed_index == i, na.rm = TRUE)
    py <- (at_risk - events) * intervals$interval_years[i] + events * intervals$interval_years[i] / 2
    ci <- binom_ci(events, at_risk)
    rci <- poisson_ci(events, py)
  } else {
    at_risk <- NA_integer_
    events <- NA_integer_
    py <- NA_real_
    ci <- c(NA_real_, NA_real_)
    rci <- c(NA_real_, NA_real_)
  }
  confirmed_rows[[i]] <- tibble(
    analysis = "confirmed_first_episode",
    analysis_label = "连续两波确认的首次受损（敏感性）",
    interval = intervals$interval[i],
    start_year = intervals$start_year[i],
    end_year = intervals$end_year[i],
    interval_years = intervals$interval_years[i],
    at_risk_n = at_risk,
    events = events,
    incidence = ifelse(is.na(events), NA_real_, events / at_risk),
    ci_lower = ci[1], ci_upper = ci[2],
    cumulative_events = NA_integer_,
    person_years_midpoint = py,
    rate_per_1000py = ifelse(is.na(events), NA_real_, 1000 * events / py),
    rate_ci_lower = ifelse(is.na(events), NA_real_, 1000 * rci[1]),
    rate_ci_upper = ifelse(is.na(events), NA_real_, 1000 * rci[2]),
    note = if (confirmable) "事件时间记为连续受损的第一波；首次确认后退出风险集" else "无法确认：需要 2022 年下一波数据"
  )
}
confirmed_tab <- bind_rows(confirmed_rows)
confirmed_tab$cumulative_events <- cumsum(ifelse(is.na(confirmed_tab$events), 0L, confirmed_tab$events))
confirmed_tab$cumulative_incidence <- confirmed_tab$cumulative_events / baseline_n
confirmed_ci <- t(vapply(seq_len(nrow(confirmed_tab)), function(i) binom_ci(confirmed_tab$cumulative_events[i], baseline_n), numeric(2)))
confirmed_tab$cumulative_ci_lower <- confirmed_ci[, 1]
confirmed_tab$cumulative_ci_upper <- confirmed_ci[, 2]

# ------------------------- C. 当前状态区间发病 -------------------------
current_rows <- vector("list", nrow(intervals))
start_vars <- c("im2011", "im2013", "im2015", "im2018")
end_vars <- c("im2013", "im2015", "im2018", "im2020")
for (i in seq_len(nrow(intervals))) {
  at_risk <- sum(wide[[start_vars[i]]] == 0, na.rm = TRUE)
  events <- sum(wide[[start_vars[i]]] == 0 & wide[[end_vars[i]]] == 1, na.rm = TRUE)
  py <- at_risk * intervals$interval_years[i]
  ci <- binom_ci(events, at_risk)
  rci <- poisson_ci(events, py)
  current_rows[[i]] <- tibble(
    analysis = "current_state",
    analysis_label = "当前状态区间新发（允许既往恢复者再次进入）",
    interval = intervals$interval[i],
    start_year = intervals$start_year[i],
    end_year = intervals$end_year[i],
    interval_years = intervals$interval_years[i],
    at_risk_n = at_risk,
    events = events,
    incidence = events / at_risk,
    ci_lower = ci[1], ci_upper = ci[2],
    cumulative_events = NA_integer_,
    person_years_midpoint = py,
    rate_per_1000py = 1000 * events / py,
    rate_ci_lower = 1000 * rci[1], rate_ci_upper = 1000 * rci[2],
    note = "分母为区间开始时当前未受损者，不是首次事件风险集"
  )
}
current_tab <- bind_rows(current_rows)

# ------------------------- D. 2020 终点 -------------------------
endpoint_ci <- binom_ci(endpoint_events, baseline_n)
endpoint_row <- tibble(
  analysis = "endpoint_2020",
  analysis_label = "2020 测评时受损（2011 正常者）",
  interval = "2011-2020",
  start_year = 2011,
  end_year = 2020,
  interval_years = 9,
  at_risk_n = baseline_n,
  events = endpoint_events,
  incidence = endpoint_events / baseline_n,
  ci_lower = endpoint_ci[1], ci_upper = endpoint_ci[2],
  cumulative_events = endpoint_events,
  cumulative_incidence = endpoint_events / baseline_n,
  cumulative_ci_lower = endpoint_ci[1], cumulative_ci_upper = endpoint_ci[2],
  person_years_midpoint = NA_real_,
  rate_per_1000py = NA_real_, rate_ci_lower = NA_real_, rate_ci_upper = NA_real_,
  note = "2020 年测评时为受损状态，减少对既往一过性事件的误读"
)

# 首次发作和连续确认的总随访结果，补充到 interval table
first_py_total <- sum(first_tab$person_years_midpoint, na.rm = TRUE)
first_rate_ci <- poisson_ci(first_events_total, first_py_total)
first_total_ci <- binom_ci(first_events_total, baseline_n)
first_total_row <- first_tab %>% slice_tail(n = 1) %>%
  mutate(interval = "2011-2020", start_year = 2011, end_year = 2020, interval_years = 9,
         at_risk_n = baseline_n, events = first_events_total,
         incidence = first_events_total / baseline_n,
         ci_lower = first_total_ci[1], ci_upper = first_total_ci[2],
         cumulative_events = first_events_total,
         cumulative_incidence = first_events_total / baseline_n,
         cumulative_ci_lower = first_total_ci[1],
         cumulative_ci_upper = first_total_ci[2],
         person_years_midpoint = first_py_total,
         rate_per_1000py = 1000 * first_events_total / first_py_total,
         rate_ci_lower = 1000 * first_rate_ci[1], rate_ci_upper = 1000 * first_rate_ci[2],
         note = "九年内至少一次受损；即使之后恢复也计为事件；事件时间区间中点近似")

confirmed_py_total <- sum(confirmed_tab$person_years_midpoint, na.rm = TRUE)
confirmed_rate_ci <- poisson_ci(confirmed_total, confirmed_py_total)
confirmed_total_ci <- binom_ci(confirmed_total, baseline_n)
confirmed_total_row <- confirmed_tab %>% slice_tail(n = 1) %>%
  mutate(interval = "2011-2020", start_year = 2011, end_year = 2020, interval_years = 9,
         at_risk_n = baseline_n, events = confirmed_total,
         incidence = confirmed_total / baseline_n,
         ci_lower = confirmed_total_ci[1], ci_upper = confirmed_total_ci[2],
         cumulative_events = confirmed_total,
         cumulative_incidence = confirmed_total / baseline_n,
         cumulative_ci_lower = confirmed_total_ci[1],
         cumulative_ci_upper = confirmed_total_ci[2],
         person_years_midpoint = confirmed_py_total,
         rate_per_1000py = 1000 * confirmed_total / confirmed_py_total,
         rate_ci_lower = 1000 * confirmed_rate_ci[1], rate_ci_upper = 1000 * confirmed_rate_ci[2],
         note = "连续两波均受损；总随访确认事件数")

results <- bind_rows(first_tab, first_total_row, confirmed_tab, confirmed_total_row, current_tab, endpoint_row) %>%
  mutate(
    incidence_pct = 100 * incidence,
    ci_lower_pct = 100 * ci_lower,
    ci_upper_pct = 100 * ci_upper,
    cumulative_incidence_pct = 100 * cumulative_incidence,
    cumulative_ci_lower_pct = 100 * cumulative_ci_lower,
    cumulative_ci_upper_pct = 100 * cumulative_ci_upper
  ) %>%
  select(analysis, analysis_label, interval, start_year, end_year, interval_years,
         at_risk_n, events, incidence, incidence_pct, ci_lower_pct, ci_upper_pct,
         cumulative_events, cumulative_incidence, cumulative_incidence_pct,
         cumulative_ci_lower_pct, cumulative_ci_upper_pct,
         person_years_midpoint, rate_per_1000py, rate_ci_lower, rate_ci_upper, note)

# ------------------------- 汇总表 -------------------------
summary_tab <- tibble(
  outcome = c("first_episode_any_wave", "confirmed_two_consecutive_waves", "endpoint_2020", "recovered_by_2020_after_first_event"),
  outcome_label = c("九年内至少一次受损", "连续两波确认受损", "2020 测评时受损", "首次事件后 2020 恢复为正常"),
  baseline_n = c(baseline_n, baseline_n, baseline_n, baseline_n),
  events = c(first_events_total, confirmed_total, endpoint_events,
             sum(base$first_event == 1 & base$endpoint_2020 == 0, na.rm = TRUE)),
  proportion = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_
)
for (i in seq_len(nrow(summary_tab))) {
  ci <- binom_ci(summary_tab$events[i], summary_tab$baseline_n[i])
  summary_tab$proportion[i] <- summary_tab$events[i] / summary_tab$baseline_n[i]
  summary_tab$ci_lower[i] <- ci[1]
  summary_tab$ci_upper[i] <- ci[2]
}
summary_tab <- summary_tab %>% mutate(
  proportion_pct = 100 * proportion,
  ci_lower_pct = 100 * ci_lower,
  ci_upper_pct = 100 * ci_upper
)

# ------------------------- 交叉核对 -------------------------
first_check <- first_events_total
endpoint_check <- endpoint_events
pre2020_first_events <- sum(base$first_index %in% 1:3, na.rm = TRUE)
recovered_pre2020 <- sum(base$first_index %in% 1:3 & base$endpoint_2020 == 0, na.rm = TRUE)
# 与既有队列的事件数核对由独立核验脚本完成，避免硬编码旧队列常量。

# ------------------------- 输出 -------------------------
write.csv(results, file.path(OUT, "interval_incidence_2026-09-20.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(summary_tab, file.path(OUT, "interval_incidence_summary_2026-09-20.csv"), row.names = FALSE, fileEncoding = "UTF-8")
saveRDS(list(base = base, wide = wide, intervals = intervals,
             first_event = first_tab, confirmed_event = confirmed_tab,
             current_state = current_tab, endpoint = endpoint_row,
             summary = summary_tab),
        file.path(OUT, "interval_incidence_events_2026-09-20.rds"))

report <- c(
  "# T04：固定队列区间新发认知受损（2026-09-20）",
  "",
  sprintf("> 数据：`cohort_change_2026-09-13.rds` ｜ 固定队列 n=%d ｜ 2011 基线正常 n=%d", length(unique(obj$cohort$KEY)), baseline_n),
  "> 结局：2011 年教育分层固定切点下的认知受损（impaired_fixed）",
  "> 随访：2011→2013→2015→2018→2020，总随访 9 年",
  "> **关键口径：** 主结局为五波唯一完全可比的 ms10（定向 5 题 + 连续减 7）教育分层固定切点，不是 2020 年单期论文中 0–31 分的 cog_B。",
  "",
  "## 一、必须先区分的三个指标",
  "",
  "| 指标 | 事件数 | 比例 | 95%CI | 含义 |",
  "|---|---:|---:|---:|---|",
  sprintf("| 九年内至少一次受损（首次发作） | %d | %s | %s | 任一随访波首次出现受损；之后恢复仍算事件 |",
          first_events_total, fmt_pct(first_events_total / baseline_n), fmt_ci_pct(binom_ci(first_events_total, baseline_n)[1], binom_ci(first_events_total, baseline_n)[2])),
  sprintf("| 连续两波确认受损 | %d | %s | %s | 降低一过性跨阈值影响 |",
          confirmed_total, fmt_pct(confirmed_total / baseline_n), fmt_ci_pct(binom_ci(confirmed_total, baseline_n)[1], binom_ci(confirmed_total, baseline_n)[2])),
  sprintf("| 2020 测评时受损 | %d | %s | %s | 2020 年测评时为受损状态 |",
          endpoint_events, fmt_pct(endpoint_events / baseline_n), fmt_ci_pct(endpoint_ci[1], endpoint_ci[2])),
  sprintf("| 首次发作后 2020 已恢复 | %d | %s | %s | 提示单次阈值跨越存在明显波动 |",
          sum(base$first_event == 1 & base$endpoint_2020 == 0, na.rm = TRUE),
          fmt_pct(mean(base$first_event == 1 & base$endpoint_2020 == 0, na.rm = TRUE)),
          fmt_ci_pct(binom_ci(sum(base$first_event == 1 & base$endpoint_2020 == 0, na.rm = TRUE), baseline_n)[1],
                     binom_ci(sum(base$first_event == 1 & base$endpoint_2020 == 0, na.rm = TRUE), baseline_n)[2])),
  "",
  sprintf("**结论：不能把 %s 直接解释为“十年持续发病率”。** %s 是九年内至少一次跨过受损阈值；%s 是 2020 年测评时处于受损状态；连续两波确认的比例更低。", fmt_pct(first_events_total / baseline_n), fmt_pct(first_events_total / baseline_n), fmt_pct(endpoint_events / baseline_n)),
  "",
  sprintf("> 在 2020 年前已发生首次事件的 %d 人中，到 2020 有 %d 人恢复为正常（%s）；该比例受不同随访时长影响，仅作描述性提示。", pre2020_first_events, recovered_pre2020, fmt_pct(recovered_pre2020 / pre2020_first_events)),
  "",
  "## 二、首次发作的区间条件发病",
  "",
  "| 区间 | 期初风险人数 | 新发事件 | 条件发病比例 | 95%CI | 累积事件 | 累积发病 |",
  "|---|---:|---:|---:|---:|---:|---:|",
  paste0("| ", first_tab$interval, " | ", first_tab$at_risk_n, " | ", first_tab$events, " | ",
         fmt_pct(first_tab$incidence), " | ", fmt_ci_pct(first_tab$ci_lower, first_tab$ci_upper), " | ",
         first_tab$cumulative_events, " | ", fmt_pct(first_tab$cumulative_incidence), " |"),
  "",
  "## 三、当前状态区间新发（敏感性口径）",
  "",
  "| 区间 | 期初当前正常人数 | 区间内新发 | 比例 | 95%CI |",
  "|---|---:|---:|---:|---:|",
  paste0("| ", current_tab$interval, " | ", current_tab$at_risk_n, " | ", current_tab$events, " | ",
         fmt_pct(current_tab$incidence), " | ", fmt_ci_pct(current_tab$ci_lower, current_tab$ci_upper), " |"),
  "",
  "> 该口径允许既往受损后恢复者再次进入风险集，不是首次事件发病率，不作为主结局。",
  "",
  "## 四、连续两波确认的首次受损",
  "",
  "| 区间 | 期初风险人数 | 确认事件 | 条件比例 | 95%CI | 累积确认事件 |",
  "|---|---:|---:|---:|---:|---:|",
  paste0("| ", confirmed_tab$interval, " | ", fmt_int(confirmed_tab$at_risk_n), " | ", fmt_int(confirmed_tab$events), " | ",
         fmt_pct(confirmed_tab$incidence), " | ", fmt_ci_pct(confirmed_tab$ci_lower, confirmed_tab$ci_upper), " | ",
         confirmed_tab$cumulative_events, " |"),
  "",
  "## 五、方法学提醒",
  "",
  "1. 事件发生时间只知道位于两次随访之间，区间人年率的点估计将事件置于区间中点，仅为近似。",
  "2. 首次发作后即退出风险集，即使后续恢复为正常，也不再次计入首次事件。",
  "3. 固定队列要求五波 ms10 均可计分，因此排除了失访、死亡和该指标缺失者，存在持续参与和幸存者偏倚。",
  "4. 本报告为未加权结果；基线权重敏感性见 T06，失访选择 IPW 见 T07。",
  "5. 2011 固定切点与逐波自身切点可能给出不同分类，逐波切点敏感性需另行计算。",
  "",
  "## 六、输出文件",
  "",
  "- `interval_incidence_2026-09-20.csv`",
  "- `interval_incidence_summary_2026-09-20.csv`",
  "- `interval_incidence_events_2026-09-20.rds`",
  "- `interval_incidence_2026-09-20.md`",
  "",
  sprintf("*生成时间：%s*", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
)

writeLines(report, file.path(OUT, "interval_incidence_2026-09-20.md"), useBytes = TRUE)

cat("\n========== T04 区间发病分析完成 ==========\n")
cat("2011 基线正常人数：", baseline_n, "\n", sep = "")
cat("九年内至少一次受损：", first_events_total, " (", fmt_pct(first_events_total / baseline_n), ")\n", sep = "")
cat("连续两波确认受损：", confirmed_total, " (", fmt_pct(confirmed_total / baseline_n), ")\n", sep = "")
cat("2020 测评时受损：", endpoint_events, " (", fmt_pct(endpoint_events / baseline_n), ")\n", sep = "")
cat("\n首次发作区间表：\n")
print(as.data.frame(first_tab %>% select(interval, at_risk_n, events, incidence, ci_lower, ci_upper, cumulative_incidence)), row.names = FALSE)
cat("\n已写：\n", file.path(OUT, "interval_incidence_2026-09-20.csv"), "\n",
    file.path(OUT, "interval_incidence_summary_2026-09-20.csv"), "\n",
    file.path(OUT, "interval_incidence_events_2026-09-20.rds"), "\n",
    file.path(OUT, "interval_incidence_2026-09-20.md"), "\n", sep = "")
