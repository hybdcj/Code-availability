# ============================================================================
# 11_interval_incidence.R —— T04: fixed cohort interval incident cognitive impairment and 95%CI
# version: v1.0(2026-09-20)
# input: output/04_longitudinal_cohort/cohort_change_2026-09-13.rds
# output: 
#   output/04_longitudinal_cohort/interval_incidence_2026-09-20.csv
#   output/04_longitudinal_cohort/interval_incidence_summary_2026-09-20.csv
#   output/04_longitudinal_cohort/interval_incidence_events_2026-09-20.rds
#   output/04_longitudinal_cohort/interval_incidence_2026-09-20.md
#
# definition: 
#   outcome: impaired_fixed(2011 year education specific fixed cutoff, same threshold at all five waves)
#   follow up: 2011→2013→2015→2018→2020, corresponding 2, 2, 3, 2 year, total follow up 9 year
#
#   analysis A(primary first episode): 
#       2011 among cognitively unimpaired participants, first occurrence at any follow up wave impaired_fixed=1.
#       exit risk set after first event(even if subsequently returning to normal, do not reenter).
#
#   analysis B(confirmed first episode, sensitivity): 
#       among unimpaired baseline participants, first impairment at two consecutive waves.event time is the first of consecutive impaired waves; 
#       2020 year lacks next wave confirmation, does not produce new confirmed event.
#
#   analysis C(current-state interval incidence): 
#       currently unimpaired at interval start, impaired at interval end; 
#       previously impaired participants who recover may reenter interval risk set, therefore not first onset incidence.
#
#   analysis D(endpoint 2020): 
#       2011 unimpaired participants to 2020 still impaired at assessment.reduces transient threshold crossing effects.
#
#   point estimate: exact binomial 95%CI(Clopper–Pearson).
#   person year rate: unknown event time, primary estimate places event at interval midpoint; approximation only, label explicitly.
# ============================================================================

invisible(try(Sys.setlocale("LC_ALL", "zh_CN.UTF-8"), silent = TRUE))
if (dir.exists("D:/R/library")) .libPaths(c("D:/R/library", .libPaths()))
suppressPackageStartupMessages({ library(dplyr); library(tidyr) })

P <- Sys.getenv("CHARLS_PROJ", "D:/BaiduSyncdisk/research/Xu_Zhonghao")
OUT <- file.path(P, "output", "04_longitudinal_cohort")
if (!dir.exists(OUT)) dir.create(OUT, recursive = TRUE)

obj <- readRDS(file.path(OUT, "cohort_change_2026-09-13.rds"))
coh <- obj$cohort

stopifnot(all(c("KEY", "year", "impaired_fixed") %in% names(coh)))
stopifnot(setequal(sort(unique(coh$year)), c(2011, 2013, 2015, 2018, 2020)))
if (anyNA(coh$impaired_fixed)) stop("impaired_fixed missing values present, audit first")

wide <- coh %>%
  select(KEY, year, impaired_fixed) %>%
  pivot_wider(names_from = year, values_from = impaired_fixed, names_prefix = "im") %>%
  arrange(KEY)

stopifnot(nrow(wide) == length(unique(coh$KEY)))

# ------------------------- general helper functions -------------------------
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
fmt_int <- function(x) ifelse(is.na(x), "cannot confirm", format(as.integer(x), big.mark = ",", scientific = FALSE, trim = TRUE))

intervals <- tibble(
  interval = c("2011-2013", "2013-2015", "2015-2018", "2018-2020"),
  start_year = c(2011, 2013, 2015, 2018),
  end_year = c(2013, 2015, 2018, 2020),
  interval_years = c(2, 2, 3, 2),
  midpoint_years = c(1, 3, 5.5, 8)
)

# ------------------------- A. first onset -------------------------
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
    analysis_label = "first onset(exit risk set after first impairment)",
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
    note = "event time at interval midpoint; exit risk set after first event"
  )
}
first_tab <- bind_rows(first_rows)
first_tab$cumulative_events <- cumsum(first_tab$events)
first_tab$cumulative_incidence <- first_tab$cumulative_events / baseline_n
first_ci <- t(vapply(seq_len(nrow(first_tab)), function(i) binom_ci(first_tab$cumulative_events[i], baseline_n), numeric(2)))
first_tab$cumulative_ci_lower <- first_ci[, 1]
first_tab$cumulative_ci_upper <- first_ci[, 2]

# ------------------------- B. first onset confirmed at two consecutive waves -------------------------
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
  confirmable <- i < nrow(intervals)   # 2018-2020 required 2022 wave required for confirmation
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
    analysis_label = "first impairment confirmed at two consecutive waves(sensitivity)",
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
    note = if (confirmable) "event time is the first of consecutive impaired waves; exit risk set after first confirmed event" else "cannot confirm: required 2022 next wave data in year"
  )
}
confirmed_tab <- bind_rows(confirmed_rows)
confirmed_tab$cumulative_events <- cumsum(ifelse(is.na(confirmed_tab$events), 0L, confirmed_tab$events))
confirmed_tab$cumulative_incidence <- confirmed_tab$cumulative_events / baseline_n
confirmed_ci <- t(vapply(seq_len(nrow(confirmed_tab)), function(i) binom_ci(confirmed_tab$cumulative_events[i], baseline_n), numeric(2)))
confirmed_tab$cumulative_ci_lower <- confirmed_ci[, 1]
confirmed_tab$cumulative_ci_upper <- confirmed_ci[, 2]

# ------------------------- C. current status interval incidence -------------------------
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
    analysis_label = "interval incidence of current impairment(previously recovered participants may reenter)",
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
    note = "denominator is currently unimpaired at interval start, not a first event risk set"
  )
}
current_tab <- bind_rows(current_rows)

# ------------------------- D. 2020 endpoint -------------------------
endpoint_ci <- binom_ci(endpoint_events, baseline_n)
endpoint_row <- tibble(
  analysis = "endpoint_2020",
  analysis_label = "2020 impaired at assessment(2011 unimpaired participants)",
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
  note = "2020 impaired at assessment in year, reduces misinterpretation of previous transient events"
)

# total follow up results for first onset and consecutive confirmation, append to interval table
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
         note = "at least one impairment within nine years; event counted despite subsequent recovery; event time approximated by interval midpoint")

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
         note = "impaired at two consecutive waves; total confirmed events during follow up")

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

# ------------------------- summary table -------------------------
summary_tab <- tibble(
  outcome = c("first_episode_any_wave", "confirmed_two_consecutive_waves", "endpoint_2020", "recovered_by_2020_after_first_event"),
  outcome_label = c("at least one impairment within nine years", "impairment confirmed at two consecutive waves", "2020 impaired at assessment", "after first event 2020 recovered to unimpaired"),
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

# ------------------------- cross check -------------------------
first_check <- first_events_total
endpoint_check <- endpoint_events
pre2020_first_events <- sum(base$first_index %in% 1:3, na.rm = TRUE)
recovered_pre2020 <- sum(base$first_index %in% 1:3 & base$endpoint_2020 == 0, na.rm = TRUE)
# independent verification script checks event counts against existing cohort, avoid hard coding historical cohort constants.

# ------------------------- output -------------------------
write.csv(results, file.path(OUT, "interval_incidence_2026-09-20.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(summary_tab, file.path(OUT, "interval_incidence_summary_2026-09-20.csv"), row.names = FALSE, fileEncoding = "UTF-8")
saveRDS(list(base = base, wide = wide, intervals = intervals,
             first_event = first_tab, confirmed_event = confirmed_tab,
             current_state = current_tab, endpoint = endpoint_row,
             summary = summary_tab),
        file.path(OUT, "interval_incidence_events_2026-09-20.rds"))

report <- c(
  "# T04: fixed cohort interval incident cognitive impairment(2026-09-20)",
  "",
  sprintf("> data: `cohort_change_2026-09-13.rds` ｜ fixed cohort n=%d ｜ 2011 unimpaired at baseline n=%d", length(unique(obj$cohort$KEY)), baseline_n),
  "> outcome: 2011 year education specific fixed cutoff cognitive impairment(impaired_fixed)",
  "> follow up: 2011→2013→2015→2018→2020, total follow up 9 year",
  "> **key definitions: ** primary outcome is the only fully comparable five wave measure ms10(orientation 5 questions + serial subtraction of 7)education specific fixed cutoff, is not 2020 in year single wave paper 0–31 points of cog_B.",
  "",
  "## one, three measures that must be distinguished",
  "",
  "| metrics | event count | proportion | 95%CI | meaning |",
  "|---|---:|---:|---:|---|",
  sprintf("| at least one impairment within nine years(first onset) | %d | %s | %s | first impairment at any follow up wave; still counts if later recovered |",
          first_events_total, fmt_pct(first_events_total / baseline_n), fmt_ci_pct(binom_ci(first_events_total, baseline_n)[1], binom_ci(first_events_total, baseline_n)[2])),
  sprintf("| impairment confirmed at two consecutive waves | %d | %s | %s | reduces transient threshold crossing effects |",
          confirmed_total, fmt_pct(confirmed_total / baseline_n), fmt_ci_pct(binom_ci(confirmed_total, baseline_n)[1], binom_ci(confirmed_total, baseline_n)[2])),
  sprintf("| 2020 impaired at assessment | %d | %s | %s | 2020 impaired at assessment in year |",
          endpoint_events, fmt_pct(endpoint_events / baseline_n), fmt_ci_pct(endpoint_ci[1], endpoint_ci[2])),
  sprintf("| after first onset 2020 recovered | %d | %s | %s | indicates substantial fluctuation in single threshold crossings |",
          sum(base$first_event == 1 & base$endpoint_2020 == 0, na.rm = TRUE),
          fmt_pct(mean(base$first_event == 1 & base$endpoint_2020 == 0, na.rm = TRUE)),
          fmt_ci_pct(binom_ci(sum(base$first_event == 1 & base$endpoint_2020 == 0, na.rm = TRUE), baseline_n)[1],
                     binom_ci(sum(base$first_event == 1 & base$endpoint_2020 == 0, na.rm = TRUE), baseline_n)[2])),
  "",
  sprintf("**conclusion: must not interpret %s directly as“ten year persistent incidence”.** %s is at least one impairment threshold crossing within nine years; %s yes 2020 impaired at assessment in year; two wave confirmed proportion is lower.", fmt_pct(first_events_total / baseline_n), fmt_pct(first_events_total / baseline_n), fmt_pct(endpoint_events / baseline_n)),
  "",
  sprintf("> at 2020 with first event before year %d among participants, to 2020 present %d participants recovered to unimpaired(%s); affected by varying follow up lengths, descriptive indication only.", pre2020_first_events, recovered_pre2020, fmt_pct(recovered_pre2020 / pre2020_first_events)),
  "",
  "## two, interval conditional incidence of first onset",
  "",
  "| interval | number at risk at interval start | incident events | conditional incidence proportion | 95%CI | cumulative events | cumulative incidence |",
  "|---|---:|---:|---:|---:|---:|---:|",
  paste0("| ", first_tab$interval, " | ", first_tab$at_risk_n, " | ", first_tab$events, " | ",
         fmt_pct(first_tab$incidence), " | ", fmt_ci_pct(first_tab$ci_lower, first_tab$ci_upper), " | ",
         first_tab$cumulative_events, " | ", fmt_pct(first_tab$cumulative_incidence), " |"),
  "",
  "## three, interval incidence of current impairment(sensitivity definition)",
  "",
  "| interval | currently unimpaired count at interval start | new events within interval | proportion | 95%CI |",
  "|---|---:|---:|---:|---:|",
  paste0("| ", current_tab$interval, " | ", current_tab$at_risk_n, " | ", current_tab$events, " | ",
         fmt_pct(current_tab$incidence), " | ", fmt_ci_pct(current_tab$ci_lower, current_tab$ci_upper), " |"),
  "",
  "> definition permits recovered participants to reenter risk set, not first event incidence, not a primary outcome.",
  "",
  "## four, first impairment confirmed at two consecutive waves",
  "",
  "| interval | number at risk at interval start | confirmed events | conditional proportion | 95%CI | cumulative confirmed events |",
  "|---|---:|---:|---:|---:|---:|",
  paste0("| ", confirmed_tab$interval, " | ", fmt_int(confirmed_tab$at_risk_n), " | ", fmt_int(confirmed_tab$events), " | ",
         fmt_pct(confirmed_tab$incidence), " | ", fmt_ci_pct(confirmed_tab$ci_lower, confirmed_tab$ci_upper), " | ",
         confirmed_tab$cumulative_events, " |"),
  "",
  "## five, methodological note",
  "",
  "1. event known only to occur between follow up visits, person year point estimate places event at interval midpoint, approximate only.",
  "2. exit risk set immediately after first onset, even with subsequent recovery, not counted again for first events.",
  "3. fixed cohort requires all five waves ms10 all scorable, therefore excludes attrition, deaths and missing measurements, subject to continued participation and survivor bias.",
  "4. report presents unweighted results; baseline weight sensitivity in T06, attrition selection IPW see T07.",
  "5. 2011 fixed and wave specific cutoffs may classify differently, wave specific cutoff sensitivity requires separate calculation.",
  "",
  "## six, output file",
  "",
  "- `interval_incidence_2026-09-20.csv`",
  "- `interval_incidence_summary_2026-09-20.csv`",
  "- `interval_incidence_events_2026-09-20.rds`",
  "- `interval_incidence_2026-09-20.md`",
  "",
  sprintf("*generation time: %s*", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
)

writeLines(report, file.path(OUT, "interval_incidence_2026-09-20.md"), useBytes = TRUE)

cat("\n========== T04 interval incidence analysis complete ==========\n")
cat("2011 number unimpaired at baseline: ", baseline_n, "\n", sep = "")
cat("at least one impairment within nine years: ", first_events_total, " (", fmt_pct(first_events_total / baseline_n), ")\n", sep = "")
cat("impairment confirmed at two consecutive waves: ", confirmed_total, " (", fmt_pct(confirmed_total / baseline_n), ")\n", sep = "")
cat("2020 impaired at assessment: ", endpoint_events, " (", fmt_pct(endpoint_events / baseline_n), ")\n", sep = "")
cat("\n first onset interval table: \n")
print(as.data.frame(first_tab %>% select(interval, at_risk_n, events, incidence, ci_lower, ci_upper, cumulative_incidence)), row.names = FALSE)
cat("\n written: \n", file.path(OUT, "interval_incidence_2026-09-20.csv"), "\n",
    file.path(OUT, "interval_incidence_summary_2026-09-20.csv"), "\n",
    file.path(OUT, "interval_incidence_events_2026-09-20.rds"), "\n",
    file.path(OUT, "interval_incidence_2026-09-20.md"), "\n", sep = "")
