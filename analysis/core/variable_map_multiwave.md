# five wave variable mapping table(for longitudinal analysis)· v0.1

> creation date: 2026-09-13 ｜ status: **in progress**(2011/2013/2015 verified, 2018 to add, 2020 existing prepared outputs)
> basis: raw data for each wave `.dta` variable labels(scripts `code/_checks/map_cognition_vars.R`, output `output/02_five_wave_data_and_quality_control/cognition_vars_by_wave.csv`)+ Harmonized D codebook
> principle: **compared with first four waves, 2020 prepared wave output(`output/01_2020_single_wave_baseline/analysis_2020.rds`)definition as reference**; all five waves must use identical 9 component definitions

---

## one, conclusion: use"reconstruct from raw modules by wave"route(route B)

verified, **Harmonized D not directly suitable for this study cognitive definition**: 

| project | Harmonized D(r1–r4) | this study 2020 definition | availability |
|---|---|---|---|
| orientation | `r*orient` ＝"summary date naming", **0–4 points** | 5 items(year/month/day/day of week/season), **0–5 points** | ✗ different definition |
| immediate/delayed recall | `r*imrc` / `r*dlrc`(0–10) | same(0–10) | ✓ |
| serial subtraction of 7 | `r*ser7`(0–5) | same(0–5) | ✓ |
| figure copying | `r*draw`(0–1) | same(0–1) | ✓ |
| CES-D | `r*cesd10`(0–30) | same(0–30, ≥10 as boundary) | ✓ |
| ADL | `r*adla_c`(0–4) / `r*adlab_c`(0–6) / `r*adlfive`(0–5) | 6 items ADL | △ select consistent version |
| IADL | `r*iadla`(0–3, **only r2–r4**) | 5 items IADL(housework/preparing meals/shopping/managing finances/taking medication) | ⚠️ 2011 raw items available this wave(`db016`–`db020`)and Harmonized single item(`r1moneya` etc), but**no summary variable**(summary includes"making telephone calls", this item at wave 2 first added) |
| sleep duration / pain / physical activity | **absent in all** | present | ✗ obtain from raw modules |
| weights | `r*wtresp` | `INDV_weight` | △ can be used for verification |

→ therefore: **cognition and other exposures taken from raw modules using 2020 reconstruct definitions**; Harmonized D cross checking only(ADL/IADL/chronic conditions/education and weights).

## two, cognition 9 components: variable names by wave

### 2.1 2011(`wave1_2011/health_status_and_functioning.dta`)✓ verified

| components | variables | notes |
|---|---|---|
| orientation-year / month / day | `dc001s1` / `dc001s2` / `dc001s3` | label "Year / Month / Day Is Correct" |
| orientation-day of week | `dc002` | Checking Day of Week |
| orientation-season | `dc003` | Checking Season |
| serial subtraction of 7(5 step) | `dc019`(100−7), `dc020`–`dc023` | and 2020 of dc007_1..5 corresponding |
| plotting | `dc025` | Drew the Picture or Failed |
| immediate recall(trial 1) | `dc006s1`–`dc006s10`(+`dc006s11` indicator flags) | word list version `dc006_version` |
| delayed recall | `dc027s1`–`dc027s10`(+`dc027s11` indicator flags) | Recall Words |
| CES-D-10 | `dc009`–`dc018`(total 10 questions) | `dc015` = "Sleep Was Restless", and 2020 of `dc022` same source |

### 2.2 2013(`wave2_2013/Health_Status_and_Functioning.dta`)✓ verified

| components | variables | notes |
|---|---|---|
| orientation 5 items | `dc001s1`,`dc001s2`,`dc001s3`,`dc002`,`dc003` | Checking Date / Day of Week / Season |
| serial subtraction of 7 | `dc019`–`dc023` | 100−7 onward |
| plotting | `dc025` | Draw Picture |
| immediate recall | `dc006_1_s1`–`s11`(trial 1), `dc006_2_*`, `dc006_3_*`(trial 2, 3) | **three learning trials**; use trial 1 and 2020 comparable |
| delayed recall | `dc027s1`–`s11` | Recall Words |
| CES-D-10 | `dc009`–`dc018` | `dc015` = My Sleep Was Restless |

### 2.3 2015(`wave3_2015/Health_Status_and_Functioning.dta`)✓ verified

| components | variables | notes |
|---|---|---|
| orientation 5 items | `dc001s1`,`dc001s2`,`dc001s3`,`dc002`,`dc003` | checking date / day of week / season |
| serial subtraction of 7 | `dc019`–`dc023` | 100 minus 7 onward |
| plotting | `dc025` | Draw picture |
| immediate recall | `dc006s1`–`dc006s12`(including indicator flags) | word list `dc006_wordlist_1_`–`_10_` |
| delayed recall | `dc027s1`–`dc027s12` | Recall Words |
| CES-D-10 | `dc009`–`dc018` | `dc015` = My Sleep Was Restless |
| (additional)number series items | `dc029_w3`–`dc035_w3` | 2015 new number series items, not used in study |

### 2.4 2018(`wave4_2018/Cognition.dta`, 219 variables)✓ verified(**important differences**)

| components | variables | notes |
|---|---|---|
| orientation 5 items | `dc001_w4`(year), `dc006_w4`(month), `dc005_w4`(day/day of week), `dc003_w4`(date), `dc002_w4`(season) | label is Checking Year / Month / Day / Date / Season, **total 5 items**, direction consistent with other waves, but**numbering order completely different**(2020 yes dc001–dc005 order corresponds to year/season/day/day of week/month) |
| serial subtraction of 7 | `dc014_w4_1` … `dc014_w4_5`(answers in `_1_`) | 5 step, same structure as 2020 |
| **immediate recall(primary definition)** | `dc028_w4_s1`–`s10`(Butter/Arm/Shore/…, 10 words)+ `s11`(None) + `s12`(Refused) | ✅ same as other four waves **10 word list**; questionnaire PROCEDURE specifies 45 assessed in full sample above age(study uses trial 1; additional trial 2/3 = `dc029_w4_s*` / `dc030_w4_s*`) |
| **delayed recall(primary definition)** | `dc047_w4_s1`–`s10`(same set 10 words)+ `s11`(None) + `s12`(Refused) | ✅ administered to full sample(n=15,688) |
| immediate/delayed recall(**not used longitudinally**) | `dc013_w4_1_s1`–`s3`(Ball/Flag/Tree, 3 words, 5 trials), `dc015_w4_s1`–`s3` | ⚠️ only for **60 years and older**administered(n=9,989 / 8,984), belongs to MMSE subset, cannot combine with other waves 10 word lists |
| figure copying | `dc024_w4`(copy overlapping pentagons: 1=for / 5=incorrect / 97=not assessed) | ✅ 2018 **present**figure copying item, and administered to full sample(n=18,155) |
| additional item block | `dc016_w4`–`dc023_w4`(Watch / Pencil / Repeat / Read / Hand / Folds / Leg / Sentence Correct) | 2018 new instruction following items, not used in this study |

> 🔧 **2026-09-13 correction**: this section previously treated `dc013_w4_1` / `dc015_w4` mistaken for 2018 of"word recall", therefore concluded"2018 only 3 words, no figure copying item"conclusion.checked questionnaire PROCEDURE then changed to primary definition 10 word list + `dc024_w4` figure copying item.
> ⚠️ **remaining comparability considerations**: number 4 round word list**content**(butter/arm/shore…)differs from first three Chinese word lists, and immediate recall used **3 learning trials**(2011/2015 only 1 times), official note: number 4/5 wave delayed recall elevated, equipercentile equating recommended before longitudinal memory analysis(Gross 2012); therefore memory trajectories descriptive only, primary outcome is"orientation 5 items + serial subtraction of 7".
### 2.5 2020(`wave5_2020/Health_Status_and_Functioning.dta`)✓ completed

see `code/variable_map_2020.md` and `output/01_2020_single_wave_baseline/analysis_2020.rds`: orientation `dc001`–`dc005`, serial subtraction of 7 `dc007_1`–`dc007_5`(answers in `_1_`), figure copying `dc009`, immediate recall `dc012_s1`–`s10`(trial 1), delayed recall `dc028_s1`–`s10`, CES-D `dc016`–`dc025`.

## three, other variables(across waves)

| variables | 2011 | 2013 | 2015 | 2018 | 2020 |
|---|---|---|---|---|---|
| individual ID | 11 positions(supplement 1 positions = 12 positions) | 12 positions | 12 positions | 12 positions | 12 positions |
| weights | `weight.dta` of `ind_weight` | `INDV_weight` | `INDV_weight` | `INDV_weight` | `INDV_weight` |
| community identifier | `communityID` | `communityID` | `communityID` | `communityID` | `communityID` |
| demographics | `demographic_background.dta` | `Demographic_Background.dta` | same as left | same as left | same as left |
| ADL/IADL | original `db*` variables | same | same | same | `db001`–`db022` |
| chronic conditions | `da002_*` / `da003_*` | same | same | same | same |
| pain | `da027` | to verify | to verify | to verify | `da027` |
| sleep duration | `da030` | to verify | to verify | to verify | `da030` |
| physical activity | `da032_*` | to verify | to verify | to verify | `da032_1_`–`3_` |

> "to verify"items handled like cognition: compare raw variable labels item by item across waves(`code/_checks/map_cognition_vars.R` pattern can be reused).

## four, cross wave comparability: three definition issues requiring decisions

## part two of four, wave specific tables(2026-09-13 generate, scripts `code/03_waves_build.R`)

output: `output/02_five_wave_data_and_quality_control/waves_cognition_table_2026-09-13.csv`(table below), `output/02_five_wave_data_and_quality_control/waves_cognition_2026-09-13.rds`(individual longitudinal data)

| year | available sample | complete cognition | mean age | female | memory score (original) | memory score (proportion) | mental status | global cognition | prevalence B | prevalence A | prevalence C | ADL limited | IADL limited | depression | hypertension | diabetes | sleep (h) | pain | physical activity | figure copying present |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 2011 | 16,831 | 12,288 | 59.1 | 51.1% | 7.42 | 0.371 | 8.28 | 15.69 | 15.9% | 15.3% | 25.8% | 14.7% | — | 34.2% | 26.3% | 6.5% | 6.40 | 31.6% | 90.2% | present |
| 2013 | 17,835 | 13,687 | 59.9 | 51.4% | 7.54 | 0.377 | 8.22 | 15.76 | 14.2% | 15.9% | 25.3% | 15.6% | 14.4% | 30.2% | 27.9% | 7.8% | 6.21 | 34.5% | 91.0% | present |
| 2015 | 19,606 | 15,859 | 60.2 | 51.3% | 7.07 | 0.353 | 8.29 | 15.36 | 17.0% | 18.0% | 27.3% | 18.0% | 14.4% | 31.2% | 33.0% | 10.1% | 6.39 | 28.4% | 90.8% | present |
| 2018 | 19,384 | 11,653 | 61.7 | 52.3% | 8.77 | 0.439 | 8.32 | 17.10 | 14.8% | 15.4% | 28.8% | 13.9% | 10.5% | 33.2% | 36.6% | 13.0% | 6.25 | 58.8% | 93.7% | present |
| 2020 | 17,730 | 14,655 | 61.5 | 53.0% | 8.81 | 0.441 | 7.42 | 16.24 | **15.6%** | **16.4%** | **28.2%** | 19.5% | 19.8% | 36.0% | 38.7% | 14.8% | 6.09 | 58.2% | 91.7% | present |

### 4.0 2018 correct year definition(2026-09-13 correction)

2018 wave `Cognition.dta`(19,744 participants, 219 variables)contains two cognitive item sets, according to round 4 round questionnaire PROCEDURE distinguish: 

| item group | variables | target respondents | nonmissing |
|---|---|---|---|
| orientation 5 items | `dc001_w4`, `dc002_w4`, `dc003_w4`, `dc005_w4`, `dc006_w4` | 45 full sample above age | 18,158 |
| serial subtraction of 7 | `dc014_w4_1_1` … `dc014_w4_5_1` | full sample | 13,614 → 10,648 |
| figure copying(overlapping pentagons) | `dc024_w4` | full sample | 18,155 |
| **10 word list immediate recall(trial 1)** | `dc028_w4_s1`–`s10`(trial 2/3 = `dc029_w4_s*` / `dc030_w4_s*`) | full sample | 16,343 |
| **10 word list delayed recall** | `dc047_w4_s1`–`s10` | full sample | 15,688 |
| place orientation, 3 word list, instruction following items, word recognition and related items | `dc007`–`dc012_w4`, `dc013_w4_1`, `dc015_w4`, `dc016`–`dc027_w4`, `dc048`–`dc067_w4` | **only 60 years and older**(MMSE / HCAP subset) | 9,989 / 7,801 |

evidence: `data/docs/04_2018/CHARLS_2018_Household_Questionnaire.pdf` number 106 page; verification script `code/_checks/verify_2018_cognition_variables.R`.

**corrected fixed cohort size**(`code/_checks/cohort_sizes_after_fix.R`)

| outcome measure | fixed cohort with data at all five waves |
|---|---|
| memory score(10 word immediate trial 1 + delayed, 0–20) | **6,338 participants** |
| "orientation + serial subtraction of 7"(ms10, 0–10) | **5,672 participants**(2011≥45 age eligible with scorable primary longitudinal outcome) |
| "orientation + serial subtraction of 7 + figure copying"(0–11) | 5,553 participants(slightly more missing figure copying, still not primary outcome) |
| orientation only(5 items) | 9,954 participants |

### 4.1 cross wave memory comparability: different word lists and learning trial counts

- 2011 / 2013 / 2015 uses Chinese 10 word list(`dc006_wordlist_*` present 4 versions; 2013 perform 3 learning trials, 2011/2015 only perform 1 times); 
- 2018 uses another word list(Butter / Arm / Shore / …), 2020 also is 10 words; both waves are 3 learning trials; 
- CHARLS official documentation(`data/docs/general/CHARLS-HCAP-User-Guidance-cn.pdf`, `data/docs/04_2018/w2018-data-alert-cn.pdf`): number 4 round word list and administration differ from first three rounds, thus elevated delayed recall, longitudinal memory should first undergo**equipercentile equating**(Gross AL, et al. 2012).

three handling methods(results in `code/08_measurement_equivalence.R` and `output/05_measurement_equivalence_ratio/measurement_equivalence_note_2026-09-13.md`): 

| scheme | approach | slope per decade |
|---|---|---|
| A. prorated(including 2018) | number of recalled words ÷ word list length(0–1) | +0.038(0.033–0.042) |
| B. within wave z | transform each wave into z points | −0.173(−0.199–−0.147) |
| C. prorated(exclude 2018) | 2018 excluded from memory analysis | +0.054(0.049–0.059) |

> A/C all are"small increase", and official note of"number 4/5 wave elevation"consistent, **cannot interpret as memory improvement**; memory trajectories descriptive and supplementary only.

### 4.2 orientation item content and order differ

all five waves use 5 items(year/month/day/day of week/season), but 2018 variable numbering order and 2020 completely different(see 2.4), must**align by item labels**, do not map by numeric codes.

### 4.3 2011 wave IADL raw items available, but no Harmonized summary variable(2026-09-19 correction)

- **2011 question was asked at this wave IADL**: raw module `db016`–`db020` = housework / preparing meals / shopping / managing finances / taking medication, coding 1–4(1 no difficulty, 2 difficulty but able, 3 difficulty requiring help, 4 unable), approximately 17,548 respondents; Harmonized D corresponding individual items also present `r1housewka`, `r1mealsa`, `r1shopa`, `r1moneya`, `r1medsa`.
- but **Harmonized D not provided 2011 wave summary variable `r1iadla`**: official definition `RwIADLA = making telephone calls (PHONEA) + managing finances (MONEYA) + taking medication (MEDSA)`, whereas"making telephone calls"starting at 2 wave first asked(codebook: *RwIADLA and RwIADLZA are only available starting at Wave 2. At Wave 1, there is no question regarding using the phone.*).
- therefore this pipeline 2011 wave `iadl` code as missing, annotated in table"2011 wave excluded", **is not"not administered"**.
- if required 2011 wave IADL: available `db019`(managing finances)+ `db020`(taking medication)perform**partial two item definition**(and 2013+ three item summary not fully comparable), or construct original five item version(also not comparable with official summary), both approaches must disclose definition differences.

> raise these points at next supervisor meeting**proactively discuss**: explain identified wave differences and proposed handling, preferable to discovering them during review.

1. complete 2018 cognitive variable mapping(including pain/sleep/physical activity); 
2. write `code/03_waves_build.R`: construct one dataset per wave"same definition"analysis dataset(9 components + same covariate set + wave specific education cutoffs); 
3. produce wave specific tables(one per wave Table 1 style table + prevalence), also count"scorable cognition at all five waves"fixed cohort participant count; 
4. before supervisor delivery, definition differences(especially 9 components in 2011 and 2013+ differences in item counts)write an explanatory paragraph.
