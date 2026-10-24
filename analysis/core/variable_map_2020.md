# CHARLS 2020(Wave 5)variable mapping table(Paper 1 using)

> basis: `data/docs/05_2020/CHARLS_2020_Codebook.pdf`, `data/docs/05_2020/CHARLS_2020_Questionnaire_Chinese.pdf`, 
> official `data/raw/harmonized_d/H_CHARLS_D_do_file.do`(1–4 wave construction code, for checking scoring definitions).
> compilation date: 2026-09-10 ｜ latest revision: **2026-09-12**(combine education variables + recomputed numbers)

## 0. important prerequisite: 2020 wave variable codes and 1–4 wave codes differ completely

official Harmonized construction code shows, 2011–2018 wave cognitive variables are `dc006`(immediate recall), `dc027`(delayed recall), 
`dc019–dc023`(serial subtraction of 7), `dc025`(figure copying); **2020 all wave codes renumbered**.
therefore **2020 wave variables must follow 2020 codebook construct independently, do not apply old mapping**——special caution
**`dc024` at 2020 wave item belongs to depression scale"felt lonely"items, not serial subtraction 7**.

## 1. cognitive domains(outcome variables)

| domain | 2020 variables | coding | scoring |
|---|---|---|---|
| orientation-year | `dc001` | 1 correct / 2 incorrect / 997 do not know / 999 refused | correct=1; incorrect and do not know=0; refused=missing |
| orientation-season | `dc002` | same as above | same as above |
| orientation-day | `dc003` | same as above | same as above |
| orientation-day of week | `dc004` | same as above | same as above |
| orientation-month | `dc005` | same as above | same as above |
| self rated memory | `dc006` | 1 excellent…5 poor | excluded from objective score |
| serial subtraction of 7 | `dc007_1`–`dc007_5`(answers in `dc007_1_1`–`dc007_5_1`) | 1 recorded answer / 997 / 999 | 0–5, **evaluate relative to respondent previous answer**(official definition); 997 counts as attempted but incorrect(score 0), 999 code as missing |
| figure copying | `dc009` | 1 correct / 2 incorrect / 3 not assessed / 997 / 999 | correct=1, incorrect=0, not assessed/do not know/refused=missing |
| immediate recall(trial 1) | `dc012_s1`–`dc012_s10`(`_s11` no words recalled, `_s12` refused, `_s13` unable to test) | each word position 1=recalled, 0=not recalled | 0–10 |
| immediate recall(trial 2, 3) | `dc013_s1`–`s10`, `dc014_s1`–`s10` | same as above | 0–10 |
| delayed recall | `dc028_s1`–`dc028_s10`(at CES-D administered afterward) | same as above | 0–10 |

**2020 wave protocol changes(important, 2026-09-13 correction)**: word recall consists of **3 learning trials(dc012/013/014)+ 1 delayed recall trials(dc028)**.
confirmed after checking every wave: **2013, 2018, 2020 all are 3 learning trials, 2011, 2015 only 1 times**(earlier"2011–2018 only 1 times"statement was incorrect).
primary analysis consistently uses **trial 1(dc012)** for comparability with earlier waves, mean of three trials used in sensitivity analysis.
additional note: number 4 round(2018)word list differs from other waves, official recommendation is equipercentile equating before longitudinal memory analysis(Gross 2012), see details in `code/variable_map_multiwave.md` number 4.1 section.

**composite**: 
- episodic memory = immediate + delayed(0–20)
- mental status = orientation(0–5)+ serial subtraction of 7(0–5)+ figure copying(0–1)= 0–11
- global cognition = episodic memory + mental status(0–31), then z standardization

**missing data handling(v3, 2026-09-10 finalized)**: 
- **"do not know"(997) → score as incorrect 0, and regard as"attempted"**(orientation, serial subtraction of 7, also applies to figure copying); 
- **"refused"(999) / "not applicable"(-1) / not tested → missing**; 
- figure copying"not assessed"(3) → missing(unable due to physical condition, informative missingness); 
- word recall: **"no words recalled" → 0 points**; **"refused recall"/"unable to test" → missing**; 
- sum available components, missing only if all components missing; **primary analysis restricted to 9 participants with all components scorable**.
- ⚠️ note: v2 previously claimed serial subtraction 7 of 997 score 0, but empty numeric field was actually coded missing, v3 corrected using indicator flags.

## 2. depression(CES-D-10)

| items | variables | notes |
|---|---|---|
| 1 bothered by small things | `dc016` | |
| 2 difficulty concentrating | `dc017` | |
| 3 felt depressed | `dc018` | |
| 4 everything was an effort | `dc019` | |
| 5 hopeful about future | `dc020` | **reverse scoring** |
| 6 felt fearful | `dc021` | |
| 7 restless sleep | `dc022` | |
| 8 felt happy | `dc023` | **reverse scoring** |
| 9 felt lonely | `dc024` | |
| 10 could not get going | `dc025` | |

coding 1–4(rarely/some of the time/moderate/most of the time)→ 0–3; 997/999 → missing.
total score 0–30; **depression = total score ≥ 10**(script requires 10 complete questions).

## 3. chronic conditions(15 classes, `da002` + `da003` joint classification)

2020 wave chronic disease module uses skip logic: previously diagnosed participants asked `da002_i`(change since last visit, value 1/2/3; 99=absent at preceding wave), 
previously unaffected participants asked `da003_i`(ever diagnosed by a doctor, 1=yes, 2=no).

**classification rule**: chronic condition present = `da002_i ∈ {1,2,3}` or `da003_i == 1`; none = `da003_i == 2`; otherwise missing.

| # | disease | variables |
|---|---|---|
| 1 | hypertension | `da002_1_` / `da003_1_` |
| 2 | dyslipidemia | `da002_2_` / `da003_2_` |
| 3 | diabetes or high blood sugar | `da002_3_` / `da003_3_` |
| 4 | cancer or malignant tumor | `_4_` |
| 5 | chronic lung disease | `_5_` |
| 6 | liver disease | `_6_` |
| 7 | heart disease | `_7_` |
| 8 | stroke | `_8_` |
| 9 | kidney disease | `_9_` |
| 10 | stomach or digestive disease | `_10_` |
| 11 | emotional or psychiatric problems | `_11_` |
| 12 | memory related disease | `_12_` |
| 13 | Parkinson disease | `_13_` |
| 14 | arthritis or rheumatism | `_14_` |
| 15 | asthma | `_15_` |

## 4. ADL / IADL(module DB)

coding: 1 no difficulty / 2 difficulty but still able / 3 difficulty requiring help / 4 unable.
**difficulty = value ≥ 2**; any difficult item classified as"ADL/IADL limited".

| ADL(6 items) | variables | IADL(6 items) | variables |
|---|---|---|---|
| dressing | `db001` | housework | `db012` |
| bathing | `db003` | preparing meals | `db014` |
| eating | `db005` | shopping | `db016` |
| getting in or out of bed | `db007` | making telephone calls | `db018` |
| toileting | `db009` | taking medication | `db020` |
| continence | `db011` | managing finances | `db022` |

(intervening `db002/db004/...` is"whether assistance received", not used in this project.)

## 5. pain / sleep / lifestyle

| variables | meaning | coding |
|---|---|---|
| `da027` | bodily pain severity(5 ordered levels) | 1 none / 2 a little / 3 some / 4 quite a bit / 5 very much; three categories for analysis(1=no pain / 2–3=mild / 4–5=severe, variables `pain_3cat`), also retain binary sensitivity definition |
| `da028_s1`–`s16` | pain location | 0/location code |
| `da030` | mean nighttime sleep in past month(hours) | continuous, -1 → missing |
| `da031` | nap duration(minutes) | continuous |
| `dc022` | CES-D items"restless sleep"(proxy for sleep quality) | 1–4; ≥3(per week 3–4 days or more)classify as restless sleep; **also include CES-D total score** |
| `da046` | ever smoked(answered only by prior wave nonsmokers) | 1 yes / 2 no |
| `zsmoke` | smoking at previous interview | 1 yes(not asked this wave da046, direct response da047) |
| `da047` | currently smoking or quit | 1 current smoker / 2 former smoker / 3 never |
| `da051` | past year drinking frequency | 1 ≥once per month / 2 less than monthly / 3 no alcohol use(`da052` monthly frequency among drinkers, not used in this project) |
| `da032_1_/2_/3_` | high/in/low intensity physical activity(per occasion ≥10 minutes) | 1 yes / 2 no |

### 5.1 variable specification and prespecified rationale(2026-09-10 set, write into SAP)

**pain**: `da027` itself is **5 level ordinal variable**, this project uses **3 categories**modeling(no pain / mild / severe).
prespecified rationale: ① dichotomization loses intensity information; ② pain cognition relationship is not monotonic, dichotomization combines opposing mild and severe effects
cancellation yields a null result.
**also retain binary version(present/none)report as sensitivity analysis alongside**(see `TableS2_pain_sleep_varform.csv`): 
binary adjusted OR 0.96(0.84–1.10, p=0.60); three categories: mild 0.89(0.76–1.03, **p=0.114, nonsignificant**), 
severe 1.17(0.99–1.38, **p=0.063, borderline**).
(numerical value is 2026-09-12 recomputed after completing education; previous version: binary 0.95, mild 0.86(p=0.039), severe 1.16(p=0.083).)
(2026-09-10 from this point all are survey weighted svyglm: INDV_weight + community clustering, adjust for age/sex/four level education factor/urban or rural; 
previous unweighted glm results——binary 0.93 p=0.15, severe 1.21 p=0.004——historical reference only.)

**sleep**: by **<6h / 6–8h / >8h** categorical modeling(using 6–8h as reference), and retain"per increase of 1 hours"continuous coding for sensitivity.
prespecified rationale: sleep duration cognition relationship is known to be **U shaped**, linear test has low power.
results(survey weighted svyglm, adjust for age/sex/four education groups/urban or rural, 2026-09-12 regenerate): 
continuous adjustment OR 1.00(0.96–1.03, p=0.80); grouping <6h **1.19(1.05–1.36, p=0.007)**, >8h 1.60(1.24–2.06, p<0.001).
(2026-09-10 previous unweighted glm result is continuous 1.01, <6h 1.18, >8h 1.81, historical reference only.)

> ⚠️ **two interpretation boundaries**
> 1. impairment rates across three pain categories are no pain 15.8% / mild **13.5%** / severe 20.2%, **minimum at"mild"**, 
>    belongs to**nonmonotonic**relationship, **must not describe as"dose response"**; can only state"severe pain associated with higher cognitive impairment odds, no such association for mild pain".
> 2. "mild pain below no pain"likely reflects**reporting bias/reverse causation**(cognitively impaired participants underreport mild or moderate pain, or proxy responses): 
>    in mild pain group ADL limitation rate(19.1%)much higher than no pain group(7.5%), yet lower cognitive impairment rate, biologically implausible.
>    should in Discussion discuss in limitations, **do not interpret as protective**.

**BMI**: 2020 wave core and all downloaded modules lack height/weight measurement; Harmonized D measured in BMI only through 2015 year(`r3mbmi`).
as decided(2026-09-10): **BMI not included in this study**(across 5 combining years introduces temporal mismatch and approximately 26% nonrandom missingness), Discussion state in limitations.

## 6. demographic and socioeconomic

| variables | meaning | coding |
|---|---|---|
| `xrage` | age | continuous |
| `xrgender` | sex | 1 male / 2 female |
| `zredu` + `ba010` | education(**2026-09-12 revision: combine complementary variables**, see explanation below; already used Harmonized `raeducl` cross verification) | 1 illiterate…11 doctoral degree |
| `ba011` | marital status | 1 married and cohabiting…6 never married |
| `ba008` | residence type | 1 urban/town center / 2 urban rural fringe / 3 rural / 4 special |
| `ba009` | household registration type | 1 agricultural / 2 nonagricultural / 3 unified resident / 4 no household registration |
| `proxy_2` | proxy interview indicator | 1 yes(exclude) |

education groups: <primary school(1–3), primary school(4), middle school(5), high school or above(6–11).

**combine education variables(2026-09-12)**: `zredu`(carried forward from previous wave, nonmissing **17,687**)and `ba010`
(newly asked this wave, nonmissing **1,708**)at 2020 wave**fully complementary**——no overlap, no gaps, total 19,395 participants.
both use identical 11 level coding(verified against birth cohort distributions), therefore combine into `edu = coalesce(zredu, ba010)`
after**no missing values**, no cross wave backfilling needed.

## 7. weights

| variables | meaning |
|---|---|
| `INDV_weight` | individual weights(**primary weights for this project**) |
| `INDV_weight_ad2` | adjusted individual weights 2 |
| `HH_weight` / `HH_weight_ad1` | household weights and adjusted version |

weighted design: `weights = INDV_weight`, `ids = communityID`(as PSU approximate).

## 8. known limitations and unresolved items(v3)

1. whether orientation included"season"(included in this project, gives 0–5; Harmonized 1–4 wave uses only 4 items).
2. figure copying"not assessed (3)"treat as missing.
3. CES-D-10 require 10 complete questions(in primary analysis sample 13,844/14,655 has scores).
4. joint chronic disease classification(da002 + da003)project specific rule, must in SAP state explicitly.
5. **primary analysis restricted to cognition 9 participants with complete components**(n=**14,655**, no missing education); participants with incomplete assessment 3,075 participants
   sensitivity analysis only.incompletely assessed participants systematically older, less healthy, less educated(see Table S1), selection bias present.
6. **cutoff stratified by education**(<within stratum mean−1SD: **6.96 / 12.38 / 14.61 / 16.33**, 2026-09-12 reestimate after education backfilling).
   ⚠️ definition by construction**substantially compresses**education and sex gradients(2026-09-12 after harmonizing education strata to weighted estimates, four groups are 16.3 / 16.5 / 13.3 / 13.3%, between groups CI overlap), therefore"stratified prevalence pattern"requires alternative definition A(overall mean−1SD)display.
7. ~~missing education excluded~~ —— **resolved(2026-09-12)**: `zredu` + `ba010` no missing education after combining, 
   primary sample from 13,592 increased to **14,655**; regression standardized residual method abandoned(method itself requires education).
8. BMI unavailable(2020 no physical measurements this wave), decision to exclude already made.
9. smoking status based on da047 as reference; few contradictory raw records(such as da046=1 and da047=3), by da047 classification.
10. **sleep quality based on CES-D items dc022("restless sleep")proxy——item itself is CES-D-10 one component of**, 
    and"depression(CES-D ≥ 10)"constitutes**shared source variables**, cannot present both as independent"multidimensional vulnerability domains"in parallel, 
    otherwise construct overlap is a concern.SAP must address explicitly in(sleep domain uses duration only/grouping, or acknowledge structural overlap and in
    sensitivity analyses mutually adjust), Discussion discuss in limitations.
