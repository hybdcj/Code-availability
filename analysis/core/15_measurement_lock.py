from pathlib import Path
import csv
from openpyxl import Workbook, load_workbook
from openpyxl.styles import Font, PatternFill, Alignment
from openpyxl.utils import get_column_letter

PROJ = Path(r'D:\BaiduSyncdisk\research\Xu_Zhonghao')
IN = PROJ / 'output' / '05_measurement_equivalence_ratio'
OUT = PROJ / 'output' / '08_plan_and_deliverables'

def read_csv(name):
    with open(IN / name, encoding='utf-8-sig', newline='') as f:
        return list(csv.DictReader(f))

common_sens = read_csv('measurement_common_score_sensitivity_2026-09-13.csv')
item_change = read_csv('measurement_item_change_2015_2018_2026-09-13.csv')
cfa_all = read_csv('measurement_cfa_allwaves_2026-09-13.csv')
cfa_pairs = read_csv('measurement_cfa_pairs_2026-09-13.csv')
memory_schemes = read_csv('measurement_memory_schemes_2026-09-13.csv')
memory_desc = read_csv('measurement_memory_wave_descriptives_2026-09-13.csv')

main = next(r for r in common_sens if r['scenario'] == 'all')
age60 = next(r for r in common_sens if r['scenario'] == 'age60_2018')
age60mem = next(r for r in common_sens if r['scenario'] == 'age60_2018_memory')
pair1518 = next(r for r in cfa_pairs if r['pair'] == '2015-2018')
delta_rmsea_1518 = f"{float(pair1518['delta_rmsea']):.4f}"
main_n = f"{int(main['n']):,}"
age60_n = f"{int(age60['n']):,}"
age60mem_n = f"{int(age60mem['n']):,}"

key_rows = [
    ('common_items', '5 orientation + 5 serial-7 items; total 0-10'),
    ('fixed_cohort_n', main['n']),
    ('item_complete_CFA_n', '5359'),
    ('mean_2011', main['m2011']), ('mean_2013', main['m2013']), ('mean_2015', main['m2015']),
    ('mean_2018', main['m2018']), ('mean_2020', main['m2020']),
    ('change_2015_2018', main['change_2015_2018']), ('ci_2015_2018', main['ci_2015_2018']),
    ('change_2011_2020', main['change_2011_2020']), ('ci_2011_2020', main['ci_2011_2020']),
    ('age60_n', age60['n']), ('age60_change_2015_2018', age60['change_2015_2018']), ('age60_ci_2015_2018', age60['ci_2015_2018']),
    ('age60_memory_n', age60mem['n']), ('age60_memory_change_2015_2018', age60mem['change_2015_2018']), ('age60_memory_ci_2015_2018', age60mem['ci_2015_2018']),
    ('delta_cfi_2015_2018', pair1518['delta_cfi']), ('delta_rmsea_2015_2018', delta_rmsea_1518),
    ('memory_n', '6338'),
    ('memory_scheme_A', 'proportional memory score, +0.0375/decade'),
    ('memory_scheme_B', 'within-wave z, -0.1728/decade'),
    ('memory_scheme_C', 'proportional excluding 2018, +0.0538/decade'),
]
with open(OUT / 'T09_measurement_key_numbers_2026-09-20.csv', 'w', encoding='utf-8', newline='') as f:
    wr = csv.writer(f)
    wr.writerow(['item', 'value'])
    wr.writerows(key_rows)

# Chinese Methods
methods_zh = f'''# T09 Methods paragraph(Chinese draft)

cognition assessed using comparable five wave 10 item common set, including orientation 5 items(year, month, date, day of week, season)and serial subtraction of 7 five steps.each correct item scores 1 points, total score 0–10 points.primary longitudinal outcome uses common total score; binary cognitive impairment based on 2011 year education specific fixed cutoff, to avoid threshold drift from wave specific recalculation.

because 2018 year cognitive modules and word lists differ, study does not combine memory words, figure copying and related items into 0–31 point total for longitudinal comparisons, instead uses consistent items across five waves 0–10 point common item set.to assess cross wave comparability of common items, report item completeness and correct response rates by wave, compare using confirmatory factor analysis 2011–2013, 2013–2015, 2015–2018 and 2018–2020 configural and scalar invariance models.measurement invariance assessed using ΔCFI and ΔRMSEA as supporting evidence, not solely chi squared tests.all waves CFA among participants complete on all five wave items n=5,359 performed in.

at 2011 year≥45 years old, five waves ms10 in fully scorable fixed cohort, total{main_n}participants.2015→2018 mean change in common item score is{main['change_2015_2018']}points(95%CI {main['ci_2015_2018'].replace(' to ', ' to ')}); retain only 2018 aged at least 60 participants above age number{age60['change_2015_2018']}points(95%CI {age60['ci_2015_2018'].replace(' to ', ' to ')}); further require 2018 year complete memory participants are{age60mem['change_2015_2018']}points(95%CI {age60mem['ci_2015_2018'].replace(' to ', ' to ')}).three samples show consistent direction, decline unlikely explained solely by 2018 year 60 module differences among younger participants.

2015→2018 changes in adjacent wave scalar invariance models are ΔCFI={pair1518['delta_cfi']}, ΔRMSEA={delta_rmsea_1518}; below conventional substantive noninvariance thresholds.but some CFA models show negative variance warnings, therefore describe results as“no noninvariance evidence sufficient to invalidate common total score”, must not describe as“complete equivalence proven”.

memory scores affected by round 4 word lists and learning trial counts differ from earlier rounds, may influence delayed recall, therefore report descriptively only; for longitudinal comparisons, first apply equipercentile equating or restrict to identical word lists.
'''
(OUT / 'T09_methods_zh_2026-09-20.md').write_text(methods_zh, encoding='utf-8')

methods_en = f'''# T09 Methods Paragraph (English Draft)

Cognitive function was assessed using a 10-item common score that was available across all five waves, comprising five orientation items (year, month, date, day of week, and season) and five serial-7 subtraction steps. Each correctly answered item contributed one point, yielding a total score of 0–10. This common score was used for longitudinal analyses. Cognitive impairment was defined using the 2011 education-specific fixed threshold to avoid threshold drift across waves.

Because the broader cognitive module and memory word lists differed in 2018, memory and drawing items were not combined with the common items into a 0–31 global score for longitudinal comparison. Measurement comparability of the common score was evaluated using item-level completeness, item-level changes, and multiple-group confirmatory factor analysis across adjacent waves and from 2011 to 2020. Configural and scalar invariance models were compared using ΔCFI and ΔRMSEA rather than relying on χ² tests alone. The all-wave CFA was conducted among 5,359 participants with complete item-level data.

The fixed cohort comprised {main_n} participants aged ≥45 years in 2011 with the common score available at all five waves. The mean 2015-to-2018 change in the common score was {main['change_2015_2018']} points (95% CI {main['ci_2015_2018'].replace(' to ', ' to ')}). The corresponding change was {age60['change_2015_2018']} (95% CI {age60['ci_2015_2018'].replace(' to ', ' to ')}) among participants aged ≥60 years in 2018 and {age60mem['change_2015_2018']} (95% CI {age60mem['ci_2015_2018'].replace(' to ', ' to ')}) among those who also had complete memory items. The consistency of these estimates suggested that the observed decline was unlikely to be solely attributable to the age-restricted 2018 module.

For the 2015-to-2018 comparison, scalar invariance was associated with ΔCFI={pair1518['delta_cfi']} and ΔRMSEA={delta_rmsea_1518}, which did not exceed conventional thresholds for substantial non-invariance. However, some CFA models produced negative-variance warnings. Therefore, the results should be interpreted as showing no evidence sufficient to invalidate the common score, rather than as definitive proof of full measurement equivalence.

Memory scores were reported descriptively only because word-list content and learning-trial procedures differed in 2018, which may affect delayed-recall scores. Percentile equating or comparisons restricted to waves using the same word list would be required before longitudinal memory conclusions are drawn.
'''
(OUT / 'T09_methods_en_2026-09-20.md').write_text(methods_en, encoding='utf-8')

lock = f'''# T09: 2018 measurement differences and common item set lock(2026-09-20)

## one, final locked definition

- primary cognitive score: orientation 5 items + serial subtraction of 7 five steps, total score 0–10; 
- primary longitudinal outcome: common item set `ms10`; 
- binary impairment: 2011 year education specific fixed cutoff; 
- 2018 wave does not use 0–31 point total for main five wave longitudinal conclusions; 
- memory score descriptive only/supplementary analysis.

## two, key results

- fixed cohort: {main_n}participants; all waves CFA item complete sample: 5,359 participants; 
- 2015→2018: {main['change_2015_2018']}points(95%CI {main['ci_2015_2018'].replace(' to ', ' to ')}); 
- 2018 aged at least 60 years old: {age60['change_2015_2018']}points(95%CI {age60['ci_2015_2018'].replace(' to ', ' to ')}); 
- 2018 aged at least 60 age eligible and complete memory: {age60mem['change_2015_2018']}points(95%CI {age60mem['ci_2015_2018'].replace(' to ', ' to ')}); 
- 2015→2018 CFA: ΔCFI={pair1518['delta_cfi']}, ΔRMSEA={delta_rmsea_1518}.

## three, permitted wording

- “no measurement noninvariance evidence sufficient to invalidate common total score”; 
- “adjacent wave scalar invariance model ΔCFI and ΔRMSEA below conventional substantive thresholds”; 
- “2018 year module differences alone do not explain primary decline”.

## four, prohibited wording

- must not state“complete measurement equivalence demonstrated”; 
- do not treat CFA hide negative variance warnings; 
- memory trajectories must not be main longitudinal conclusions; 
- do not treat 0–31 point total and 0–10 point common set are not interchangeable.

## five, corresponding files

- `T09_methods_zh_2026-09-20.md`
- `T09_methods_en_2026-09-20.md`
- `T09_measurement_equivalence_supplement_2026-09-20.xlsx`
- `T09_measurement_key_numbers_2026-09-20.csv`
'''
(OUT / 'T09_measurement_lock_2026-09-20.md').write_text(lock, encoding='utf-8')

# Clean valid Excel supplement.
wb = Workbook()
ws = wb.active; ws.title = 'notes'
ws.append(['project', 'content'])
notes = [
    ('analysis participants', '2011 year≥45 years old, five waves ms10 fully scorable fixed cohort'),
    ('primary outcome', 'orientation 5 questions+serial subtraction of 7, total score 0–10; binary impairment uses 2011 year education specific fixed cutoff'),
    ('common item set', '10 items included in all five waves; item completeness, item change and CFA used to assess comparability'),
    ('2018 wave handling', '2018 year broader cognitive module and word list changed; do not use 0–31 total score as primary longitudinal outcome'),
    ('measurement invariance', 'CFA no substantive noninvariance detected, but some models have negative variance warnings, cannot claim complete equivalence proven'),
    ('memory score handling', 'memory scores reported descriptively only; raw longitudinal memory affected by word lists and learning trials'),
    ('scope of conclusions', 'CHARLS not a surgical cohort; supports preoperative related cognitive measurement only, no inference about postoperative outcomes'),
]
for row in notes: ws.append(row)

def convert_cell(v):
    if v is None or v == '':
        return v
    try:
        return float(v)
    except (TypeError, ValueError):
        return v

def display_cell(header, v):
    if header in {'p', 'scalar_p'} and v is not None:
        try:
            if float(v) == 0:
                return '<0.001'
        except (TypeError, ValueError):
            pass
    return convert_cell(v)

def add_csv_sheet(name, rows):
    s = wb.create_sheet(name)
    if not rows:
        return
    headers = list(rows[0].keys())
    s.append(headers)
    for r in rows:
        s.append([display_cell(h, r.get(h, '')) for h in headers])

add_csv_sheet('sample sensitivity', common_sens)
add_csv_sheet('item specific change', item_change)
add_csv_sheet('CFA_all waves', cfa_all)
add_csv_sheet('CFA_adjacent waves', cfa_pairs)
add_csv_sheet('memory scoring', memory_schemes)
add_csv_sheet('memory description', memory_desc)

ws = wb.create_sheet('common item set')
ws.append(['domain', 'project', 'scoring', 'total score'])
items = [
    ('orientation', 'year', 'correct=1, remaining=0', 'orientation 0–5'),
    ('orientation', 'month', 'correct=1, remaining=0', 'orientation 0–5'),
    ('orientation', 'date', 'correct=1, remaining=0', 'orientation 0–5'),
    ('orientation', 'day of week', 'correct=1, remaining=0', 'orientation 0–5'),
    ('orientation', 'season', 'correct=1, remaining=0', 'orientation 0–5'),
    ('serial subtraction of 7', 'step one=93', 'answer 93=1, otherwise 0', 'serial subtraction of 7: 0–5'),
    ('serial subtraction of 7', 'step two', 'subtract from previous response 7 correct=1, otherwise 0', 'serial subtraction of 7: 0–5'),
    ('serial subtraction of 7', 'step three', 'subtract from previous response 7 correct=1, otherwise 0', 'serial subtraction of 7: 0–5'),
    ('serial subtraction of 7', 'step four', 'subtract from previous response 7 correct=1, otherwise 0', 'serial subtraction of 7: 0–5'),
    ('serial subtraction of 7', 'step five', 'subtract from previous response 7 correct=1, otherwise 0', 'serial subtraction of 7: 0–5'),
]
for row in items: ws.append(row)

header_fill = PatternFill('solid', fgColor='1F4E78')
for s in wb.worksheets:
    for cell in s[1]:
        cell.fill = header_fill
        cell.font = Font(color='FFFFFF', bold=True)
        cell.alignment = Alignment(horizontal='center', vertical='center')
    s.freeze_panes = 'A2'
    for col in range(1, s.max_column + 1):
        max_len = max(len(str(s.cell(r, col).value or '')) for r in range(1, min(s.max_row, 100) + 1))
        s.column_dimensions[get_column_letter(col)].width = min(max(max_len + 2, 12), 55)

xlsx = OUT / 'T09_measurement_equivalence_supplement_2026-09-20.xlsx'
wb.save(xlsx)
# Reload must succeed.
check = load_workbook(xlsx, read_only=True, data_only=True)
required = {'notes','common item set','sample sensitivity','item specific change','CFA_all waves','CFA_adjacent waves','memory scoring','memory description'}
assert required.issubset(set(check.sheetnames))
assert check['sample sensitivity'].max_row > 1
print('T09 Excel supplement regenerated and validated')
