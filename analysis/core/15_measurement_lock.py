from pathlib import Path
import csv
from openpyxl import Workbook, load_workbook
from openpyxl.styles import Font, PatternFill, Alignment
from openpyxl.utils import get_column_letter

PROJ = Path(r'D:\BaiduSyncdisk\科研\徐中豪')
IN = PROJ / 'output' / '05_测量等价比'
OUT = PROJ / 'output' / '08_方案与交付'

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
methods_zh = f'''# T09 Methods 段落（中文草稿）

认知功能使用五波均可比的10项共同题集评估，包括定向5项（年、月、日期、星期、季节）和连续减7五步。各项目正确计1分，总分0–10分。主要纵向结局使用该共同总分；二分类认知受损按2011年教育分层固定切点定义，以避免各波重新计算切点造成阈值漂移。

由于2018年认知模块和记忆词表存在差异，本研究未将记忆词、画图等项目合并为0–31分总分进行纵向比较，而是使用五波项目一致的0–10分共同题集。为评估该共同题集的跨波可比性，研究逐波报告逐题完整性和正确率，并以验证性因子分析比较2011–2013、2013–2015、2015–2018和2018–2020的形态等值与标量等值模型。测量等值性以ΔCFI和ΔRMSEA辅助判断，而不是仅依据卡方检验。全波CFA在五波逐题完整者n=5,359中进行。

在2011年≥45岁、五波ms10均可计分的固定队列中，共{main_n}人。2015→2018共同题集平均变化为{main['change_2015_2018']}分（95%CI {main['ci_2015_2018'].replace(' to ', ' 至 ')}）；仅保留2018年满60岁者为{age60['change_2015_2018']}分（95%CI {age60['ci_2015_2018'].replace(' to ', ' 至 ')}）；进一步要求2018年记忆题完整者为{age60mem['change_2015_2018']}分（95%CI {age60mem['ci_2015_2018'].replace(' to ', ' 至 ')}）。三种样本的结果方向一致，提示下降不太可能仅由2018年60岁以下受试者的模块差异造成。

2015→2018相邻波的标量等值模型变化为ΔCFI={pair1518['delta_cfi']}、ΔRMSEA={delta_rmsea_1518}；未超过常用实质性非等值阈值。但部分CFA模型出现负方差警告，因此结果应表述为“未发现足以推翻共同题总分的非等值证据”，不能写成“已证明完全等值”。

记忆分因第4轮词表内容和学习试次数与前几轮不同，可能影响延迟回忆成绩，因此仅作描述性报告；如需纵向比较，应先按等百分位等分或限定同一词表波次。
'''
(OUT / 'T09_methods_zh_2026-09-20.md').write_text(methods_zh, encoding='utf-8')

methods_en = f'''# T09 Methods Paragraph (English Draft)

Cognitive function was assessed using a 10-item common score that was available across all five waves, comprising five orientation items (year, month, date, day of week, and season) and five serial-7 subtraction steps. Each correctly answered item contributed one point, yielding a total score of 0–10. This common score was used for longitudinal analyses. Cognitive impairment was defined using the 2011 education-specific fixed threshold to avoid threshold drift across waves.

Because the broader cognitive module and memory word lists differed in 2018, memory and drawing items were not combined with the common items into a 0–31 global score for longitudinal comparison. Measurement comparability of the common score was evaluated using item-level completeness, item-level changes, and multiple-group confirmatory factor analysis across adjacent waves and from 2011 to 2020. Configural and scalar invariance models were compared using ΔCFI and ΔRMSEA rather than relying on χ² tests alone. The all-wave CFA was conducted among 5,359 participants with complete item-level data.

The fixed cohort comprised {main_n} participants aged ≥45 years in 2011 with the common score available at all five waves. The mean 2015-to-2018 change in the common score was {main['change_2015_2018']} points (95% CI {main['ci_2015_2018'].replace(' 至 ', ' to ')}). The corresponding change was {age60['change_2015_2018']} (95% CI {age60['ci_2015_2018'].replace(' 至 ', ' to ')}) among participants aged ≥60 years in 2018 and {age60mem['change_2015_2018']} (95% CI {age60mem['ci_2015_2018'].replace(' 至 ', ' to ')}) among those who also had complete memory items. The consistency of these estimates suggested that the observed decline was unlikely to be solely attributable to the age-restricted 2018 module.

For the 2015-to-2018 comparison, scalar invariance was associated with ΔCFI={pair1518['delta_cfi']} and ΔRMSEA={delta_rmsea_1518}, which did not exceed conventional thresholds for substantial non-invariance. However, some CFA models produced negative-variance warnings. Therefore, the results should be interpreted as showing no evidence sufficient to invalidate the common score, rather than as definitive proof of full measurement equivalence.

Memory scores were reported descriptively only because word-list content and learning-trial procedures differed in 2018, which may affect delayed-recall scores. Percentile equating or comparisons restricted to waves using the same word list would be required before longitudinal memory conclusions are drawn.
'''
(OUT / 'T09_methods_en_2026-09-20.md').write_text(methods_en, encoding='utf-8')

lock = f'''# T09：2018 测量差异说明与共同题集锁定（2026-09-20）

## 一、最终锁定口径

- 主认知分数：定向5项 + 连续减7五步，总分0–10；
- 主纵向结局：共同题集 `ms10`；
- 二分类受损：2011年教育分层固定切点；
- 2018期不将0–31分总分用于五波纵向主结论；
- 记忆分只作描述性/补充分析。

## 二、关键结果

- 固定队列：{main_n}人；全波CFA逐题完整样本：5,359人；
- 2015→2018：{main['change_2015_2018']}分（95%CI {main['ci_2015_2018'].replace(' to ', ' 至 ')}）；
- 2018年满60岁：{age60['change_2015_2018']}分（95%CI {age60['ci_2015_2018'].replace(' to ', ' 至 ')}）；
- 2018年满60岁且记忆完整：{age60mem['change_2015_2018']}分（95%CI {age60mem['ci_2015_2018'].replace(' to ', ' 至 ')}）；
- 2015→2018 CFA：ΔCFI={pair1518['delta_cfi']}，ΔRMSEA={delta_rmsea_1518}。

## 三、允许的写法

- “未发现足以推翻共同题总分的测量非等值证据”；
- “相邻波标量等值模型的ΔCFI和ΔRMSEA未超过常用实质性阈值”；
- “2018年模块差异未单独解释主要下降”。

## 四、禁止的写法

- 不得写“已证明完全测量等值”；
- 不得把CFA的负方差警告隐藏；
- 不得把记忆轨迹作为主要纵向结论；
- 不得把0–31分总分与0–10分共同题集混用。

## 五、对应文件

- `T09_methods_zh_2026-09-20.md`
- `T09_methods_en_2026-09-20.md`
- `T09_measurement_equivalence_supplement_2026-09-20.xlsx`
- `T09_measurement_key_numbers_2026-09-20.csv`
'''
(OUT / 'T09_measurement_lock_2026-09-20.md').write_text(lock, encoding='utf-8')

# Clean valid Excel supplement.
wb = Workbook()
ws = wb.active; ws.title = '说明'
ws.append(['项目', '内容'])
notes = [
    ('分析对象', '2011年≥45岁、五波ms10均可计分的固定队列'),
    ('主结局', '定向5题+连续减7，总分0–10；二分类受损使用2011年教育分层固定切点'),
    ('共同题集', '10个项目在五波均纳入；逐题完整性、逐题变化和CFA用于评估可比性'),
    ('2018期处理', '2018年较广认知模块和记忆词表有变化；不使用0–31总分作为纵向主结局'),
    ('测量等值性', 'CFA未发现实质性非等值证据，但部分模型有负方差警告，不能声称已证明完全等值'),
    ('记忆分处理', '记忆分仅描述性报告；原始纵向记忆分受词表内容和学习试次变化影响'),
    ('结论边界', 'CHARLS不是手术队列；本材料只支持术前相关认知测量，不推断术后结局'),
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

add_csv_sheet('样本敏感性', common_sens)
add_csv_sheet('逐题变化', item_change)
add_csv_sheet('CFA_全波', cfa_all)
add_csv_sheet('CFA_相邻波', cfa_pairs)
add_csv_sheet('记忆计分', memory_schemes)
add_csv_sheet('记忆描述', memory_desc)

ws = wb.create_sheet('共同题集')
ws.append(['维度', '项目', '计分', '总分'])
items = [
    ('定向', '年', '正确=1，其余=0', '定向0–5'),
    ('定向', '月', '正确=1，其余=0', '定向0–5'),
    ('定向', '日期', '正确=1，其余=0', '定向0–5'),
    ('定向', '星期', '正确=1，其余=0', '定向0–5'),
    ('定向', '季节', '正确=1，其余=0', '定向0–5'),
    ('连续减7', '第一步=93', '答案93=1，否则0', '连续减7：0–5'),
    ('连续减7', '第二步', '按前一步减7正确=1，否则0', '连续减7：0–5'),
    ('连续减7', '第三步', '按前一步减7正确=1，否则0', '连续减7：0–5'),
    ('连续减7', '第四步', '按前一步减7正确=1，否则0', '连续减7：0–5'),
    ('连续减7', '第五步', '按前一步减7正确=1，否则0', '连续减7：0–5'),
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
required = {'说明','共同题集','样本敏感性','逐题变化','CFA_全波','CFA_相邻波','记忆计分','记忆描述'}
assert required.issubset(set(check.sheetnames))
assert check['样本敏感性'].max_row > 1
print('T09 Excel supplement regenerated and validated')
