"""Wave6 QC v1. No outcome threshold selection, AUC calculation or model fitting.
Run export_wave6_qc_inputs_20261113.R first from the project root.
"""
from pathlib import Path
import json, hashlib
import numpy as np
import pandas as pd
import pyreadstat

P=Path(__file__).resolve().parents[1]
O=P/'output/wave6_qc_20261113'; O.mkdir(exist_ok=True)
inputs=[]; dictionary=[]
def read(rel):
    f=P/rel;inputs.append(f);d,m=pyreadstat.read_dta(str(f))
    assert not d.ID.duplicated().any(), rel
    assert d.ID.str.fullmatch(r'\d{12}').all(), rel
    return d,m
def save(df,name): df.to_csv(O/name,index=False,encoding='utf-8-sig')
def recall(d,prefix):
    cols=[prefix+'_s'+str(i) for i in range(1,11)]
    x=d[cols].apply(pd.to_numeric); cnt=x.gt(0).sum(axis=1)
    none=d[prefix+'_s11'].eq(11);ref=d[prefix+'_s12'].eq(12)
    unable=d.get(prefix+'_s13',pd.Series(np.nan,index=d.index)).eq(13)
    badcode=pd.concat([~x[c].isin([0,i+1]) & x[c].notna() for i,c in enumerate(cols)],axis=1).any(axis=1)
    conflict=(none & cnt.gt(0)) | (ref & (none|cnt.gt(0))) | (unable & (none|cnt.gt(0)))
    observed=x.notna().all(axis=1)&(none|cnt.gt(0))
    valid=observed & ~ref & ~unable & ~conflict & ~badcode
    reason=np.select([ref,unable,conflict|badcode,~observed],['refused','unable','conflict','not_observed'],default='valid')
    return cnt.where(valid).astype(float),reason,conflict|badcode

parts=[];annual=[]
for year in [2021,2022,2023]:
    rel=f'data/raw/wave6_2021_2023/{year}/Cognition_Depression.dta'
    d,m=read(rel); im,dl=('dc039','dc061') if year==2021 else ('dc151','dc153')
    ir,ireason,ic=recall(d,im);dr,dreason,dc=recall(d,dl)
    proxy=pd.to_numeric(d['dcproxy' if year==2021 else 'proxy'])
    select=pd.to_numeric(d.get('dcselect',pd.Series(np.nan,index=d.index)))
    # 2022/23 proxy indicates eligibility; dcselect=7 confirms direct normal interview.
    allowed=~proxy.eq(1) | select.eq(7)
    adapted=proxy.eq(1)&select.isin([1,2])
    proxy_unresolved=proxy.eq(1)&~select.isin([1,2,7])
    mem=ir+dr
    r=pd.DataFrame({'ID':d.ID,'packyear':year,'community_w6':d.communityID,
        'immediate':ir,'delayed':dr,'memory_items':mem,'memory':mem.where(allowed),
        'memory_adapted':mem.where(allowed|adapted),'immediate_status':ireason,'delayed_status':dreason,
        'proxy':proxy,'dcselect':select,'adapted_test':adapted,'proxy_unresolved':proxy_unresolved,
        'item_conflict':ic|dc,'wordlist':d.xwordlist_1_.map({'苹果':'A','葡萄':'B'}),
        'official_sum':pd.to_numeric(d.xscorewr)+pd.to_numeric(d.xscoredr)})
    for raw,new in [('xscorewr',ir),('xscoredr',dr)]:
        ok=new.notna()&d[raw].notna();assert (new[ok]==pd.to_numeric(d.loc[ok,raw])).all()
    r['reason']=np.select([proxy_unresolved,adapted,ic|dc,pd.Series(ireason).eq('refused')|pd.Series(dreason).eq('refused'),pd.Series(ireason).eq('unable'),mem.isna()],['proxy_unresolved','adapted_test','item_conflict','refused','unable','not_observed'],default='valid')
    annual.append({'year':year,'rows':len(r),'official_complete':int(r.official_sum.notna().sum()),'item_valid':int(mem.notna().sum()),'primary_valid':int(r.memory.notna().sum()),'adapted_valid':int((adapted&mem.notna()).sum()),'unresolved_proxy_with_score':int((proxy_unresolved&mem.notna()).sum()),'item_conflict':int((ic|dc).sum())})
    demo,dm=read(f'data/raw/wave6_2021_2023/{year}/Demographics.dta')
    for c in ['xiwyear','xiwmonth']:
        r[c]=pd.to_numeric(r.ID.map(demo.set_index('ID')[c]),errors='coerce') if c in demo else np.nan
    for c in ['ID','dcproxy','proxy','dcselect',im+'_s1',im+'_s11',im+'_s12',im+'_s13',dl+'_s1',dl+'_s11',dl+'_s12','xscorewr','xscoredr','xscorewrmean','xwordlist_1_']:
        if c in d:dictionary.append({'year':year,'file':rel,'variable':c,'label':m.column_names_to_labels[c],'values':json.dumps(m.variable_value_labels.get(c,{}),ensure_ascii=False)})
    parts.append(r)
a=pd.concat(parts,ignore_index=True)
s,sm=read('data/raw/wave6_2021_2023/2023/Sample_Infor.dta')
for c in ['res_wave','iyear','imonth','died']:s[c]=pd.to_numeric(s[c],errors='coerce')
a=a.merge(s[['ID','res_wave','iyear','imonth','died']],on='ID',how='left',validate='many_to_one')
a['year_match']=a.packyear.eq(a.res_wave)
a['sample_month_valid']=a.imonth.between(1,12)&a.iyear.eq(a.packyear)
a['demo_date_conflict']=a.xiwyear.notna() & (a.xiwyear.ne(a.iyear)|a.xiwmonth.ne(a.imonth))
a['records_per_id']=a.groupby('ID').ID.transform('size')
a['record_decision']=np.select([~a.year_match,a.died.eq(1),a.memory.isna()],['exclude_off_response_year','exclude_dead','selected_but_cognition_unusable'],default='selected_valid')
sel=a[a.year_match].copy();assert not sel.ID.duplicated().any()
save(a,'cognition_records_internal.csv')
save(a[a.records_per_id.gt(1)|~a.year_match],'record_conflicts_internal.csv')
save(sel,'selected_wave6_internal.csv')
save(pd.DataFrame(annual),'annual_qc.csv');save(pd.DataFrame(dictionary),'variable_dictionary.csv')
dups=a[a.records_per_id.gt(1)].groupby('ID').agg(records=('ID','size'),valid_records=('memory','count'),different_valid_scores=('memory','nunique'),matched_year_records=('year_match','sum'))
save(dups.reset_index(),'duplicate_id_audit_internal.csv')

# Audit the existing 2018 memory construction against original item responses.
bpath=O/'baseline2018_internal.csv';inputs.append(bpath)
b=pd.read_csv(bpath,dtype={'IDn':str,'communityID':str}).rename(columns={'IDn':'ID'})
assert not b.ID.duplicated().any()
c18,cm=read('data/raw/wave4_2018/Cognition.dta');i18,_,_=recall(c18,'dc028_w4');d18,_,_=recall(c18,'dc047_w4')
raw18=pd.DataFrame({'ID':c18.ID,'im18_checked':i18,'dl18_checked':d18})
b=b.merge(raw18,on='ID',how='left',validate='one_to_one')
disagree=~((b.imrc.eq(b.im18_checked)|(b.imrc.isna()&b.im18_checked.isna()))&(b.dlrc.eq(b.dl18_checked)|(b.dlrc.isna()&b.dl18_checked.isna())))
save(b[disagree],'baseline_memory_disagreements_internal.csv')
b['baseline_memory']=b.im18_checked+b.dl18_checked
baseinfo,_=read('data/raw/wave4_2018/Sample_Infor.dta')
for c in ['iyear','imonth','died']:baseinfo[c]=pd.to_numeric(baseinfo[c],errors='coerce')
root=baseinfo[['ID','iyear','imonth','died']].rename(columns={'iyear':'year2018','imonth':'month2018','died':'died2018'}).merge(b,on='ID',how='left',validate='one_to_one')
flow=[]
def step(label,mask):flow.append({'step':label,'n':int(mask.sum())})
mask=pd.Series(True,index=root.index);step('2018 Sample_Infor total',mask)
mask&=root.died2018.eq(0);step('2018 alive',mask)
mask&=root.age.ge(45);step('age >=45',mask)
mask&=root.baseline_memory.notna();step('baseline memory observed',mask)
mask&=root.baseline_memory.gt(4);step('baseline memory >4',mask)
root['score_predictors_complete']=root[['age','edu_group','cesd10']].notna().all(axis=1)
root['four_predictors_complete']=root[['age','edu_group','cesd10','adl']].notna().all(axis=1)
step('three score predictors complete',mask&root.score_predictors_complete)
mask&=root.four_predictors_complete;step('four predictors complete: inception risk set',mask)
risk=root[mask].copy()
save(root,'baseline_all_internal.csv')
risk=risk.merge(s[['ID','res_wave','iyear','imonth','died']],on='ID',how='left',validate='one_to_one')
risk=risk.merge(sel[['ID','packyear','memory','memory_adapted','reason','wordlist','sample_month_valid','demo_date_conflict']],on='ID',how='left',validate='one_to_one')
for name in ['development_ids_internal.csv','source_ids_internal.csv']:inputs.append(O/name)
dev=set(pd.read_csv(O/'development_ids_internal.csv',dtype=str).KEY)
src=set(pd.read_csv(O/'source_ids_internal.csv',dtype=str).KEY)
risk['development_overlap']=risk.ID.isin(dev);risk['source_overlap']=risk.ID.isin(src)
risk['status']=np.select([risk.died.eq(1),risk.died.eq(0)&risk.memory.notna(),risk.died.eq(0)&risk.packyear.isna(),risk.died.eq(0)],['dead_at_recorded_wave','alive_memory_observed','alive_no_matched_cognition','alive_cognition_unusable'],default='unmatched_status_unknown')
risk['months_interval']=(risk.iyear-risk.year2018)*12+risk.imonth-risk.month2018
dateok=risk.month2018.between(1,12)&risk.imonth.between(1,12)&risk.iyear.eq(risk.packyear)&risk.months_interval.gt(0)
risk['date_usable']=dateok&~risk.demo_date_conflict.fillna(False).astype(bool)
risk['approx4y_42_54m']=risk.months_interval.between(42,54)&risk.date_usable
risk['tight4y_45_51m']=risk.months_interval.between(45,51)&risk.date_usable
risk['primary_candidate']=risk.status.eq('alive_memory_observed')&risk.date_usable
risk['interval_group']=pd.cut(risk.months_interval,[0,42,54,100],right=False,labels=['<42 months','42 to <54 months','>=54 months']).astype(str)
save(risk,'riskset2018_audit_internal.csv')
step('alive valid memory (official response year)',risk.status.eq('alive_memory_observed'))
step('date consistent primary candidates',risk.primary_candidate)
step('exclude actual development participants',risk.primary_candidate&~risk.development_overlap)
step('exclude entire development source',risk.primary_candidate&~risk.source_overlap)
save(pd.DataFrame(flow),'sample_flow.csv')
save(risk.groupby('status',dropna=False).size().rename('n').reset_index(),'followup_status.csv')
save(risk.groupby(['status','reason'],dropna=False).size().rename('n').reset_index(),'missingness_reasons.csv')
duration=[]
for name,mk in [('all',risk.primary_candidate),('exclude_development',risk.primary_candidate&~risk.development_overlap),('exclude_source',risk.primary_candidate&~risk.source_overlap)]:
    for yr in [2021,2022,2023]:
        z=risk[mk&risk.packyear.eq(yr)]
        duration.append({'population':name,'year':yr,'n':len(z),'min_months':z.months_interval.min(),'median_months':z.months_interval.median(),'max_months':z.months_interval.max(),'42_54_months':int(z.approx4y_42_54m.sum()),'45_51_months':int(z.tight4y_45_51m.sum())})
save(pd.DataFrame(duration),'followup_duration.csv')
chars=[]
for st,z in risk.groupby('status'):
    chars.append({'status':st,'n':len(z),'age_mean':z.age.mean(),'age_sd':z.age.std(),'female_n_observed':int(z.female.notna().sum()),'female_fraction':z.female.mean(),'rural_n_observed':int(z.rural.notna().sum()),'rural_fraction':z.rural.mean(),'depression_fraction':z.cesd10.ge(10).mean(),'adl_fraction':z.adl.mean(),'baseline_memory_mean':z.baseline_memory.mean(),'less_than_primary_fraction':z.edu_group.eq('Less than primary').mean()})
save(pd.DataFrame(chars),'baseline_by_followup_status.csv')

# Dedup sensitivity uses package year only; month unavailable for off-response-year records.
sens=[]
valid=a[a.memory.notna()&a.died.eq(0)].sort_values(['ID','packyear'])
for name,z in [('official_response_year',sel[sel.memory.notna()&sel.died.eq(0)]),('first_valid_package',valid.drop_duplicates('ID',keep='first')),('last_valid_package',valid.drop_duplicates('ID',keep='last'))]:
    q=risk[['ID','development_overlap','source_overlap']].merge(z[['ID','memory','packyear']],on='ID',validate='one_to_one')
    ref=sel.set_index('ID')
    matchscore=q.ID.map(ref.memory);matchyear=q.ID.map(ref.packyear)
    sens.append({'rule':name,'all_wave6_alive_valid':len(z),'riskset_matched':len(q),'exclude_development':int((~q.development_overlap).sum()),'exclude_source':int((~q.source_overlap).sum()),'changed_year_vs_official':int((matchyear.notna()&q.packyear.ne(matchyear)).sum()),'changed_score_vs_official':int((matchscore.notna()&q.memory.ne(matchscore)).sum())})
save(pd.DataFrame(sens),'dedup_sensitivity.csv')
save(sel[sel.demo_date_conflict], 'date_conflicts_internal.csv')
save(risk[risk.primary_candidate].groupby(['packyear','wordlist']).size().rename('n').reset_index(),'wordlist_by_year.csv')
save(a.groupby(['packyear','proxy','dcselect'],dropna=False).agg(n=('ID','size'),item_valid=('memory_items','count'),primary_valid=('memory','count')).reset_index(),'proxy_audit.csv')
summary={'annual':annual,'cognition_rows':len(a),'unique_ids':a.ID.nunique(),'duplicate_extra_records':len(a)-a.ID.nunique(),'duplicated_people':len(dups),'duplicate_people_valid_scores_differ':int(dups.different_valid_scores.gt(1).sum()),'off_response_year_records':int((~a.year_match).sum()),'no_matching_year_people':int((~a.ID.isin(sel.ID)).groupby(a.ID).any().sum()),'raw2018_memory_disagreements':int(disagree.sum()),'flow':flow,'followup_status':risk.status.value_counts().to_dict(),'duration':duration,'dedup':sens,'w6_primary_alive_valid':int((sel.memory.notna()&sel.died.eq(0)).sum()),'w6_adapted_alive_valid':int((sel.memory_adapted.notna()&sel.died.eq(0)).sum()),'primary_date_missing_or_conflict':int((risk.status.eq('alive_memory_observed')&~risk.date_usable).sum()),'demo_date_conflicts_selected':int(sel.demo_date_conflict.sum()),'months_distribution':risk.loc[risk.primary_candidate,'months_interval'].value_counts().sort_index().to_dict(),'matched_death_with_cognitive_score':int((sel.died.eq(1)&sel.memory.notna()).sum()),'score_only_baseline_extra':int((root.died2018.eq(0)&root.age.ge(45)&root.baseline_memory.gt(4)&root.score_predictors_complete&~root.four_predictors_complete).sum())}
(O/'qc_summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding='utf-8')
save(pd.DataFrame([{'path':str(f),'sha256':hashlib.sha256(f.read_bytes()).hexdigest()} for f in dict.fromkeys(inputs)]),'input_manifest.csv')
assert risk.ID.nunique()==len(risk)
assert risk.status.value_counts().sum()==len(risk)
assert risk.loc[risk.primary_candidate,'memory'].between(0,20).all()
print(json.dumps(summary,ensure_ascii=False,indent=2))
