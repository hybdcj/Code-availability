from pathlib import Path
import hashlib,json
import numpy as np
import pandas as pd
from scipy.special import expit
from statsmodels.stats.proportion import proportion_confint

P=Path(__file__).resolve().parents[1];O=P/'output/wave6_events_20261114'
qc=P/'output/wave6_qc_20261113/riskset2018_audit_internal.csv'
cf=P/'output/12_模型定稿/A14_最终模型系数.csv'
r=pd.read_csv(qc,dtype={'ID':str,'communityID':str})
dev=pd.read_csv(O/'development_internal.csv',dtype={'KEY':str,'communityID':str})
beta=pd.read_csv(cf).set_index('term').beta
r['lp']=beta['(Intercept)']+beta['age']*r.age+beta['dep']*r.cesd10.ge(10)+beta['adl']*r.adl
for edu in ['Primary school','Middle school','High school or above']:
    r['lp']+=beta['edu_group'+edu]*r.edu_group.eq(edu)
r['p_frozen']=expit(r.lp)
ref=pd.read_csv(O/'prediction_reference_internal.csv',dtype={'ID':str})
ck=r[['ID','p_frozen']].merge(ref,on='ID',validate='one_to_one')
assert np.max(abs(ck.p_frozen-ck.p_R))<1e-12
r['event']=r.memory.le(4).astype(float).where(r.memory.notna())
assert r.ID.is_unique and len(r)==11108
assert r.loc[r.status.eq('dead_at_recorded_wave'),'event'].isna().all()
assert dev.incident_2015.sum()==658 and len(dev)==4441
main=r[r.primary_candidate & ~r.development_overlap].copy()
assert len(main)==5901 and not set(main.ID)&set(dev.KEY)
def stats(d,name):
    assert d.event.notna().all() and d.communityID.notna().all()
    n=len(d);e=int(d.event.sum());lo,hi=proportion_confint(e,n,alpha=.05,method='wilson')
    cl=d.groupby('communityID').event.agg(['count','sum']);k=len(cl)
    # Identical samples get identical uncertainty intervals regardless of table order.
    sample_hash=hashlib.sha256(cl.to_csv().encode()).hexdigest()
    rng=np.random.default_rng(20261114+int(sample_hash[:8],16))
    ix=rng.integers(0,k,size=(2000,k));boots=cl['sum'].to_numpy()[ix].sum(axis=1)/cl['count'].to_numpy()[ix].sum(axis=1)
    blo,bhi=np.quantile(boots,[.025,.975])
    return dict(population=name,n=n,events=e,nonevents=n-e,rate=e/n,wilson_low=lo,wilson_high=hi,cluster_low=blo,cluster_high=bhi,clusters=k,predicted4y_mean=d.p_frozen.mean(),predicted4y_median=d.p_frozen.median(),memory_mean=d.memory.mean() if 'memory' in d else np.nan)
rows=[stats(dev.rename(columns={'incident_2015':'event'}),'development_2011_2015')]
for y in [2021,2022,2023]:rows.append(stats(main[main.packyear.eq(y)],f'nonoverlap_{y}'))
rows.append(stats(main,'nonoverlap_all_DESCRIPTIVE_MIXED_HORIZONS'))
rows.append(stats(main[main.approx4y_42_54m],'nonoverlap_42_54m'))
rows.append(stats(main[main.tight4y_45_51m],'nonoverlap_45_51m'))
strict=r[r.primary_candidate&~r.source_overlap]
for y in [2021,2022,2023]:rows.append(stats(strict[strict.packyear.eq(y)],f'exclude_source_{y}'))
restore=r[r.status.eq('alive_memory_observed')&~r.development_overlap]
for y in [2021,2022,2023]:rows.append(stats(restore[restore.packyear.eq(y)],f'date_conflicts_retained_{y}'))
rates=pd.DataFrame(rows)
rates['rate_difference_vs_development_pp']=100*(rates.rate-658/4441)
rates.to_csv(O/'event_rates.csv',index=False,encoding='utf-8-sig')
comp=[]
for scope,z in [('all_baseline',r),('exclude_development',r[~r.development_overlap])]:
    for st,d in z.groupby('status'):
        comp.append(dict(scope=scope,status=st,n=len(d),age_mean=d.age.mean(),depression_fraction=d.cesd10.ge(10).mean(),adl_fraction=d.adl.mean(),predicted4y_mean=d.p_frozen.mean(),predicted4y_median=d.p_frozen.median(),predicted4y_p25=d.p_frozen.quantile(.25),predicted4y_p75=d.p_frozen.quantile(.75)))
pd.DataFrame(comp).to_csv(O/'baseline_predicted_risk_by_followup.csv',index=False,encoding='utf-8-sig')
word=[]
for (year,wl),d in main.groupby(['packyear','wordlist']):
    word.append(dict(year=int(year),wordlist=wl,n=len(d),events=int(d.event.sum()),rate=d.event.mean(),predicted4y_mean=d.p_frozen.mean()))
pd.DataFrame(word).to_csv(O/'wordlist_event_counts.csv',index=False,encoding='utf-8-sig')
r[['ID','event','p_frozen','lp','status','primary_candidate','development_overlap','source_overlap','packyear','months_interval']].to_csv(O/'event_predictions_internal.csv',index=False)
summary=dict(rates=rows,baseline_risk=comp,wordlist=word,max_probability_difference_R=float(abs(ck.p_frozen-ck.p_R).max()))
(O/'summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding='utf-8')
files=[qc,cf,O/'development_internal.csv',O/'prediction_reference_internal.csv',O/'analysis_protocol.md',Path(__file__)]
pd.DataFrame([dict(file=str(f),sha256=hashlib.sha256(f.read_bytes()).hexdigest()) for f in files]).to_csv(O/'input_manifest.csv',index=False)
assert sum(x['events'] for x in rows[1:4])==int(main.event.sum())
assert sum(x['n'] for x in rows[1:4])==5901
print(rates.round(4).to_string(index=False))
print(pd.DataFrame(comp).round(4).to_string(index=False))
