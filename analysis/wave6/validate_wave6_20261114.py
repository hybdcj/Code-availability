from pathlib import Path
import json,hashlib
import numpy as np
import pandas as pd
from scipy.special import expit
from scipy.optimize import brentq
import statsmodels.api as sm
from patsy import dmatrix,build_design_matrices
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

P=Path(__file__).resolve().parents[1];O=P/'output/wave6_validation_20261114';O.mkdir(exist_ok=True)
qf=P/'output/wave6_qc_20261113/riskset2018_audit_internal.csv'
pf=P/'output/wave6_events_20261114/event_predictions_internal.csv'
d=pd.read_csv(qf,dtype={'ID':str,'communityID':str})
pr=pd.read_csv(pf,dtype={'ID':str})
d=d.merge(pr[['ID','lp','p_frozen','event']],on='ID',validate='one_to_one')
d['score']=2*np.floor((d.age-45)/5)+d.edu_group.map({'Less than primary':0,'Primary school':-2,'Middle school':-7,'High school or above':-10})+3*d.cesd10.ge(10)
d['stratum']=pd.cut(d.score,[-np.inf,-3,4,np.inf],labels=['Low','Middle','High'])
assert d.score.notna().all() and np.max(abs(expit(d.lp)-d.p_frozen))<1e-12
main=d[d.primary_candidate&~d.development_overlap]
sets={f'year{y}':main[main.packyear.eq(y)].copy() for y in [2021,2022,2023]}
sets['tight4y']=main[main.tight4y_45_51m].copy()
sets['exclude_source4y']=d[d.primary_candidate&~d.source_overlap&d.approx4y_42_54m].copy()
sets['restore_dates4y']=d[d.status.eq('alive_memory_observed')&~d.development_overlap&d.months_interval.between(42,54)].copy()
assert [len(sets[f'year{y}']) for y in [2021,2022,2023]]==[4041,1188,672]

def logitfit(X,y,w):
    b=np.zeros(X.shape[1])
    for it in range(100):
        p=expit(X@b); H=X.T@((w*p*(1-p))[:,None]*X);g=X.T@(w*(y-p))
        step=np.linalg.solve(H,g);step=step/max(1,np.max(abs(step))/3)
        b+=step
        if np.max(abs(step))<1e-9:return b
    raise RuntimeError('logistic fit failed to converge')
def auc_prepare(x):
    _,inv=np.unique(x,return_inverse=True)
    def auc(y,w):
        pos=np.bincount(inv,weights=w*y);neg=np.bincount(inv,weights=w*(1-y))
        return np.sum(pos*(np.cumsum(neg)-.5*neg))/(pos.sum()*neg.sum())
    return auc
rows=[];strata=[];curves=[];bins=[];boots_all=[];checks=[]
for name,z in sets.items():
    z=z.reset_index(drop=True);y=z.event.to_numpy();lp=z.lp.to_numpy();p=expit(lp);n=len(z)
    assert z.communityID.notna().all() and np.isfinite(y).all()
    X=np.column_stack([np.ones(n),lp]);g,cluster=np.unique(z.communityID,return_inverse=True);K=len(g)
    au=auc_prepare(lp);aus=auc_prepare(z.score.to_numpy())
    strmask=[z.stratum.eq(s).to_numpy() for s in ['Low','Middle','High']]
    def metrics(w):
        a=brentq(lambda a:np.sum(w*(y-expit(lp+a))),-30,30,xtol=1e-12)
        b=logitfit(X,y,w);sw=w.sum();rate=np.sum(w*y)/sw
        A=au(y,w);S=aus(y,w)
        return np.array([A,S,A-S,a,b[1],b[0],np.sum(w*(p-y)**2)/sw,np.sum(w*(expit(lp+a)-y)**2)/sw,rate*(1-rate)]+[np.sum(w[m]*y[m])/np.sum(w[m]) if np.sum(w[m]) else np.nan for m in strmask])
    keys=['auc_continuous','auc_score','auc_difference','offset_alpha','calibration_slope','joint_gamma','brier_frozen','brier_intercept_updated_apparent','brier_null_apparent','rate_Low','rate_Middle','rate_High']
    pt=metrics(np.ones(n));f=sm.GLM(y,X,family=sm.families.Binomial()).fit();fo=sm.GLM(y,np.ones((n,1)),offset=lp,family=sm.families.Binomial()).fit()
    assert np.max(abs(f.params-pt[[5,4]]))<1e-7 and abs(fo.params[0]-pt[3])<1e-7
    checks.append(dict(sample=name,max_slope_fit_difference=float(np.max(abs(f.params-pt[[5,4]]))),offset_difference=float(abs(fo.params[0]-pt[3]))))
    design=dmatrix('bs(lp, df=4, degree=3, include_intercept=False)',{'lp':lp},return_type='dataframe')
    sx=np.asarray(design);grid=np.linspace(np.quantile(lp,.01),np.quantile(lp,.99),100)
    gx=np.asarray(build_design_matrices([design.design_info],{'lp':grid})[0]);smooth=expit(gx@logitfit(sx,y,np.ones(n)))
    fingerprint=hashlib.sha256(','.join(z.ID).encode()).hexdigest();rng=np.random.default_rng(20261114+int(fingerprint[:8],16))
    boot=[];curveboot=[];failed=0;curvefailed=0
    for b in range(2000):
        counts=np.bincount(rng.integers(0,K,K),minlength=K);w=counts[cluster].astype(float)
        try:boot.append(metrics(w))
        except (RuntimeError,np.linalg.LinAlgError,ValueError):failed+=1;continue
        if name=='year2022':
            try:curveboot.append(expit(gx@logitfit(sx,y,w)))
            except (RuntimeError,np.linalg.LinAlgError):curvefailed+=1
    bt=np.array(boot);lo,hi=np.nanquantile(bt,[.025,.975],axis=0)
    row=dict(sample=name,n=n,events=int(y.sum()),clusters=K,bootstrap_ok=len(bt),bootstrap_failed=failed,curve_bootstrap_failed=curvefailed)
    for j,k in enumerate(keys):row[k]=float(pt[j]);row[k+'_low']=float(lo[j]);row[k+'_high']=float(hi[j])
    row['brier_change_apparent']=float(pt[7]-pt[6]);delta=bt[:,7]-bt[:,6];row['brier_change_apparent_low'],row['brier_change_apparent_high']=np.quantile(delta,[.025,.975])
    if name in ['year2021','year2023']:
        for key in list(row):
            if key.startswith(('brier_intercept_updated','brier_change_apparent')):row[key]=np.nan
    rows.append(row)
    for j,s in enumerate(['Low','Middle','High']):
        m=strmask[j];strata.append(dict(sample=name,stratum=s,n=int(m.sum()),events=int(y[m].sum()),rate=float(pt[9+j]),low=float(lo[9+j]),high=float(hi[9+j])))
    c=pd.DataFrame({'sample':name,'predicted':expit(grid),'observed_spline':smooth})
    if curveboot:c['lower'],c['upper']=np.quantile(curveboot,[.025,.975],axis=0)
    curves.append(c)
    z['decile']=pd.qcut(z.p_frozen,10,duplicates='drop')
    bi=z.groupby('decile',observed=True).agg(predicted=('p_frozen','mean'),observed=('event','mean'),n=('ID','size'));bi['sample']=name;bins.append(bi.reset_index(drop=True))
    btframe=pd.DataFrame(bt,columns=keys);btframe['sample']=name;boots_all.append(btframe)
    if name in ['year2021','year2023']:btframe['brier_intercept_updated_apparent']=np.nan
    if name=='year2022':
        z['p_updated_apparent']=expit(lp+pt[3]);z[['ID','communityID','event','lp','p_frozen','score','stratum','p_updated_apparent']].to_csv(O/'fouryear_predictions_internal.csv',index=False)
    print(name,'n=',n,'auc=',round(pt[0],4),'alpha=',round(pt[3],4),'slope=',round(pt[4],4),'bootstrap=',len(bt),flush=True)
res=pd.DataFrame(rows);res.to_csv(O/'performance.csv',index=False,encoding='utf-8-sig')
pd.DataFrame(strata).to_csv(O/'score_strata.csv',index=False,encoding='utf-8-sig')
cv=pd.concat(curves);bn=pd.concat(bins);cv.to_csv(O/'calibration_curves.csv',index=False);bn.to_csv(O/'calibration_bins.csv',index=False)
pd.concat(boots_all).to_csv(O/'bootstrap_metrics.csv',index=False)
pd.DataFrame(checks).to_csv(O/'independent_fit_checks.csv',index=False)
fig,axes=plt.subplots(1,3,figsize=(12.5,4.2),sharex=True,sharey=True)
for ax,year in zip(axes,[2021,2022,2023]):
    name=f'year{year}';c=cv[cv['sample'].eq(name)];b=bn[bn['sample'].eq(name)]
    ax.plot([0,.7],[0,.7],ls='--',color='gray',label='Ideal')
    ax.plot(c.predicted,c.observed_spline,color='#14667a',lw=2,label='Spline calibration')
    if year==2022:ax.fill_between(c.predicted,c.lower,c.upper,color='#14667a',alpha=.15,label='95% pointwise cluster CI')
    ax.scatter(b.predicted,b.observed,s=28,color='#b96b33',label='Grouped observations',zorder=3)
    ax.set_title(f'2018 to {year} (~{year-2018} years)');ax.set_xlabel('Frozen 4-year predicted probability');ax.set_xlim(0,.65);ax.set_ylim(0,.75);ax.grid(alpha=.15)
axes[0].set_ylabel('Observed low-memory proportion');axes[1].legend(fontsize=7,loc='upper left')
fig.suptitle('Temporal transport of frozen model: different follow-up durations',fontsize=12)
fig.tight_layout();fig.savefig(O/'calibration_by_year.png',dpi=200);fig.savefig(O/'calibration_by_year.pdf');plt.close(fig)
manifest=[qf,pf,O/'analysis_protocol.md',Path(__file__)]
pd.DataFrame([dict(path=str(f),sha256=hashlib.sha256(f.read_bytes()).hexdigest()) for f in manifest]).to_csv(O/'input_manifest.csv',index=False)
(O/'summary.json').write_text(json.dumps(res.replace({np.nan:None}).to_dict('records'),ensure_ascii=False,indent=2,allow_nan=False),encoding='utf-8')
print(res[['sample','auc_continuous','auc_score','offset_alpha','calibration_slope','brier_frozen','brier_intercept_updated_apparent']].to_string(index=False),flush=True)
