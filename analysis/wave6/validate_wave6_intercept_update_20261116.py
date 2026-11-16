"""Cluster bootstrap optimism correction for a fixed intercept-only update."""
from pathlib import Path
import hashlib,json
import numpy as np
import pandas as pd
from scipy.special import expit
from scipy.optimize import brentq

P=Path(__file__).resolve().parents[1]
O=P/'output/wave6_intercept_validation_20261116'
O.mkdir(exist_ok=True)
F=P/'output/wave6_validation_20261114/fouryear_predictions_internal.csv'
PERF=P/'output/wave6_validation_20261114/performance.csv'
DEV=P/'output/wave6_qc_20261113/development_ids_internal.csv'
sources=[F,PERF,DEV,O/'analysis_protocol.md',Path(__file__)]
before={str(f):hashlib.sha256(f.read_bytes()).hexdigest() for f in sources}
d=pd.read_csv(F,dtype={'ID':str,'communityID':str})
assert len(d)==1188 and d.event.sum()==269 and d.ID.is_unique
assert d[['communityID','event','lp','p_frozen']].notna().all().all()
assert set(d.event.unique())=={0,1}
assert not set(d.ID)&set(pd.read_csv(DEV,dtype=str).KEY)
y=d.event.to_numpy(float); lp=d.lp.to_numpy(float); p=expit(lp)
assert np.max(abs(p-d.p_frozen))<1e-12
clusters,ix=np.unique(d.communityID,return_inverse=True); K=len(clusters)
assert K==212
def fit(w):
    a=brentq(lambda a:np.sum(w*(y-expit(lp+a))),-30,30,xtol=1e-12)
    assert abs(np.sum(w*(y-expit(lp+a))))<1e-7
    return a
one=np.ones(len(d));alpha=fit(one);pa=expit(lp+alpha)
frozen=float(np.mean((p-y)**2));apparent=float(np.mean((pa-y)**2))
saved=pd.read_csv(PERF).set_index('sample').loc['year2022']
assert abs(alpha-saved.offset_alpha)<1e-10
assert abs(frozen-saved.brier_frozen)<1e-12
assert abs(apparent-saved.brier_intercept_updated_apparent)<1e-12
assert np.max(abs(pa-d.p_updated_apparent))<1e-10
rng=np.random.default_rng(20261116); B=10000
counts=np.empty((B,K),dtype=np.int16);records=[];failures=[]
for b in range(B):
    count=np.bincount(rng.integers(K,size=K),minlength=K);counts[b]=count
    w=count[ix].astype(float)
    try:
        a=fit(w);loss=(expit(lp+a)-y)**2
        train=float(np.average(loss,weights=w));test=float(loss.mean())
        records.append(dict(replicate=b+1,n_boot=int(w.sum()),events_boot=int(w@y),
            alpha=a,brier_train=train,brier_test=test,correction=test-train,
            frozen_brier_train=float(np.average((p-y)**2,weights=w))))
    except (ValueError,AssertionError) as e:
        failures.append(dict(replicate=b+1,error=str(e)))
    if (b+1)%2000==0:print('completed',b+1,flush=True)
r=pd.DataFrame(records)
assert not failures, failures
penalty=float(r.correction.mean());mcse=float(r.correction.std(ddof=1)/np.sqrt(len(r)))
corrected=apparent+penalty
summary=dict(n=len(d),events=int(y.sum()),clusters=K,bootstrap_requested=B,
    bootstrap_success=len(r),bootstrap_failed=len(failures),seed=20261116,
    alpha=alpha,brier_original_frozen=frozen,brier_updated_apparent=apparent,
    optimism_loss_increment=penalty,optimism_increment_mcse=mcse,
    brier_updated_corrected=corrected,delta_apparent=apparent-frozen,
    delta_corrected=corrected-frozen,
    correction_as_fraction_apparent_gain=penalty/(frozen-apparent),
    frozen_bootstrap_mean_shift=float(r.frozen_brier_train.mean()-frozen),
    n_boot_min=int(r.n_boot.min()),n_boot_max=int(r.n_boot.max()),
    sampling_ci_computed=False,validation_type='cluster bootstrap internal validation')
r.to_csv(O/'bootstrap_results.csv',index=False)
np.savez_compressed(O/'bootstrap_counts_internal.npz',counts=counts,cluster_ids=clusters)
pd.DataFrame(counts[:50],columns=[f'g{i+1}' for i in range(K)]).to_csv(O/'first50_counts_internal.csv',index=False)
pd.DataFrame({'communityID':clusters,'cluster_index':np.arange(1,K+1)}).to_csv(O/'cluster_index_internal.csv',index=False)
pd.DataFrame([summary]).to_csv(O/'performance.csv',index=False,encoding='utf-8-sig')
(O/'summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2,allow_nan=False),encoding='utf-8')
assert before=={str(f):hashlib.sha256(f.read_bytes()).hexdigest() for f in sources}
pd.DataFrame([{'path':f,'sha256':h} for f,h in before.items()]).to_csv(O/'input_manifest.csv',index=False)
print(json.dumps(summary,ensure_ascii=False,indent=2),flush=True)
