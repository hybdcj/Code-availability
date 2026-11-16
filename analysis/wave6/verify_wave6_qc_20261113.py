"""Independent conservation, provenance and date checks on QC outputs."""
from pathlib import Path
import hashlib
import pandas as pd
P=Path(__file__).resolve().parents[1];O=P/'output/wave6_qc_20261113'
def get(n):return pd.read_csv(O/n,dtype={'ID':str})
a=get('cognition_records_internal.csv');s=get('selected_wave6_internal.csv');r=get('riskset2018_audit_internal.csv')
assert a.shape[0]==18870 and a.ID.nunique()==17968
assert s.ID.is_unique and r.ID.is_unique
assert not set(s.ID)-set(a.ID)
assert (s.packyear==s.res_wave).all()
assert a.loc[a.immediate_status.eq('refused')|a.delayed_status.eq('refused'),'memory'].isna().all()
assert a.loc[a.memory.notna(),'memory'].eq(a.loc[a.memory.notna(),'immediate']+a.loc[a.memory.notna(),'delayed']).all()
assert (a.loc[a.memory.notna(),'proxy'].ne(1)|a.loc[a.memory.notna(),'dcselect'].eq(7)).all()
assert r.baseline_memory.gt(4).all() and r.age.ge(45).all()
assert r[['age','edu_group','cesd10','adl']].notna().all().all()
assert set(r.edu_group)=={'Less than primary','Primary school','Middle school','High school or above'}
assert r.cesd10.between(0,30).all() and r.adl.isin([0,1]).all()
assert r.loc[r.primary_candidate,'died'].eq(0).all()
assert r.loc[r.primary_candidate,'months_interval'].eq((r.loc[r.primary_candidate,'iyear']-r.loc[r.primary_candidate,'year2018'])*12+r.loc[r.primary_candidate,'imonth']-r.loc[r.primary_candidate,'month2018']).all()
assert r.groupby('status').size().sum()==11108
assert get('baseline_memory_disagreements_internal.csv').empty
for row in get('input_manifest.csv').itertuples():
    assert hashlib.sha256(Path(row.path).read_bytes()).hexdigest()==row.sha256
print('PASS: unique IDs, refusal handling, proxy eligibility, baseline risk set, month arithmetic, status conservation, raw2018 reconstruction, input hashes')
