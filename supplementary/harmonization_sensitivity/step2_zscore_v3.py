"""Replicate run_zscores_all.sh's per-node z-score step on UNHARMONIZED input.
Identical algorithm: OLS on age+sex within combat_reference=='yes', z = (y - yhat)/resid_sd,
min 10 reference subjects per tract-node-scalar. Only the input file differs."""
import os
import pandas as pd, numpy as np
from numpy.linalg import lstsq
BASE=os.environ["PROJECT_DIR"] + "/hcp_tractometry"
scalars=['dti_fa','dti_md']   # only FA/MD needed downstream

print("loading unharmonized profiles...")
df=pd.read_csv(f'{BASE}/harmonization_sensitivity/combined_profiles_harmonized_v3.csv',low_memory=False)
# match run_zscores_all.sh numeric coercion
df['age_numeric']=pd.to_numeric(df['age'].replace('90 or older','90'),errors='coerce')
df['sex_numeric']=(df['sex']=='F').astype(int)
df.rename(columns={'tractID':'tract'},inplace=True)
print(f"  {df['subject'].nunique()} subjects, {len(df)} rows")

ref=df[df['combat_reference']=='yes']
tracts=sorted(df['tract'].unique()); nodes=sorted(df['nodeID'].unique())
out=[]
for tract in tracts:
    dt=df[df['tract']==tract]; rt=ref[ref['tract']==tract]
    for node in nodes:
        dn=dt[dt['nodeID']==node]; rn=rt[rt['nodeID']==node]
        for scalar in scalars:
            r=rn.dropna(subset=[scalar,'age_numeric','sex_numeric'])
            if len(r)<10: continue
            X=np.column_stack([r['age_numeric'].values,r['sex_numeric'].values,np.ones(len(r))])
            c,_,_,_=lstsq(X,r[scalar].values,rcond=None)
            sd=np.std(r[scalar].values-X@c)
            s=dn[['subject','cohort','site','group_1','group_2',scalar,'age_numeric','sex_numeric']].copy()
            z=(s[scalar].values-np.column_stack([s.age_numeric.values,s.sex_numeric.values,np.ones(len(s))])@c)/sd
            out.append(pd.DataFrame({'subject':s.subject.values,'tract':tract,'nodeID':node,
                'scalar':scalar,'z_score':z,'cohort':s.cohort.values,'group_1':s.group_1.values,
                'group_2':s.group_2.values}))
    print(f"  {tract} done")
z=pd.concat(out,ignore_index=True)
z.to_csv('node_zscores_v3.csv',index=False)
print(f"saved {z.shape} -> node_zscores_v3.csv")
