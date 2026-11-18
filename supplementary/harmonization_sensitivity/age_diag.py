import os
import pandas as pd, numpy as np
from scipy.stats import spearmanr
BASE=os.environ["PROJECT_DIR"] + "/hcp_tractometry"
bins=[(22,35,'22-35'),(36,50,'36-50'),(51,65,'51-65'),(66,200,'66+')]

def load(path,cols=None):
    return pd.read_csv(path,usecols=cols,low_memory=False)

def diag(zpath,label,agecol_from_profiles=True):
    z=load(zpath,['subject','tract','scalar','z_score','cohort'])
    z=z[z.scalar.isin(['dti_fa','dti_md'])]
    prof=load(f'{BASE}/combined_profiles.csv',['subject','age','cohort','combat_reference'])
    pm=prof.groupby('subject').first()
    pm['age_n']=pd.to_numeric(pm['age'].replace('90 or older','90'),errors='coerce')
    sub=z.groupby(['subject','scalar'],as_index=False).z_score.mean()
    sub=sub.merge(pm[['age_n','combat_reference','cohort']],left_on='subject',right_index=True)
    ref=sub[sub.combat_reference=='yes']
    out=[]
    for lo,hi,lbl in bins:
        b=ref[(ref.age_n>=lo)&(ref.age_n<=hi)]
        for s in ['dti_fa','dti_md']:
            v=b[b.scalar==s].z_score
            out.append({'version':label,'age_bin':lbl,'scalar':s,'n':len(v),
                        'median_z':round(v.median(),4) if len(v) else np.nan})
    return pd.DataFrame(out)

res=[]
res.append(diag(f'{BASE}/node_zscores.csv','v1_combat_as_run'))
res.append(diag('node_zscores_unharmonized.csv','v2_no_combat'))
pd.concat(res).to_csv('age_diag_2a_partial.csv',index=False)
print(pd.concat(res).pivot_table(index=['version','age_bin'],columns='scalar',values='median_z').to_string())
