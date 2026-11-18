import os
import pandas as pd, numpy as np
from scipy.stats import spearmanr
BASE=os.environ["PROJECT_DIR"] + "/hcp_tractometry"
bins=[(22,35,'22-35'),(36,50,'36-50'),(51,65,'51-65'),(66,200,'66+')]
prof=pd.read_csv(f'{BASE}/combined_profiles.csv',usecols=['subject','age','cohort','combat_reference'],low_memory=False)
pm=prof.groupby('subject').first()
pm['age_n']=pd.to_numeric(pm['age'].replace('90 or older','90'),errors='coerce')
VERS=[('ComBat as run',f'{BASE}/node_zscores.csv'),
      ('No ComBat','node_zscores_unharmonized.csv'),
      ('ComBat + age/sex, ref-only','node_zscores_v3.csv')]
out=[]
for label,path in VERS:
    z=pd.read_csv(path,usecols=['subject','scalar','z_score','cohort'],low_memory=False)
    z=z[z.scalar.isin(['dti_fa','dti_md'])]
    s=z.groupby(['subject','scalar'],as_index=False).z_score.mean()
    s=s.merge(pm[['age_n','combat_reference']],left_on='subject',right_index=True)
    ref=s[s.combat_reference=='yes']
    for lo,hi,lbl in bins:
        b=ref[(ref.age_n>=lo)&(ref.age_n<=hi)]
        for sc in ['dti_fa','dti_md']:
            v=b[b.scalar==sc].z_score
            out.append({'version_label':label,'age_bin':lbl,'scalar':sc,'n':len(v),
                        'median_z':round(v.median(),4) if len(v) else np.nan})
d=pd.DataFrame(out); d.to_csv('age_diag_2a.csv',index=False)
print(d.pivot_table(index=['version_label','age_bin'],columns='scalar',values='median_z').to_string())
