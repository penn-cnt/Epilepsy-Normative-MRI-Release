"""Version 3: neuroCombat, batch=site, age (continuous) + sex (categorical) PROTECTED,
estimates fit on the 2,455 reference subjects only; remaining clinical subjects
transformed via neuroCombatFromTraining. F_L imputed only so ComBat can run, then restored to NaN."""
import os
import pandas as pd, numpy as np
from neuroCombat import neuroCombat, neuroCombatFromTraining
BASE=os.environ["PROJECT_DIR"] + "/hcp_tractometry"
OUT=os.environ["PROJECT_DIR"] + "/hcp_tractometry/harmonization_sensitivity"
scalars=['dti_fa','dti_md']

df=pd.read_csv(f'{BASE}/combined_profiles.csv',low_memory=False)
tracts=sorted(df.tractID.unique()); nodes=sorted(df.nodeID.unique())
meta=df.groupby('subject').first()[['cohort','site','age','sex','group_1','group_2','combat_reference']].reset_index()
meta['age_numeric']=pd.to_numeric(meta['age'].replace('90 or older','90'),errors='coerce')
meta['sex_numeric']=(meta['sex']=='F').astype(int)
subject_order=meta['subject'].values
ref_mask=(meta['combat_reference']=='yes').values
print(f"subjects {len(meta)} | reference {ref_mask.sum()} | non-reference {(~ref_mask).sum()}")
print("sites:",sorted(meta.site.unique()))

harm={}
for scalar in scalars:
    feats,names=[],[]
    for tract in tracts:
        piv=df[df.tractID==tract].pivot_table(index='subject',columns='nodeID',values=scalar).reindex(subject_order)
        for node in nodes:
            if node in piv.columns:
                feats.append(piv[node].values); names.append(f"{tract}_node{node}_{scalar}")
    data=np.array(feats,dtype=float)
    nan_mask=np.isnan(data)                      # remember F_L gaps
    print(f"{scalar}: {data.shape[0]} features x {data.shape[1]} subjects | NaN cells {nan_mask.sum()}")
    for i in range(data.shape[0]):               # impute ONLY so ComBat can run
        if nan_mask[i].any(): data[i,nan_mask[i]]=np.nanmean(data[i])

    covars=pd.DataFrame({'batch':meta['site'].values,
                         'age':meta['age_numeric'].values,
                         'sex':meta['sex_numeric'].values})
    # --- fit on reference subjects only ---
    res=neuroCombat(dat=data[:,ref_mask], covars=covars[ref_mask].reset_index(drop=True),
                    batch_col='batch', continuous_cols=['age'], categorical_cols=['sex'])
    out=np.empty_like(data); out[:,ref_mask]=res['data']
    # --- apply to held-out clinical subjects ---
    ho=neuroCombatFromTraining(dat=data[:,~ref_mask], batch=meta.loc[~ref_mask,'site'].values,
                               estimates=res['estimates'])
    out[:,~ref_mask]=ho['data']
    out[nan_mask]=np.nan                          # restore F_L missingness
    for i,n in enumerate(names): harm[n]=out[i]

rows=[]
for i,sub in enumerate(subject_order):
    sm=meta.iloc[i]
    for tract in tracts:
        for node in nodes:
            r={'subject':sub,'cohort':sm.cohort,'site':sm.site,'age':sm.age,'sex':sm.sex,
               'group_1':sm.group_1,'group_2':sm.group_2,'combat_reference':sm.combat_reference,
               'tractID':tract,'nodeID':node}
            for s in scalars:
                k=f"{tract}_node{node}_{s}"; r[s]=harm[k][i] if k in harm else np.nan
            rows.append(r)
pd.DataFrame(rows).to_csv(f'{OUT}/combined_profiles_harmonized_v3.csv',index=False)
print("saved combined_profiles_harmonized_v3.csv")
