import pandas as pd, numpy as np, sys
from scipy.stats import spearmanr
sys.path.insert(0,'.')
from common import BASE, META, mwu_auc, bootstrap_auc_ci, fold_auc, morph_cohort
morph=morph_cohort(); subs=set(morph.SubjectID)
m=pd.read_csv(META).rename(columns={'record_id':'SubjectID'}).drop_duplicates('SubjectID',keep='first')
ages=m.assign(age_n=pd.to_numeric(m.age_MRI,errors='coerce')).set_index('SubjectID')['age_n']
rows=[];aucs=[]
for lab,path in [('ComBat as run',f'{BASE}/node_zscores.csv'),
                 ('No ComBat','node_zscores_unharmonized.csv'),
                 ('ComBat + age/sex, ref-only','node_zscores_v3.csv')]:
    z=pd.read_csv(path,usecols=['subject','scalar','z_score','cohort'],low_memory=False)
    z=z[(z.cohort=='Clinical')&(z.scalar=='dti_md')&(z.subject.isin(subs))]
    s=z.groupby('subject',as_index=False).z_score.mean().rename(columns={'z_score':'md_z'})
    s['age']=s.subject.map(ages)
    s=s.merge(morph.rename(columns={'SubjectID':'subject'})[['subject','group_1','group_2']],on='subject')
    r,p=spearmanr(s.age.dropna(),s.loc[s.age.notna(),'md_z'])
    rows.append({'version':lab,'spearman_age_vs_MDz':round(r,3),'p':f"{p:.2g}",'n':int(s.age.notna().sum())})
    d=s.dropna(subset=['age','md_z']).copy()
    c=np.polyfit(d.age,d.md_z,1)
    for how,dd,col in [('raw',d,'md_z'),('age-residualized',d.assign(res=d.md_z-np.polyval(c,d.age)),'res')]:
        epi=dd[dd.group_2.isin(['drug-resistant','drug-responsive'])]
        a=epi[epi.group_1=='focal'][col].dropna().values; b=epi[epi.group_1=='generalized'][col].dropna().values
        auc,pp=mwu_auc(a,b); lo,hi=bootstrap_auc_ci(a,b); A,L,H=fold_auc(auc,lo,hi)
        aucs.append({'version':lab,'method':how,'AUC':round(A,3),'CI':f"({L:.3f}-{H:.3f})",'p':round(pp,4),'n_focal':len(a),'n_gen':len(b)})
    gg=d[(d.group_1=='generalized')&(d.group_2.isin(['drug-resistant','drug-responsive']))]
    lo_,hi_=gg.age.quantile([0.05,0.95]); mm=d[(d.age>=lo_)&(d.age<=hi_)]
    epi=mm[mm.group_2.isin(['drug-resistant','drug-responsive'])]
    a=epi[epi.group_1=='focal'].md_z.dropna().values; b=epi[epi.group_1=='generalized'].md_z.dropna().values
    auc,pp=mwu_auc(a,b); lo,hi=bootstrap_auc_ci(a,b); A,L,H=fold_auc(auc,lo,hi)
    aucs.append({'version':lab,'method':f'age-matched({lo_:.0f}-{hi_:.0f}y)','AUC':round(A,3),'CI':f"({L:.3f}-{H:.3f})",'p':round(pp,4),'n_focal':len(a),'n_gen':len(b)})
pd.DataFrame(rows).to_csv('age_diag_2b_spearman.csv',index=False)
pd.DataFrame(aucs).to_csv('age_diag_2c_focal_gen_md.csv',index=False)
print(pd.DataFrame(rows).to_string(index=False)); print(); print(pd.DataFrame(aucs).to_string(index=False))
