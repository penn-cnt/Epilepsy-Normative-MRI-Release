"""2b: Spearman age vs tract-mean MD z within the 355 patients.
   2c: focal vs generalized MD AUC with age accounted for (residualized + age-matched)."""
import pandas as pd, numpy as np, sys
from scipy.stats import spearmanr, mannwhitneyu
sys.path.insert(0,'.')
from common import BASE, META, bilateral_tracts, tract_order, mwu_auc, bootstrap_auc_ci, fold_auc, morph_cohort
from statsmodels.stats.multitest import multipletests

morph=morph_cohort(); subs=set(morph.SubjectID)
m=pd.read_csv(META).rename(columns={'record_id':'SubjectID'}).drop_duplicates('SubjectID',keep='first')
m['age_n']=pd.to_numeric(m.age_MRI,errors='coerce')
ages=m.set_index('SubjectID')['age_n']

def tractmean_md(zpath,label):
    z=pd.read_csv(zpath,usecols=['subject','tract','scalar','z_score','cohort'],low_memory=False)
    z=z[(z.cohort=='Clinical')&(z.scalar=='dti_md')&(z.subject.isin(subs))]
    s=z.groupby('subject',as_index=False).z_score.mean().rename(columns={'z_score':'md_z'})
    s['age']=s.subject.map(ages)
    s=s.merge(morph.rename(columns={'SubjectID':'subject'})[['subject','group_1','group_2']],on='subject')
    r,p=spearmanr(s.age.dropna(),s.loc[s.age.notna(),'md_z'])
    print(f"2b {label}: Spearman(age, tract-mean MD z) r={r:+.3f} p={p:.2g} n={s.age.notna().sum()}")
    return s,r,p

def auc_fg(s,label,how):
    epi=s[s.group_2.isin(['drug-resistant','drug-responsive'])].copy()
    f=epi[epi.group_1=='focal']; g=epi[epi.group_1=='generalized']
    a,b=f.val.dropna().values,g.val.dropna().values
    auc,p=mwu_auc(a,b); lo,hi=bootstrap_auc_ci(a,b)
    A,L,H=fold_auc(auc,lo,hi)
    print(f"2c {label} [{how}]: AUC={A:.3f} ({L:.3f}-{H:.3f}) p={p:.4f} nF={len(a)} nG={len(b)}")
    return {'version':label,'method':how,'AUC':round(A,3),'CI':f"({L:.3f}-{H:.3f})",'p':round(p,4),'n_focal':len(a),'n_gen':len(b)}

rows=[]
for zp,lab in [(f'{BASE}/node_zscores.csv','v1_combat_as_run'),('node_zscores_unharmonized.csv','v2_no_combat')]:
    s,r,p=tractmean_md(zp,lab)
    rows.append({'version':lab,'spearman_r':round(r,3),'p':f"{p:.2g}",'n':int(s.age.notna().sum())})
    # raw
    s['val']=s.md_z; auc_fg(s,lab,'raw')
    # residualize on age within patients
    d=s.dropna(subset=['age','md_z']).copy()
    c=np.polyfit(d.age,d.md_z,1); d['val']=d.md_z-np.polyval(c,d.age)
    auc_fg(d,lab,'age-residualized')
    # age-matched subset: restrict both groups to overlapping age range of generalized
    gg=d[(d.group_2.isin(['drug-resistant','drug-responsive']))&(d.group_1=='generalized')]
    lo,hi=gg.age.quantile([0.05,0.95])
    mm=d[(d.age>=lo)&(d.age<=hi)].copy(); mm['val']=mm.md_z
    auc_fg(mm,lab,f'age-matched({lo:.0f}-{hi:.0f}y)')
pd.DataFrame(rows).to_csv('age_diag_2b_partial.csv',index=False)
