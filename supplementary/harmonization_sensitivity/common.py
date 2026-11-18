import os
import pandas as pd, numpy as np
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests
RNG_SEED=42
BASE=os.environ["PROJECT_DIR"] + "/hcp_tractometry"
META=os.environ["BIDS_ROOT"] + "/metadata/metadata_drug_responsiveness.csv"

bilateral_tracts={'Cingulum':['C_L','C_R'],'Fornix':['F_L','F_R'],'Uncinate':['UF_L','UF_R'],
 'Arcuate':['AF_L','AF_R'],'IFOF':['IFOF_L','IFOF_R'],'CC Forceps Major':['CC_ForcepsMajor'],
 'CC Forceps Minor':['CC_ForcepsMinor'],'CC Mid':['CCMid'],'CST':['CST_L','CST_R']}
tract_order=['Cingulum','Fornix','Uncinate','Arcuate','IFOF','CC Forceps Major','CC Forceps Minor','CC Mid','CST']
comparisons=[('Epilepsy','FS','Epilepsy','Functional Seizures'),
             ('Focal','Generalized','Focal','Generalized'),
             ('DR','DS','Drug-Resistant','Drug-Responsive')]

def morph_cohort():
    cb_dir=f'{BASE}/centilebrain_zscores'
    ct=pd.concat([pd.read_csv(f'{cb_dir}/zscore_CorticalThickness_male.csv'),
                  pd.read_csv(f'{cb_dir}/zscore_CorticalThickness_female.csv')],ignore_index=True)
    sv=pd.concat([pd.read_csv(f'{cb_dir}/zscore_SubcorticalVolume_male.csv'),
                  pd.read_csv(f'{cb_dir}/zscore_SubcorticalVolume_female.csv')],ignore_index=True)
    both=set(ct['SubjectID'])&set(sv['SubjectID'])
    ct=ct[ct['SubjectID'].isin(both)]
    meta=pd.read_csv(META).rename(columns={'record_id':'SubjectID'})
    m=ct[['SubjectID']].merge(meta[['SubjectID','group_1','group_2']],on='SubjectID',how='left')
    m=m.drop_duplicates(subset='SubjectID',keep='first')
    tract_subs=set(pd.read_csv(f'{BASE}/node_zscores.csv',usecols=['subject'])['subject'].unique())
    return m[m['SubjectID'].isin(tract_subs)].copy()

def build_wide(zdf, morph, bilat=None):
    """zdf: long per-node z-scores with subject,tract,scalar,z_score (cohort/Clinical filtered upstream)"""
    bilat = bilat or bilateral_tracts
    clin=zdf[zdf['scalar'].isin(['dti_fa','dti_md'])].copy()
    clin=clin[clin['subject'].isin(set(morph['SubjectID']))]
    clin['feature']=clin['tract']+'_'+clin['scalar']
    wide=clin.pivot_table(index='subject',columns='feature',values='z_score').reset_index()
    wide=wide.merge(morph.rename(columns={'SubjectID':'subject'})[['subject','group_1','group_2']],
                    on='subject',how='left').drop_duplicates(subset='subject',keep='first')
    for name,tracts in bilat.items():
        for scalar in ['dti_fa','dti_md']:
            col=f'zb_{name}_{scalar}'.replace(' ','_')
            fc=[f'{t}_{scalar}' for t in tracts]
            miss=[c for c in fc if c not in wide.columns]
            assert not miss, f'{name}: missing {miss}'
            wide[col]=wide[fc].mean(axis=1)
    return wide

def groups(wide):
    epi=wide[wide['group_2'].isin(['drug-resistant','drug-responsive'])].copy()
    return {'Epilepsy':epi,'FS':wide[wide['group_1']=='pnee'].copy(),
            'Focal':epi[epi['group_1']=='focal'].copy(),
            'Generalized':epi[epi['group_1']=='generalized'].copy(),
            'DR':epi[epi['group_2']=='drug-resistant'].copy(),
            'DS':epi[epi['group_2']=='drug-responsive'].copy()}

def mwu_auc(a,b):
    U,p=mannwhitneyu(a,b,alternative='two-sided'); return U/(len(a)*len(b)),p

def bootstrap_auc_ci(a,b,n_boot=5000,seed=RNG_SEED):
    rng=np.random.default_rng(seed); na,nb=len(a),len(b); aucs=np.empty(n_boot)
    for k in range(n_boot):
        U,_=mannwhitneyu(rng.choice(a,na,True),rng.choice(b,nb,True),alternative='two-sided')
        aucs[k]=U/(na*nb)
    return np.percentile(aucs,2.5),np.percentile(aucs,97.5)

def auc_family(ga,gb,cols):
    aucs,los,his,ps=[],[],[],[]
    for col in cols:
        a=ga[col].dropna().values; b=gb[col].dropna().values
        auc,p=mwu_auc(a,b); lo,hi=bootstrap_auc_ci(a,b)
        aucs.append(auc); los.append(lo); his.append(hi); ps.append(p)
    return {'auc':np.array(aucs),'ci_lo':np.array(los),'ci_hi':np.array(his),'p':np.array(ps)}

def tract_auc_family(ga,gb,torder):
    fa_cols=[f'zb_{t}_dti_fa'.replace(' ','_') for t in torder]
    md_cols=[f'zb_{t}_dti_md'.replace(' ','_') for t in torder]
    fa=auc_family(ga,gb,fa_cols); md=auc_family(ga,gb,md_cols)
    allp=np.concatenate([fa['p'],md['p']])
    q=multipletests(allp,method='fdr_bh')[1]
    n=len(torder); fa['q'],md['q']=q[:n],q[n:]
    return fa,md

def fold_auc(auc,lo,hi):
    return (auc,lo,hi) if auc>=0.5 else (1-auc,1-hi,1-lo)
