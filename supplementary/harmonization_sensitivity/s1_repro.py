import os
import pandas as pd, numpy as np
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests
RNG_SEED=42
base=os.environ["PROJECT_DIR"] + "/hcp_tractometry"
cb_dir=f'{base}/centilebrain_zscores'
meta_path=os.environ["BIDS_ROOT"] + "/metadata/metadata_drug_responsiveness.csv"

# ---- cell 3: build morph cohort (defines the subject set) ----
ct=pd.concat([pd.read_csv(f'{cb_dir}/zscore_CorticalThickness_male.csv'),
              pd.read_csv(f'{cb_dir}/zscore_CorticalThickness_female.csv')],ignore_index=True)
sv=pd.concat([pd.read_csv(f'{cb_dir}/zscore_SubcorticalVolume_male.csv'),
              pd.read_csv(f'{cb_dir}/zscore_SubcorticalVolume_female.csv')],ignore_index=True)
ct=ct.rename(columns=lambda c:c.replace('_entorhil_','_entorhinal_').replace('_supramargil_','_supramarginal_'))
both=set(ct['SubjectID'])&set(sv['SubjectID'])
ct=ct[ct['SubjectID'].isin(both)].copy(); sv=sv[sv['SubjectID'].isin(both)].copy()
sv_z=[c for c in sv.columns if c not in ['SITE','SubjectID','Vendor','FreeSurfer_Version','age']]
morph=ct.merge(sv[['SubjectID']+sv_z],on='SubjectID',how='inner')
meta=pd.read_csv(meta_path).rename(columns={'record_id':'SubjectID'})
morph=morph.merge(meta[['SubjectID','group_1','group_2']],on='SubjectID',how='left')
morph=morph.drop_duplicates(subset='SubjectID',keep='first')
tract_subs=set(pd.read_csv(f'{base}/node_zscores.csv',usecols=['subject'])['subject'].unique())
morph=morph[morph['SubjectID'].isin(tract_subs)].copy()
print(f"morph cohort: {len(morph)}")
np.save('morph_subjects.npy',morph['SubjectID'].values)
