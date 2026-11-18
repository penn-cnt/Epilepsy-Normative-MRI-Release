#!/bin/bash
#SBATCH --job-name=zscores_fnode
#SBATCH --time=2:00:00
#SBATCH --mem=48G
#SBATCH --cpus-per-task=4
#SBATCH --output=${PROJECT_DIR}/hcp_tractometry/code/zscores_fnode_%j.o
#SBATCH --error=${PROJECT_DIR}/hcp_tractometry/code/zscores_fnode_%j.e

# Site configuration comes from the repo-root .env (see .env.example); run `source env.sh` before sbatch.
: "${PROJECT_DIR:?PROJECT_DIR is not set: source env.sh at the repository root first}"

source ~/.bashrc
conda activate "${CONDA_ENV}"

python3 << 'PYEOF'
import os
import pandas as pd
import numpy as np
from numpy.linalg import lstsq
import os

base = os.environ["PROJECT_DIR"] + '/hcp_tractometry'
input_path = f'{base}/combined_profiles_harmonized_pernode.csv'

print("Loading harmonized data (per-node combat)...")
df = pd.read_csv(input_path, low_memory=False)
print(f"Loaded {df['subject'].nunique()} subjects")

# Swap in real HCP-YA ages in memory for accurate regression (not saved to output)
ages_path = os.environ["HCPYA_SUBJECTS_CSV"]
if os.path.exists(ages_path):
    ages_csv = pd.read_csv(ages_path)
    if 'real_age' in ages_csv.columns:
        age_map = dict(zip('sub-' + ages_csv['hcp_id'].astype(str), ages_csv['real_age']))
        hcpya_mask = df['cohort'] == 'HCP-YA'
        df.loc[hcpya_mask, 'age'] = df.loc[hcpya_mask, 'subject'].map(age_map).astype(str)
        print(f"Using real HCP-YA ages for regression ({hcpya_mask.sum()} rows updated in memory)")
    else:
        print("WARNING: real_age column not found, using placeholder ages")
else:
    print("WARNING: HCP-YA ages file not found, using placeholder ages")

# Save original placeholder ages for output
orig_ages = pd.read_csv(input_path, low_memory=False, usecols=['subject', 'cohort', 'age'])
placeholder_ages = orig_ages[orig_ages['cohort'] == 'HCP-YA'].groupby('subject')['age'].first().to_dict()

df['age_numeric'] = pd.to_numeric(df['age'].replace('90 or older', '90'), errors='coerce')
df['sex_numeric'] = (df['sex'] == 'F').astype(int)

scalars = ['dti_fa', 'dti_md', 'dti_rd', 'dti_ad']
all_tracts = sorted(df['tractID'].unique())
all_nodes = sorted(df['nodeID'].unique())
ref = df[df['combat_reference'] == 'yes']

print(f"\nTracts: {len(all_tracts)}, Nodes: {len(all_nodes)}, Scalars: {len(scalars)}")
print(f"Reference subjects: {ref['subject'].nunique()}")

# Compute per-node z-scores
print("\n=== Computing per-node z-scores (from per-node combat) ===")
zscore_rows = []
total = len(all_tracts) * len(all_nodes) * len(scalars)
done = 0

for tract in all_tracts:
    for node in all_nodes:
        for scalar in scalars:
            done += 1
            ref_tn = ref[(ref['tractID'] == tract) & (ref['nodeID'] == node)].dropna(subset=[scalar, 'age_numeric', 'sex_numeric'])
            if len(ref_tn) < 10:
                continue
            X_ref = np.column_stack([ref_tn['age_numeric'].values, ref_tn['sex_numeric'].values, np.ones(len(ref_tn))])
            y_ref = ref_tn[scalar].values
            coeffs, _, _, _ = lstsq(X_ref, y_ref, rcond=None)
            residual_sd = np.std(y_ref - X_ref @ coeffs)

            subset = df[(df['tractID'] == tract) & (df['nodeID'] == node)]
            for _, row in subset.iterrows():
                if pd.isna(row[scalar]) or pd.isna(row['age_numeric']):
                    z = np.nan
                else:
                    x = np.array([row['age_numeric'], row['sex_numeric'], 1.0])
                    z = (row[scalar] - x @ coeffs) / residual_sd
                zscore_rows.append({
                    'subject': row['subject'],
                    'tract': row['tractID'],
                    'nodeID': node,
                    'scalar': scalar,
                    'raw_value': row[scalar],
                    'z_score': z,
                    'cohort': row['cohort'],
                    'site': row['site'],
                    'age': row['age'],
                    'sex': row['sex'],
                    'group_1': row['group_1'],
                    'group_2': row['group_2'],
                })
    print(f"  {tract} done ({done}/{total} tract-node-scalar combos)")

df_z = pd.DataFrame(zscore_rows)

# Restore placeholder ages for HCP-YA subjects
df_z.loc[df_z['subject'].map(lambda s: s in placeholder_ages), 'age'] = \
    df_z['subject'].map(lambda s: placeholder_ages.get(s, s))

out_path = f'{base}/node_zscores_pernode.csv'
df_z.to_csv(out_path, index=False)
print(f"\nSaved {df_z.shape} to {out_path}")

# Sanity check
ref_subs = ref['subject'].unique()
ref_z = df_z[df_z['subject'].isin(ref_subs)]
print(f"\nReference z-score sanity check:")
print(f"  Mean z: {ref_z['z_score'].mean():.3f}")
print(f"  SD z: {ref_z['z_score'].std():.3f}")

print("\nDone!")
PYEOF
