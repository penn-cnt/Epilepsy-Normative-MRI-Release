#!/bin/bash
#SBATCH --job-name=combat_wb
#SBATCH --time=4:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --output=${PROJECT_DIR}/hcp_tractometry/code/combat_wb_%j.o
#SBATCH --error=${PROJECT_DIR}/hcp_tractometry/code/combat_wb_%j.e

# Site configuration comes from the repo-root .env (see .env.example); run `source env.sh` before sbatch.
: "${PROJECT_DIR:?PROJECT_DIR is not set: source env.sh at the repository root first}"

source ~/.bashrc
conda activate "${CONDA_ENV}"

python3 << 'PYEOF'
import os
import pandas as pd
import numpy as np
from neuroCombat import neuroCombat
import os

base = os.environ["PROJECT_DIR"] + '/hcp_tractometry'
# Use temp file with real ages if it exists (from update_ages.py), otherwise use original
import os
temp_path = f'{base}/combined_profiles_realages_tmp.csv'
input_path = temp_path if os.path.exists(temp_path) else f'{base}/combined_profiles.csv'
output_path = f'{base}/combined_profiles_harmonized.csv'

# Load combined profiles
print("Loading data...")
df = pd.read_csv(input_path, low_memory=False)
print(f"Loaded {df['subject'].nunique()} subjects, {len(df)} rows")

tracts = sorted(df['tractID'].unique())
scalars = ['dti_fa', 'dti_md', 'dti_rd', 'dti_ad']
nodes = sorted(df['nodeID'].unique())
print(f"Tracts: {len(tracts)}, Nodes: {len(nodes)}, Scalars: {len(scalars)}")

# Build subject-level metadata
meta = df.groupby('subject').first()[['cohort', 'site', 'age', 'sex', 'group_1', 'group_2', 'combat_reference']].reset_index()
subject_order = meta['subject'].values

# Convert age to numeric
meta['age_numeric'] = pd.to_numeric(meta['age'].replace('90 or older', '90'), errors='coerce')
# Convert sex to numeric
meta['sex_numeric'] = (meta['sex'] == 'F').astype(int)
# Batch encoding
batch = pd.Categorical(meta['site']).codes + 1

print(f"\nSite encoding:")
for site in sorted(meta['site'].unique()):
    n = (meta['site'] == site).sum()
    n_ref = ((meta['site'] == site) & (meta['combat_reference'] == 'yes')).sum()
    print(f"  {site}: {n} total, {n_ref} reference")

covars = pd.DataFrame({
    'batch': batch,
    'age': meta['age_numeric'].values,
    'sex': meta['sex_numeric'].values
})

ref_mask = (meta['combat_reference'] == 'yes').values
print(f"\nReference subjects: {ref_mask.sum()}")
print(f"Non-reference subjects: {(~ref_mask).sum()}")

# Pivot and run ComBat per scalar
print("\nPivoting to wide format...")
all_harmonized = {}

for scalar in scalars:
    print(f"\nProcessing {scalar}...")
    features = []
    feature_names = []

    for tract in tracts:
        tract_data = df[df['tractID'] == tract][['subject', 'nodeID', scalar]]
        tract_pivot = tract_data.pivot(index='subject', columns='nodeID', values=scalar)
        tract_pivot = tract_pivot.reindex(subject_order)

        for node in nodes:
            if node in tract_pivot.columns:
                features.append(tract_pivot[node].values)
                feature_names.append(f"{tract}_node{node}_{scalar}")

    data = np.array(features)  # features x subjects

    # Fill NaN with feature mean for ComBat
    nan_features = np.isnan(data).any(axis=1).sum()
    if nan_features > 0:
        print(f"  {nan_features} features have NaN - filling with feature mean")
        for i in range(data.shape[0]):
            nan_cols = np.isnan(data[i])
            if nan_cols.any():
                data[i, nan_cols] = np.nanmean(data[i])

    print(f"  Running ComBat on {data.shape[0]} features x {data.shape[1]} subjects...")
    combat_result = neuroCombat(
        dat=data,
        covars=covars,
        batch_col='batch',
    )

    harmonized = combat_result['data']
    for idx, fname in enumerate(feature_names):
        all_harmonized[fname] = harmonized[idx]

print(f"\nRebuilding long format...")
rows_out = []
for i, sub in enumerate(subject_order):
    sub_meta = meta[meta['subject'] == sub].iloc[0]
    for tract in tracts:
        for node in nodes:
            row = {
                'subject': sub,
                'cohort': sub_meta['cohort'],
                'site': sub_meta['site'],
                'age': sub_meta['age'],
                'sex': sub_meta['sex'],
                'group_1': sub_meta['group_1'],
                'group_2': sub_meta['group_2'],
                'combat_reference': sub_meta['combat_reference'],
                'tractID': tract,
                'nodeID': node,
            }
            for scalar in scalars:
                fname = f"{tract}_node{node}_{scalar}"
                if fname in all_harmonized:
                    row[scalar] = all_harmonized[fname][i]
                else:
                    row[scalar] = np.nan
            rows_out.append(row)

df_out = pd.DataFrame(rows_out)

# If we used real ages, restore placeholder ages in the output so restricted data isn't saved
if os.path.exists(temp_path):
    print("\nRestoring placeholder ages in output (restricted ages not saved)...")
    orig = pd.read_csv(f'{base}/combined_profiles.csv', low_memory=False, usecols=['subject', 'cohort', 'age'])
    orig_ages = orig[orig['cohort'] == 'HCP-YA'].groupby('subject')['age'].first().to_dict()
    hcpya_mask = df_out['cohort'] == 'HCP-YA'
    df_out.loc[hcpya_mask, 'age'] = df_out.loc[hcpya_mask, 'subject'].map(orig_ages)
    print(f"  Restored placeholder ages for {hcpya_mask.sum()} HCP-YA rows")
    # Clean up temp file
    os.remove(temp_path)
    print(f"  Deleted temporary file: {temp_path}")

df_out.to_csv(output_path, index=False)
size = os.path.getsize(output_path) / (1024**2)
print(f"\nSaved to {output_path} ({size:.1f} MB)")
print(f"Shape: {df_out.shape}")
print("Done!")
PYEOF
