#!/bin/bash
#SBATCH --job-name=combat_tract
#SBATCH --time=4:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --output=${PROJECT_DIR}/hcp_tractometry/code/combat_tract_%j.o
#SBATCH --error=${PROJECT_DIR}/hcp_tractometry/code/combat_tract_%j.e

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
temp_path = f'{base}/combined_profiles_realages_tmp.csv'
input_path = temp_path if os.path.exists(temp_path) else f'{base}/combined_profiles.csv'
output_path = f'{base}/combined_profiles_harmonized_pertract.csv'

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
meta['sex_numeric'] = (meta['sex'] == 'F').astype(int)
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

# Run ComBat per tract (each tract harmonized independently)
# Following Chandio et al. (2024) bundle-wise approach:
# each tract's 100 nodes are the features for one ComBat call
print("\nRunning ComBat per tract (bundle-wise)...")
harmonized_values = {}

for scalar in scalars:
    print(f"\nProcessing {scalar}...")
    for tract in tracts:
        # Build feature matrix: 100 nodes x N subjects for this tract
        features = []
        feature_names = []

        for node in nodes:
            tract_node_data = df[(df['tractID'] == tract) & (df['nodeID'] == node)][['subject', scalar]]
            pivot = tract_node_data.set_index('subject')[scalar].reindex(subject_order)
            features.append(pivot.values)
            feature_names.append((tract, node, scalar))

        data = np.array(features)  # 100 nodes x N subjects

        # Fill NaN with feature mean
        for i in range(data.shape[0]):
            nan_cols = np.isnan(data[i])
            if nan_cols.any():
                data[i, nan_cols] = np.nanmean(data[i])

        # Run ComBat on 100 features (nodes) for this tract
        combat_result = neuroCombat(
            dat=data,
            covars=covars,
            batch_col='batch',
        )

        harmonized = combat_result['data']
        for idx, key in enumerate(feature_names):
            harmonized_values[key] = harmonized[idx]

        print(f"  {tract}: {data.shape[0]} nodes x {data.shape[1]} subjects")

print(f"\nTotal harmonized features: {len(harmonized_values)}")

# Rebuild long format
print("Rebuilding long format...")
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
                key = (tract, node, scalar)
                if key in harmonized_values:
                    row[scalar] = harmonized_values[key][i]
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
    # Don't delete temp file - per-node combat may still need it
    print(f"  Temp file preserved for other combat runs: {temp_path}")

df_out.to_csv(output_path, index=False)
size = os.path.getsize(output_path) / (1024**2)
print(f"\nSaved to {output_path} ({size:.1f} MB)")
print(f"Shape: {df_out.shape}")
print("Done!")
PYEOF
