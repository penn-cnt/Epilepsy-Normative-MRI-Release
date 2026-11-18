#!/bin/bash
#SBATCH --job-name=zscores_wb
#SBATCH --time=1:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --output=${PROJECT_DIR}/hcp_tractometry/code/zscores_wb_%j.o
#SBATCH --error=${PROJECT_DIR}/hcp_tractometry/code/zscores_wb_%j.e

# Site configuration comes from the repo-root .env (see .env.example); run `source env.sh` before sbatch.
: "${PROJECT_DIR:?PROJECT_DIR is not set: source env.sh at the repository root first}"

source ~/.bashrc
conda activate "${CONDA_ENV}"

python3 << 'PYEOF'
import os
import pandas as pd
import numpy as np
from numpy.linalg import lstsq

input_path = os.environ["PROJECT_DIR"] + '/hcp_tractometry/combined_profiles_harmonized.csv'
output_path = os.environ["PROJECT_DIR"] + '/hcp_tractometry/tract_zscores.csv'

print("Loading harmonized data...")
df = pd.read_csv(input_path, low_memory=False)
print(f"Loaded {df['subject'].nunique()} subjects")

# Convert age to numeric
df['age_numeric'] = pd.to_numeric(df['age'].replace('90 or older', '90'), errors='coerce')
df['sex_numeric'] = (df['sex'] == 'F').astype(int)

scalars = ['dti_fa', 'dti_md', 'dti_rd', 'dti_ad']

# Define bilateral tract pairs and midline tracts
bilateral_pairs = {
    'AF': ('AF_L', 'AF_R'),
    'CST': ('CST_L', 'CST_R'),
    'C': ('C_L', 'C_R'),
    'F': ('F_L', 'F_R'),
    'IFOF': ('IFOF_L', 'IFOF_R'),
    'UF': ('UF_L', 'UF_R'),
}
midline_tracts = ['MCP', 'CCMid', 'CC_ForcepsMajor', 'CC_ForcepsMinor']

# Step 1: Compute mean across nodes for each tract per subject
print("Computing mean tract values...")
tract_means = df.groupby(['subject', 'tractID'])[scalars].mean().reset_index()

# Add metadata back
meta_cols = ['subject', 'cohort', 'site', 'age', 'age_numeric', 'sex', 'sex_numeric', 'group_1', 'group_2', 'combat_reference']
subject_meta = df.groupby('subject')[meta_cols[1:]].first().reset_index()
tract_means = tract_means.merge(subject_meta, on='subject')

# Step 2: Average L+R for bilateral tracts
print("Averaging bilateral tracts...")
rows = []
for sub in tract_means['subject'].unique():
    sub_data = tract_means[tract_means['subject'] == sub]
    sub_meta = subject_meta[subject_meta['subject'] == sub].iloc[0]

    # Bilateral: average L and R
    for tract_name, (left, right) in bilateral_pairs.items():
        left_data = sub_data[sub_data['tractID'] == left]
        right_data = sub_data[sub_data['tractID'] == right]

        row = {
            'subject': sub,
            'tract': tract_name,
            'cohort': sub_meta['cohort'],
            'site': sub_meta['site'],
            'age': sub_meta['age'],
            'age_numeric': sub_meta['age_numeric'],
            'sex': sub_meta['sex'],
            'sex_numeric': sub_meta['sex_numeric'],
            'group_1': sub_meta['group_1'],
            'group_2': sub_meta['group_2'],
            'combat_reference': sub_meta['combat_reference'],
        }

        for scalar in scalars:
            vals = []
            if len(left_data) > 0:
                vals.append(left_data[scalar].values[0])
            if len(right_data) > 0:
                vals.append(right_data[scalar].values[0])
            row[scalar] = np.mean(vals) if vals else np.nan

        rows.append(row)

    # Midline: use as-is
    for tract_name in midline_tracts:
        t_data = sub_data[sub_data['tractID'] == tract_name]
        row = {
            'subject': sub,
            'tract': tract_name,
            'cohort': sub_meta['cohort'],
            'site': sub_meta['site'],
            'age': sub_meta['age'],
            'age_numeric': sub_meta['age_numeric'],
            'sex': sub_meta['sex'],
            'sex_numeric': sub_meta['sex_numeric'],
            'group_1': sub_meta['group_1'],
            'group_2': sub_meta['group_2'],
            'combat_reference': sub_meta['combat_reference'],
        }
        for scalar in scalars:
            if len(t_data) > 0:
                row[scalar] = t_data[scalar].values[0]
            else:
                row[scalar] = np.nan
        rows.append(row)

df_tracts = pd.DataFrame(rows)
print(f"Tract-level data: {df_tracts['subject'].nunique()} subjects x {df_tracts['tract'].nunique()} tracts")

# Step 3: Age/sex-adjusted z-scores
print("Computing age/sex-adjusted z-scores...")
ref = df_tracts[df_tracts['combat_reference'] == 'yes']
tracts = sorted(df_tracts['tract'].unique())

zscore_rows = []
for tract in tracts:
    for scalar in scalars:
        # Fit age+sex model on reference subjects using numpy lstsq
        ref_tract = ref[ref['tract'] == tract].dropna(subset=[scalar, 'age_numeric', 'sex_numeric'])
        X_ref = np.column_stack([ref_tract['age_numeric'].values, ref_tract['sex_numeric'].values, np.ones(len(ref_tract))])
        y_ref = ref_tract[scalar].values

        coeffs, _, _, _ = lstsq(X_ref, y_ref, rcond=None)

        # Get residual SD from reference
        y_pred_ref = X_ref @ coeffs
        residual_sd = np.std(y_ref - y_pred_ref)

        # Z-score all subjects
        for _, row in df_tracts[df_tracts['tract'] == tract].iterrows():
            if pd.isna(row[scalar]) or pd.isna(row['age_numeric']):
                z = np.nan
            else:
                x = np.array([row['age_numeric'], row['sex_numeric'], 1.0])
                predicted = x @ coeffs
                z = (row[scalar] - predicted) / residual_sd

            zscore_rows.append({
                'subject': row['subject'],
                'tract': tract,
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

df_z = pd.DataFrame(zscore_rows)
df_z.to_csv(output_path, index=False)

print(f"\nSaved to {output_path}")
print(f"Shape: {df_z.shape}")

# Summary stats for epilepsy subjects
epi = df_z[df_z['cohort'] == 'Clinical']
epi_not_sync = epi[epi['group_1'] != 'syncope']
print(f"\nEpilepsy z-score summary (n={epi_not_sync['subject'].nunique()}):")
print(f"  Mean z: {epi_not_sync['z_score'].mean():.3f}")
print(f"  SD z: {epi_not_sync['z_score'].std():.3f}")
print(f"  Range: [{epi_not_sync['z_score'].min():.2f}, {epi_not_sync['z_score'].max():.2f}]")

# Sanity check: reference subjects should have mean ~0, SD ~1
ref_z = df_z[df_z['subject'].isin(ref['subject'].unique())]
print(f"\nReference z-score sanity check (n={ref_z['subject'].nunique()}):")
print(f"  Mean z: {ref_z['z_score'].mean():.3f}")
print(f"  SD z: {ref_z['z_score'].std():.3f}")

print("\nDone!")
PYEOF
