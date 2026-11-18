#!/bin/bash
#SBATCH --job-name=zscores_all
#SBATCH --time=2:00:00
#SBATCH --mem=48G
#SBATCH --cpus-per-task=4
#SBATCH --output=${PROJECT_DIR}/hcp_tractometry/code/zscores_all_%j.o
#SBATCH --error=${PROJECT_DIR}/hcp_tractometry/code/zscores_all_%j.e

# Site configuration comes from the repo-root .env (see .env.example); run `source env.sh` before sbatch.
: "${PROJECT_DIR:?PROJECT_DIR is not set: source env.sh at the repository root first}"

source ~/.bashrc
conda activate "${CONDA_ENV}"

python3 << 'PYEOF'
import os
import pandas as pd
import numpy as np
from numpy.linalg import lstsq

base = os.environ["PROJECT_DIR"] + '/hcp_tractometry'
input_path = f'{base}/combined_profiles_harmonized.csv'

print("Loading harmonized data...")
df = pd.read_csv(input_path, low_memory=False)
print(f"Loaded {df['subject'].nunique()} subjects")

# Swap in real HCP-YA ages in memory for accurate regression (not saved to output)
import os
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
meta_cols = ['cohort', 'site', 'age', 'sex', 'group_1', 'group_2', 'combat_reference']
subject_meta = df.groupby('subject')[meta_cols + ['age_numeric', 'sex_numeric']].first().reset_index()

bilateral_pairs = {
    'AF': ('AF_L', 'AF_R'), 'CST': ('CST_L', 'CST_R'),
    'C': ('C_L', 'C_R'), 'F': ('F_L', 'F_R'),
    'IFOF': ('IFOF_L', 'IFOF_R'), 'UF': ('UF_L', 'UF_R'),
}
midline_tracts = ['MCP', 'CCMid', 'CC_ForcepsMajor', 'CC_ForcepsMinor']

def compute_zscores(df_input, ref_mask_col='combat_reference', ref_val='yes'):
    """Compute age/sex-adjusted z-scores for each tract-scalar combo."""
    ref = df_input[df_input[ref_mask_col] == ref_val]
    tract_col = 'tract' if 'tract' in df_input.columns else 'tractID'
    tracts = sorted(df_input[tract_col].unique())

    zscore_rows = []
    for tract in tracts:
        for scalar in scalars:
            ref_t = ref[ref[tract_col] == tract].dropna(subset=[scalar, 'age_numeric', 'sex_numeric'])
            if len(ref_t) < 10:
                continue
            X_ref = np.column_stack([ref_t['age_numeric'].values, ref_t['sex_numeric'].values, np.ones(len(ref_t))])
            y_ref = ref_t[scalar].values
            coeffs, _, _, _ = lstsq(X_ref, y_ref, rcond=None)
            residual_sd = np.std(y_ref - X_ref @ coeffs)

            for _, row in df_input[df_input[tract_col] == tract].iterrows():
                if pd.isna(row[scalar]) or pd.isna(row['age_numeric']):
                    z = np.nan
                else:
                    x = np.array([row['age_numeric'], row['sex_numeric'], 1.0])
                    z = (row[scalar] - x @ coeffs) / residual_sd

                out = {
                    'subject': row['subject'],
                    'tract': row[tract_col],
                    'scalar': scalar,
                    'raw_value': row[scalar],
                    'z_score': z,
                    'cohort': row['cohort'],
                    'site': row['site'],
                    'age': row['age'],
                    'sex': row['sex'],
                    'group_1': row['group_1'],
                    'group_2': row['group_2'],
                }
                if 'nodeID' in row.index:
                    out['nodeID'] = row['nodeID']
                zscore_rows.append(out)

    return pd.DataFrame(zscore_rows)


# =====================================================
# 1. BILATERAL (10 tracts) - existing behavior
# =====================================================
print("\n=== Computing 10-tract bilateral z-scores ===")
tract_means = df.groupby(['subject', 'tractID'])[scalars].mean().reset_index()
tract_means = tract_means.merge(subject_meta, on='subject')

rows = []
for sub in tract_means['subject'].unique():
    sub_data = tract_means[tract_means['subject'] == sub]
    sm = subject_meta[subject_meta['subject'] == sub].iloc[0]
    for tract_name, (left, right) in bilateral_pairs.items():
        row = {'subject': sub, 'tract': tract_name}
        for c in meta_cols + ['age_numeric', 'sex_numeric']:
            row[c] = sm[c]
        for scalar in scalars:
            vals = []
            ld = sub_data[sub_data['tractID'] == left]
            rd = sub_data[sub_data['tractID'] == right]
            if len(ld) > 0: vals.append(ld[scalar].values[0])
            if len(rd) > 0: vals.append(rd[scalar].values[0])
            row[scalar] = np.mean(vals) if vals else np.nan
        rows.append(row)
    for tract_name in midline_tracts:
        td = sub_data[sub_data['tractID'] == tract_name]
        row = {'subject': sub, 'tract': tract_name}
        for c in meta_cols + ['age_numeric', 'sex_numeric']:
            row[c] = sm[c]
        for scalar in scalars:
            row[scalar] = td[scalar].values[0] if len(td) > 0 else np.nan
        rows.append(row)

df_bilateral = pd.DataFrame(rows)
df_z_bilateral = compute_zscores(df_bilateral)
df_z_bilateral.loc[df_z_bilateral['subject'].map(lambda s: s in placeholder_ages), 'age'] = \
    df_z_bilateral['subject'].map(lambda s: placeholder_ages.get(s, s))
out1 = f'{base}/tract_zscores.csv'
df_z_bilateral.to_csv(out1, index=False)
print(f"Saved {df_z_bilateral.shape} to {out1}")


# =====================================================
# 2. LEFT/RIGHT (16 tracts) - no averaging
# =====================================================
print("\n=== Computing 16-tract L/R z-scores ===")
tract_means_16 = df.groupby(['subject', 'tractID'])[scalars].mean().reset_index()
tract_means_16 = tract_means_16.merge(subject_meta, on='subject')
tract_means_16.rename(columns={'tractID': 'tract'}, inplace=True)

df_z_16 = compute_zscores(tract_means_16)
df_z_16.loc[df_z_16['subject'].map(lambda s: s in placeholder_ages), 'age'] = \
    df_z_16['subject'].map(lambda s: placeholder_ages.get(s, s))
out2 = f'{base}/tract_zscores_16tracts.csv'
df_z_16.to_csv(out2, index=False)
print(f"Saved {df_z_16.shape} to {out2}")


# =====================================================
# 3. PER-NODE (16 tracts x 100 nodes)
# =====================================================
print("\n=== Computing per-node z-scores ===")
node_data = df.copy()
node_data.rename(columns={'tractID': 'tract'}, inplace=True)

# z-score per tract per node
all_tracts = sorted(node_data['tract'].unique())
all_nodes = sorted(node_data['nodeID'].unique())
ref_node = node_data[node_data['combat_reference'] == 'yes']

zscore_rows = []
total = len(all_tracts) * len(all_nodes) * len(scalars)
done = 0
for tract in all_tracts:
    for node in all_nodes:
        for scalar in scalars:
            done += 1
            ref_tn = ref_node[(ref_node['tract'] == tract) & (ref_node['nodeID'] == node)].dropna(subset=[scalar, 'age_numeric', 'sex_numeric'])
            if len(ref_tn) < 10:
                continue
            X_ref = np.column_stack([ref_tn['age_numeric'].values, ref_tn['sex_numeric'].values, np.ones(len(ref_tn))])
            y_ref = ref_tn[scalar].values
            coeffs, _, _, _ = lstsq(X_ref, y_ref, rcond=None)
            residual_sd = np.std(y_ref - X_ref @ coeffs)

            subset = node_data[(node_data['tract'] == tract) & (node_data['nodeID'] == node)]
            for _, row in subset.iterrows():
                if pd.isna(row[scalar]) or pd.isna(row['age_numeric']):
                    z = np.nan
                else:
                    x = np.array([row['age_numeric'], row['sex_numeric'], 1.0])
                    z = (row[scalar] - x @ coeffs) / residual_sd
                zscore_rows.append({
                    'subject': row['subject'],
                    'tract': tract,
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

df_z_node = pd.DataFrame(zscore_rows)
df_z_node.loc[df_z_node['subject'].map(lambda s: s in placeholder_ages), 'age'] = \
    df_z_node['subject'].map(lambda s: placeholder_ages.get(s, s))
out3 = f'{base}/node_zscores.csv'
df_z_node.to_csv(out3, index=False)
print(f"Saved {df_z_node.shape} to {out3}")


# =====================================================
# Sanity checks
# =====================================================
ref_subs = subject_meta[subject_meta['combat_reference'] == 'yes']['subject'].unique()
for name, df_z in [('bilateral', df_z_bilateral), ('16-tract', df_z_16), ('per-node', df_z_node)]:
    ref_z = df_z[df_z['subject'].isin(ref_subs)]
    print(f"\n{name} reference check: mean={ref_z['z_score'].mean():.3f}, sd={ref_z['z_score'].std():.3f}")

print("\nDone!")
PYEOF
