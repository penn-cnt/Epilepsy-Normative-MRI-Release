"""
Create a temporary combined_profiles with real HCP-YA ages for ComBat.
Real ages are used only during processing and NOT saved to the final outputs.
The temporary file is deleted after ComBat completes.

Expects hcpya_subjects.csv to have columns: hcp_id, placeholder_age, real_age
"""
import pandas as pd
import os

base = os.environ["PROJECT_DIR"] + "/hcp_tractometry"
ages_path = os.environ["HCPYA_SUBJECTS_CSV"]
profiles_path = f'{base}/combined_profiles.csv'
temp_path = f'{base}/combined_profiles_realages_tmp.csv'

# Load real ages
ages = pd.read_csv(ages_path)
print(f"Loaded {len(ages)} HCP-YA subjects from hcpya_subjects.csv")

if 'real_age' not in ages.columns:
    raise ValueError("hcpya_subjects.csv must have a 'real_age' column")

missing = ages['real_age'].isna().sum()
if missing > 0:
    print(f"WARNING: {missing} subjects have no real_age value")

# Build mapping: sub-{hcp_id} -> real_age
age_map = dict(zip('sub-' + ages['hcp_id'].astype(str), ages['real_age']))

# Load combined profiles
print(f"Loading {profiles_path}...")
df = pd.read_csv(profiles_path, low_memory=False)
print(f"Loaded {len(df)} rows, {df['subject'].nunique()} subjects")

# Update ages for HCP-YA subjects in memory
updated = 0
for sub, real_age in age_map.items():
    mask = df['subject'] == sub
    if mask.any() and pd.notna(real_age):
        df.loc[mask, 'age'] = str(real_age)
        updated += 1

print(f"Updated ages for {updated}/{len(age_map)} HCP-YA subjects")

# Verify
hcpya_ages = df[df['cohort'] == 'HCP-YA'].groupby('subject')['age'].first()
print(f"\nHCP-YA age distribution after update:")
print(f"  Min: {pd.to_numeric(hcpya_ages, errors='coerce').min()}")
print(f"  Max: {pd.to_numeric(hcpya_ages, errors='coerce').max()}")
print(f"  Mean: {pd.to_numeric(hcpya_ages, errors='coerce').mean():.1f}")

# Save to TEMPORARY file (not overwriting the original)
df.to_csv(temp_path, index=False)
size = os.path.getsize(temp_path) / (1024**2)
print(f"\nSaved temporary file to {temp_path} ({size:.1f} MB)")
print("NOTE: Original combined_profiles.csv is unchanged.")
print("The temporary file will be deleted after the pipeline completes.")
print("Done!")
