"""
Generate CentileBrain input CSVs from FreeSurfer stats files.
Outputs:
  - cortical_thickness.csv (subject, age, sex, LThickness, RThickness, 68 ROI columns)
  - subcortical_volumes.csv (subject, age, sex, ICV, 14 subcortical structure columns)
"""
import os
import pandas as pd
import re
from pathlib import Path

# paths
base = Path(os.environ["PROJECT_DIR"]) / "outputs"
output_dir = Path(os.environ["PROJECT_DIR"]) / "outputs"
meta_path = Path(os.environ["BIDS_ROOT"]) / "metadata" / "metadata_drug_responsiveness.csv"

freesurfer_dir = base / 'freesurfer'
subjects = sorted([d.name for d in freesurfer_dir.glob('sub-*') if d.is_dir()])
print(f'Found {len(subjects)} FreeSurfer subjects')

# load metadata for age/sex (deduplicate by subject)
meta = pd.read_csv(meta_path)
meta = meta.rename(columns={'record_id': 'subject'})
meta = meta.drop_duplicates(subset='subject', keep='first')
meta_lookup = meta.set_index('subject')[['age_MRI', 'sex']].to_dict('index')

# subcortical structure mapping (aseg name -> CentileBrain column)
subcort_map = {
    'Left-Thalamus': 'Lthal', 'Right-Thalamus': 'Rthal',
    'Left-Caudate': 'Lcaud', 'Right-Caudate': 'Rcaud',
    'Left-Putamen': 'Lput', 'Right-Putamen': 'Rput',
    'Left-Pallidum': 'Lpal', 'Right-Pallidum': 'Rpal',
    'Left-Hippocampus': 'Lhippo', 'Right-Hippocampus': 'Rhippo',
    'Left-Amygdala': 'Lamyg', 'Right-Amygdala': 'Ramyg',
    'Left-Accumbens-area': 'Laccumb', 'Right-Accumbens-area': 'Raccumb',
}

ct_rows = []
sv_rows = []
skipped = []

for sub in subjects:
    stats_dir = freesurfer_dir / sub / 'stats'
    lh_file = stats_dir / 'lh.aparc.stats'
    rh_file = stats_dir / 'rh.aparc.stats'
    aseg_file = stats_dir / 'aseg.stats'

    if not (lh_file.exists() and rh_file.exists() and aseg_file.exists()):
        skipped.append(sub)
        continue

    info = meta_lookup.get(sub, {})
    age = info.get('age_MRI', '')
    sex = info.get('sex', '')
    if pd.isna(age): age = ''
    if pd.isna(sex): sex = ''

    # --- Cortical thickness ---
    ct_row = {'subject': sub, 'age': age, 'sex': sex}

    lh_text = lh_file.read_text()
    m = re.search(r'# Measure Cortex, MeanThickness, .+?, ([\d.]+),', lh_text)
    ct_row['LThickness'] = float(m.group(1)) if m else ''

    df_lh = pd.read_csv(lh_file, comment='#', sep=r'\s+', header=None,
                        names=['StructName','NumVert','SurfArea','GrayVol',
                               'ThickAvg','ThickStd','MeanCurv','GausCurv','FoldInd','CurvInd'])
    for _, r in df_lh.iterrows():
        ct_row[f"L_{r['StructName']}_thickavg"] = r['ThickAvg']

    rh_text = rh_file.read_text()
    m = re.search(r'# Measure Cortex, MeanThickness, .+?, ([\d.]+),', rh_text)
    ct_row['RThickness'] = float(m.group(1)) if m else ''

    df_rh = pd.read_csv(rh_file, comment='#', sep=r'\s+', header=None,
                        names=['StructName','NumVert','SurfArea','GrayVol',
                               'ThickAvg','ThickStd','MeanCurv','GausCurv','FoldInd','CurvInd'])
    for _, r in df_rh.iterrows():
        ct_row[f"R_{r['StructName']}_thickavg"] = r['ThickAvg']

    ct_rows.append(ct_row)

    # --- Subcortical volumes ---
    sv_row = {'subject': sub, 'age': age, 'sex': sex}

    aseg_text = aseg_file.read_text()
    m = re.search(r'EstimatedTotalIntraCranialVol.+?, ([\d.]+),', aseg_text)
    sv_row['ICV'] = float(m.group(1)) if m else ''

    df_aseg = pd.read_csv(aseg_file, comment='#', sep=r'\s+', header=None,
                          names=['Index','SegId','NVoxels','Volume_mm3',
                                 'StructName','normMean','normStdDev','normMin','normMax','normRange'])
    for _, r in df_aseg.iterrows():
        if r['StructName'] in subcort_map:
            sv_row[subcort_map[r['StructName']]] = r['Volume_mm3']

    sv_rows.append(sv_row)

# --- Build and save cortical thickness CSV ---
ct_df = pd.DataFrame(ct_rows)
ct_col_order = ['subject', 'age', 'sex', 'LThickness', 'RThickness'] + \
    [c for c in ct_df.columns if c.endswith('_thickavg') and c.startswith('L_')] + \
    [c for c in ct_df.columns if c.endswith('_thickavg') and c.startswith('R_')]
ct_df = ct_df[[c for c in ct_col_order if c in ct_df.columns]]
ct_df.to_csv(output_dir / 'cortical_thickness.csv', index=False)
print(f'cortical_thickness.csv: {ct_df.shape[0]} subjects, {ct_df.shape[1]} columns')

# --- Build and save subcortical volumes CSV ---
sv_col_order = ['subject', 'age', 'sex', 'ICV',
                'Lthal', 'Rthal', 'Lcaud', 'Rcaud', 'Lput', 'Rput',
                'Lpal', 'Rpal', 'Lhippo', 'Rhippo', 'Lamyg', 'Ramyg', 'Laccumb', 'Raccumb']
sv_df = pd.DataFrame(sv_rows)[sv_col_order]
sv_df.to_csv(output_dir / 'subcortical_volumes.csv', index=False)
print(f'subcortical_volumes.csv: {sv_df.shape[0]} subjects, {sv_df.shape[1]} columns')

# summary
n_with_age = ct_df[ct_df['age'] != ''].shape[0]
print(f'\n{n_with_age} subjects have age/sex metadata (ready for CentileBrain)')
print(f'{ct_df.shape[0] - n_with_age} subjects missing age/sex')
if skipped:
    print(f'\nSkipped {len(skipped)} subjects (missing stats files):')
    for s in skipped:
        print(f'  {s}')
