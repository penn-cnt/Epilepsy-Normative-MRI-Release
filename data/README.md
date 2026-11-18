# Input data

`analysis/composite_figures.ipynb`, `analysis/metadata_summary.ipynb` and `analysis/abnormality_burden_summary.ipynb` read these files relative to `analysis/`. The name in parentheses in the first column is the file's path in the Pennsieve dataset.

| File | What |
|---|---|
| `metadata_drug_responsiveness.csv` (dataset `participants/metadata_drug_responsiveness.csv`) | One row per subject: sex, ages, `group_1` (focal / generalized / mixed / unclassified / pnee / syncope / one_sz_only), `group_2` (drug-resistant / drug-responsive / ...), lesion flags (`mri_lesion`, the twelve lesion-type columns from `lesion_mts` to `lesion_other`, `lesion_left`, `lesion_right`, `lesion_bilateral`, used by the burden notebook). |
| `centilebrain_zscores/zscore_CorticalThickness_{male,female}.csv` (dataset `morphometry/centilebrain_zscores/`) | CentileBrain cortical-thickness z-scores, 68 Desikan-Killiany regions, 355 subjects. |
| `centilebrain_zscores/zscore_SubcorticalVolume_{male,female}.csv` (dataset `morphometry/centilebrain_zscores/`) | CentileBrain subcortical-volume z-scores, 14 structures. |
| `node_zscores.csv` (dataset `tractometry/node_zscores_clinical.csv`) | Node-level tractometry z-scores, long format (`subject, cohort, tract, scalar, nodeID, z_score, ...`). The notebooks use only `cohort == "Clinical"` and `scalar in {dti_fa, dti_md}`. |
