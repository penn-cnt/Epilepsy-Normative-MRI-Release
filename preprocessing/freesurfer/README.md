# FreeSurfer pipeline

- FreeSurfer 7.3.2 standalone Apptainer image (`freesurfer-7.3.2.sif`), build `freesurfer-linux-centos8_x86_64-7.3.2-20220804-6354275`.
- Input: raw BIDS `sub-*/ses-clinical01/anat/sub-*_ses-clinical01_acq-mprage_T1w.nii.gz`. T1w only, no FLAIR or T2 pial refinement.
- Command (`recon_all.sh`): `recon-all -all -i <T1w> -subjid <sub> -sd <SUBJECTS_DIR> -openmp 8`.
- Median runtime about 2 hours on 8 CPUs.

| Script | Role |
|---|---|
| `recon_all.sh` | Per-subject recon-all (the command above) |
| `recon_all_array.sh` | SLURM array worker that picks one subject from the subject list |
| `run_recon_all.sh` | Submits the cohort as SLURM arrays, split into waves under the queue cap |

Site-specific paths are read from `.env` (see the top-level README). The per-subject `scripts/` provenance files and `stats/` tables are part of the data release, not this repository.
