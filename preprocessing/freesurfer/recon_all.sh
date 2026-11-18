#!/bin/bash
#SBATCH --job-name=recon_all
#SBATCH --time=16:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=8

# Site configuration comes from the repo-root .env (see .env.example); run `source env.sh` before sbatch.
: "${PROJECT_DIR:?PROJECT_DIR is not set: source env.sh at the repository root first}"

sub=${1}
ses=clinical01

# Directories
base_dir=${BIDS_ROOT}
bids_dir=${base_dir}/outputs/BIDS_Penn
fs_dir=${PROJECT_DIR}/outputs/freesurfer

# Standalone FreeSurfer 7.3.2 container.
fs_img=${FS_CONTAINER}

# FreeSurfer license.
fs_license=${FS_LICENSE}

# Input T1w
t1w=${bids_dir}/${sub}/ses-${ses}/anat/${sub}_ses-${ses}_acq-mprage_T1w.nii.gz

if [ ! -f "$t1w" ]; then
    echo "ERROR: T1w not found for ${sub}: ${t1w}"
    exit 1
fi

if [ ! -f "$fs_license" ]; then
    echo "ERROR: FreeSurfer license not found: ${fs_license}"
    exit 1
fi

if [ ! -f "$fs_img" ]; then
    echo "ERROR: container not found: ${fs_img}"
    exit 1
fi

if [ -d "${fs_dir}/${sub}/stats" ]; then
    echo "Subject ${sub} already has FreeSurfer stats, skipping"
    exit 0
fi

echo "Running recon-all for ${sub}"
echo "Input: ${t1w}"
echo "Output: ${fs_dir}/${sub}"

export APPTAINER_BINDPATH=/mnt:/mnt

apptainer exec --cleanenv --containall \
    -B /mnt:/mnt \
    --env FS_LICENSE=${fs_license} \
    ${fs_img} \
    /usr/local/freesurfer/bin/recon-all \
    -all \
    -i ${t1w} \
    -subjid ${sub} \
    -sd ${fs_dir} \
    -openmp 8

exit_code=$?
if [ $exit_code -eq 0 ]; then
    echo "SUCCESS: recon-all completed for ${sub}"
else
    echo "FAILED: recon-all exited with code ${exit_code} for ${sub}"
fi
