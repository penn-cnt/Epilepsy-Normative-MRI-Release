#!/bin/bash
#SBATCH --job-name=pyafq_hcpya
#SBATCH --time=12:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4

# Site configuration comes from the repo-root .env (see .env.example); run `source env.sh` before sbatch.
: "${PROJECT_DIR:?PROJECT_DIR is not set: source env.sh at the repository root first}"

# Per-subject pyAFQ processing for HCP-YA
# Usage: sbatch pyafq_hcpya.sh <sub-ID>
# Example: sbatch pyafq_hcpya.sh sub-100206

sub=$1
hcp_id=${sub#sub-}

hcp_raw_dir=${HCPYA_RAW_DIR}
clone_dir=${PROJECT_DIR}/hcp_tractometry/hcpya_datalad
bids_dir=${PROJECT_DIR}/hcp_tractometry/bids_hcpya
output_dir=${PROJECT_DIR}/hcp_tractometry/pyafq_hcpya

source ~/.bashrc
conda activate "${CONDA_ENV}"

# Git safe.directory for accessing the shared HCP datalad repositories
export GIT_CONFIG_COUNT=1
export GIT_CONFIG_KEY_0=safe.directory
export GIT_CONFIG_VALUE_0="*"

# Step 1: Clone datalad repo if not already cloned
if [ ! -d "${clone_dir}/${hcp_id}" ]; then
    echo "Cloning datalad repo for ${sub}..."
    datalad clone ${hcp_raw_dir}/${hcp_id} ${clone_dir}/${hcp_id}
    if [ $? -ne 0 ]; then
        echo "ERROR: datalad clone failed for ${sub}"
        exit 1
    fi
fi

# Step 2: Fetch diffusion data
echo "Fetching data for ${sub} via datalad..."
cd ${clone_dir}/${hcp_id}/T1w
datalad get Diffusion/data.nii.gz Diffusion/bvals Diffusion/bvecs
if [ $? -ne 0 ]; then
    echo "ERROR: datalad get failed for ${sub}"
    exit 1
fi
echo "Data fetched successfully"

# Step 3: Update BIDS symlinks to point to cloned repo
echo "Updating BIDS symlinks..."
dwi_dir=${bids_dir}/derivatives/hcp/${sub}/dwi
mkdir -p ${dwi_dir}
clone_diffusion=${clone_dir}/${hcp_id}/T1w/Diffusion
ln -sf ${clone_diffusion}/data.nii.gz ${dwi_dir}/${sub}_dwi.nii.gz
ln -sf ${clone_diffusion}/bvals ${dwi_dir}/${sub}_dwi.bval
ln -sf ${clone_diffusion}/bvecs ${dwi_dir}/${sub}_dwi.bvec

# Step 4: Run pyAFQ
echo "Running pyAFQ for ${sub}..."
python - <<PYTHON
from AFQ.api.group import GroupAFQ
import AFQ.data.fetch as afd

myafq = GroupAFQ(
    bids_path='${bids_dir}',
    preproc_pipeline='hcp',
    output_dir='${output_dir}',
    participant_labels=['${hcp_id}'],
    bundle_info=afd.read_hcp_atlas(n_bundles=16),
    scalars=["dti_fa", "dti_md", "dti_rd", "dti_ad"],
    bids_validate=False,
    b0_threshold=50,
    min_bval=500,
    max_bval=1500,
)

myafq.export("profiles")
myafq.export("bundles")
PYTHON

pyafq_exit=$?

# Step 5: Drop the data to free space
echo "Dropping data for ${sub} via datalad..."
cd ${clone_dir}/${hcp_id}/T1w
datalad drop Diffusion/data.nii.gz
echo "Data dropped"

if [ ${pyafq_exit} -ne 0 ]; then
    echo "ERROR: pyAFQ failed for ${sub}"
    exit 1
fi

echo "pyAFQ complete for ${sub}"
