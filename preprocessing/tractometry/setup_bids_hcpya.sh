#!/bin/bash

# Site configuration comes from the repo-root .env (see .env.example); run `source env.sh` before sbatch.
: "${PROJECT_DIR:?PROJECT_DIR is not set: source env.sh at the repository root first}"
# setup_bids_hcpya.sh
# Creates a BIDS-compatible directory structure for HCP-YA data
# by symlinking to the HCP1200 preprocessed diffusion files.
#
# pyAFQ requires BIDS-formatted input. HCP-YA data is in HCP's own format,
# so this script creates a BIDS wrapper with symlinks.
#
# Usage: bash setup_bids_hcpya.sh

hcp_dir=${HCPYA_RAW_DIR}
bids_dir=${PROJECT_DIR}/hcp_tractometry/bids_hcpya

sub_counter=0
skip_counter=0

for sub_dir in ${hcp_dir}/[0-9]*; do
    hcp_id=$(basename ${sub_dir})
    sub="sub-${hcp_id}"
    dwi_dir=${bids_dir}/${sub}/dwi

    # Skip if already set up
    if [ -d "${dwi_dir}" ]; then
        skip_counter=$((skip_counter + 1))
        continue
    fi

    # Check that diffusion files exist
    hcp_diffusion=${sub_dir}/T1w/Diffusion
    if [ ! -L "${hcp_diffusion}/data.nii.gz" ] && [ ! -f "${hcp_diffusion}/data.nii.gz" ]; then
        echo "WARNING: No diffusion data for ${hcp_id}, skipping"
        continue
    fi

    # Create BIDS subject directory
    mkdir -p ${dwi_dir}

    # Create symlinks with BIDS naming
    ln -s ${hcp_diffusion}/data.nii.gz ${dwi_dir}/${sub}_dwi.nii.gz
    ln -s ${hcp_diffusion}/bvals ${dwi_dir}/${sub}_dwi.bval
    ln -s ${hcp_diffusion}/bvecs ${dwi_dir}/${sub}_dwi.bvec

    sub_counter=$((sub_counter + 1))
done

echo "Setup complete: ${sub_counter} subjects created, ${skip_counter} already existed"
