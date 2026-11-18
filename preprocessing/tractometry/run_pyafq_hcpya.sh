#!/bin/bash

# Site configuration comes from the repo-root .env (see .env.example); run `source env.sh` before sbatch.
: "${PROJECT_DIR:?PROJECT_DIR is not set: source env.sh at the repository root first}"
# run_pyafq_hcpya.sh
# Batch submission script for HCP-YA pyAFQ processing
# Usage: bash run_pyafq_hcpya.sh

n_subs=1064  # Full run

bids_dir=${PROJECT_DIR}/hcp_tractometry/bids_hcpya
pyafq_dir=${PROJECT_DIR}/hcp_tractometry/pyafq_hcpya
code_dir=${PROJECT_DIR}/hcp_tractometry/code
logs_dir=${PROJECT_DIR}/hcp_tractometry/logs

mkdir -p ${pyafq_dir}
mkdir -p ${logs_dir}

sub_counter=0

for sub_dir in ${bids_dir}/derivatives/hcp/sub-*; do
    sub=$(basename ${sub_dir})
    echo ${sub}

    # Skip if pyafq already processed
    if [ -d "${pyafq_dir}/${sub}" ]; then
        echo "  Already processed, skipping"
        continue
    fi

    echo "  Submitting pyAFQ job"

    logs_output=${logs_dir}/${sub}_pyafq_job-%j.o
    logs_error=${logs_dir}/${sub}_pyafq_job-%j.e
    sbatch --output=${logs_output} --error=${logs_error} ${code_dir}/pyafq_hcpya.sh ${sub}

    sub_counter=$((sub_counter + 1))

    if [[ $sub_counter -ge $n_subs ]]; then
        echo "Reached limit of ${n_subs} subjects"
        break
    fi
done

echo "Submitted ${sub_counter} jobs"
