#!/bin/bash

# Site configuration comes from the repo-root .env (see .env.example); run `source env.sh` before sbatch.
: "${PROJECT_DIR:?PROJECT_DIR is not set: source env.sh at the repository root first}"
# Pipeline to re-run harmonization with real HCP-YA ages
# Usage: bash run_pipeline_with_real_ages.sh
#
# Prerequisites:
#   Marc updates hcpya_subjects.csv with a 3rd column 'real_age'
#   (columns: hcp_id, placeholder_age, real_age)

BASE=${PROJECT_DIR}/hcp_tractometry
CODE=${BASE}/code
AGES_CSV=${HCPYA_SUBJECTS_CSV}

# Check that real ages exist
if ! head -1 "${AGES_CSV}" | grep -q "real_age"; then
    echo "ERROR: hcpya_subjects.csv does not have a 'real_age' column."
    echo "Marc needs to add a 'real_age' column with the actual HCP-YA ages."
    echo "Expected format: hcp_id,placeholder_age,real_age"
    exit 1
fi

echo "=== Step 1: Update HCP-YA ages in combined_profiles.csv ==="
JOB1=$(sbatch --parsable \
    --job-name=update_ages \
    --time=0:30:00 \
    --mem=16G \
    --cpus-per-task=2 \
    --output=${CODE}/update_ages_%j.o \
    --error=${CODE}/update_ages_%j.e \
    --wrap="source ~/.bashrc && conda activate "${CONDA_ENV}" && python3 ${CODE}/update_ages.py")
echo "Submitted update_ages job: ${JOB1}"

echo "=== Step 2: Per-tract ComBat harmonization (runs after Step 1) ==="
echo "  (bundle-wise approach per Chandio et al. 2024)"
JOB2=$(sbatch --parsable \
    --dependency=afterok:${JOB1} \
    ${CODE}/run_combat_pertract.sh)
echo "Submitted combat_pertract job: ${JOB2} (depends on ${JOB1})"

echo "=== Step 3: Per-node z-scores (runs after Step 2) ==="
JOB3=$(sbatch --parsable \
    --dependency=afterok:${JOB2} \
    ${CODE}/run_zscores_pernode.sh)
echo "Submitted zscores_pernode job: ${JOB3} (depends on ${JOB2})"

echo ""
echo "Pipeline submitted! Monitor with: squeue -u \$(whoami)"
echo "  Step 1 (update ages):        ${JOB1}"
echo "  Step 2 (ComBat per-tract):   ${JOB2}"
echo "  Step 3 (z-scores per-node):  ${JOB3}"
