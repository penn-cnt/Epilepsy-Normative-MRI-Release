#!/bin/bash
#SBATCH --cpus-per-task=8
#SBATCH --mem=64GB
#SBATCH --job-name=qsiprep

# Site configuration comes from the repo-root .env (see .env.example); run `source env.sh` before sbatch.
: "${PROJECT_DIR:?PROJECT_DIR is not set: source env.sh at the repository root first}"

SES=${1}

# Subject list passed via environment variable
read -ra SUBS <<< "${SUBJECTS}"
SUBJ=${SUBS[${SLURM_ARRAY_TASK_ID}]}

if [ -z "$SUBJ" ]; then
    echo "ERROR: No subject for array task ID ${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

echo "Array task ${SLURM_ARRAY_TASK_ID}: Processing ${SUBJ}"

# Specify temp directory for apptainer
export APPTAINER_TMPDIR=${BIDS_ROOT}/software/apptainer_tmp

# Define key directories
BASE=${BIDS_ROOT}
DATA=${BASE}/outputs/BIDS_Penn
QSIPREP_DIR=${BASE}/outputs/derivatives/ses-${SES}/qsiprep
QSIPREP_IMG=${BASE}/code/qsiprep/qsiprep-1.0.0.sif
FSLIC=${BASE}/jobs/qsiprep/license.txt
WORK_DIR=${BASE}/jobs/qsiprep/work; mkdir -p $WORK_DIR

# Skip if already processed
if [ -d "${QSIPREP_DIR}/${SUBJ}/ses-${SES}" ]; then
    echo "Subject ${SUBJ} already processed, skipping"
    exit 0
fi

# Skip if no DWI
if [ ! -d "${DATA}/${SUBJ}/ses-${SES}/dwi" ]; then
    echo "Subject ${SUBJ} does not have a dwi directory, skipping"
    exit 0
fi

cmd="apptainer run \
-B $BASE:$BASE \
-B $FSLIC:$BASE/jobs/qsiprep/license.txt \
$QSIPREP_IMG \
${DATA} \
${QSIPREP_DIR} \
participant \
--participant-label $SUBJ \
--session-id $SES \
--work-dir $WORK_DIR \
--skip-bids-validation \
--output-resolution 1.25 \
--unringing-method rpg \
--denoise-method dwidenoise \
--fs-license-file $FSLIC"

echo $cmd
$cmd

# Change the permission of the output subject directory to cntgroup
chgrp -R cntgroup ${QSIPREP_DIR}/${SUBJ}
