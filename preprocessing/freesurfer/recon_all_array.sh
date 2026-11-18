#!/bin/bash
#SBATCH --job-name=recon_all
#SBATCH --time=16:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=8

# Site configuration comes from the repo-root .env (see .env.example); run `source env.sh` before sbatch.
: "${PROJECT_DIR:?PROJECT_DIR is not set: source env.sh at the repository root first}"

# SLURM array worker: resolves this task's subject from the subject list,
# then defers to recon_all.sh.
#
# Submitted by run_recon_all.sh -- do not sbatch directly.

code_dir=${PROJECT_DIR}/code/freesurfer
sub_list=${PROJECT_DIR}/outputs/freesurfer/subjects.txt

sub=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${sub_list}")

if [ -z "$sub" ]; then
    echo "ERROR: no subject at line ${SLURM_ARRAY_TASK_ID} of ${sub_list}"
    exit 1
fi

echo "Array task ${SLURM_ARRAY_TASK_ID} -> ${sub}"
exec bash ${code_dir}/recon_all.sh "${sub}"
