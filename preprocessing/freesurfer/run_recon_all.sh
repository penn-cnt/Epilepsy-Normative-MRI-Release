#!/bin/bash

# Site configuration comes from the repo-root .env (see .env.example); run `source env.sh` before sbatch.
: "${PROJECT_DIR:?PROJECT_DIR is not set: source env.sh at the repository root first}"
# Submit FreeSurfer 7.3.2 recon-all for the full cohort as a throttled SLURM array.
#
# Reads one subject ID per line from outputs/freesurfer/subjects.txt.
# recon_all.sh skips any subject that already has stats/, so resubmitting is safe.
#
# Usage: bash run_recon_all.sh [--dry-run]

code_dir=${PROJECT_DIR}/code/freesurfer
fs_dir=${PROJECT_DIR}/outputs/freesurfer
logs_dir=${fs_dir}/logs
sub_list=${fs_dir}/subjects.txt

concurrency=16          # 16 jobs x 8 cpus = 128 cpus = full compute node
wave_size=200           # the cluster QOS caps submitted jobs at 250
dry_run=0
[ "${1}" = "--dry-run" ] && dry_run=1

n_subs=$(grep -c '^sub-RID' "${sub_list}")
echo "Subject list: ${sub_list} (${n_subs} subjects)"

if [ "$dry_run" -eq 1 ]; then
    echo "DRY RUN - nothing will be created or submitted"
    echo "Would create:  ${fs_dir} and ${logs_dir}"
    echo "Would submit:  waves of ${wave_size} tasks (QOS cap 250), max ${concurrency} concurrent, chained"
    exit 0
fi

mkdir -p "${logs_dir}"

# fsaverage is needed by recon-all's surface registration steps
if [ ! -e "${fs_dir}/fsaverage" ]; then
    ln -s /usr/local/freesurfer/subjects/fsaverage "${fs_dir}/fsaverage"
fi

# The cluster partition enforces a QOS with MaxSubmitJobsPerUser=250,
# so a single 355-task array is rejected. Split into waves under the cap and
# chain each wave to start after the previous one finishes.
echo ""
echo "Submitting ${n_subs} tasks in waves of ${wave_size}, max ${concurrency} concurrent"

prev_jobid=""
start=1
while [ ${start} -le ${n_subs} ]; do
    end=$(( start + wave_size - 1 ))
    [ ${end} -gt ${n_subs} ] && end=${n_subs}

    dep=""
    [ -n "${prev_jobid}" ] && dep="--dependency=afterany:${prev_jobid}"

    jobid=$(sbatch --parsable ${dep} \
        --array=${start}-${end}%${concurrency} \
        --output=${logs_dir}/%A_%a_recon-all.o \
        --error=${logs_dir}/%A_%a_recon-all.e \
        ${code_dir}/recon_all_array.sh)

    if [ -z "${jobid}" ]; then
        echo "ERROR: submission failed for tasks ${start}-${end}"
        exit 1
    fi

    if [ -n "${prev_jobid}" ]; then
        echo "  tasks ${start}-${end}: job ${jobid} (starts after ${prev_jobid})"
    else
        echo "  tasks ${start}-${end}: job ${jobid}"
    fi

    prev_jobid="${jobid}"
    start=$(( end + 1 ))
done

echo ""
echo "Monitor with:  squeue -u "${SLURM_USER}""
echo "Logs in:       ${logs_dir}"
