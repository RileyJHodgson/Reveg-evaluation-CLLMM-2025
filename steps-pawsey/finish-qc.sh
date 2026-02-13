#!/bin/bash --login
#SBATCH --account=pawsey1018
#SBATCH --partition=work
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=800M  
#SBATCH --time=00:10:00

set -euo pipefail

sample=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $SAMP_MAN)
[[ -z "$sample" ]] && { echo "No sample for index $SLURM_ARRAY_TASK_ID"; exit 1; }

out_dir="${OUT_DIR:-./results}/${sample}"

mkdir -p "$out_dir"

echo "$(date -Is) finish.sh ran for sample=${sample} SLURM_JOB_ID=${SLURM_JOB_ID:-NA}"  > "${out_dir}/FINISH_RAN.txt"

# storage directories
sample=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $SAMP_MAN)
[[ -z "$sample" ]] && { echo "No sample for index $SLURM_ARRAY_TASK_ID"; exit 1; }

out_root="${OUT_DIR:-$submit_dir/results}"
out_root="$(realpath "$out_root")"

qc_dir="${out_root}/${sample}/qc"
qc_dir="$(realpath "$qc_dir")"

# Remove spent files
rm -r $qc_dir

echo "Cleaned qc directory with fastqs for sample $sample"
