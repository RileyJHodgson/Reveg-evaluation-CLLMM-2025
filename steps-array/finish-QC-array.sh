#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --time=00:02:00
#SBATCH --mem=100M
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1

set -euo pipefail

sample=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $SAMP_MAN)
[[ -z "$sample" ]] && { echo "No sample for index $SLURM_ARRAY_TASK_ID"; exit 1; }

out_dir="${OUT_DIR:-./results}/${sample}"

mkdir -p "$out_dir"

echo "$(date -Is) finish.sh ran for sample=${sample} SLURM_JOB_ID=${SLURM_JOB_ID:-NA}"  > "${out_dir}/FINISH_RAN.txt"


# storage directories
sample=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $SAMP_MAN)
[[ -z "$sample" ]] && { echo "No sample for index $SLURM_ARRAY_TASK_ID"; exit 1; }

#out_root="${OUT_DIR:-$submit_dir/results}"
#out_root="$(realpath "$out_root")"

#qc_dir="${out_root}/${sample}/qc"
#qc_dir="$(realpath "$qc_dir")"

# Remove spent files
#rm ${qc_dir}/${sample}_R1.good.fastq.gz # these delete the cleaned reads
#rm ${qc_dir}/${sample}_R2.good.fastq.gz # these delete the raw reads
#rm ${qc_dir}/${sample}_R1.single.fastq.gz # these delete the cleaned reads
#rm ${qc_dir}/${sample}_R2.single.fastq.gz # these delete the raw reads

#echo "Cleaned qc fastqs for sample $sample removed"

