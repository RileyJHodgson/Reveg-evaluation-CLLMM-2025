#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --time=4-0
#SBATCH --mem=64000M
#SBATCH --cpus-per-task=16
#SBATCH --nodes=1

set -euo pipefail

set +u
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate superfocusonly
set -u

#source "/home/user/hodg0248/miniconda3/etc/profile.d/conda.sh"

sample=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $SAMP_MAN)
[[ -z "$sample" ]] && { echo "No sample for index $SLURM_ARRAY_TASK_ID"; exit 1; }

# storage directories
submit_dir="${SLURM_SUBMIT_DIR:-$PWD}"

out_root="${OUT_DIR:-$submit_dir/results}"
mkdir -p "$out_root"
out_root="$(realpath "$out_root")"

qc_dir="${out_root}/${sample}/qc"
fun_dir="${out_root}/${sample}/func"
mkdir -p "$fun_dir"

# specify superfocus path. Will always be hard coded within a specific conda environment
database="/home/hodg0248/miniconda3/envs/superfocusonly/lib/python3.10/site-packages/superfocus_app"

# Temporary directory on compute node
tmpdir="${SLURM_TMPDIR:-/scratch/user/${USER}/${SLURM_JOB_ID}}"
mkdir -p "$tmpdir"
trap 'rm -rf "$tmpdir"' EXIT
echo "Using tmpdir: $tmpdir"

# Transfer files
cp "${qc_dir}/${sample}_R1.good.fastq.gz" "$tmpdir/"

# change directory to tmpdir
cd "$tmpdir"
gunzip "${sample}_R1.good.fastq.gz"

superfocus --query "${sample}_R1.good.fastq" \
    --output_directory "super_out_${sample}" \
    --alternate_directory "${database}" \
    --output_prefix "norm_${sample}_" \
    --aligner diamond --database DB_100 --normalise_output 1

# Copy relevant file back to storage
cp "super_out_${sample}/norm_${sample}_all_levels_and_function.xls" "${fun_dir}/"
echo "superfocus complete for sample: $sample"

rm -r "super_out_${sample}"

# End
