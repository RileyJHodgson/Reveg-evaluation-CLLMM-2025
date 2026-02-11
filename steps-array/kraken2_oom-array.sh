#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --time=1-0
#SBATCH --mem=8000M
#SBATCH --cpus-per-task=2
#SBATCH --nodes=1

set -euo pipefail

set +u
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate k2DB_env
set -u

sample=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $SAMP_MAN)
[[ -z "$sample" ]] && { echo "No sample for index $SLURM_ARRAY_TASK_ID"; exit 1; }

# storage directories
submit_dir="${SLURM_SUBMIT_DIR:-$PWD}"

out_root="${OUT_DIR:-$submit_dir/results}"
mkdir -p "$out_root"
out_root="$(realpath "$out_root")"

tax_dir="${out_root}/${sample}/tax"
mkdir -p "$tax_dir"

# database
db_root="$submit_dir/resources"
db_root="$(realpath "$db_root")"

# /scratch/user/hodg0248/goyder/resources/k2_Oomycota_20250113
database_oom="${db_root}/k2_Oomycota_20250113"
database_oom="$(realpath "$database_oom")"

# Temporary directory on compute node
tmpdir="${SLURM_TMPDIR:-/scratch/user/${USER}/${SLURM_JOB_ID}}"
mkdir -p "$tmpdir"
trap 'rm -rf "$tmpdir"' EXIT

echo "Using tmpdir: $tmpdir"

# Transfer files
cp "${qc_dir}/${sample}_R1.good.fastq.gz" "$tmpdir/"
cp "${qc_dir}/${sample}_R2.good.fastq.gz" "$tmpdir/"

# tmpdir
cd $tmpdir
echo "Using tmpdir: $tmpdir"

# kraken2
kraken2 --db ${database_oom} \
    --threads 16 --memory-mapping --report-zero-counts --use-names \
    --paired "${sample}_R1.good.fastq.gz" "${sample}_R2.good.fastq.gz" \
    --output "${sample}.k2_oom_output" \
    --report "${sample}.k2_oom_report"

# Copy results back to storage
# cp "${sample}.k2_oom_output" "${tax_dir}/"
cp "${sample}.k2_oom_report" "${tax_dir}/"

echo "Kraken2 (oom) complete for sample $sample"

# End
