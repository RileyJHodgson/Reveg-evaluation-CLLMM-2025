#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --time=4-0
#SBATCH --mem=64000M
#SBATCH --cpus-per-task=16
#SBATCH --nodes=1

set -euo pipefail

set +u
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate k2DB_env
set -u

sample=$(sed -n "${SLURM_ARRAY_TASK_ID}" $SAMP_MAN)
[[ -z "$sample" ]] && { echo "No sample for index $SLURM_ARRAY_TASK_ID"; exit 1; }

# storage directories
submit_dir="${SLURM_SUBMIT_DIR:-$PWD}"

out_root="${OUT_DIR:-$submit_dir/results}"
mkdir -p "$out_root"
out_root="$(realpath "$out_root")"

qc_dir="${out_root}/${sample}/qc"
tax_dir="${out_root}/${sample}/tax"
mkdir -p "$qc_dir" "$tax_dir"

# database paths. k2_standard_20251015 is default but can be set in workflow .sh file
database="${TAX_DB:-$submit_dir/resources/k2_standard_20251015}"
database="$(realpath "$database")"

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
kraken2 --db ${database} \
    --threads 16 --memory-mapping --report-zero-counts --use-names \
    --paired "${sample}_R1.good.fastq.gz" "${sample}_R2.good.fastq.gz" \
    --output "${sample}.k2_output" \
    --report "${sample}.k2_report"

# Copy results back to storage
# cp "${sample}.k2_output" "${tax_dir}/"
cp "${sample}.k2_report" "${tax_dir}/"

echo "Kraken2 complete for sample $sample"

# End
