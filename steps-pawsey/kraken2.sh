#!/bin/bash --login
#SBATCH --account=pawsey1018
#SBATCH --partition=work
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=4000M  
#SBATCH --time=24:00:00

set -euo pipefail

set +u
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate /scratch/pawsey1018/rhodgson/envs/k2DB_env
set -u

sample=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $SAMP_MAN)
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

# prep temporary directory
base_tmp="${SLURM_TMPDIR:-/scratch/pawsey1018/rhodgson/tmp}"
bucket=$((SLURM_JOB_ID % 100))
tmpdir="$base_tmp/bucket_${bucket}/${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
mkdir -p "$tmpdir"
trap 'rm -rf "$tmpdir"; echo "Cleaned tmpdir: $tmpdir"' EXIT
cd "$tmpdir"

# Transfer files
cp "${qc_dir}/${sample}_R1.good.fastq.gz" "$tmpdir/"
cp "${qc_dir}/${sample}_R2.good.fastq.gz" "$tmpdir/"

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
