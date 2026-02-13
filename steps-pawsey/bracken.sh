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
conda activate /scratch/pawsey1018/rhodgson/envs/kraken-suite
set -u

sample=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $SAMP_MAN)
#[[ -z "$sample" ]] && { echo "No sample for index $SLURM_ARRAY_TASK_ID"; exit 1; }

# storage directories
submit_dir="${SLURM_SUBMIT_DIR:-$PWD}"

out_root="${OUT_DIR:-$submit_dir/results}"
mkdir -p "$out_root"
out_root="$(realpath "$out_root")"

tax_dir="${out_root}/${sample}/tax"
mkdir -p "$tax_dir"

# database
database="${TAX_DB:-$submit_dir/resources/k2_standard_202510155}"
database="$(realpath "$database")"

# prep temporary directory
base_tmp="${SLURM_TMPDIR:-/scratch/pawsey1018/rhodgson/tmp}"
bucket=$((SLURM_JOB_ID % 100))
tmpdir="$base_tmp/bucket_${bucket}/${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
mkdir -p "$tmpdir"
trap 'rm -rf "$tmpdir"; echo "Cleaned tmpdir: $tmpdir"' EXIT
cd "$tmpdir"

# Transfer files
cp "${tax_dir}/${sample}.k2_report" "$tmpdir/"

# brcaken
bracken -r 100 -l S -t 16 \
    -d "${database}" \
    -i "${sample}.k2_report" \
    -o "${sample}.bracken_output" \
    -w "${sample}.bracken_report"

# Copy results back to storage
# cp "${sample}.bracken_output" "${tax_dir}/"
cp "${sample}.bracken_report" "${tax_dir}/"

echo "bracken complete for sample $sample"

rm "${tax_dir}/${sample}.k2_report" # these delete the k2 files
echo "k2 reports for sample $sample removed"

# End
