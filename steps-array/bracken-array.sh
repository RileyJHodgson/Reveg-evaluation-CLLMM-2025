#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --time=3-0
#SBATCH --mem=64000M
#SBATCH --cpus-per-task=16
#SBATCH --nodes=1

set -euo pipefail

set +u
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate kraken-suite
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

# Temporary directory on compute node
tmpdir="${SLURM_TMPDIR:-/scratch/user/${USER}/${SLURM_JOB_ID}}"
mkdir -p "$tmpdir"
trap 'rm -rf "$tmpdir"' EXIT
echo "Using tmpdir: $tmpdir"

# Transfer files
cp "${tax_dir}/${sample}.k2_report" "$tmpdir/"

# tmpdir
echo "Using tmpdir: $tmpdir"
cd "$tmpdir"

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
