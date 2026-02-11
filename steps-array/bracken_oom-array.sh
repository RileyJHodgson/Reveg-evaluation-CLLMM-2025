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
[[ -z "$sample" ]] && { echo "No sample for index $SLURM_ARRAY_TASK_ID"; exit 1; }


# storage directories
submit_dir="${SLURM_SUBMIT_DIR:-$PWD}"

out_root="${OUT_DIR:-$submit_dir/results}"
mkdir -p "$out_root"
out_root="$(realpath "$out_root")"

tax_dir="${out_root}/${sample}/tax"
mkdir -p "$qc_dir" "$tax_dir"

# database
db_root="${TAX_DB:-$submit_dir/resources}"
db_root="$(realpath "$db_root")"
database_oom="${tax_root}/k2_Oomycota_20250113"
database_oom="$(realpath "$database_oom")"

# Temporary directory on compute node
tmpdir="${SLURM_TMPDIR:-/scratch/user/${USER}/${SLURM_JOB_ID}}"
mkdir -p "$tmpdir"
trap 'rm -rf "$tmpdir"' EXIT
echo "Using tmpdir: $tmpdir"

# Transfer files
cp "${tax_dir}/${sample}.k2_oom_report" "$tmpdir/"

# tmpdir
echo "Using tmpdir: $tmpdir"
cd "$tmpdir"

# braken
bracken -r 100 -l S -t 16 \
    -d "${database_oom}" \
    -i "${sample}.k2_oom_output" \
    -o "${sample}.bracken_oom_output" \
    -w "${sample}.bracken_oom_report"

# Copy useful results back to storage
# cp "${sample}.bracken_oom_output" "${tax_dir}/"
cp "${sample}.bracken_oom_report" "${tax_dir}/"

echo "bracken (oom) complete for sample $sample"

# Clean up tmpdir
rm -rf "$tmpdir"
# rm "${tax_dir}/${sample}.k2_oom_output" # these delete the k2 files
# rm "${tax_dir}/${sample}.k2_oom_report" # these delete the k2 files

# End
