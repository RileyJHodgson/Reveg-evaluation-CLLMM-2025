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

# prep temporary directory
base_tmp="${SLURM_TMPDIR:-/scratch/pawsey1018/rhodgson/tmp}"
bucket=$((SLURM_JOB_ID % 100))
tmpdir="$base_tmp/bucket_${bucket}/${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
mkdir -p "$tmpdir"
trap 'rm -rf "$tmpdir"; echo "Cleaned tmpdir: $tmpdir"' EXIT
cd "$tmpdir"

# Transfer files
cp "${tax_dir}/${sample}.k2_oom_report" "$tmpdir/"

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
# rm "${tax_dir}/${sample}.k2_oom_output" # these delete the k2 files
rm "${tax_dir}/${sample}.k2_oom_report" # these delete the k2 files

# End
