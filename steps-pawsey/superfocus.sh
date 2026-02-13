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
conda activate /scratch/pawsey1018/rhodgson/envs/superfocusonly
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
database="/scratch/pawsey1018/rhodgson/envs/superfocusonly/lib/python3.11/site-packages/superfocus_app/"

# prep temporary directory
base_tmp="${SLURM_TMPDIR:-/scratch/pawsey1018/rhodgson/tmp}"
bucket=$((SLURM_JOB_ID % 100))
tmpdir="$base_tmp/bucket_${bucket}/${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
mkdir -p "$tmpdir"
trap 'rm -rf "$tmpdir"; echo "Cleaned tmpdir: $tmpdir"' EXIT
cd "$tmpdir"

# Transfer files
cp "${qc_dir}/${sample}_R1.good.fastq.gz" "$tmpdir/"

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
