#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --time=00:02:00
#SBATCH --mem=100M
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1

set -euo pipefail

sample=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $SAMP_MAN)
[[ -z "$sample" ]] && { echo "No sample for index $SLURM_ARRAY_TASK_ID"; exit 1; }

out_root="${OUT_DIR:-$submit_dir/results}"
out_root="$(realpath "$out_root")"

fun_dir="${out_root}/${sample}/func"
fun_dir="$(realpath "$fun_dir")"

# Remove spent files
rm -r $fun_dir

echo "Cleaned directories containing functions files for sample $sample"
