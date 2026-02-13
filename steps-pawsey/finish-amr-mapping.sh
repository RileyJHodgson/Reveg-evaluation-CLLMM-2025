#!/bin/bash --login
#SBATCH --account=pawsey1018
#SBATCH --partition=work
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=800M  
#SBATCH --time=00:10:00

set -euo pipefail

sample=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $SAMP_MAN)
[[ -z "$sample" ]] && { echo "No sample for index $SLURM_ARRAY_TASK_ID"; exit 1; }

out_root="${OUT_DIR:-$submit_dir/results}"
out_root="$(realpath "$out_root")"

amr_dir="${out_root}/${sample}/amr"
amr_dir="$(realpath "$amr_dir")"

bti_dir="${out_root}/${sample}/indexes"
bti_dir="$(realpath "$bti_dir")"

btm_dir="${out_root}/${sample}/mapped_reads"
btm_dir="$(realpath "$btm_dir")"

spf_dir="${out_root}/${sample}/singlem"
spf_dir="$(realpath "$spf_dir")"

# Remove spent files
rm -r $amr_dir
rm -r $bti_dir
rm -r $btm_dir
rm -r $spf_dir

echo "Cleaned directories containing amr, bowtie2 indexes, bowtie2 mapped and singlem files for sample $sample"
