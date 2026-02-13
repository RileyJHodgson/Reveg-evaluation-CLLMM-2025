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

contig_dir="${out_root}/${sample}/contigs"
contig_dir="$(realpath "$contig_dir")"

# Remove spent files
rm -r $contig_dir
echo "Cleaned directory containing contigs files for sample $sample"

