#!/bin/bash --login
#SBATCH --account=pawsey1018
#SBATCH --partition=work
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=4000M  
#SBATCH --time=12:00:00

set -euo pipefail

set +u
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate /scratch/pawsey1018/rhodgson/envs/amrfinder
set -u

sample=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $SAMP_MAN)
[[ -z "$sample" ]] && { echo "No sample for index $SLURM_ARRAY_TASK_ID"; exit 1; }

# storage directories
submit_dir="${SLURM_SUBMIT_DIR:-$PWD}"

out_root="${OUT_DIR:-$submit_dir/results}"
mkdir -p "$out_root"
out_root="$(realpath "$out_root")"

contig_dir="${out_root}/${sample}/contigs"
amr_dir="${out_root}/${sample}/amr"
mkdir -p "$amr_dir"

# prep temporary directory
base_tmp="${SLURM_TMPDIR:-/scratch/pawsey1018/rhodgson/tmp}"
bucket=$((SLURM_JOB_ID % 100))
tmpdir="$base_tmp/bucket_${bucket}/${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
mkdir -p "$tmpdir"
trap 'rm -rf "$tmpdir"; echo "Cleaned tmpdir: $tmpdir"' EXIT
cd "$tmpdir"

# Transfer files
cp "$contig_dir/final.contigs.fa" "$tmpdir"

# AMRfinderplus
amrfinder --plus --threads 16 \
    --name "${sample}" \
    --nucleotide "final.contigs.fa" \
    --output "amrfp_${sample}.tsv"

# Copy results back to storage
cp "amrfp_${sample}.tsv" "$amr_dir/"

echo "amrfinder complete for sample $sample"

# End
