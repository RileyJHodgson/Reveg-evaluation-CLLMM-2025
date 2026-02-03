#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --time=4-0
#SBATCH --mem=64000
#SBATCH --cpus-per-task=16
#SBATCH --nodes=1

set -euo pipefail

set +u
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate hostremoval
set -u

sample=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $SAMP_MAN)
[[ -z "$sample" ]] && { echo "No sample for index $SLURM_ARRAY_TASK_ID"; exit 1; }

# storage directories
submit_dir="${SLURM_SUBMIT_DIR:-$PWD}"

out_root="${OUT_DIR:-$submit_dir/results}"
mkdir -p "$out_root"
out_root="$(realpath "$out_root")"

contig_dir="${out_root}/${sample}/contigs"
index_dir="${out_root}/${sample}/indexes"
mkdir -p "$contig_dir" "index_dir"

# Temporary directory on compute node
tmpdir="${SLURM_TMPDIR:-/scratch/user/${USER}/${SLURM_JOB_ID}}"
mkdir -p "$tmpdir"
trap 'rm -rf "$tmpdir"' EXIT

# Transfer files
cp "$contig_dir/final.contigs.fa" "$tmpdir"

# tmpdir
echo "Using tmpdir: $tmpdir"
cd "$tmpdir" || exit 1

mkdir -p bt2_index

# bowtie-build
bowtie2-build final.contigs.fa bt2_index/ref_${sample}

# Copy results back to storage
cp -r bt2_index $index_dir

# Clean up tmpdir
rm -rf $tmpdir

echo "bowtie2-build complete for sample $sample"

# End
