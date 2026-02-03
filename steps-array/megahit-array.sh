#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --time=4-0
#SBATCH --mem=64000M
#SBATCH --cpus-per-task=16
#SBATCH --nodes=1

set -euo pipefail

set +u
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate mega-contig
set -u

sample=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $SAMP_MAN)
[[ -z "$sample" ]] && { echo "No sample for index $SLURM_ARRAY_TASK_ID"; exit 1; }

# storage directories
submit_dir="${SLURM_SUBMIT_DIR:-$PWD}"

out_root="${OUT_DIR:-$submit_dir/results}"
mkdir -p "$out_root"
out_root="$(realpath "$out_root")"

qc_dir="${out_root}/${sample}/qc"
contig_dir="${out_root}/${sample}/contigs"
mkdir -p "$qc_dir" "$contig_dir"

# Temporary directory on compute node
tmpdir="${SLURM_TMPDIR:-/scratch/user/${USER}/${SLURM_JOB_ID}}"
mkdir -p "$tmpdir"
trap 'rm -rf "$tmpdir"' EXIT
echo "Using tmpdir: $tmpdir"

# Transfer files
cp "${qc_dir}/${sample}_R1.good.fastq.gz" "$tmpdir/"
cp "${qc_dir}/${sample}_R2.good.fastq.gz" "$tmpdir/"
cp "${qc_dir}/${sample}_R1.single.fastq.gz" "$tmpdir/"
cp "${qc_dir}/${sample}_R2.single.fastq.gz" "$tmpdir/"

# tmpdir
echo "Using tmpdir: $tmpdir"
cd "$tmpdir"

# megahit
megahit \
  --presets meta-large \
  -1 "${sample}_R1.good.fastq.gz" \
  -2 "${sample}_R2.good.fastq.gz" \
  -r "${sample}_R1.single.fastq.gz","${sample}_R2.single.fastq.gz" \
  -o "${sample}_contigs"

# Copy results back to storage
cp "${sample}_contigs/final.contigs.fa" "$contig_dir/"

# Clean up tmpdir
rm -rf $tmpdir

echo "contigs assembly complete for sample $sample"

# End
