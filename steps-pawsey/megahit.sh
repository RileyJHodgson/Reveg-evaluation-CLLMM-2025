#!/bin/bash --login
#SBATCH --account=pawsey1018
#SBATCH --partition=long
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=4000M  
#SBATCH --time=96:00:00

set -euo pipefail

set +u
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate /scratch/pawsey1018/rhodgson/envs/mega-contig
set -u

sample=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $SAMP_MAN)
[[ -z "$sample" ]] && { echo "No sample for index $SLURM_ARRAY_TASK_ID"; exit 1; }

# storage directories
# storage directories
submit_dir="${SLURM_SUBMIT_DIR:-$PWD}"

out_root="${OUT_DIR:-$submit_dir/results}"
mkdir -p "$out_root"
out_root="$(realpath "$out_root")"

qc_dir="${out_root}/${sample}/qc"
contig_dir="${out_root}/${sample}/contigs"
mkdir -p "$qc_dir" "$contig_dir"

# prep temporary directory
base_tmp="${SLURM_TMPDIR:-/scratch/pawsey1018/rhodgson/tmp}"
bucket=$((SLURM_JOB_ID % 100))
tmpdir="$base_tmp/bucket_${bucket}/${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
mkdir -p "$tmpdir"
trap 'rm -rf "$tmpdir"; echo "Cleaned tmpdir: $tmpdir"' EXIT
cd "$tmpdir"

# Transfer files
cp "${qc_dir}/${sample}_R1.good.fastq.gz" "$tmpdir/"
cp "${qc_dir}/${sample}_R2.good.fastq.gz" "$tmpdir/"
cp "${qc_dir}/${sample}_R1.single.fastq.gz" "$tmpdir/"
cp "${qc_dir}/${sample}_R2.single.fastq.gz" "$tmpdir/"

# megahit
megahit \
  --presets meta-large \
  -1 "${sample}_R1.good.fastq.gz" \
  -2 "${sample}_R2.good.fastq.gz" \
  -r "${sample}_R1.single.fastq.gz","${sample}_R2.single.fastq.gz" \
  -o "${sample}_contigs"

# Copy results back to storage
cp "${sample}_contigs/final.contigs.fa" "$contig_dir/"

echo "contigs assembly complete for sample $sample"

# End
