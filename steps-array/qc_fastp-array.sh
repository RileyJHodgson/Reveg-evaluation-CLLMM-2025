#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --time=4-0
#SBATCH --mem=16000M
#SBATCH --cpus-per-task=4
#SBATCH --nodes=1

set -euo pipefail

set +u
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate hostremoval
set -u

# storage directories. As a lot of file paths are relative to where script submitted, we need to make sure the paths are carefully specified
# NOTE: "${SLURM_SUBMIT_DIR:-$PWD}" # here we say use the variable $SLURM_SUBMIT_DIR, but if it doesn't exist, use $PWD


sample=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $SAMP_MAN)
[[ -z "$sample" ]] && { echo "No sample for index $SLURM_ARRAY_TASK_ID"; exit 1; }

echo "QC for sample: $sample"

submit_dir="${SLURM_SUBMIT_DIR:-$PWD}"
# sample="${1:?No sample ID provided}"

# specify input paths
raw_reads_dir="${DATA_DIR:-$submit_dir/data}"
raw_reads_dir="$(realpath "$raw_reads_dir")"

# specify output paths
out_root="${OUT_DIR:-$submit_dir/results}"
mkdir -p "$out_root"
out_root="$(realpath "$out_root")"
qc_dir="${out_root}/${sample}/qc"
mkdir -p "$qc_dir"
qc_dir="$(realpath "$qc_dir")"

# Temporary directory on compute node
tmpdir="${SLURM_TMPDIR:-/scratch/user/${USER}/${SLURM_JOB_ID}}"
mkdir -p "$tmpdir"
trap 'rm -rf "$tmpdir"' EXIT

echo "Using tmpdir: $tmpdir"

# Transfer files
cp "${raw_reads_dir}/${sample}_R1.fastq.gz" "$tmpdir/"
cp "${raw_reads_dir}/${sample}_R2.fastq.gz" "$tmpdir/"

# tmpdir
echo "Using tmpdir: $tmpdir"
cd "$tmpdir"

# fastp
fastp --length_required 60 \
    --average_qual 25 \
    --n_base_limit 1 \
    --dedup \
    --trim_front1 10 \
    --trim_tail1 5 \
    --trim_front2 10 \
    --trim_tail2 5 \
    --cut_mean_quality 30 \
    --cut_window_size 10 \
    --cut_front \
    --cut_tail \
    --out1 ${sample}_R1.good.fastq.gz \
    --out2 ${sample}_R2.good.fastq.gz \
    --unpaired1 ${sample}_R1.single.fastq.gz \
    --unpaired2 ${sample}_R2.single.fastq.gz \
    --in1 ${sample}_R1.fastq.gz \
    --in2 ${sample}_R2.fastq.gz \
    --html ${sample}_fastp.html \
    --json ${sample}_fastp.json

# Copy results back to storage
cp "${sample}_R1.good.fastq.gz" "$qc_dir/"
cp "${sample}_R2.good.fastq.gz" "$qc_dir/"
cp "${sample}_R1.single.fastq.gz" "$qc_dir/"
cp "${sample}_R2.single.fastq.gz" "$qc_dir/"
cp "${sample}_fastp.html" "$qc_dir/"
cp "${sample}_fastp.json" "$qc_dir/"

echo "QC complete for sample $sample"

# rm ${raw_reads_dir}/${sample}_R1.fastq.gz #these delete the raw reads
# rm ${raw_reads_dir}/${sample}_R2.fastq.gz #these delete the raw reads
# echo "Raw fastqs for sample $sample removed"

# End
