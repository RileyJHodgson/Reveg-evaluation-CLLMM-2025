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
conda activate /scratch/pawsey1018/rhodgson/envs/hostremoval
set -u

sample=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $SAMP_MAN)
[[ -z "$sample" ]] && { echo "No sample for index $SLURM_ARRAY_TASK_ID"; exit 1; }

# storage directories
submit_dir="${SLURM_SUBMIT_DIR:-$PWD}"

out_root="${OUT_DIR:-$submit_dir/results}"
mkdir -p "$out_root"
out_root="$(realpath "$out_root")"

qc_dir="${out_root}/${sample}/qc"
index_dir="${out_root}/${sample}/indexes"
map_dir="${out_root}/${sample}/mapped_reads"
mkdir -p "$map_dir"

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
cp -r "$index_dir" "$tmpdir/"

bowtie2 --local -q \
    -x "./indexes/ref_${sample}" \
    -1 "${sample}_R1.good.fastq.gz" -2 "${sample}_R2.good.fastq.gz" \
    -U "${sample}_R1.single.fastq.gz","${sample}_R2.single.fastq.gz" \
    -p 16 2> "bowtie2_mapLog_${sample}.log" | \
              samtools view -bh | \
              samtools sort -o "${sample}.sorted.bam" && \
              samtools index "${sample}.sorted.bam" && \
              samtools idxstats "${sample}.sorted.bam" > "idxStats_${sample}.txt" && \
              samtools coverage "${sample}.sorted.bam" > "covStats_${sample}.txt"

# Copy results back to storage
cp "bowtie2_mapLog_${sample}.log" "$map_dir"
cp "covStats_${sample}.txt" "$map_dir"

echo "bowtie2 complete for sample $sample"

# End
