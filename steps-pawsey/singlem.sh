#!/bin/bash --login
#SBATCH --account=pawsey1018
#SBATCH --partition=work
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4000M  
#SBATCH --time=24:00:00

set -euo pipefail

set +u
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate /scratch/pawsey1018/rhodgson/envs/singlem
set -u

# storage directories. As a lot of file paths are relative to where script submitted, we need to make sure the paths are carefully specified
# NOTE: "${SLURM_SUBMIT_DIR:-$PWD}" # here we say use the variable $SLURM_SUBMIT_DIR, but if it doesn't exist, use $PWD

sample=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $SAMP_MAN)
[[ -z "$sample" ]] && { echo "No sample for index $SLURM_ARRAY_TASK_ID"; exit 1; }

echo "QC for sample: $sample"

submit_dir="${SLURM_SUBMIT_DIR:-$PWD}"
# sample="${1:?No sample ID provided}"

# specify input/output paths
out_root="${OUT_DIR:-$submit_dir/results}"
out_root="$(realpath "$out_root")"

# specify input paths
qc_dir="${out_root}/${sample}/qc"
qc_dir="$(realpath "$qc_dir")"

# specify output paths
spf_dir="${out_root}/${sample}/singlem"
mkdir -p "$spf_dir"
spf_dir="$(realpath "$spf_dir")"

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

# singlem
# run this once ever. It will donload S5.4.0.GTDB_r226.metapackage_20250331.smpkg.zb and privide the export command below
# singlem data --output-directory /scratch/user/hodg0248/goyder/resources/singlem_data

# call this command from workflow.sh script so easier to change in other projects
# export SINGLEM_METAPACKAGE_PATH='/scratch/user/hodg0248/goyder/resources/singlem_data/S5.4.0.GTDB_r226.metapackage_20250331.smpkg.zb'

# Step one (paired ends) -----------------------------------------------------------------
singlem pipe \
  -1 "${sample}_R1.good.fastq.gz" \
  -2 "${sample}_R2.good.fastq.gz" \
  --threads 4 \
  -p "${sample}.singlem.pe_profile.tsv" \
  --otu-table "${sample}.singlem.pe_otu.tsv"

singlem prokaryotic_fraction \
  -1 "${sample}_R1.good.fastq.gz" \
  -2 "${sample}_R2.good.fastq.gz" \
  -p "${sample}.singlem.pe_profile.tsv" > "${sample}.singlem.pe_spf.tsv"

cp "${sample}.singlem.pe_spf.tsv" "$spf_dir/"
echo "singlem complete for paired ends sample $sample : {sample}.singlem.pe_spf.tsv"

# Step two (single ends: 1 and 2) -----------------------------------------------------------------
singlem pipe \
  -1 "${sample}_R1.single.fastq.gz" \
  --threads 4 \
  -p "${sample}.singlem.s1_profile.tsv" \
  --otu-table "${sample}.singlem.s1_otu.tsv"

singlem prokaryotic_fraction \
  -1 "${sample}_R1.single.fastq.gz" \
  -p "${sample}.singlem.s1_profile.tsv" > "${sample}.singlem.s1_spf.tsv"

cp "${sample}.singlem.s1_spf.tsv" "$spf_dir/"
echo "singlem complete for single end R1 reads sample $sample : {sample}.singlem.s1_spf.tsv"

singlem pipe \
  -1 "${sample}_R2.single.fastq.gz" \
  --threads 4 \
  -p "${sample}.singlem.s2_profile.tsv" \
  --otu-table "${sample}.singlem.s2_otu.tsv"

singlem prokaryotic_fraction \
  -1 "${sample}_R2.single.fastq.gz" \
  -p "${sample}.singlem.s2_profile.tsv" > "${sample}.singlem.s2_spf.tsv"

cp "${sample}.singlem.s2_spf.tsv" "$spf_dir/"
echo "singlem complete for single end R2 reads sample $sample : {sample}.singlem.s2_spf.tsv"

# Step three (R1 reads only -------------------------------------------------------------
singlem pipe \
  -1 "${sample}_R1.good.fastq.gz" \
  --threads 4 \
  -p "${sample}.singlem.r1_profile.tsv" \
  --otu-table "${sample}.singlem.r1_otu.tsv"

singlem prokaryotic_fraction \
  -1 "${sample}_R1.good.fastq.gz" \
  -p "${sample}.singlem.r1_profile.tsv" > "${sample}.singlem.r1_spf.tsv"

cp "${sample}.singlem.r1_spf.tsv" "$spf_dir/"
echo "singlem complete for R1 reads sample $sample : {sample}.singlem.r1_spf.tsv"

# END -----------------------------------------------------------------------------------
