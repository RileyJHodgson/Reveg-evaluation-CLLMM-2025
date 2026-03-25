#!/bin/bash --login
#SBATCH --account=pawsey1018
#SBATCH --partition=long
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem-per-cpu=4000M
#SBATCH --time=96:00:00

set -euo pipefail

# Activate Kraken/Bracken env
set +u
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate k2DB_env
set -u

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

# Sample
sample=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SAMP_MAN")
[[ -z "$sample" ]] && { echo "No sample for index $SLURM_ARRAY_TASK_ID"; exit 1; }

echo "Processing sample: $sample"

# Storage directories
submit_dir="${SLURM_SUBMIT_DIR:-$PWD}"
out_root="${OUT_DIR:-$submit_dir/results}"
out_root="$(realpath "$out_root")"

qc_dir="${out_root}/${sample}/qc"
tax_dir="${out_root}/${sample}/tax"
singlem_dir="${out_root}/${sample}/singlem"

mkdir -p "$tax_dir"

# Check if pipeline already completed for this sample
DONE_FILE="${out_root}/${sample}/tax/out_mpa_reports/${sample}.mpa_breport.CPM_prok.tsv"

if [[ -s "$DONE_FILE" ]]; then
    echo "Kraken2-Bracken-MPA-CPM already completed for $sample — skipping."
    exit 0
fi

database="${TAX_DB:-$submit_dir/resources/k2_standard_20251015}"
database="$(realpath "$database")"

# Scratch setup
base_tmp="${SLURM_TMPDIR:-/scratch/pawsey1018/rhodgson/tmp}"
bucket=$((SLURM_JOB_ID % 100))
tmpdir="$base_tmp/bucket_${bucket}/${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}"

mkdir -p "$tmpdir"
trap 'rm -rf "$tmpdir"; echo "Cleaned tmpdir: $tmpdir"' EXIT
cd "$tmpdir"

echo "Using tmpdir: $tmpdir"

# Copy FASTQ to scratch
cp -n "${qc_dir}/${sample}_R1.good.fastq.gz" .
cp -n "${qc_dir}/${sample}_R2.good.fastq.gz" .

# Kraken2
echo "Running Kraken2..."
nproc

kraken2 --db "$database" \
    --threads $OMP_NUM_THREADS \
    --memory-mapping \
    --report-zero-counts \
    --use-names \
    --paired "${sample}_R1.good.fastq.gz" "${sample}_R2.good.fastq.gz" \
    --output "${sample}.k2_output" \
    --report "${sample}.k2_report"

echo "Kraken2 complete"

# Bracken
echo "Running Bracken..."
nproc

bracken -r 100 -l S -t $OMP_NUM_THREADS \
    -d "$database" \
    -i "${sample}.k2_report" \
    -o "${sample}.bracken_output" \
    -w "${sample}.bracken_report"

echo "Bracken complete"

# Switch to python env
set +u
conda deactivate
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate python_base
set -u

# Convert to MPA
KR_TOOL_DIR="$submit_dir/resources/KrakenTools-1.2"
echo "Running python kreport2mpa.py ..."
nproc

python "$KR_TOOL_DIR/kreport2mpa.py" \
    -r "${sample}.bracken_report" \
    -o "${sample}.mpa_breport" \
    --display-header

echo "MPA conversion complete"

# CPM normalisation
qc_json="$qc_dir/${sample}_fastp.json"
spf_file="$singlem_dir/${sample}.singlem.spf.tsv"

echo "Running python normalise_mpa_cpm_prok.py ..."

python steps-pawsey/normalise_mpa_cpm_prok.py \
  --qc-json "$qc_json" \
  --spf-file "$spf_file" \
  --mpa-file "${sample}.mpa_breport" \
  --out-file "${sample}.mpa_breport.CPM_prok.tsv"

echo "CPM normalisation complete"

# Save ONLY final result
mkdir -p "$tax_dir/out_mpa_reports"

cp "${sample}.mpa_breport.CPM_prok.tsv" "$tax_dir/out_mpa_reports/"

echo "Final output saved"
echo "Pipeline complete for sample: $sample"
