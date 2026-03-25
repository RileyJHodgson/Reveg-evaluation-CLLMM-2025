#!/bin/bash --login
#SBATCH --account=pawsey1018
#SBATCH --partition=work
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=4000M
#SBATCH --time=24:00:00

set -euo pipefail

# Initialise conda (no global activation)
set +u
source "$(conda info --base)/etc/profile.d/conda.sh"
set -u

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

# Get sample from manifest
sample=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SAMP_MAN")
[[ -z "$sample" ]] && { echo "No sample for index $SLURM_ARRAY_TASK_ID"; exit 1; }

submit_dir="${SLURM_SUBMIT_DIR:-$PWD}"
out_root="${OUT_DIR:-$submit_dir/results}"
mkdir -p "$out_root"
out_root="$(realpath "$out_root")"

qc_dir="${out_root}/${sample}/qc"
singlem_dir="${out_root}/${sample}/singlem"
fun_dir="${out_root}/${sample}/func"
mkdir -p "$fun_dir"

# Required input files
qc_json="${qc_dir}/${sample}_fastp.json"
spf_file="${singlem_dir}/${sample}.singlem.spf.tsv"

OUT1="${fun_dir}/norm_${sample}_all_levels_and_function-cpm_prok.xls"
if [[ -s "$OUT1" ]]; then
    echo "QC already completed for $sample — skipping."
    exit 0
fi

# Superfocus database path
database="/scratch/pawsey1018/rhodgson/conda_envs/super_env/lib/python3.14/site-packages/superfocus_app/"

# Temporary working directory
base_tmp="${SLURM_TMPDIR:-/scratch/pawsey1018/rhodgson/tmp}"
bucket=$((SLURM_JOB_ID % 100))
tmpdir="$base_tmp/bucket_${bucket}/${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
mkdir -p "$tmpdir"
trap 'rm -rf "$tmpdir"; echo "Cleaned tmpdir: $tmpdir"' EXIT
cd "$tmpdir"

echo "Working in temporary directory: $tmpdir"

# Copy + unzip FASTQ
cp -n "${qc_dir}/${sample}_R1.good.fastq.gz" .
gunzip "${sample}_R1.good.fastq.gz"

# Run Superfocus (superfocusonly env)
echo "Running Superfocus for sample: $sample"
nproc

conda run -n super_env superfocus \
  --query "${sample}_R1.good.fastq" \
  --output_directory "super_out" \
  --alternate_directory "${database}" \
  --output_prefix "norm_${sample}_" \
  --threads $OMP_NUM_THREADS \
  --aligner diamond \
  --database DB_100 \
  --normalise_output 1

sf_raw="super_out/norm_${sample}_all_levels_and_function.xls"
sf_cpm="${fun_dir}/norm_${sample}_all_levels_and_function-cpm_prok.xls"

# Run CPM_prok normalisation (base_python env)
echo "Running CPM_prok normalisation for sample: $sample"
nproc

conda run -n base_python python \
  "$submit_dir/steps-pawsey/normalise_sf_cpm_prok.py" \
  --qc-json "$qc_json" \
  --spf-file "$spf_file" \
  --sf-file "$sf_raw" \
  --out-file "$sf_cpm"

echo "Superfocus + CPM_prok complete for sample: $sample"
echo "Final output: $sf_cpm"
nproc

# tmpdir auto-cleaned via trap
