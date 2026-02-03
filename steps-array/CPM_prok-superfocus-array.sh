#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --time=1-0
#SBATCH --mem=12000M
#SBATCH --cpus-per-task=2
#SBATCH --nodes=1

set -euo pipefail

submit_dir="${SLURM_SUBMIT_DIR:-$PWD}"

sample=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $SAMP_MAN)
[[ -z "$sample" ]] && { echo "No sample for index $SLURM_ARRAY_TASK_ID"; exit 1; }

out_root="${OUT_DIR:-$submit_dir/results}"
out_root="$(realpath "$out_root")"

# input directories
qc_dir="$(realpath "$out_root/$sample/qc")"
singlem_dir="$(realpath "$out_root/$sample/singlem")"
sf_dir="$(realpath "$out_root/$sample/func")"

# input files
qc_json="$qc_dir/${sample}_fastp.json"
spf_file="$singlem_dir/${sample}.singlem.r1_spf.tsv"
sf_file="$sf_dir/norm_${sample}_all_levels_and_function.xls"
out_file="$sf_dir/norm_${sample}_all_levels_and_function-cpm_prok.xls"

python steps-array/normalise_sf_cpm_prok.py \
  --qc-json "$qc_json" \
  --spf-file "$spf_file" \
  --sf-file "$sf_file" \
  --out-file "$out_file"

echo "CPM_prok normalisation complete for sample $sample"
echo "Output: $out_file"
