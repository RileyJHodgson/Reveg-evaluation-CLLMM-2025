#!/bin/bash --login
#SBATCH --account=pawsey1018
#SBATCH --partition=work
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=4000M  
#SBATCH --time=24:00:00

set -euo pipefail

submit_dir="${SLURM_SUBMIT_DIR:-$PWD}"

sample=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $SAMP_MAN)
[[ -z "$sample" ]] && { echo "No sample for index $SLURM_ARRAY_TASK_ID"; exit 1;}

out_root="${OUT_DIR:-$submit_dir/results}"
out_root="$(realpath "$out_root")"

# input directories
qc_dir="$(realpath "$out_root/$sample/qc")"
singlem_dir="$(realpath "$out_root/$sample/singlem")"
mpa_dir="$(realpath "$out_root/$sample/tax/out_mpa_reports")"

# input files
qc_json="$qc_dir/${sample}_fastp.json"
spf_file="$singlem_dir/${sample}.singlem.pe_spf.tsv"
mpa_file="$mpa_dir/${sample}.mpa_breport"
out_file="$mpa_dir/${sample}.mpa_breport.CPM_prok.tsv"

# use python script to extract and normalise values
python steps-array/normalise_mpa_cpm_prok.py \
  --qc-json "$qc_json" \
  --spf-file "$spf_file" \
  --mpa-file "$mpa_file" \
  --out-file "$out_file"

echo "CPM_prok normalisation complete for sample $sample"
echo "Output: $out_file"

rm $mpa_file
echo "removed the original mpa file for sample: $sample"

