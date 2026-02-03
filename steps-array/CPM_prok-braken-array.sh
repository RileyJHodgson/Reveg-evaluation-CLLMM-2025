#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --time=1-0
#SBATCH --mem=12000M
#SBATCH --cpus-per-task=2
#SBATCH --nodes=1

set -euo pipefail

# --- conda ---
#set +u
#source "$(conda info --base)/etc/profile.d/conda.sh"
#conda activate singlem
#set -u

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
spf_file="$singlem_dir/${sample}.singlem.spf.tsv"
mpa_file="$mpa_dir/${sample}.mpa_breport"
out_file="$mpa_dir/${sample}.mpa_breport.CPM_prok.tsv"

# use python script to extract and normalise values
python steps-array/normalise_mpa_cpm_prok.py \
  --qc-json "$qc_json" \
  --spf-file "$spf_file" \
  --mpa-file "$mpa_file" \
  --out-file "$out_file"

# --- extract values ---
# total_reads=$(jq '.summary.after_filtering.total_reads' "$qc_json")
#total_reads=$(python3 -c "import json; print(json.load(open('$qc_json'))['summary']['after_filtering']['total_reads'])")
#p rok_frac=$(awk 'NR==2 { print $4 }' "$spf_file")

# denominator: total prokaryotic reads
#denom=$(awk -v tr="$total_reads" -v pf="$prok_frac" \
#  'BEGIN { print tr * (pf / 100.0) }')

# --- normalise ---
#awk -v denom="$denom" '
#BEGIN { OFS="\t" }
#NR==1 { print $0, "CPM_prok"; next }
#{
#  $2 = ($2 / denom) * 1e6
#  print
#}
#' "$mpa_file" > "$out_file"

echo "CPM_prok normalisation complete for sample $sample"
echo "Output: $out_file"
