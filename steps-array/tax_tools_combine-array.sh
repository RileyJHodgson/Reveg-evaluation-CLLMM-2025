#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --time=0-1
#SBATCH --mem=4000M
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1

set -euo pipefail

set +u
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate kraken-suite
set -u

# Stable base
submit_dir="${SLURM_SUBMIT_DIR:-$PWD}"

# Resolve output root (same convention as the pipeline)
out_root="${OUT_DIR:-$submit_dir/results}"
mkdir -p "$out_root"
out_root="$(realpath "$out_root")"

# KrakenTools location (same convention as your other step)
KR_TOOL_DIR="$submit_dir/resources/KrakenTools-1.2"
KR_TOOL_DIR="$(realpath "$KR_TOOL_DIR")"


# Collect inputs safely
shopt -s nullglob
# mpa_files=( "$out_root"/*/tax/out_mpa_reports/*.mpa_breport ) # for standard mpa output
mpa_files=( "$out_root"/*/tax/out_mpa_reports/*.mpa_breport.CPM_prok.tsv ) # for prok_frac normalised output
shopt -u nullglob

# Combine step
python "${KR_TOOL_DIR}/combine_mpa.py" \
  --input "${mpa_files[@]}" \
  --output "${out_root}/combined_standard.mpa"

echo "Combined MPA written to: ${out_root}/combined_standard.mpa"

# Clean up old files
rm "${mpa_files[@]}"
