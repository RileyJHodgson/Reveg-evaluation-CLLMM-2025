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

export PYTHONPATH=/home/hodg0248/miniconda3/envs/kraken-suite/lib/python3.10/site-packages:$PYTHONPATH

sample=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $SAMP_MAN)
[[ -z "$sample" ]] && { echo "No sample for index $SLURM_ARRAY_TASK_ID"; exit 1; }

# storage directories
submit_dir="${SLURM_SUBMIT_DIR:-$PWD}"

# Specify root path (based on OUT_DIR predefined, or from $submit_dir if not )
# Then define full path for out_root so we don't get lost when changing relative paths
out_root="${OUT_DIR:-$submit_dir/results}"
mkdir -p "$out_root"
out_root="$(realpath "$out_root")"

# Paths to kraken tools scripts
KR_TOOL_DIR="$submit_dir/resources/KrakenTools-1.2"
KR_TOOL_DIR="$(realpath "$KR_TOOL_DIR")"

echo "Looking for  Kraken tools scripts in: $KR_TOOL_DIR"

# Paths to output scripts
tax_dir="${out_root}/${sample}/tax"

mpa_dir="${tax_dir}/out_mpa_reports"
krona_dir="${tax_dir}/out_krona_reports"
mkdir -p "$mpa_dir" "$krona_dir"

mpa_dir="$(realpath "$mpa_dir")"
krona_dir="$(realpath "$krona_dir")"

# run scripts
python  "$KR_TOOL_DIR"/kreport2mpa.py \
    -r "$tax_dir"/"${sample}.bracken_report" \
    -o "$mpa_dir"/"${sample}.mpa_breport" \
    --display-header
python "$KR_TOOL_DIR"/kreport2krona.py \
        --report "$tax_dir"/"${sample}.bracken_report" \
        --output "$krona_dir"/"${sample}.krona" \
        --no-intermediate-ranks

echo "completed kreport2mpa.py and kreport2krona.py jobs"

# Clean up old files
# rm "$tax_dir"/"${sample}.bracken_report"
# rm "$tax_dir"/"${sample}.bracken_output"
