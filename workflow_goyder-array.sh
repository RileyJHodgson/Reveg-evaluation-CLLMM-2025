#!/bin/bash
#SBATCH --job-name=goyder
#SBATCH --time=00:10:00
#SBATCH --mem=4000M
#SBATCH --cpus-per-task=1
#SBATCH --output=workflow_%j.out
#SBATCH --error=workflow_%j.err

# FULL example: sbatch workflow_goyder-array.sh --manifest ./resources/sample_manifest_full-unique.txt --raw-fastqs ./merged_fastqs --output ./complete_results
# TEST example: sbatch workflow_goyder-array.sh --manifest ./resources/sample-manifest.txt --raw-fastqs ./fastqs --output ./results_array

# The k2 database can be manually set using: --database ./resources/k2_standard_20251015

set -euo pipefail

if [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then
  cd "$SLURM_SUBMIT_DIR"
fi

echo "Working directory: $(pwd)"

set +u
source /home/hodg0248/miniconda3/etc/profile.d/conda.sh
#source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate base_python
set -u

# to address python pandas issues...
export PYTHONPATH=/home/hodg0248/miniconda3/envs/base_python/lib/python3.10/site-packages:$PYTHONPATH

# Draft pipeline for amr and eskape query and assignment
# Riley 2026-1-13

# Default values. NOTE may need to change TAX_DB if using a different version...
DATA_DIR="./data"
OUT_DIR="./results"
TAX_DB="./resources/k2_standard_20251015"
SAMP_MAN="./sample-manifest.txt"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --raw-fastqs|-f)
      [[ -n "${2:-}" && "${2:0:1}" != "-" ]] || { echo "Error: $1 requires a path"; exit 1; }
      DATA_DIR="$2"
      shift 2
      ;;
    --output|-o)
      [[ -n "${2:-}" && "${2:0:1}" != "-" ]] || { echo "Error: $1 requires a path"; exit 1; }
      OUT_DIR="$2"
      shift 2
      ;;
    --database|-d)
      [[ -n "${2:-}" && "${2:0:1}" != "-" ]] || { echo "Error: $1 requires a path"; exit 1; }
      TAX_DB="$2"
      shift 2
      ;;
    --manifest|-m)
      [[ -n "${2:-}" && "${2:0:1}" != "-" ]] || { echo "Error: $1 requires a path to sample-manifest.txt"; exit 1; }
      SAMP_MAN="$2"
      shift 2
      ;;
    --) # end of flags
      shift
      break
      ;;
    *)
      echo "Unknown flag: $1"
      exit 1
      ;;
  esac
done

# Export so step-scripts can see them
export DATA_DIR
export OUT_DIR
export TAX_DB
export SAMP_MAN
export SINGLEM_METAPACKAGE_PATH='/scratch/user/hodg0248/goyder/resources/singlem_data/S5.4.0.GTDB_r226.metapackage_20250331.smpkg.zb'

echo "Sample manifest:    $SAMP_MAN"
echo "Data dir:           $DATA_DIR"
echo "Output dir:         $OUT_DIR"
echo "Taxonomy DB:        $TAX_DB"

# Check if taxonomy database exists
if [[ -d "$TAX_DB" ]]; then
  echo "Taxonomy database directory found: $TAX_DB"
else
  echo "ERROR: Required taxonomy database directory missing at $TAX_DB. Exiting."
  exit 1
fi

# Create directory for log files (.out,.err)
logs_dir="${OUT_DIR}/logs"
mkdir -p "$logs_dir"

# Array steps ----------------------------------------------------------------------------
NSAMPLES=$(wc -l < "$SAMP_MAN")

# *** Phase 1 *** ------------------------------------------------------------------------
# --- Step 1.1 QC on raw reads ---
jid_QC=$(sbatch \
  --export=DATA_DIR,OUT_DIR,SAMP_MAN \
  --parsable \
  --array=1-"$NSAMPLES" \
  --output="${logs_dir}/fp_qc_%A_%a.out" \
  --error="${logs_dir}/fp_qc_%A_%a.err" \
  steps-array/qc_fastp-array.sh)

# --- Step 1.2 Estimate prokaryotic fractions ---
jid_SPF=$(sbatch \
  --export=DATA_DIR,OUT_DIR,SAMP_MAN,SINGLEM_METAPACKAGE_PATH \
  --parsable \
  --array=1-"$NSAMPLES" \
  --output="${logs_dir}/spf_%A_%a.out" \
  --error="${logs_dir}/spf_%A_%a.err" \
  --dependency=afterok:$jid_QC \
  steps-array/singlem-array.sh)

# *** Phase 2 *** ------------------------------------------------------------------------
# --- Step 2.1 Standard db assignment ---
jid_kraken=$(sbatch --export=DATA_DIR,OUT_DIR,SAMP_MAN,TAX_DB \
  --parsable \
  --array=1-"$NSAMPLES" \
  --output="${logs_dir}/k2_%A_%a.out" \
  --error="${logs_dir}/k2_%A_%a.err" \
  --dependency=afterok:$jid_QC \
  steps-array/kraken2-array.sh)

# --- Step 2.2 bracken estimate abundances ---
jid_bracken=$(sbatch --export=DATA_DIR,OUT_DIR,SAMP_MAN,TAX_DB \
  --parsable \
  --array=1-"$NSAMPLES" \
  --output="${logs_dir}/br_%A_%a.out" \
  --error="${logs_dir}/br_%A_%a.err" \
  --dependency=afterok:$jid_kraken \
  steps-array/bracken-array.sh)

# --- Step 2.3 reformat taxonomy reports to mpa ---
jid_ktools=$(sbatch --export=ALL \
  --parsable \
  --array=1-"$NSAMPLES" \
  --output="${logs_dir}/kt_%A_%a.out" \
  --error="${logs_dir}/kt_%A_%a.err" \
  --dependency=afterok:$jid_bracken \
  steps-array/tax_tools-array.sh)

# --- Step 2.4 Noramlise taxonomy output to counts per million prokaryotes ---
jid_cpm=$(sbatch --export=ALL \
  --parsable \
  --array=1-"$NSAMPLES" \
  --output="${logs_dir}/cpm_tax_%A_%a.out" \
  --error="${logs_dir}/cpm_tax_%A_%a.err" \
  --dependency=afterok:$jid_ktools \
  steps-array/CPM_prok-braken-array.sh)

# --- Step 2.5 Combine Kraken2-bracken standard db (non-array job) ---
jid_combine_pfp=$(sbatch --export=ALL \
  --parsable  \
  --output="${logs_dir}/combine_mpa_%j.out"  \
  --error="${logs_dir}/combine_mpa_%j.err" \
  --dependency=afterok:$jid_cpm \
  steps-array/tax_tools_combine-array.sh)

# --- Step 3.1 FUN superfocus ---
jid_superfocus=$(sbatch --export=DATA_DIR,OUT_DIR,SAMP_MAN \
  --parsable \
  --array=1-"$NSAMPLES" \
  --output="${logs_dir}/sf_%A_%a.out" \
  --error="${logs_dir}/sf_%A_%a.err" \
  --dependency=afterok:$jid_QC \
  steps-array/superfocus-array.sh)

# --- Step 3.2 Noramlise superfocus output to counts per million prokaryotes ---
jid_CPM_sf=$(sbatch --export=ALL \
  --parsable \
  --array=1-"$NSAMPLES" \
  --output="${logs_dir}/prok_sf_%A_%a.out" \
  --error="${logs_dir}/prok_sf_%A_%a.err" \
  --dependency=afterok:$jid_superfocus \
  steps-array/CPM_prok-superfocus-array.sh)

# --- Step 3.2 Combine superfocus outputs (non-array job) ---
jid_combine_sf=$(sbatch --export=ALL \
  --parsable \
  --output="${logs_dir}/combine_sf_%j.out" \
  --error="${logs_dir}/combine_sf_%j.err" \
  --dependency=afterok:$jid_CPM_sf \
  steps-array/combine_superfocus_output-array.py \
  --recursive \
  --output-dir "$OUT_DIR")

# *** Phase 3 *** ------------------------------------------------------------------------
# --- Step 4 Contigs megahit ---
jid_megahit=$(sbatch --export=ALL \
  --parsable \
  --array=1-"$NSAMPLES" \
  --output="${logs_dir}/mh_%A_%a.out" \
  --error="${logs_dir}/mh_%A_%a.err" \
  --dependency=afterok:$jid_QC \
  steps-array/megahit-array.sh)

# --- Step 5 Index bowtie2 ---
jid_bt2_index=$(sbatch --export=ALL \
  --parsable \
  --array=1-"$NSAMPLES" \
  --output="${logs_dir}/bti_%A_%a.out" \
  --error="${logs_dir}/bti_%A_%a.err" \
  --dependency=afterok:$jid_megahit \
  steps-array/bt2_index-array.sh)

# --- Step 6 Map bowtie2 ---
jid_bt2_map=$(sbatch --export=DATA_DIR,OUT_DIR,SAMP_MAN \
  --parsable \
  --array=1-"$NSAMPLES" \
  --output="${logs_dir}/btm_%A_%a.out" \
  --error="${logs_dir}/btm_%A_%a.err" \
  --dependency=afterok:$jid_bt2_index \
  steps-array/bt2_map-array.sh)

# --- Step 7 AMR amrfinderplus ---
jid_amr=$(sbatch --export=ALL \
  --parsable \
  --array=1-"$NSAMPLES" \
  --output="${logs_dir}/amr_%A_%a.out" \
  --error="${logs_dir}/amr_%A_%a.err" \
  --dependency=afterok:$jid_megahit \
  steps-array/amrfinder-array.sh)

# --- Step 8 combine AMR genes and mapping files ---
jid_amr_map_comb=$(sbatch --export=ALL \
  --parsable \
  --output="${logs_dir}/fpkm_%j.out" \
  --error="${logs_dir}/fpkm_%j.err" \
  --time=1-0 \
  --mem=8000M \
  --cpus-per-task=2 \
  --dependency=afterok:$jid_bt2_map:$jid_amr \
  steps-array/combine-amr-metrics.py)

# *** Phase 4 *** ------------------------------------------------------------------------
# --- Finalise outputs/clean-up files ---

# remove fastp outputs once completed jid_kraken, jid_kraken_oom, jid_superfocus and jid_bt2_map complete
jid_finish_qc=$(sbatch --export=ALL \
  --parsable \
  --array=1-"$NSAMPLES" \
  --output="${logs_dir}/rm_qc_%A_%a.out" \
  --error="${logs_dir}/rm_qc_%A_%a.err" \
  --dependency=afterok:$jid_kraken:$jid_superfocus:$jid_bt2_map:$jid_SPF \
  steps-array/finish-QC-array.sh)

# remove functional output directories after completing: jid_combine_sf
jid_finish_qc=$(sbatch --export=ALL \
  --parsable \
  --array=1-"$NSAMPLES" \
  --output="${logs_dir}/rm_sf_%A_%a.out" \
  --error="${logs_dir}/rm_sf_%A_%a.err" \
  --dependency=afterok:$jid_combine_sf \
  steps-array/finish-func-array.sh)

# # remove the contigs files once jid_megahit, jid_bt2_index and jid_amr complete
jid_finish_contigs=$(sbatch --export=ALL \
  --parsable \
  --array=1-"$NSAMPLES" \
  --output="${logs_dir}/rm_contigs_%A_%a.out" \
  --error="${logs_dir}/rm_contigs_%A_%a.err" \
  --dependency=afterok:$jid_megahit:$jid_bt2_index:$jid_amr \
  steps-array/finish-contigs-array.sh)

# # remove the amr gene and large bt2 mapping files once jid_amr_map_comb complete
jid_finish_map_amr_waste=$(sbatch --export=ALL \
  --parsable \
  --array=1-"$NSAMPLES" \
  --output="${logs_dir}/rm_amr_%A_%a.out" \
  --error="${logs_dir}/rm_amr_%A_%a.err" \
  --dependency=afterok:$jid_amr_map_comb \
  steps-array/finish-amr-mapping-array.sh)

# Summary
echo ""
echo ""
echo "Summary:"
echo "  - Taxonomy output files will be found in $OUT_DIR/combined_standard.mpa, a sample x species table with kraken2 standard database assign taxonomy in counts per million prokaryotic reads."
echo "  - Functional output files will be found in $OUT_DIR/functions_counts.csv, a table with all SEED subsystem assigned functions normalised to produce counts per million prokaryotic reads."
echo "  - AMR gene count output will be found in $OUT_DIR/amr-reads-summary.csv. The main metric calculated is FPKM_prokaryotes, see below."
echo ""
echo ""
echo "Details:"
echo "  Taxonomy, functional and AMR outputs have all been normalised for prokaryotic fraction using singlem, and sampling depth via total abundances"
echo "  - counts_per_million_prokaryotes = read_counts / ((total_reads_persample x singlem_prokaryotic_fraction) ) x 10^6"
echo "  - fragments_per_kilobase_per_million_prokaryotes = (10^9 x ARG_reads)/(ARG_length x total_reads_persample x singlem_prokaryotic_fraction)"
echo ""
echo "  This accounts for:..."
echo "          1 *Sequencing depth* (library size) via the CPM-style scaling"
echo "          2 *Between-sample variation in prokaryotic content* via the SingleM prokaryotic fraction in the denominator"
echo "   I.E., Estimated prokaryote-normalised read abundance: 'reads assigned to a taxon per million estimated prokaryotic reads in that sample'."
echo ""
echo ""
echo "All interim files will deleted upon conclusion of the pipeline"

# *** END *** ------------------------------------------------------------------------
