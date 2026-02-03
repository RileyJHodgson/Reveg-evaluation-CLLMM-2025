# Reveg-evaluation-CLLMM-2025
Scripts for analysis of revegetation through CLLMM region 2025 surveys

Scripts for analysis of: 
- Vegetation communities (Veg-analysis-CLLMM.R),
- Bird communities (Bird-analysis-CLLMM.R)
- Soil physicochemical conditions (Soil-analysis-CLLMM.R)
- Home made functions used in these scripts (Permute-LMEM-Toolkit.R)
- Shotgun metagenomics bioinformatics workflows intended to be run on SLURM (workflow_goyder-array.sh)
  - Individual step scripts used by the workflow are located in steps-array/.
  - These step scripts are designed to be executed via SLURM job arrays and job dependencies, not run directly
  - Only the top-level workflow script (workflow_goyder-array.sh) is submitted directly to SLURM; step scripts are called internally by the workflow


Scripts are a work in progress

To be added:
metagenomic analysis (tax and functions .R)
