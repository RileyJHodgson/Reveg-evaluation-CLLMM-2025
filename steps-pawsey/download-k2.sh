#!/bin/bash --login
#SBATCH --account=pawsey1018
#SBATCH --job-name=k2_download
#SBATCH --partition=copy
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=48:00:00

cd /scratch/pawsey1018/rhodgson/goyder/resources/k2_pluspf_20251015

wget https://genome-idx.s3.amazonaws.com/kraken/k2_pluspf_20251015.tar.gz

# extract database files
tar -xvf k2_pluspf_20251015.tar.gz

# tidy up?
rm k2_pluspf_20251015.tar.gz