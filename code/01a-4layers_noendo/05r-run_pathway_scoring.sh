#!/bin/bash

# Simple SLURM sbatch example
#SBATCH --job-name=pthwy_scoring
#SBATCH --ntasks=1
#SBATCH --time=3-00
#SBATCH --mem-per-cpu=48G
#SBATCH --partition=ncpu

set -e -o pipefail

cd /nemo/lab/gandhis/home/users/grantpm/
ml purge 
. load_panpipes.sh

cd MRIxST/code/01a-4layers_noendo
echo "Start running python"

python3 "05-pathway_scoring.py"

echo "All done."
