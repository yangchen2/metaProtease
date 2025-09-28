#!/bin/bash -l
#SBATCH --job-name=tree
#SBATCH --output=../outputs/reference-map/slurm_out/tree_%j.out
#SBATCH --error=../outputs/reference-map/slurm_out/tree_%j.err
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --mail-user=yac027@ucsd.edu
#SBATCH --mail-type=END,FAIL

# Activate environment
source activate metaProtease

echo 'Checking for input FASTA'
# Define input FASTA (update path as needed)
FASTA="/home/yac027/Gallo_projs/gallo_AD_peds/protease_mining/fasta/strep_proteases_renamed.fasta"

# Run your alignment script
echo 'Running MAFFT'
cd ../
./metaProtease/tree.sh "$FASTA"

echo 'Finished'
