#!/bin/bash -l
#SBATCH --job-name=refmap_index
#SBATCH --output=../outputs/reference-map/slurm_out/%x_%j.out
#SBATCH --error=../outputs/reference-map/slurm_out/%x_%j.err
#SBATCH --cpus-per-task=4
#SBATCH --time=01:00:00
#SBATCH --mem=16G
#SBATCH --partition=highmem
#SBATCH --mail-user=yac027@ucsd.edu
#SBATCH --mail-type=END,FAIL

FASTA=/home/yac027/Gallo_projs/gallo_AD_peds/protease_mining/fasta/strep_proteases.fasta
OUTDIR=../outputs/reference-map # Needs to exist all the way down to slurm_out

source activate metaProtease

echo "Building Bowtie2 index..."
reference-map index \
  --fasta "$FASTA" \
  --outdir "$OUTDIR" \
  --threads $SLURM_CPUS_PER_TASK

