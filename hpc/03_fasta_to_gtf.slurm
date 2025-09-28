#!/bin/bash -l
#SBATCH --job-name=create_gtf
#SBATCH --output=../outputs/reference-map/slurm_out/create_gtf_%A.out
#SBATCH --error=../outputs/reference-map/slurm_out/create_gtf_%A.err
#SBATCH --time=01:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1

source activate metaProtease

cd ../metaProtease

python create_gtf.py \
  -i /home/yac027/Gallo_projs/gallo_AD_peds/protease_mining/fasta/strep_proteases_renamed.fasta \
  -o /home/yac027/metaProtease/outputs/reference-map/featureCounts \
  -l INFO

echo 'Finished'
