#!/bin/bash -l
#SBATCH --job-name=woltka_stratify-species-gene
#SBATCH --output=/home/yac027/metaProtease/outputs/db-map/slurm_out/woltka-stratify-species-gene_%j.out
#SBATCH --error=/home/yac027/metaProtease/outputs/db-map/slurm_out/woltka-stratify-species-gene_%j.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --mail-user=yac027@ucsd.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --partition=short

source activate metaProtease

GENE_BIOM="/home/yac027/metaProtease/outputs/db-map/tables/gene.biom"
SPECIES_MAP="/ddn_scratch/yac027/00_references/genome-to-species.map"
OUTPUT_BIOM="/home/yac027/metaProtease/outputs/db-map/tables/gene_by_species.biom"

mkdir -p "$(dirname "$OUTPUT_BIOM")"

echo "[INFO] Stratifying gene table by species..."
woltka collapse \
  --input "$GENE_BIOM" \
  --map "$SPECIES_MAP" \
  --output "$OUTPUT_BIOM" \
  -f 1 -e

echo "[DONE] Woltka species|gene stratification complete."
echo "Result: $OUTPUT_BIOM"

