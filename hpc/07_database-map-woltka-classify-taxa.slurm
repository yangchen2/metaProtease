#!/bin/bash -l
#SBATCH --job-name=woltka_classify-taxa
#SBATCH --output=/home/yac027/metaProtease/outputs/db-map/slurm_out/woltka-classify-taxa_%j.out
#SBATCH --error=/home/yac027/metaProtease/outputs/db-map/slurm_out/woltka-classify-taxa_%j.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --mail-user=yac027@ucsd.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --partition=short

# ===========================
# 1. Activate conda environment
# ===========================
source activate metaProtease

# ===========================
# 2. Define input/output paths
# ===========================
SAM_DIR="/home/yac027/metaProtease/outputs/db-map/sams"
OUTPUT_PATH="/home/yac027/metaProtease/outputs/db-map/tables"
SLURM_OUT="/home/yac027/metaProtease/outputs/db-map/slurm_out"

mkdir -p "$OUTPUT_PATH" "$SLURM_OUT"

echo "[INFO] Running Woltka classify"
echo "       Input directory:  $SAM_DIR"
echo "       Output directory: $OUTPUT_PATH"

# ===========================
# 3. Wolr2 taxonomy files
# ===========================
TAXONOMY_DIR="/projects/wol/qiyun/wol2/taxonomy"
MAP_FILE="${TAXONOMY_DIR}/taxid.map"
NODES_FILE="${TAXONOMY_DIR}/nodes.dmp"
NAMES_FILE="${TAXONOMY_DIR}/names.dmp"

# ===========================
# 4. Run classification
# ===========================
woltka classify \
  --input "$SAM_DIR" \
  --output "$OUTPUT_PATH" \
  --map "$MAP_FILE" \
  --nodes "$NODES_FILE" \
  --names "$NAMES_FILE" \
  --rank phylum,genus,species

echo "[DONE] Finished Woltka classify"

