#!/bin/bash -l
#SBATCH --job-name=woltka_classify-gene
#SBATCH --output=/home/yac027/metaProtease/outputs/db-map/slurm_out/woltka-classify-gene_%j.out
#SBATCH --error=/home/yac027/metaProtease/outputs/db-map/slurm_out/woltka-classify-gene_%j.err
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
OUTPUT_BIOM="/home/yac027/metaProtease/outputs/db-map/tables/gene.biom"
COORDS_FILE="/projects/wol/qiyun/wol2/proteins/coords.txt.xz"

# Create output directory if it doesn’t exist
mkdir -p "$(dirname "$OUTPUT_BIOM")"

echo "[INFO] Starting Woltka classify (gene-level)"
echo "       Input:  $SAM_DIR"
echo "       Coords: $COORDS_FILE"
echo "       Output: $OUTPUT_BIOM"

# ===========================
# 3. Run Woltka classify (gene-level)
# ===========================
woltka classify \
  --input "$SAM_DIR" \
  --coords "$COORDS_FILE" \
  --output "$OUTPUT_BIOM"

# ===========================
# 4. Summary
# ===========================
echo "[DONE] Woltka gene-level classification complete."
echo "Result: $OUTPUT_BIOM"

