#!/bin/bash -l
#SBATCH --job-name=refcount_filtsam
#SBATCH --output=../outputs/reference-map/slurm_out/refcount_filtsam_%A.out
#SBATCH --error=../outputs/reference-map/slurm_out/refcount_filtsam_%A.err
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#SBATCH --mail-user=yac027@ucsd.edu
#SBATCH --mail-type=END,FAIL

source activate metaProtease

SAM_DIR=/home/yac027/metaProtease/outputs/reference-map/sams
INPUTS=$(ls ${SAM_DIR}/*.sam 2>/dev/null)

if [ -z "$INPUTS" ]; then
  echo "No SAM files found in ${SAM_DIR}"
  exit 1
fi

for sam in $INPUTS; do
  reference-count \
    -i "$sam" \
    -o /home/yac027/metaProtease/outputs/reference-map \
    -t $SLURM_CPUS_PER_TASK
done

echo "Finished"

