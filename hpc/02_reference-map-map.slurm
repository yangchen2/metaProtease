#!/bin/bash -l
#SBATCH --job-name=refmap_map
#SBATCH --output=../outputs/reference-map/slurm_out/%x_%A_%a.out
#SBATCH --error=../outputs/reference-map/slurm_out/%x_%A_%a.err
#SBATCH --cpus-per-task=4
#SBATCH --time=04:00:00
#SBATCH --mem=16G
#SBATCH --partition=highmem
#SBATCH --array=0-86
#SBATCH --mail-user=yac027@ucsd.edu
#SBATCH --mail-type=END,FAIL

FASTQ_DIR=/qmounts/qiita_data/per_sample_FASTQ/193232
INDEX_DIR=/home/yac027/metaProtease/outputs/reference-map/index
OUTDIR=../outputs/reference-map

FASTQS=($(ls ${FASTQ_DIR}/*_R1*.fastq.gz | sort))
R1=${FASTQS[$SLURM_ARRAY_TASK_ID]}
R2=${R1/_R1/_R2}

if [ ! -f "$R1" ] || [ ! -f "$R2" ]; then
    echo "Error: FASTQ pair not found for task $SLURM_ARRAY_TASK_ID"
    exit 1
fi

source activate metaProtease

echo "[$SLURM_ARRAY_TASK_ID] Mapping $R1 and $R2"
reference-map map \
  --r1 "$R1" \
  --r2 "$R2" \
  --index-dir "$INDEX_DIR" \
  --outdir "$OUTDIR" \
  --threads $SLURM_CPUS_PER_TASK

echo "Finished!"
