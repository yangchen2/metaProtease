#!/bin/bash -l
#SBATCH --job-name=db_map
#SBATCH --output=../outputs/db-map/slurm_out/%x_%A_%a.out
#SBATCH --error=../outputs/db-map/slurm_out/%x_%A_%a.err
#SBATCH --cpus-per-task=4
#SBATCH --time=12:00:00
#SBATCH --mem=64G
#SBATCH --partition=short
#SBATCH --array=0-86
#SBATCH --mail-user=yac027@ucsd.edu
#SBATCH --mail-type=END,FAIL

FASTQ_DIR=/qmounts/qiita_data/per_sample_FASTQ/193232
INDEX_DIR=/scratch/qp-woltka/WoLr2/WoLr2
DB_NAME=wolr2
OUT_DIR=../outputs/db-map

FASTQS=($(ls ${FASTQ_DIR}/*_R1*.fastq.gz | sort))
R1=${FASTQS[$SLURM_ARRAY_TASK_ID]}
R2=${R1/_R1/_R2}

if [ ! -f "$R1" ] || [ ! -f "$R2" ]; then
    echo "Error: FASTQ pair not found for task $SLURM_ARRAY_TASK_ID"
    exit 1
fi

source activate metaProtease

echo "[$SLURM_ARRAY_TASK_ID] Mapping $R1 and $R2"
metaProtease db-map map \
    --db-name "$DB_NAME" \
    --index-path "$INDEX_DIR" \
    --r1 "$R1" \
    --r2 "$R2" \
    --outdir "$OUT_DIR" \
    --threads 16 \
    --force

echo 'Finished'
