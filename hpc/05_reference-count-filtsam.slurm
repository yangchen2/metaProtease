#!/bin/bash -l
#SBATCH --job-name=refcount_featurecounts
#SBATCH --output=../outputs/reference-map/slurm_out/refcount_featurecounts_%A.out
#SBATCH --error=../outputs/reference-map/slurm_out/refcount_featurecounts_%A.err
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#SBATCH --mail-user=yac027@ucsd.edu
#SBATCH --mail-type=END,FAIL

source activate metaProtease

GTF=/home/yac027/metaProtease/outputs/reference-map/featureCounts/strep_proteases_renamed.gtf
BAM_DIR=/home/yac027/metaProtease/outputs/reference-map/bams
OUTDIR=/home/yac027/metaProtease/outputs/reference-map/featureCounts

mkdir -p $OUTDIR

# Step 1: check BAMs and relabel empties
for bam in ${BAM_DIR}/*_filt.bam; do
    n=$(samtools view -c -f 0x2 "$bam")
    if [ "$n" -eq 0 ]; then
        new="${bam/_filt.bam/_filt_empty.bam}"
        echo "Renaming $bam → $new (no paired reads)"
        mv "$bam" "$new"
    fi
done

# Step 2: collect only BAMs with valid paired reads
BAMS=$(ls ${BAM_DIR}/*_filt.bam 2>/dev/null)

if [ -z "$BAMS" ]; then
  echo "No valid filtered BAM files found in ${BAM_DIR}"
  exit 1
fi

# Step 3: run featureCounts in paired-end mode
featureCounts \
  -a $GTF \
  -o ${OUTDIR}/counts.txt \
  -t gene \
  -g gene_id \
  -Q 30 \
  -T $SLURM_CPUS_PER_TASK \
  -p \
  $BAMS

# Step 4: convert to TSV for easier downstream parsing
awk 'BEGIN{OFS="\t"} {print $0}' ${OUTDIR}/counts.txt > ${OUTDIR}/counts.tsv

echo "Finished"

