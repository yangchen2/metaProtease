#!/bin/bash -l
#SBATCH --job-name=tree
#SBATCH --output=../outputs/reference-map/slurm_out/tree_%j.out
#SBATCH --error=../outputs/reference-map/slurm_out/tree_%j.err
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --mail-user=yac027@ucsd.edu
#SBATCH --mail-type=END,FAIL

set -euo pipefail

echo "[INFO] Starting tree pipeline"
echo "[INFO] Environment: metaProtease"

# Define input FASTAs
DNA_FASTA="/home/yac027/Gallo_projs/gallo_AD_peds/protease_mining/fasta/strep_proteases_renamed.fasta"
PROT_FASTA="/home/yac027/Gallo_projs/gallo_AD_peds/protease_mining/fasta/strep_proteases_renamed.faa"
OUTDIR="../outputs/reference-map/tree"

# Define base names
DNA_BASENAME=$(basename "$DNA_FASTA" .fasta)
PROT_BASENAME=$(basename "$PROT_FASTA" .faa)

# Expected outputs
DNA_ALIGN="$OUTDIR/${DNA_BASENAME}.aln"
DNA_TREE="$OUTDIR/${DNA_BASENAME}.orig.nwk"
PROT_ALIGN="$OUTDIR/${PROT_BASENAME}_protein.aln"
PROT_TREE="$OUTDIR/${PROT_BASENAME}_protein.orig.nwk"

echo "[INFO] Checking for existing files in $OUTDIR"

# Step 1: Run DNA section if missing
if [[ ! -f "$DNA_ALIGN" || ! -f "$DNA_TREE" ]]; then
    echo "[INFO] Missing DNA alignment or tree — running DNA pipeline."
    cd ../
    ./metaProtease/tree.sh "$DNA_FASTA" dna
    cd -
else
    echo "[INFO] DNA alignment and tree already exist — skipping DNA pipeline."
fi

# Step 2: Run protein section if missing
if [[ ! -f "$PROT_ALIGN" || ! -f "$PROT_TREE" ]]; then
    echo "[INFO] Missing protein alignment or tree — running protein pipeline."
    cd ../
    ./metaProtease/tree.sh "$PROT_FASTA" protein
    cd -
else
    echo "[INFO] Protein alignment and tree already exist — skipping protein pipeline."
fi

echo "[DONE] All tree construction tasks complete."

