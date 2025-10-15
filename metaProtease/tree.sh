#!/bin/bash -l

set -euo pipefail

# === Input validation ===
if [[ $# -lt 2 ]]; then
    echo "Usage: $0 <reference_multi_fasta> <mode: dna|protein>"
    exit 1
fi

FASTA=$1
MODE=$2  # dna or protein

# === Set up paths ===
BASENAME=$(basename "$FASTA")
BASENAME="${BASENAME%.*}"

OUTDIR="outputs/reference-map/tree"
mkdir -p "$OUTDIR"

# Adjust output names based on mode
if [[ "$MODE" == "protein" ]]; then
    ALIGN="$OUTDIR/${BASENAME}_protein.aln"
    CLEAN_ALIGN="$OUTDIR/${BASENAME}_protein.clean.aln"
    MAPFILE="$OUTDIR/${BASENAME}_protein.header_map.tsv"
    TREE_SHORT="$OUTDIR/${BASENAME}_protein.nwk"
    TREE_FULL="$OUTDIR/${BASENAME}_protein.orig.nwk"
else
    ALIGN="$OUTDIR/${BASENAME}.aln"
    CLEAN_ALIGN="$OUTDIR/${BASENAME}.clean.aln"
    MAPFILE="$OUTDIR/${BASENAME}.header_map.tsv"
    TREE_SHORT="$OUTDIR/${BASENAME}.nwk"
    TREE_FULL="$OUTDIR/${BASENAME}.orig.nwk"
fi

echo "[INFO] Running $MODE tree building pipeline"
echo "[INFO] Input FASTA: $FASTA"
echo "[INFO] Output directory: $OUTDIR"

# === Step 1: MAFFT alignment ===
if [[ -f "$ALIGN" ]]; then
    echo "[INFO] Alignment already exists: $ALIGN — skipping."
else
    echo "[INFO] Aligning sequences with MAFFT..."
    mafft --auto "$FASTA" > "$ALIGN"
    echo "[INFO] Alignment complete: $ALIGN"
fi

# === Step 2: Shorten headers & create mapping ===
if [[ -f "$CLEAN_ALIGN" && -f "$MAPFILE" ]]; then
    echo "[INFO] Shortened alignment and mapping already exist — skipping."
else
    echo "[INFO] Creating compact headers and mapping file..."
    awk '/^>/{i++; short=sprintf("seq%03d", i); print ">" short; print short "\t" substr($0,2) >> "'"$MAPFILE"'" ; next} {print}' \
        "$ALIGN" > "$CLEAN_ALIGN"
    echo "[INFO] Shortened alignment written to $CLEAN_ALIGN"
    echo "[INFO] Mapping file written to $MAPFILE"
fi

# === Step 3: Build tree with FastTree ===
if [[ -f "$TREE_SHORT" ]]; then
    echo "[INFO] Tree already exists: $TREE_SHORT — skipping."
else
    echo "[INFO] Building tree with FastTree ($MODE mode)..."
    if [[ "$MODE" == "protein" ]]; then
        FastTree -wag "$CLEAN_ALIGN" > "$TREE_SHORT"
    else
        FastTree -nt -gtr -gamma "$CLEAN_ALIGN" > "$TREE_SHORT"
    fi
    echo "[INFO] Tree complete: $TREE_SHORT"
fi

# === Step 4: Map short IDs back to full headers ===
if [[ -f "$TREE_FULL" ]]; then
    echo "[INFO] Full-header tree already exists: $TREE_FULL — skipping."
else
    echo "[INFO] Replacing short IDs with original headers..."
    cp "$TREE_SHORT" "$TREE_FULL"
    while IFS=$'\t' read -r short full; do
        sed -i "s/\b${short}\b/${full}/g" "$TREE_FULL"
    done < "$MAPFILE"
    echo "[INFO] Tree with original headers written to $TREE_FULL"
fi

echo "[DONE] $MODE tree pipeline completed successfully."

