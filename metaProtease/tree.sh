#!/bin/bash -l

set -euo pipefail

echo 'Checking for FASTA file input.'
if [[ $# -lt 1 ]]; then
    echo "Usage: $0 <reference_multi_fasta>"
    exit 1
fi

# Setting variables and directories
FASTA=$1
BASENAME=$(basename "$FASTA" .fasta)
BASENAME=$(basename "$BASENAME" .fa)
OUTDIR="outputs/reference-map/tree"
ALIGN="$OUTDIR/${BASENAME}.aln"
CLEAN_ALIGN="$OUTDIR/${BASENAME}.clean.aln"
MAPFILE="$OUTDIR/${BASENAME}.header_map.tsv"
TREE_SHORT="$OUTDIR/${BASENAME}.nwk"
TREE_FULL="$OUTDIR/${BASENAME}.orig.nwk"

mkdir -p "$OUTDIR"

# Step 1: Alignment
if [[ -f "$ALIGN" ]]; then
    echo "[INFO] Alignment already exists: $ALIGN — skipping."
else
    echo "[INFO] Aligning $FASTA → $ALIGN"
    mafft --auto "$FASTA" > "$ALIGN"
    echo "[INFO] Alignment complete. Output: $ALIGN"
fi

# Step 2: Shorten headers & create mapping
if [[ -f "$CLEAN_ALIGN" && -f "$MAPFILE" ]]; then
    echo "[INFO] Shortened alignment and mapping already exist — skipping."
else
    echo "[INFO] Creating compact headers and mapping file"
    awk '/^>/{i++; short=sprintf("seq%03d", i); print ">" short; print short "\t" substr($0,2) >> "'"$MAPFILE"'" ; next} {print}' \
        "$ALIGN" > "$CLEAN_ALIGN"
    echo "[INFO] Shortened alignment written to $CLEAN_ALIGN"
    echo "[INFO] Mapping file written to $MAPFILE"
fi

# Step 3: Tree building with FastTree
if [[ -f "$TREE_SHORT" ]]; then
    echo "[INFO] Tree already exists: $TREE_SHORT — skipping."
else
    echo "[INFO] Making tree from cleaned alignment"
    FastTree -nt -gtr -gamma "$CLEAN_ALIGN" > "$TREE_SHORT"
    echo "[INFO] Tree complete. Output: $TREE_SHORT"
fi

# Step 4: Map short IDs back to original headers
if [[ -f "$TREE_FULL" ]]; then
    echo "[INFO] Tree with full headers already exists: $TREE_FULL — skipping."
else
    echo "[INFO] Replacing short IDs with original headers"
    cp "$TREE_SHORT" "$TREE_FULL"
    while IFS=$'\t' read -r short full; do
        sed -i "s/\b${short}\b/${full}/g" "$TREE_FULL"
    done < "$MAPFILE"
    echo "[INFO] Tree with original headers written to $TREE_FULL"
fi

