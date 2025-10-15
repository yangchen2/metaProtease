#!/bin/bash -l
#SBATCH --job-name=woltka_annotate_stream
#SBATCH --output=/home/yac027/metaProtease/outputs/db-map/slurm_out/woltka-annotate_%j.out
#SBATCH --error=/home/yac027/metaProtease/outputs/db-map/slurm_out/woltka-annotate_%j.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --mail-user=yac027@ucsd.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --partition=short

# Activate environment
source activate metaProtease

# Define paths
SPECIES_GENE_BIOM="/home/yac027/metaProtease/outputs/db-map/tables/gene_by_species.biom"
EGGNOG_ANNOT="/projects/wol/qiyun/wol2/function/eggnog/raw/emapper.annotations.xz"
OUT_TSV="/home/yac027/metaProtease/outputs/db-map/tables/gene_by_species_annotated.tsv.gz"

echo "[INFO] Starting lightweight annotation ..."
date

# Python annotation (streaming)
python <<'PYCODE'
import pandas as pd
from biom import load_table
from tqdm import tqdm

biom_path = "/home/yac027/metaProtease/outputs/db-map/tables/gene_by_species.biom"
anno_path = "/projects/wol/qiyun/wol2/function/eggnog/raw/emapper.annotations.xz"
out_tsv = "/home/yac027/metaProtease/outputs/db-map/tables/gene_by_species_annotated.tsv.gz"

# === Extract only feature IDs ===
print("[INFO] Extracting feature IDs from BIOM ...")
table = load_table(biom_path)
features = pd.Series(table.ids(axis="observation"), name="feature")

# Split into species and gene cluster IDs
split = features.str.extract(r'(?P<species>.+)\s+\|?\s*(?P<gene_cluster>G\d{9})')
split["gene_cluster"] = split["gene_cluster"].fillna(features.str.extract(r'(G\d{9})')[0])
split["feature"] = features
print(f"[INFO] Extracted {len(split)} feature IDs")

# === Stream EggNOG annotations and join ===
print("[INFO] Reading EggNOG annotations in chunks ...")
chunksize = 500000
annotated_list = []

for chunk in tqdm(pd.read_csv(
        anno_path, sep='\t', comment='#', header=None,
        usecols=[0, 8, 9],
        names=['gene_cluster_full', 'description', 'gene_name'],
        chunksize=chunksize)):

    chunk["gene_cluster"] = chunk["gene_cluster_full"].str.extract(r'(G\d{9})')
    merged = split.merge(chunk[["gene_cluster","gene_name","description"]],
                         on="gene_cluster", how="inner")
    if not merged.empty:
        annotated_list.append(merged[["feature","species","gene_cluster","gene_name","description"]])

if not annotated_list:
    raise RuntimeError("No matches found between BIOM and EggNOG annotations!")

result = pd.concat(annotated_list, ignore_index=True)

# === Write output (compressed) ===
print(f"[INFO] Writing annotated table to {out_tsv} ...")
result.to_csv(out_tsv, sep='\t', index=False, compression='gzip')
print(f"[DONE] Saved {len(result)} annotated entries to {out_tsv}")
PYCODE

date
echo "[ALL DONE] Annotation complete."
echo "Output: $OUT_TSV"

