#!/usr/bin/env python3
"""
Translate DNA FASTA (CDS) to protein FASTA.
Assumes bacterial genetic code (NCBI table 11).
Stops translation at the first stop codon.
"""

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os, sys

if len(sys.argv) < 2:
    sys.exit("Usage: python translate_to_protein.py <input.fasta>")

input_fasta = sys.argv[1]
output_faa = os.path.splitext(input_fasta)[0] + ".faa"

records = []
for record in SeqIO.parse(input_fasta, "fasta"):
    seq = record.seq.upper()
    protein_seq = seq.translate(table=11, to_stop=True)
    prot_record = SeqRecord(
        protein_seq,
        id=record.id,
        description=record.description
    )
    records.append(prot_record)

SeqIO.write(records, output_faa, "fasta")

print(f"[DONE] Translated {len(records)} sequences → {output_faa}")

