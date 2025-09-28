#!/usr/bin/env python3

import argparse
import os
import sys
import logging
from Bio import SeqIO


def setup_logging(outdir, log_level=logging.INFO):
    """Set up logging to console + file."""
    os.makedirs(outdir, exist_ok=True)
    log_file = os.path.join(outdir, "fasta_to_gtf.log")

    logging.basicConfig(
        level=log_level,
        format="[%(levelname)s] %(message)s",
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler(sys.stdout),
        ],
    )
    logger = logging.getLogger(__name__)
    logger.info(f"Logging to {log_file}")
    return logger


def fasta_to_gtf(fasta, outdir, logger):
    """
    Convert a FASTA file into a simple GTF annotation.
    Each FASTA record becomes a 'gene' entry.
    """
    # Ensure output directory exists
    os.makedirs(outdir, exist_ok=True)

    # Build output filename
    fasta_base = os.path.splitext(os.path.basename(fasta))[0]
    gtf_file = os.path.join(outdir, f"{fasta_base}.gtf")

    with open(gtf_file, "w") as out:
        for record in SeqIO.parse(fasta, "fasta"):
            seqid = record.id.split()[0]  # first token of header
            length = len(record.seq)
            gtf_line = f"{seqid}\tcustom\tgene\t1\t{length}\t.\t+\t.\tgene_id \"{seqid}\";\n"
            out.write(gtf_line)

    logger.info(f"GTF written to {gtf_file}")
    return gtf_file


def main():
    parser = argparse.ArgumentParser(
        description="Convert a FASTA file to a GTF annotation file for featureCounts."
    )
    parser.add_argument("-i", "--fasta", required=True, help="Input FASTA file (protease reference)")
    parser.add_argument("-o", "--outdir", required=True, help="Output directory (GTF will be saved here)")
    parser.add_argument("-l", "--log-level", default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR"],
                        help="Logging level (default: INFO)")

    args = parser.parse_args()

    # Set up logger
    log_level = getattr(logging, args.log_level.upper(), logging.INFO)
    logger = setup_logging(args.outdir, log_level)

    logger.info(f"Input FASTA: {args.fasta}")
    logger.info(f"Output directory: {args.outdir}")

    fasta_to_gtf(args.fasta, args.outdir, logger)


if __name__ == "__main__":
    main()

