#!/usr/bin/env python3
import argparse
import sys
import os
import logging
from . import __version__
from .reference_map import ReferenceMapper


def run_reference_index(args):
    """Build Bowtie2 index for reference mapping."""
    mapper = ReferenceMapper(
        outdir=args.outdir,
        threads=args.threads,
        force=args.force,
        log_level=getattr(logging, args.log_level.upper(), logging.INFO),
    )
    mapper.check_executable("bowtie2-build")
    mapper.build_index(args.fasta)


def run_reference_map(args):
    """Run the full reference-based read mapping and assembly workflow."""
    mapper = ReferenceMapper(
        outdir=args.outdir,
        threads=args.threads,
        force=args.force,
        log_level=getattr(logging, args.log_level.upper(), logging.INFO),
        score_min=args.score_min,
    )

    mapper.check_executable("bowtie2")
    mapper.check_executable("samtools")
    mapper.check_executable("bedtools")
    mapper.check_executable("metaspades.py")

    index_prefix = args.index_path

    mapper.logger.info(f"=== Starting reference mapping workflow ===")
    mapper.logger.info(f"Index prefix: {index_prefix}")
    mapper.logger.info(f"FASTQs: {args.r1}, {args.r2}")

    # Step 2: Map reads to reference (Bowtie2)
    sam = mapper.map_reads(index_prefix, args.r1, args.r2)

    # Step 3: Convert SAM → BAM (sorted, mapped-only)
    bam_mapped = mapper.sam_to_bam(sam, index=args.index)

    # Step 4: Convert mapped BAM → paired FASTQs
    r1_mapped, r2_mapped = mapper.bam_to_fastq(bam_mapped)

    # Step 5: Assemble mapped reads using metaSPAdes
    assembly_dir = mapper.assembled_mapped_reads(r1_mapped, r2_mapped)

    # Step 6: Map assembled scaffolds back to reference
    scaffolds_file = os.path.join(assembly_dir, "scaffolds.fasta")
    contigs_file = os.path.join(assembly_dir, "contigs.fasta")

    mapper.logger.info(f"=== Step 6: Checking for scaffolds in {assembly_dir} ===")
    if os.path.exists(scaffolds_file):
        mapper.map_scaffolds2refs(index_prefix, scaffolds_file)
    elif os.path.exists(contigs_file):
        mapper.logger.warning(f"scaffolds.fasta not found — using contigs.fasta for {assembly_dir}")
        mapper.map_scaffolds2refs(index_prefix, contigs_file)
    else:
        mapper.logger.warning(f"No scaffolds.fasta or contigs.fasta found for {assembly_dir} — skipping remapping.")


def main():
    parser = argparse.ArgumentParser(
        description="metaProtease: Identify protease genes from metagenome data."
    )
    parser.add_argument("--version", action="version", version=f"metaProtease {__version__}")

    subparsers = parser.add_subparsers(dest="command")

    # === reference-map group ===
    ref_parser = subparsers.add_parser("reference-map", help="Run reference-based workflows")
    ref_subparsers = ref_parser.add_subparsers(dest="subcommand")

    # --- reference-map index ---
    ref_index_parser = ref_subparsers.add_parser("index", help="Build Bowtie2 index from FASTA")
    ref_index_parser.add_argument("--fasta", required=True, help="Path to input FASTA file for index building")
    ref_index_parser.add_argument("--outdir", default=None, help="Base output directory")
    ref_index_parser.add_argument("--threads", type=int, default=16, help="Number of threads for index building")
    ref_index_parser.add_argument("--force", action="store_true", help="Force rebuild even if index exists")
    ref_index_parser.add_argument(
        "--log-level", default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        help="Set logging verbosity"
    )
    ref_index_parser.set_defaults(func=run_reference_index)

    # --- reference-map map ---
    ref_map_parser = ref_subparsers.add_parser("map", help="Map reads to a reference database using Bowtie2")
    ref_map_parser.add_argument("--index-path", required=True, help="Path to Bowtie2 index prefix")
    ref_map_parser.add_argument("--r1", required=True, help="Path to R1 FASTQ file")
    ref_map_parser.add_argument("--r2", required=True, help="Path to R2 FASTQ file")
    ref_map_parser.add_argument("--outdir", default=None, help="Base output directory")
    ref_map_parser.add_argument("--threads", type=int, default=16, help="Number of threads for mapping")
    ref_map_parser.add_argument("--force", action="store_true", help="Overwrite existing SAM/BAM files if set")
    ref_map_parser.add_argument("--index", action="store_true", help="Also create BAM index (.bai)")
    ref_map_parser.add_argument(
        "--score-min", default="L,0,-0.6",
        help="Bowtie2 --score-min parameter (default: L,0,-0.6)"
    )
    ref_map_parser.add_argument(
        "--log-level", default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        help="Set logging verbosity"
    )
    ref_map_parser.set_defaults(func=run_reference_map)

    # === Parse and dispatch ===
    args = parser.parse_args()

    if hasattr(args, "func"):
        args.func(args)
    else:
        parser.print_help()


if __name__ == "__main__":
    main()

