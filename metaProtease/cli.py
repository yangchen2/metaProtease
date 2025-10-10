#!/usr/bin/env python3
import argparse
import sys
from . import __version__
from .database_map import DatabaseMapper


def run_reference_map(args):
    """Run the reference-based read mapping workflow."""
    mapper = DatabaseMapper(
        db_name=args.db_name,
        index_path=args.index_path,
        outdir=args.outdir,
        threads=args.threads,
        force=args.force
    )
    mapper.check_executable("bowtie2")
    mapper.map_reads(args.r1, args.r2)


def main():
    parser = argparse.ArgumentParser(
        description="metaProtease: A pipeline for the identification of protease genes from metagenome data"
    )
    parser.add_argument("--version", action="version", version=f"metaProtease {__version__}")

    subparsers = parser.add_subparsers(dest="command")

    # === reference-map group ===
    ref_parser = subparsers.add_parser("reference-map", help="Run reference-based workflows")
    ref_subparsers = ref_parser.add_subparsers(dest="subcommand")

    # --- reference-map map ---
    ref_map_parser = ref_subparsers.add_parser("map", help="Map reads to a reference database using Bowtie2")
    ref_map_parser.add_argument("--db-name", required=True, help="Database name (e.g., wolr2)")
    ref_map_parser.add_argument("--index-path", required=True, help="Path to Bowtie2 index prefix")
    ref_map_parser.add_argument("--r1", required=True, help="Path to R1 FASTQ file")
    ref_map_parser.add_argument("--r2", required=True, help="Path to R2 FASTQ file")
    ref_map_parser.add_argument("--outdir", default=None, help="Base output directory (default: outputs/database-map/{db_name})")
    ref_map_parser.add_argument("--threads", type=int, default=16, help="Number of threads for mapping (default: 16)")
    ref_map_parser.add_argument("--force", action="store_true", help="Overwrite existing SAM files if set")
    ref_map_parser.set_defaults(func=run_reference_map)

    # === db-map group ===
    db_parser = subparsers.add_parser("db-map", help="Run database-based workflows")
    db_subparsers = db_parser.add_subparsers(dest="subcommand")

    # --- db-map map ---
    db_map_parser = db_subparsers.add_parser("map", help="Map reads to a reference database using Bowtie2")
    db_map_parser.add_argument("--db-name", required=True, help="Database name (e.g., wolr2)")
    db_map_parser.add_argument("--index-path", required=True, help="Path to Bowtie2 index prefix")
    db_map_parser.add_argument("--r1", required=True, help="Path to R1 FASTQ file")
    db_map_parser.add_argument("--r2", required=True, help="Path to R2 FASTQ file")
    db_map_parser.add_argument("--outdir", default=None, help="Base output directory (default: outputs/database-map/{db_name})")
    db_map_parser.add_argument("--threads", type=int, default=16, help="Number of threads for mapping (default: 16)")
    db_map_parser.add_argument("--force", action="store_true", help="Overwrite existing SAM files if set")
    db_map_parser.set_defaults(func=run_reference_map)

    # === Parse and dispatch ===
    args = parser.parse_args()

    if hasattr(args, "func"):
        args.func(args)
    else:
        parser.print_help()


if __name__ == "__main__":
    main()

