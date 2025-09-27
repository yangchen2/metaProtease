#!/usr/bin/env python3

import argparse
import os
import sys
import subprocess
import shutil
import logging


class ReferenceMapper:
    def __init__(self, outdir=None, threads=4, log_level=logging.INFO, score_min="L,0,-0.1"):
        self.threads = threads
        self.score_min = score_min

        # Default output structure
        self.base_outdir = outdir or os.path.join("outputs", "reference-map")
        self.index_dir = os.path.join(self.base_outdir, "index")
        self.log_dir = os.path.join(self.base_outdir, "logs")
        self.sams_dir = os.path.join(self.base_outdir, "sams")
        self.slurm_out_dir = os.path.join(self.base_outdir, "slurm_out")

        os.makedirs(self.index_dir, exist_ok=True)
        os.makedirs(self.log_dir, exist_ok=True)
        os.makedirs(self.sams_dir, exist_ok=True)
        os.makedirs(self.slurm_out_dir, exist_ok=True)

        # Configure logging to both console and file (no datetime)
        log_file = os.path.join(self.log_dir, "reference-map.log")
        logging.basicConfig(
            level=log_level,
            format="[%(levelname)s] %(message)s",
            handlers=[
                logging.FileHandler(log_file),
                logging.StreamHandler(sys.stdout),
            ],
        )
        self.logger = logging.getLogger(__name__)
        self.logger.info(f"Logging to {log_file}")

    @staticmethod
    def check_executable(program):
        """Check if a required program is in PATH."""
        if shutil.which(program) is None:
            logging.error(f"{program} not found in PATH.")
            sys.exit(1)

    @staticmethod
    def run_command(cmd):
        """Run a system command with error handling."""
        try:
            subprocess.run(cmd, check=True)
        except subprocess.CalledProcessError:
            logging.error(f"Command failed: {' '.join(cmd)}")
            sys.exit(1)

    def build_index(self, fasta):
        """Build Bowtie2 index from FASTA file with fixed prefix 'reference'."""
        index_prefix = os.path.join(self.index_dir, "reference")

        self.logger.info(f"Building Bowtie2 index: {index_prefix}")
        cmd = ["bowtie2-build", fasta, index_prefix]
        self.run_command(cmd)

        self.logger.info(f"Bowtie2 index created at {index_prefix}.*")
        return index_prefix

    def map_reads(self, index_prefix, r1, r2):
        """Map paired-end FASTQs to SAM using Bowtie2."""
        # Auto-derive sample name from R1 filename
        sample_name = os.path.basename(r1)
        if "_R1" in sample_name:
            sample_name = sample_name.split("_R1")[0]
        elif ".fastq" in sample_name:
            sample_name = sample_name.replace(".fastq.gz", "").replace(".fastq", "")

        sam_file = os.path.join(self.sams_dir, f"{sample_name}.sam")

        cmd = [
            "bowtie2",
            "--score-min", self.score_min,
            "-x", index_prefix,
            "-1", r1,
            "-2", r2,
            "-S", sam_file,
            "-p", str(self.threads)
        ]

        self.logger.info(f"Mapping sample {sample_name} → {sam_file}")
        self.run_command(cmd)

        self.logger.info(f"SAM file created: {sam_file}")
        return sam_file


def main():
    parser = argparse.ArgumentParser(
        description="Reference mapping pipeline (class-based)."
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    # Subcommand: index
    index_parser = subparsers.add_parser("index", help="Build Bowtie2 index")
    index_parser.add_argument("--fasta", "-f", required=True, help="Input multi-FASTA file")
    index_parser.add_argument(
        "--outdir",
        "-o",
        help="Base output directory (default: outputs/reference-map/)",
    )
    index_parser.add_argument("--threads", "-t", type=int, default=4, help="Number of threads")
    index_parser.add_argument(
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        help="Set logging level (default: INFO)",
    )

    # Subcommand: map
    map_parser = subparsers.add_parser("map", help="Map paired-end FASTQs to SAM")
    map_parser.add_argument("--r1", required=True, help="FASTQ R1 file")
    map_parser.add_argument("--r2", required=True, help="FASTQ R2 file")
    map_parser.add_argument(
        "--index-dir", "-x", required=True, help="Directory containing Bowtie2 index"
    )
    map_parser.add_argument(
        "--outdir",
        "-o",
        help="Base output directory (default: outputs/reference-map/)",
    )
    map_parser.add_argument("--threads", "-t", type=int, default=4, help="Number of threads")
    map_parser.add_argument(
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        help="Set logging level (default: INFO)",
    )
    map_parser.add_argument(
        "--score-min",
        default="L,0,-0.1",
        help="Bowtie2 --score-min parameter (default: L,0,-0.1 for ~90% identity)"
    )

    args = parser.parse_args()

    # Convert log level string → logging constant
    log_level = getattr(logging, args.log_level.upper(), logging.INFO)

    if args.command == "index":
        ReferenceMapper.check_executable("bowtie2-build")
        mapper = ReferenceMapper(outdir=args.outdir, threads=args.threads, log_level=log_level)
        mapper.build_index(args.fasta)

    elif args.command == "map":
        ReferenceMapper.check_executable("bowtie2")
        index_prefix = os.path.join(args.index_dir, "reference")
        mapper = ReferenceMapper(
            outdir=args.outdir,
            threads=args.threads,
            log_level=log_level,
            score_min=args.score_min
        )
        mapper.map_reads(index_prefix, args.r1, args.r2)


if __name__ == "__main__":
    main()

