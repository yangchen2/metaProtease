#!/usr/bin/env python

import argparse
import os
import sys
import subprocess
import shutil
import logging


class DatabaseMapper:
    def __init__(self, db_name, index_path, outdir=None, threads=16, log_level=logging.INFO, score_min="L,0,-0.1", force=False):
        self.db_name = db_name.lower()
        self.index_path = os.path.abspath(index_path)
        self.threads = threads
        self.score_min = score_min
        self.force = force

        # === Output directory structure ===
        self.base_outdir = outdir or os.path.join("outputs", "database-map", self.db_name)
        self.log_dir = os.path.join(self.base_outdir, "logs")
        self.sams_dir = os.path.join(self.base_outdir, "sams")
        self.bams_dir = os.path.join(self.base_outdir, "bams")
        self.slurm_out_dir = os.path.join(self.base_outdir, "slurm_out")
        self.index_dir = os.path.dirname(self.index_path)

        # Create directories if not exist
        for d in [self.log_dir, self.sams_dir, self.bams_dir, self.slurm_out_dir]:
            os.makedirs(d, exist_ok=True)

        # === Logging setup ===
        log_file = os.path.join(self.log_dir, f"reference-map_{self.db_name}.log")
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
        self.logger.info(f"Database: {self.db_name}")
        self.logger.info(f"Index path: {self.index_path}")

    # --------------------------------------------------------
    # === Utility functions ===
    # --------------------------------------------------------

    @staticmethod
    def check_executable(program):
        """Check if a required program is in PATH."""
        if shutil.which(program) is None:
            logging.error(f"{program} not found in PATH.")
            sys.exit(1)

    @staticmethod
    def run_command(cmd, **kwargs):
        """Run a system command with error handling."""
        try:
            subprocess.run(cmd, check=True, **kwargs)
        except subprocess.CalledProcessError:
            logging.error("Command failed: " + (" ".join(cmd) if isinstance(cmd, list) else cmd))
            sys.exit(1)

    # --------------------------------------------------------
    # === Bowtie2 mapping function ===
    # --------------------------------------------------------

    def map_reads(self, r1, r2):
        """Map paired-end FASTQs to SAM using Bowtie2."""
        sample_name = os.path.basename(r1)
        if "_R1" in sample_name:
            sample_name = sample_name.split("_R1")[0]
        elif ".fastq" in sample_name:
            sample_name = sample_name.replace(".fastq.gz", "").replace(".fastq", "")

        sam_file = os.path.join(self.sams_dir, f"{sample_name}.sam")

        if os.path.exists(sam_file) and not self.force:
            self.logger.info(f"SAM file already exists: {sam_file} — skipping mapping.")
            return sam_file

        self.logger.info(f"Mapping sample {sample_name} → {sam_file}")

        cmd = [
            "bowtie2",
            "--no-exact-upfront",
            "--no-1mm-upfront",
            "--score-min", self.score_min,
            "-x", self.index_path,
            "-1", r1,
            "-2", r2,
            "-S", sam_file,
            "-p", str(self.threads),
            "--no-unal",
            "--very-sensitive",
            "--seed", "42",
            "-k", "16"
        ]

        self.run_command(cmd)
        self.logger.info(f"Finished mapping {sample_name}")
        return sam_file


# --------------------------------------------------------
# === Argument Parser ===
# --------------------------------------------------------

def parse_args():
    parser = argparse.ArgumentParser(
        description="Map reads to a specified database index (e.g., WoLr2, MetaPhlAn)."
    )
    parser.add_argument(
        "-d", "--db-name",
        required=True,
        help="Name of the database (e.g., wolr2, metaphlan)."
    )
    parser.add_argument(
        "-i", "--index-dir",
        required=True,
        help="Path to the Bowtie2 index prefix."
    )
    parser.add_argument(
        "-r1", "--read1",
        required=True,
        help="Path to R1 FASTQ file."
    )
    parser.add_argument(
        "-r2", "--read2",
        required=True,
        help="Path to R2 FASTQ file."
    )
    parser.add_argument(
        "-o", "--outdir",
        default=None,
        help="Base output directory (default: outputs/database-map/{db_name})."
    )
    parser.add_argument(
        "-t", "--threads",
        type=int,
        default=16,
        help="Number of threads for mapping (default: 16)."
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Overwrite existing SAM files if set."
    )
    return parser.parse_args()


# --------------------------------------------------------
# === Main ===
# --------------------------------------------------------

if __name__ == "__main__":
    args = parse_args()

    mapper = DatabaseMapper(
        db_name=args.db_name,
        index_path=args.index_path,
        outdir=args.outdir,
        threads=args.threads,
        force=args.force
    )

    mapper.check_executable("bowtie2")
    mapper.map_reads(args.read1, args.read2)

