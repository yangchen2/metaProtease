#!/usr/bin/env python3

import argparse
import os
import sys
import subprocess
import shutil
import logging


class ReferenceMapper:
    def __init__(self, outdir=None, threads=4, log_level=logging.INFO, score_min="L,0,-0.1", force=False):
        self.threads = threads
        self.score_min = score_min
        self.force = force  # NEW flag

        # Default output structure
        self.base_outdir = outdir or os.path.join("outputs", "reference-map")
        self.index_dir = os.path.join(self.base_outdir, "index")
        self.log_dir = os.path.join(self.base_outdir, "logs")
        self.sams_dir = os.path.join(self.base_outdir, "sams")
        self.bams_dir = os.path.join(self.base_outdir, "bams")
        self.slurm_out_dir = os.path.join(self.base_outdir, "slurm_out")

        os.makedirs(self.index_dir, exist_ok=True)
        os.makedirs(self.log_dir, exist_ok=True)
        os.makedirs(self.sams_dir, exist_ok=True)
        os.makedirs(self.bams_dir, exist_ok=True)
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
    def run_command(cmd, **kwargs):
        """Run a system command with error handling."""
        try:
            subprocess.run(cmd, check=True, **kwargs)
        except subprocess.CalledProcessError:
            logging.error("Command failed: " + (" ".join(cmd) if isinstance(cmd, list) else cmd))
            sys.exit(1)

    def build_index(self, fasta):
        """Build Bowtie2 index from FASTA file with fixed prefix 'reference'."""
        index_prefix = os.path.join(self.index_dir, "reference")
        expected_files = [f"{index_prefix}.{ext}" for ext in
                          ["1.bt2", "2.bt2", "3.bt2", "4.bt2", "rev.1.bt2", "rev.2.bt2"]]

        if all(os.path.exists(f) for f in expected_files) and not self.force:
            self.logger.info(f"Bowtie2 index already exists at {index_prefix}.* — skipping.")
            return index_prefix

        self.logger.info(f"Building Bowtie2 index: {index_prefix}")
        cmd = ["bowtie2-build", fasta, index_prefix]
        self.run_command(cmd)

        self.logger.info(f"Bowtie2 index created at {index_prefix}.*")
        return index_prefix

    def map_reads(self, index_prefix, r1, r2):
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
            "--score-min", self.score_min,
            "-x", index_prefix,
            "-1", r1,
            "-2", r2,
            "-S", sam_file,
            "-p", str(self.threads)
        ]
        self.run_command(cmd)

        self.logger.info(f"SAM file created: {sam_file}")
        return sam_file

    def sam_to_sorted_bam(self, sam_file, index=False):
        """Convert SAM file to a sorted BAM file using Samtools."""
        sample_name = os.path.splitext(os.path.basename(sam_file))[0]
        bam_file = os.path.join(self.bams_dir, f"{sample_name}.bam")
        bai_file = bam_file + ".bai"

        if os.path.exists(bam_file) and not self.force:
            self.logger.info(f"Sorted BAM already exists: {bam_file} — skipping conversion.")
            if index and os.path.exists(bai_file) and not self.force:
                self.logger.info(f"BAM index already exists: {bai_file} — skipping indexing.")
            elif index:
                self.logger.info(f"Indexing BAM file: {bam_file}")
                cmd_index = ["samtools", "index", bam_file]
                self.run_command(cmd_index)
            return bam_file

        self.logger.info(f"Converting {sam_file} → {bam_file} (sorted BAM)")
        cmd = f"samtools view -@ {self.threads} -bS {sam_file} | samtools sort -@ {self.threads} -o {bam_file}"
        self.run_command(cmd, shell=True)
        self.logger.info(f"Sorted BAM file created: {bam_file}")

        if index:
            self.logger.info(f"Indexing BAM file: {bam_file}")
            cmd_index = ["samtools", "index", bam_file]
            self.run_command(cmd_index)

        return bam_file


def main():
    parser = argparse.ArgumentParser(
        description="Reference mapping pipeline (class-based)."
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    # Subcommand: index
    index_parser = subparsers.add_parser("index", help="Build Bowtie2 index")
    index_parser.add_argument("--fasta", "-f", required=True, help="Input multi-FASTA file")
    index_parser.add_argument("--outdir", "-o", help="Base output directory")
    index_parser.add_argument("--threads", "-t", type=int, default=4, help="Number of threads")
    index_parser.add_argument(
        "--log-level", default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        help="Set logging level"
    )
    index_parser.add_argument(
        "--force", action="store_true",
        help="Force rebuild index even if files exist"
    )

    # Subcommand: map
    map_parser = subparsers.add_parser("map", help="Map paired-end FASTQs to SAM and convert to BAM")
    map_parser.add_argument("--r1", required=True, help="FASTQ R1 file")
    map_parser.add_argument("--r2", required=True, help="FASTQ R2 file")
    map_parser.add_argument("--index-dir", "-x", required=True, help="Directory containing Bowtie2 index")
    map_parser.add_argument("--outdir", "-o", help="Base output directory")
    map_parser.add_argument("--threads", "-t", type=int, default=4, help="Number of threads")
    map_parser.add_argument(
        "--log-level", default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        help="Set logging level"
    )
    map_parser.add_argument(
        "--score-min", default="L,0,-0.1",
        help="Bowtie2 --score-min parameter"
    )
    map_parser.add_argument(
        "--index", action="store_true",
        help="Also create BAM index (.bai). Default = off"
    )
    map_parser.add_argument(
        "--force", action="store_true",
        help="Force re-run mapping and BAM conversion even if files exist"
    )

    args = parser.parse_args()
    log_level = getattr(logging, args.log_level.upper(), logging.INFO)

    if args.command == "index":
        ReferenceMapper.check_executable("bowtie2-build")
        mapper = ReferenceMapper(
            outdir=args.outdir, threads=args.threads, log_level=log_level, force=args.force
        )
        mapper.build_index(args.fasta)

    elif args.command == "map":
        ReferenceMapper.check_executable("bowtie2")
        ReferenceMapper.check_executable("samtools")
        index_prefix = os.path.join(args.index_dir, "reference")
        mapper = ReferenceMapper(
            outdir=args.outdir, threads=args.threads, log_level=log_level,
            score_min=args.score_min, force=args.force
        )
        sam = mapper.map_reads(index_prefix, args.r1, args.r2)
        mapper.sam_to_sorted_bam(sam, index=args.index)


if __name__ == "__main__":
    main()

