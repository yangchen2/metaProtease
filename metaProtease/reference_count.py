#!/usr/bin/env python3

import os
import sys
import argparse
import shutil
import logging
import subprocess


class ReferenceCount:
    def __init__(self, outdir="outputs/reference-map", threads=4, log_level=logging.INFO):
        self.threads = threads

        # Output directories
        self.base_outdir = outdir
        self.sam_dir = os.path.join(self.base_outdir, "sams")
        self.bam_dir = os.path.join(self.base_outdir, "bams")
        self.log_dir = os.path.join(self.base_outdir, "logs")

        os.makedirs(self.sam_dir, exist_ok=True)
        os.makedirs(self.bam_dir, exist_ok=True)
        os.makedirs(self.log_dir, exist_ok=True)

        # Logging
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

    def filter_sam(self, sam):
        """Filter SAM file to only properly paired primary alignments."""
        base = os.path.basename(sam)
        if base.endswith("_filt.sam"):
            sam_name = base
        else:
            sam_name = base.replace(".sam", "_filt.sam")

        out_sam = os.path.join(self.sam_dir, sam_name)
        cmd = ["samtools", "view", "-h", "-f", "0x2", "-F", "0x100", "-o", out_sam, sam]
        self.run_command(cmd)

        self.logger.info(f"Filtered SAM saved to {out_sam}")
        return out_sam

    def build_bam(self, bam):
        """Filter BAM file to only properly paired primary alignments and return sorted BAM."""
        base = os.path.basename(bam)
        if base.endswith("_filt.bam"):
            bam_name = base
        else:
            bam_name = base.replace(".bam", "_filt.bam")

        out_bam = os.path.join(self.bam_dir, bam_name)

        tmp_filt = out_bam + ".tmp.filter"
        cmd = ["samtools", "view", "-b", "-f", "0x2", "-F", "0x100", "-o", tmp_filt, bam]
        self.run_command(cmd)

        tmp_sort = out_bam + ".tmp.sort"
        self.run_command(["samtools", "sort", "-o", tmp_sort, tmp_filt])
        os.replace(tmp_sort, out_bam)

        if os.path.exists(tmp_filt):
            os.remove(tmp_filt)

        self.logger.info(f"Filtered + sorted BAM saved to {out_bam}")
        return out_bam


def main():
    parser = argparse.ArgumentParser(
        description="Filter alignments to properly paired primaries (SLURM array-compatible)."
    )
    parser.add_argument("-i", "--input", required=True, help="Input SAM or BAM file")
    parser.add_argument("-o", "--outdir", default="outputs/reference-map", help="Base output directory")
    parser.add_argument("-t", "--threads", type=int, default=4, help="Number of threads")

    args = parser.parse_args()

    rc = ReferenceCount(outdir=args.outdir, threads=args.threads)

    # Step 1: filter input
    if args.input.endswith(".sam"):
        filt = rc.filter_sam(args.input)
        bam = filt.replace(".sam", ".bam")
        rc.run_command(["samtools", "view", "-bS", filt, "-o", bam])
        bam = rc.build_bam(bam)
    elif args.input.endswith(".bam"):
        bam = rc.build_bam(args.input)
    else:
        rc.logger.error("Input must be a SAM or BAM file")
        sys.exit(1)


if __name__ == "__main__":
    main()

