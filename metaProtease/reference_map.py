#!/usr/bin/env python3
"""
===========================================================
ReferenceMapper Pipeline
Author: Yang Chen (yac027@ucsd.edu)
Date: 2025-10-17

Description:
------------
This class-based script automates reference-based mapping,
assembly, and remapping for identifying protease gene fragments
in metagenomic reads.

Steps:
--------------
1. Build Bowtie2 index of protease reference genes.
2. Map paired-end metagenomic reads to the reference (Bowtie2 → SAM).
3. Convert SAM to sorted BAM (samtools) and filter to mapped reads only.
4. Convert mapped reads back to paired gzipped FASTQs (bedtools bamtofastq).
5. Assemble mapped reads using metaSPAdes (per-sample assemblies).
6. Map resulting contigs/scaffolds back to the reference to evaluate
   coverage and confirm partial gene reconstruction.

Directories created automatically:
---------------------------------
outputs/reference-map/
├── index/                   (Bowtie2 index files)
├── sams/                    (SAM alignment files)
│   └── scaffold_mappings/   (SAMs from scaffolds remapping)
├── bams/                    (Sorted BAMs)
│   └── bams_mapped/         (Mapped-only BAMs)
├── fastq_mapped/            (FASTQs extracted from mapped reads)
├── metaspades/              (metaSPAdes assemblies)
├── logs/                    (Log files)
└── slurm_out/               (Cluster job logs)

===========================================================
"""

import argparse
import os
import sys
import subprocess
import shutil
import logging
import re


class ReferenceMapper:
    def __init__(self, outdir=None, threads=4, log_level=logging.INFO, score_min="L,0,-0.6", force=False):
        self.threads = threads
        self.score_min = score_min
        self.force = force

        # === Output directory structure ===
        self.base_outdir = outdir or os.path.join("outputs", "reference-map")
        self.index_dir = os.path.join(self.base_outdir, "index")
        self.log_dir = os.path.join(self.base_outdir, "logs")
        self.sams_dir = os.path.join(self.base_outdir, "sams")
        self.bams_dir = os.path.join(self.base_outdir, "bams")
        self.bams_mapped_dir = os.path.join(self.bams_dir, "bams_mapped")
        self.fastqs_mapped_dir = os.path.join(self.base_outdir, "fastq_mapped")
        self.metaspades_dir = os.path.join(self.base_outdir, "metaspades")
        self.mapping_output_dir = os.path.join(self.sams_dir, "scaffold_mappings")
        self.slurm_out_dir = os.path.join(self.base_outdir, "slurm_out")

        # === Create directories if missing ===
        for d in [
            self.index_dir, self.log_dir, self.sams_dir, self.bams_dir,
            self.bams_mapped_dir, self.fastqs_mapped_dir, self.metaspades_dir,
            self.mapping_output_dir, self.slurm_out_dir
        ]:
            os.makedirs(d, exist_ok=True)

        # === Configure logging ===
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

    # ----------------------------------------------------
    # Utility methods
    # ----------------------------------------------------
    @staticmethod
    def check_executable(program):
        """Check if a required program is in PATH."""
        if shutil.which(program) is None:
            logging.error(f"{program} not found in PATH.")
            sys.exit(1)

    @staticmethod
    def run_command(cmd, **kwargs):
        """Run a system command with stderr capture."""
        result = subprocess.run(cmd, text=True, capture_output=True, **kwargs)
        if result.returncode != 0:
            logging.error("Command failed: " + (" ".join(cmd) if isinstance(cmd, list) else cmd))
            logging.error(f"STDERR:\n{result.stderr.strip()}")
            sys.exit(1)

    # ----------------------------------------------------
    # Step 1: Build Bowtie2 index
    # ----------------------------------------------------
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

    # ----------------------------------------------------
    # Step 2: Map reads → SAM
    # ----------------------------------------------------
    def map_reads(self, index_prefix, r1, r2):
        """Map paired-end FASTQs to SAM using Bowtie2."""
        sample_name = os.path.basename(r1).split("_R1")[0]
        sam_file = os.path.join(self.sams_dir, f"{sample_name}.sam")

        self.logger.info(f"=== Processing sample {sample_name} ===")

        if os.path.exists(sam_file) and not self.force:
            self.logger.info(f"SAM file already exists: {sam_file} — skipping mapping.")
            return sam_file

        self.logger.info(f"Mapping sample {sample_name} → {sam_file}")
        cmd = [
            "bowtie2",
            "--no-exact-upfront", "--no-1mm-upfront",
            "--score-min", self.score_min,
            "-x", index_prefix,
            "-1", r1, "-2", r2,
            "-S", sam_file,
            "-p", str(self.threads),
            "--no-unal",
            "--very-sensitive",
            "--seed", "42",
            "-k", "16"
        ]
        self.run_command(cmd)
        self.logger.info(f"SAM file created: {sam_file}")
        return sam_file

    # ----------------------------------------------------
    # Step 3: SAM → BAM → mapped-only BAM
    # ----------------------------------------------------
    def sam_to_bam(self, sam_file, index=False):
        """Convert SAM file to a sorted BAM file using Samtools."""
        sample_name = os.path.splitext(os.path.basename(sam_file))[0]
        bam_file = os.path.join(self.bams_dir, f"{sample_name}.bam")
        bai_file = bam_file + ".bai"
        bam_file_mapped = os.path.join(self.bams_mapped_dir, f"{sample_name}_mapped.bam")

        # Case 1: Both exist
        if os.path.exists(bam_file) and os.path.exists(bam_file_mapped) and not self.force:
            self.logger.info(f"Both sorted BAM and mapped BAM already exist for {sample_name}. Skipping all steps.")
            if index and not os.path.exists(bai_file):
                self.run_command(["samtools", "index", bam_file])
            return bam_file_mapped

        # Convert SAM → BAM → mapped-only BAM
        self.logger.info(f"Converting {sam_file} → {bam_file} (sorted BAM)")
        cmd_convert = ["samtools", "sort", "-@", str(self.threads), "-o", bam_file, sam_file]
        self.run_command(cmd_convert)

        self.logger.info(f"Filtering mapped reads → {bam_file_mapped}")
        cmd_mapped = ["samtools", "view", "-@", str(self.threads), "-b", "-F", "4", "-o", bam_file_mapped, bam_file]
        self.run_command(cmd_mapped)

        if os.path.getsize(bam_file_mapped) == 0:
            self.logger.warning(f"No mapped reads found in {sample_name} — mapped BAM is empty.")

        if index:
            self.run_command(["samtools", "index", bam_file])

        return bam_file_mapped

    # ----------------------------------------------------
    # Step 4: BAM → gzipped FASTQ (mapped reads only)
    # ----------------------------------------------------
    def bam_to_fastq(self, bam_file_mapped):
        """Convert mapped reads back to gzipped FASTQs for metaSPAdes assembly."""
        sample_name = os.path.splitext(os.path.basename(bam_file_mapped))[0].replace("_mapped", "")
        mapped_R1_file = os.path.join(self.fastqs_mapped_dir, f"{sample_name}_mapped_R1.fastq.gz")
        mapped_R2_file = os.path.join(self.fastqs_mapped_dir, f"{sample_name}_mapped_R2.fastq.gz")

        if not os.path.exists(bam_file_mapped):
            raise FileNotFoundError(f"BAM file not found: {bam_file_mapped}")

        if os.path.exists(mapped_R1_file) and os.path.exists(mapped_R2_file) and not self.force:
            self.logger.info(f"Gzipped FASTQ files already exist for {sample_name} — skipping conversion.")
            return mapped_R1_file, mapped_R2_file

        cmd = (
            f"bedtools bamtofastq -i {bam_file_mapped} "
            f"-fq {mapped_R1_file.replace('.gz','')} "
            f"-fq2 {mapped_R2_file.replace('.gz','')} && "
            f"gzip -f {mapped_R1_file.replace('.gz','')} {mapped_R2_file.replace('.gz','')}"
        )

        self.logger.info(f"Converting BAM → gzipped FASTQ for {sample_name} ...")
        self.run_command(cmd, shell=True)
        return mapped_R1_file, mapped_R2_file

    # ----------------------------------------------------
    # Step 5: Assemble mapped reads (metaSPAdes)
    # ----------------------------------------------------
    def assembled_mapped_reads(self, mapped_R1_file, mapped_R2_file):
        """Assemble mapped reads with metaSPAdes."""
        sample_name = os.path.basename(mapped_R1_file)
        sample_name = re.sub(r"_mapped_R1(?:\.fastq(?:\.gz)?)?$", "", sample_name)
        metaspades_persample_outdir = os.path.join(self.metaspades_dir, sample_name)

        if not (os.path.exists(mapped_R1_file) and os.path.exists(mapped_R2_file)):
            raise FileNotFoundError(f"Missing pre-mapped FASTQs for {sample_name}")

        if os.path.exists(metaspades_persample_outdir) and not self.force:
            self.logger.info(f"metaSPAdes output for {sample_name} already exists — skipping assembly.")
            return metaspades_persample_outdir

        cmd = (
            f"metaspades.py -1 {mapped_R1_file} -2 {mapped_R2_file} "
            f"-o {metaspades_persample_outdir} --only-assembler -t {self.threads}"
        )
        self.logger.info(f"Running metaSPAdes for {sample_name} ...")
        self.run_command(cmd, shell=True)
        return metaspades_persample_outdir

    # ----------------------------------------------------
    # Step 6: Map assembled scaffolds back to reference
    # ----------------------------------------------------
    def map_scaffolds2refs(self, index_prefix, scaffolds_file):
        """Map assembled scaffolds back to Bowtie2 reference."""
        sample_name = os.path.basename(os.path.abspath(os.path.dirname(scaffolds_file)))
        contigs_file = os.path.join(os.path.dirname(scaffolds_file), "contigs.fasta")
        sam_out = os.path.join(self.mapping_output_dir, f"{sample_name}.sam")

        # Convert to absolute paths
        scaffolds_file = os.path.abspath(scaffolds_file)
        sam_out = os.path.abspath(sam_out)
        index_prefix = os.path.abspath(index_prefix)

        if not os.path.exists(scaffolds_file):
            if os.path.exists(contigs_file):
                self.logger.warning(f"No scaffolds.fasta found; using contigs.fasta instead for {sample_name}.")
                scaffolds_file = os.path.abspath(contigs_file)
            else:
                raise FileNotFoundError(f"No scaffolds or contigs file found for {sample_name}")

        if os.path.getsize(scaffolds_file) == 0:
            self.logger.warning(f"{scaffolds_file} is empty — skipping mapping.")
            return None

        if os.path.exists(sam_out) and not self.force:
            self.logger.info(f"SAM already exists for {sample_name} — skipping remapping.")
            return sam_out

        cmd = [
            "bowtie2",
            "--no-exact-upfront", "--no-1mm-upfront",
            "--score-min", "L,0,2.0", # threshold for contigs or scaffolds mapping to gene
            "--mp", "1,1", "--np", "1",     # ← optional soft penalties
            "--rdg", "0,1", "--rfg", "0,1", # ← optional gap open/ext penalties
            "-x", index_prefix,
            "-f", "-U", scaffolds_file,
            "--local",
            "--very-sensitive",
            "--seed", "42",
            "-k", "16",
            "-p", str(self.threads),
            "-S", sam_out
        ]

        self.logger.info(f"Mapping scaffolds for {sample_name} back to reference ...")
        result = subprocess.run(cmd, text=True, capture_output=True)
        if result.returncode != 0:
            self.logger.error(f"Bowtie2 failed for {sample_name}")
            self.logger.error(f"STDERR:\n{result.stderr.strip()}")
            raise RuntimeError(f"Bowtie2 mapping failed for {sample_name}")
        else:
            self.logger.info(f"Mapping complete for {sample_name}. SAM saved at: {sam_out}")

        return sam_out


# ============================================================
# CLI Entry Point
# ============================================================
def main():
    parser = argparse.ArgumentParser(description="Reference mapping pipeline (metaProtease).")
    subparsers = parser.add_subparsers(dest="command", required=True)

    # --- Build index ---
    idx = subparsers.add_parser("index", help="Build Bowtie2 index")
    idx.add_argument("--fasta", "-f", required=True)
    idx.add_argument("--outdir", "-o", help="Output directory")
    idx.add_argument("--threads", "-t", type=int, default=8)
    idx.add_argument("--force", action="store_true")
    idx.set_defaults(func="index")

    # --- Map full pipeline ---
    map_cmd = subparsers.add_parser("map", help="Run full mapping + assembly pipeline")
    map_cmd.add_argument("--r1", required=True)
    map_cmd.add_argument("--r2", required=True)
    map_cmd.add_argument("--index-path", "-x", required=True)
    map_cmd.add_argument("--outdir", "-o", help="Output directory")
    map_cmd.add_argument("--threads", "-t", type=int, default=8)
    map_cmd.add_argument("--force", action="store_true")

    args = parser.parse_args()

    if args.command == "index":
        mapper = ReferenceMapper(outdir=args.outdir, threads=args.threads, force=args.force)
        mapper.build_index(args.fasta)

    elif args.command == "map":
        mapper = ReferenceMapper(outdir=args.outdir, threads=args.threads, force=args.force)
        index_prefix = args.index_path
        sam = mapper.map_reads(index_prefix, args.r1, args.r2)
        bam_mapped = mapper.sam_to_bam(sam)
        r1_mapped, r2_mapped = mapper.bam_to_fastq(bam_mapped)
        assembly_dir = mapper.assembled_mapped_reads(r1_mapped, r2_mapped)
        scaffolds_file = os.path.join(assembly_dir, "scaffolds.fasta")
        mapper.map_scaffolds2refs(index_prefix, scaffolds_file)


if __name__ == "__main__":
    main()

