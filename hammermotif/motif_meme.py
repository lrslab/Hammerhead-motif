#!/usr/bin/env python
#-*- coding: utf-8 -*-
# @Time    : 14/3/2025 9:18 pm
# @Author  : Runsheng
# @File    : motif_meme.py

"""
This module handles running MEME and managing the MEME process.
"""

import os
import logging
from pathlib import Path
from .utils import read_bed_file, read_fasta_file, create_background_model, extract_sequences
from .meme_runner import run_meme, parse_meme_output
from .stat import assess_motif_significance

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)


class MotifFinder:
    def __init__(self,
                 bed_file: str,
                 fasta_file: str,
                 output_dir: str,
                 window_size: int = 20,
                 motif_min_width: int = 4,
                 motif_max_width: int = 8,
                 meme_path: str = "meme",
                 num_motifs: int = 3,
                 pvalue_threshold: float = 0.05):
        self.bed_file = bed_file
        self.fasta_file = fasta_file
        self.output_dir = Path(output_dir)
        self.window_size = window_size
        self.motif_min_width = motif_min_width
        self.motif_max_width = motif_max_width
        self.meme_path = meme_path
        self.num_motifs = num_motifs
        self.pvalue_threshold = pvalue_threshold

        # Create output directory
        self.output_dir.mkdir(parents=True, exist_ok=True)

    def run(self):
        logging.info("Starting motif finding pipeline")

        # Read input files
        sites = read_bed_file(self.bed_file)
        genome = read_fasta_file(self.fasta_file)

        # Extract sequences around modification sites
        sequences_file = self.output_dir / "modification_sites.fasta"
        extract_sequences(sites, genome, sequences_file,
                          window_size=self.window_size)

        # Create background model from the genome
        bg_model_file = self.output_dir / "background.txt"
        create_background_model(genome, bg_model_file)

        # Run MEME with optimized parameters for methylation motifs
        meme_output_dir = self.output_dir / "meme_output"
        meme_params = {
            "mod": "zoops",  # Zero or One Occurrence Per Sequence
            "nmotifs": self.num_motifs,
            "minw": self.motif_min_width,
            "maxw": self.motif_max_width,
            "objfun": "classic",  # Classic objective function
            "markov_order": 3,  # 3rd order Markov background
            "dna": True,
            "revcomp": True,  # Search both strands
            "pal": True,  # Look for palindromic motifs (common in methylation)
            "maxsize": 100000,  # Maximum dataset size
            "bfile": str(bg_model_file)
        }

        run_meme(str(sequences_file),
                 str(meme_output_dir),
                 self.meme_path,
                 **meme_params)

        # Parse MEME output and assess significance
        motifs = parse_meme_output(meme_output_dir / "meme.txt")
        significant_motifs = assess_motif_significance(
            motifs,
            self.pvalue_threshold
        )

        # Generate report
        self._generate_report(significant_motifs)

        logging.info(f"Motif finding complete. Results in {self.output_dir}")
        return significant_motifs

    def _generate_report(self, motifs):
        report_file = self.output_dir / "motif_report.txt"
        with open(report_file, 'w') as f:
            f.write("Bacterial DNA Methylation Motif Report\n")
            f.write("=====================================\n\n")

            for i, motif in enumerate(motifs, 1):
                f.write(f"Motif {i}:\n")
                f.write(f"Consensus: {motif['consensus']}\n")
                f.write(f"E-value: {motif['evalue']:.2e}\n")
                f.write(f"Sites: {motif['nsites']}\n")
                f.write(f"Information Content: {motif['ic']:.2f}\n")
                f.write("\nPosition Weight Matrix:\n")
                f.write(str(motif['pwm']))
                f.write("\n\n")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Find methylation motifs in bacterial DNA using MEME."
    )
    parser.add_argument("bed_file",
                        help="Input BED file with modification sites")
    parser.add_argument("fasta_file",
                        help="Input genome FASTA file")
    parser.add_argument("output_dir",
                        help="Output directory for results")
    parser.add_argument("--window-size",
                        type=int,
                        default=20,
                        help="Window size around modification sites (default: 20)")
    parser.add_argument("--motif-min-width",
                        type=int,
                        default=4,
                        help="Minimum motif width (default: 4)")
    parser.add_argument("--motif-max-width",
                        type=int,
                        default=8,
                        help="Maximum motif width (default: 8)")
    parser.add_argument("--meme-path",
                        default="meme",
                        help="Path to MEME binary")
    parser.add_argument("--num-motifs",
                        type=int,
                        default=3,
                        help="Number of motifs to find (default: 3)")
    parser.add_argument("--pvalue-threshold",
                        type=float,
                        default=0.05,
                        help="P-value threshold for significance (default: 0.05)")

    args = parser.parse_args()

    finder = MotifFinder(
        args.bed_file,
        args.fasta_file,
        args.output_dir,
        args.window_size,
        args.motif_min_width,
        args.motif_max_width,
        args.meme_path,
        args.num_motifs,
        args.pvalue_threshold
    )

    finder.run()