#!/usr/bin/env python3
"""
Comprehensive Bacterial DNA Methylation Motif Caller

This script provides a command-line interface for calling bacterial DNA methylation motifs
from a single genome FASTA file and a BED file containing modification sites.

Usage:
    conda activate py3
    python methylation_motif_caller.py <genome_file> <bed_file> [output_directory]
"""

import os
import sys
import logging
from pathlib import Path
import argparse

# Import the MethylationMotifCaller class from the hammermotif package
from hammermotif.caller import MethylationMotifCaller

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

def main():
    """Main entry point for single genome analysis."""
    print("Bacterial DNA Methylation Motif Caller")
    print("======================================\n")
    
    
    
    # Parse command line arguments
    parser = argparse.ArgumentParser(
        description='Bacterial DNA Methylation Motif Caller for a single genome.'
    )
    parser.add_argument(
        'genome_file', 
        help='Path to the genome FASTA file.'
    )
    parser.add_argument(
        'bed_file', 
        help='Path to the BED file with modification sites.'
    )
    parser.add_argument(
        '--output-dir', 
        default=None, 
        help='Optional: Output directory for results. Defaults to results_<genome_name>.'
    )
    
    args = parser.parse_args()
    
    genome_file = Path(args.genome_file)
    bed_file = Path(args.bed_file)
    output_dir = args.output_dir

    if not genome_file.exists():
        logger.error("Genome file not found: {}".format(genome_file))
        sys.exit(1)
    if not bed_file.exists():
        logger.error("BED file not found: {}".format(bed_file))
        sys.exit(1)

    logger.info(f"Starting analysis for genome: {genome_file.name}")
    try:
        caller = MethylationMotifCaller(
            str(genome_file),
            str(bed_file),
            output_dir
        )
        discovered_motifs = caller.run_complete_analysis()
        logger.info(f"Analysis complete. Discovered motifs: {discovered_motifs}")
    except Exception as e:
        logger.error(f"An error occurred during analysis: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()