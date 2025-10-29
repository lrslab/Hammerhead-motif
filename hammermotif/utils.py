#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 14/3/2025 9:15 pm
# @Author  : Runsheng
# @File    : utils.py

import os
import logging
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import GC
from itertools import product
from collections import defaultdict
import numpy as np


def read_bed_file(bed_file):
    """Read modification sites from BED file."""
    sites = []
    with open(bed_file, 'r') as file:
        for line in file:
            if line.startswith("#"):
                continue
            parts = line.strip().split()
            chrom = parts[0]
            start = int(parts[1])
            end = int(parts[2])
            score = float(parts[4]) if len(parts) > 4 else 1.0
            strand = parts[5] if len(parts) > 5 else '+'
            sites.append({
                'chrom': chrom,
                'start': start,
                'end': end,
                'score': score,
                'strand': strand
            })
    logging.info(f"Read {len(sites)} modification sites from {bed_file}")
    return sites


def read_fasta_file(fasta_file):
    """
    Read genome from FASTA file using BioPython.
    Returns a dictionary of chromosome IDs to Seq objects for efficient manipulation.
    """
    try:
        genome = {record.id: record.seq
                  for record in SeqIO.parse(fasta_file, "fasta")}
        logging.info(f"Successfully read genome from {fasta_file}")
        return genome
    except Exception as e:
        logging.error(f"Error reading FASTA file {fasta_file}: {str(e)}")
        raise


def extract_sequences(sites, genome, output_file, window_size=20):
    """
    Extract sequences around modification sites with configurable window size.

    Parameters:
    -----------
    sites : list
        List of dictionaries containing modification site information
    genome : dict
        Dictionary of chromosome sequences
    output_file : str
        Path to output FASTA file
    window_size : int
        Size of the window around modification sites (half on each side)
    """
    with open(output_file, 'w') as f:
        for i, site in enumerate(sites):
            chrom = site['chrom']
            center = (site['start'] + site['end']) // 2
            start = max(0, center - window_size)
            end = min(len(genome[chrom]), center + window_size + 1)

            sequence = genome[chrom][start:end]

            # Handle strand orientation
            if site['strand'] == '-':
                sequence = str(Seq(sequence).reverse_complement())

            # Write sequence to FASTA file
            f.write(f">{chrom}_{start}_{end}_{site['strand']}\n")
            f.write(f"{sequence}\n")

    logging.info(f"Extracted {len(sites)} sequences to {output_file}")


def create_background_model(genome, output_file, order=3):
    """
    Create Markov background model from genome sequences.
    Uses Biopython Seq objects for efficient sequence manipulation.
    """
    # Initialize counts for k-mers
    kmer_counts = defaultdict(int)
    bases = ['A', 'C', 'G', 'T']

    # Count k-mers in genome using Biopython Seq objects
    for seq in genome.values():
        # Convert to string only once for iteration
        seq_str = str(seq)
        for i in range(len(seq_str) - order):
            kmer = seq_str[i:i + order + 1]
            if set(kmer).issubset(set(bases)):  # Only count valid DNA sequences
                kmer_counts[kmer] += 1

    # Calculate probabilities and write background model
    with open(output_file, 'w') as f:
        f.write(f"# order-{order} Markov background model\n")

        # Write transition probabilities
        for prefix in [''.join(p) for p in product(bases, repeat=order)]:
            total = sum(kmer_counts[prefix + b] for b in bases)
            if total > 0:
                for b in bases:
                    prob = kmer_counts[prefix + b] / total
                    f.write(f"{prefix}\t{b}\t{prob:.6f}\n")

    logging.info(f"Created background model at {output_file}")


def calculate_gc_content(sequence):
    """
    Calculate GC content of a sequence using Biopython's GC function.
    Much faster than manual calculation for large sequences.

    Parameters:
    -----------
    sequence : str or Bio.Seq.Seq
        DNA sequence

    Returns:
    --------
    float
        GC content percentage
    """
    if isinstance(sequence, str):
        sequence = Seq(sequence)
    return GC(sequence)


def reverse_complement(sequence):
    """
    Get reverse complement of a DNA sequence using Biopython.
    Much faster than manual calculation, especially for large sequences.

    Parameters:
    -----------
    sequence : str or Bio.Seq.Seq
        DNA sequence

    Returns:
    --------
    str
        Reverse complement sequence
    """
    if isinstance(sequence, str):
        sequence = Seq(sequence)
    return str(sequence.reverse_complement())


def validate_dna_sequence(sequence):
    """
    Validate that a sequence contains only valid DNA bases.

    Parameters:
    -----------
    sequence : str or Bio.Seq.Seq
        DNA sequence to validate

    Returns:
    --------
    bool
        True if sequence is valid, False otherwise
    """
    if isinstance(sequence, Seq):
        sequence = str(sequence)
    return all(base in 'ACGTNacgtn' for base in sequence)


def chunk_sequence(sequence, chunk_size=1000):
    """
    Generator to process large sequences in chunks.
    Useful for memory-efficient processing of large genomes.

    Parameters:
    -----------
    sequence : Bio.Seq.Seq or str
        Input sequence
    chunk_size : int
        Size of chunks to yield

    Yields:
    -------
    Bio.Seq.Seq
        Sequence chunk
    """
    if isinstance(sequence, str):
        sequence = Seq(sequence)

    for i in range(0, len(sequence), chunk_size):
        yield sequence[i:i + chunk_size]


# Example usage of memory-efficient sequence processing
def process_large_genome(genome_dict):
    """
    Example of processing a large genome efficiently.

    Parameters:
    -----------
    genome_dict : dict
        Dictionary of chromosome IDs to sequences
    """
    for chrom, seq in genome_dict.items():
        gc_contents = []
        for chunk in chunk_sequence(seq):
            gc_contents.append(calculate_gc_content(chunk))
        avg_gc = np.mean(gc_contents)
        logging.info(f"Chromosome {chrom} average GC content: {avg_gc:.2f}%")


def calculate_gc_content(seq):
    """Calculate GC content of a sequence."""
    gc_count = seq.count('G') + seq.count('C')
    return gc_count / len(seq) if len(seq) > 0 else 0


def is_low_complexity(seq, threshold=0.7):
    """Check if sequence is low complexity (e.g., mostly one or two bases)."""
    if len(seq) < 3:
        return True
    
    # Check for single base repeats
    for base in 'ACGT':
        if seq.count(base) / len(seq) > threshold:
            return True
    
    # Check for dinucleotide repeats
    for i in range(len(seq) - 1):
        dinuc = seq[i:i+2]
        if seq.count(dinuc) * 2 / len(seq) > threshold:
            return True
    
    return False


def compute_chi2_fast(obs, total_obs, obs_bg, total_bg, gc_bias_correction=1.0):
    """
    Fast chi-square calculation with error handling and GC bias correction.
    """
    if total_obs <= 0 or total_bg <= 0 or obs > total_obs or obs_bg > total_bg:
        return 0, 1
    
    # Expected frequency with GC bias correction
    expected = (obs + obs_bg) * total_obs / (total_obs + total_bg) * gc_bias_correction
    
    if expected <= 0:
        return 0, 1
    
    # Simple chi-square approximation for speed
    chi2 = (obs - expected) ** 2 / expected
    
    # Very rough p-value approximation (not exact but fast)
    if chi2 > 10.83:  # ~0.001 significance
        p = 0.0001
    elif chi2 > 6.64:  # ~0.01 significance
        p = 0.01
    elif chi2 > 3.84:  # ~0.05 significance
        p = 0.05
    else:
        p = 0.5
    
    return chi2, p


if __name__ == "__main__":
    # Set up logging
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )

    # Example usage
    try:
        # Read a FASTA file
        genome = read_fasta_file("example.fasta")

        # Process each sequence
        for chrom, seq in genome.items():
            gc_content = calculate_gc_content(seq)
            logging.info(f"Chromosome {chrom} GC content: {gc_content:.2f}%")

            # Example of efficient reverse complement
            rev_comp = reverse_complement(seq)
            logging.info(f"Processed reverse complement for {chrom}")

    except Exception as e:
        logging.error(f"Error processing genome: {str(e)}")