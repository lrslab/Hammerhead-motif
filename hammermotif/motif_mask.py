# !/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 13/2/2025 5:45
# @Author  : Runsheng
# @File    : motif_mask.py

# !/usr/bin/env python3
import os
import re
from collections import Counter
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from scipy.stats import chi2_contingency


##############################################
# 1. Data Parsing and Sequence Extraction
##############################################

def parse_bed_file(bed_file_path, flank_size=7):
    """
    Parse a BED file and return a list of genomic regions with added flanking regions.
    Each region is a dictionary with keys: chrom, start, end, strand.
    """
    regions = []
    with open(bed_file_path) as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.strip().split()
            chrom = parts[0]
            start = int(parts[1])
            end = int(parts[2])
            strand = parts[3] if len(parts) > 3 else '+'
            # Extend region by flank_size on both ends.
            start = max(0, start - flank_size)
            end += flank_size
            regions.append({"chrom": chrom, "start": start, "end": end, "strand": strand})
    return regions


def extract_sequences(regions, genome):
    """
    Given a list of regions and a genome (dict from SeqIO.to_dict()),
    extract the corresponding sequences (as uppercase strings).
    """
    sequences = []
    for region in regions:
        chrom = region["chrom"]
        start = region["start"]
        end = region["end"]
        strand = region["strand"]
        if chrom not in genome:
            print(f"Warning: {chrom} not found in genome.")
            continue
        chrom_len = len(genome[chrom].seq)
        if end > chrom_len:
            end = chrom_len
        seq = genome[chrom].seq[start:end]
        if strand == '-':
            seq = seq.reverse_complement()
        sequences.append(str(seq).upper())
    return sequences


##############################################
# 2. K-mer Counting and Chi-Square Statistics
##############################################

def count_kmers_in_seqs(seqs, k):
    """
    Count all k-mers of length k in a list of sequences.
    Only k-mers with letters A, C, G, T are counted.
    Returns a Counter mapping each k-mer to its count.
    """
    counts = Counter()
    for s in seqs:
        for i in range(len(s) - k + 1):
            kmer = s[i:i + k]
            if all(base in "ACGT" for base in kmer):
                counts[kmer] += 1
    return counts


def total_possible_positions(seqs, k):
    """
    Compute the total number of possible positions for a motif of length k in a list of sequences.
    """
    total = 0
    for s in seqs:
        total += max(0, len(s) - k + 1)
    return total


def compute_chi2(obs, total_obs, obs_bg, total_bg):
    """
    Build a 2x2 contingency table and compute the chi-square statistic (and p-value)
    using chi2_contingency.
    """
    table = [[obs, obs_bg], [total_obs - obs, total_bg - obs_bg]]
    chi2, p, _, _ = chi2_contingency(table)
    return chi2, p


##############################################
# 3. Greedy Motif Extraction
##############################################

def greedy_motif_extraction(seqs, bg_seqs, k, chi2_threshold=500, max_motifs=20):
    """
    Iteratively extract motifs from seqs using a greedy algorithm.
    For each candidate k-mer:
      - Count occurrences in seqs (obs) and in bg_seqs (obs_bg)
      - Compute chi-square statistic.
    The candidate with the highest chi2 (above threshold) is selected,
    and all sequences containing it are masked (removed) for further rounds.
    Returns a list of tuples: (motif, chi2_value).
    """
    motifs_found = []
    remaining = seqs.copy()
    iteration = 1
    while remaining and len(motifs_found) < max_motifs:
        print(f"Iteration {iteration}: {len(remaining)} sequences remaining")
        kmer_counts = count_kmers_in_seqs(remaining, k)
        if not kmer_counts:
            break
        total_obs = total_possible_positions(remaining, k)
        bg_counts = count_kmers_in_seqs(bg_seqs, k)
        total_bg = total_possible_positions(bg_seqs, k)
        best_motif = None
        best_chi2 = 0
        for motif, count in kmer_counts.items():
            bg_count = bg_counts.get(motif, 0)
            chi2, p = compute_chi2(count, total_obs, bg_count, total_bg)
            if chi2 > best_chi2 and chi2 >= chi2_threshold:
                best_chi2 = chi2
                best_motif = motif
        if best_motif is None:
            print("No candidate motif exceeded the chi2 threshold.")
            break
        motifs_found.append((best_motif, best_chi2))
        print(f"Selected motif: {best_motif} (chi2={best_chi2:.2f})")
        # Remove sequences that contain the best motif.
        pattern = re.compile(best_motif)
        remaining = [s for s in remaining if not pattern.search(s)]
        iteration += 1
    return motifs_found


##############################################
# 4. Masking Function
##############################################

def mask_motifs_in_seqs(seqs, motifs):
    """
    For each sequence in seqs, mask (replace with 'N') all occurrences of any motif in motifs.
    Returns a new list of sequences.
    """
    masked_seqs = []
    for s in seqs:
        masked = s
        for motif in motifs:
            # Replace all occurrences with N's of equal length.
            masked = re.sub(motif, "N" * len(motif), masked)
        masked_seqs.append(masked)
    return masked_seqs


##############################################
# 5. Motif Merging via Masking Iterative Approach
##############################################

def iterative_motif_extraction(mod_seqs, ref_seqs, k_small=4, k_large=7,
                               chi2_threshold_small=500, chi2_threshold_large=500,
                               max_motifs_small=5, max_motifs_large=10):
    """
    First, extract significant motifs using a small k-mer (k_small).
    Then, mask their occurrences in mod_seqs.
    Finally, run the greedy motif extraction on the masked sequences with a larger k (k_large).
    Returns a tuple (small_motifs, large_motifs) of the extracted motif lists.
    """
    print("=== Small k-mer extraction ===")
    small_results = greedy_motif_extraction(mod_seqs, ref_seqs, k_small,
                                            chi2_threshold=chi2_threshold_small,
                                            max_motifs=max_motifs_small)
    small_motifs = [m for m, chi2_val in small_results]
    print("\nSmall k-mer motifs:")
    for m, chi2_val in small_results:
        print(f"{m} (chi2={chi2_val:.2f})")

    # Mask small motifs in the modified sequences.
    masked_mod_seqs = mask_motifs_in_seqs(mod_seqs, small_motifs)
    print(f"\nMasked sequences; first 3 examples:")
    for s in masked_mod_seqs[:3]:
        print(s)

    print("\n=== Large k-mer extraction on masked sequences ===")
    large_results = greedy_motif_extraction(masked_mod_seqs, ref_seqs, k_large,
                                            chi2_threshold=chi2_threshold_large,
                                            max_motifs=max_motifs_large)
    large_motifs = [m for m, chi2_val in large_results]
    print("\nLarge k-mer motifs:")
    for m, chi2_val in large_results:
        print(f"{m} (chi2={chi2_val:.2f})")

    return small_results, large_results


##############################################
# 6. Main Pipeline
##############################################

def main(bed_file, fasta_file, k_small=4, k_large=7,
         chi2_threshold_small=500, chi2_threshold_large=500,
         max_motifs_small=5, max_motifs_large=10, flank=7):
    """
    Complete pipeline:
      1. Load reference genome.
      2. Parse BED file and extract modified sequences.
      3. Use the entire genome as background.
      4. Run small k-mer extraction.
      5. Mask found small motifs.
      6. Run large k-mer extraction on masked sequences.
      7. Output both sets of motifs.
    """
    # Load genome as dictionary of SeqRecords.
    genome = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))
    # Parse BED file.
    regions = parse_bed_file(bed_file, flank_size=flank)
    # Extract modified sequences.
    mod_seqs = extract_sequences(regions, genome)
    if not mod_seqs:
        print("No modified sequences extracted from BED.")
        return
    print(f"Extracted {len(mod_seqs)} modified sequences from BED.")
    # Background: use entire genome.
    ref_seqs = [str(rec.seq).upper() for rec in genome.values()]
    # Iterative motif extraction:
    small_motifs, large_motifs = iterative_motif_extraction(mod_seqs, ref_seqs,
                                                            k_small=k_small,
                                                            k_large=k_large,
                                                            chi2_threshold_small=chi2_threshold_small,
                                                            chi2_threshold_large=chi2_threshold_large,
                                                            max_motifs_small=max_motifs_small,
                                                            max_motifs_large=max_motifs_large)
    print("\nFinal extracted motifs:")
    print("Small k-mer motifs:")
    for m, chi2_val in small_motifs:
        print(f"{m} (chi2={chi2_val:.2f})")
    print("\nLarge k-mer motifs:")
    for m, chi2_val in large_motifs:
        print(f"{m} (chi2={chi2_val:.2f})")


if __name__ == "__main__":
    # Adjust file paths as needed.
    bed_file = "methylation_sites.bed"  # BED file with modified site coordinates.
    fasta_file = "genome.fasta"  # Reference genome FASTA.
    # Run the pipeline with small k=4 and then large k=7.
    main(bed_file, fasta_file, k_small=4, k_large=7,
         chi2_threshold_small=500, chi2_threshold_large=500,
         max_motifs_small=5, max_motifs_large=10, flank=7)