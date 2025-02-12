#!/usr/bin/env python3
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
            # Extend the region by flank_size on both ends.
            start = max(0, start - flank_size)
            end += flank_size
            regions.append({"chrom": chrom, "start": start, "end": end, "strand": strand})
    return regions


def extract_sequences(regions, genome):
    """
    Given a list of regions and a genome (dictionary from SeqIO.to_dict()),
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
        # Ensure we do not exceed chromosome length.
        chrom_len = len(genome[chrom].seq)
        if end > chrom_len:
            end = chrom_len
        seq = genome[chrom].seq[start:end]
        if strand == '-':
            seq = seq.reverse_complement()
        sequences.append(str(seq).upper())
    return sequences


##############################################
# 2. K-mer Counting and Statistics
##############################################

def count_kmers_in_seqs(seqs, k):
    """
    Count all k-mers of length k in a list of sequences.
    Only k-mers with letters A, C, G, T are counted.
    Returns a Counter object mapping each k-mer to its count.
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
    Build a 2x2 contingency table:
          | Motif   | Other sites |
    -------------------------------
    Modified|   obs    | total_obs - obs |
    BG      |  obs_bg  | total_bg - obs_bg |

    Compute the chi-square statistic (and p-value) for the table using chi2_contingency.
    """
    table = [[obs, obs_bg], [total_obs - obs, total_bg - obs_bg]]
    chi2, p, _, _ = chi2_contingency(table)
    return chi2, p


##############################################
# 3. Greedy Motif Extraction (Final Motif Finding)
##############################################

def greedy_motif_extraction(mod_seqs, ref_seqs, k, chi2_threshold=500, max_motifs=20):
    """
    Iteratively extract motifs from the modified sequence set using a greedy algorithm.

    For each candidate k-mer (extracted from mod_seqs):
      - Count its occurrences in mod_seqs (obs) and in ref_seqs (obs_bg).
      - Compute the total possible positions in mod_seqs and ref_seqs.
      - Calculate a chi-square statistic comparing obs vs. obs_bg.

    The best candidate (with highest chi2 above chi2_threshold) is chosen as a motif.
    Then, all sequences in mod_seqs that contain this motif are removed,
    and the process repeats until no candidate motif passes the threshold
    or max_motifs have been extracted.

    Returns a list of tuples: (motif, chi2_value).
    """
    motifs_found = []
    remaining_seqs = mod_seqs.copy()
    iteration = 1
    while remaining_seqs and len(motifs_found) < max_motifs:
        print(f"Iteration {iteration}: {len(remaining_seqs)} sequences remaining")
        kmer_counts = count_kmers_in_seqs(remaining_seqs, k)
        if not kmer_counts:
            break
        total_obs = total_possible_positions(remaining_seqs, k)
        bg_counts = count_kmers_in_seqs(ref_seqs, k)
        total_bg = total_possible_positions(ref_seqs, k)

        best_motif = None
        best_chi2 = 0
        for motif, count in kmer_counts.items():
            bg_count = bg_counts.get(motif, 0)
            chi2, p = compute_chi2(count, total_obs, bg_count, total_bg)
            # Choose the candidate with highest chi2 above the threshold.
            if chi2 > best_chi2 and chi2 >= chi2_threshold:
                best_chi2 = chi2
                best_motif = motif
        if best_motif is None:
            print("No candidate motif exceeded the chi2 threshold.")
            break
        motifs_found.append((best_motif, best_chi2))
        print(f"Selected motif: {best_motif} (chi2={best_chi2:.2f})")
        # Remove all sequences that contain the best motif.
        pattern = re.compile(best_motif)
        new_remaining = [s for s in remaining_seqs if not pattern.search(s)]
        remaining_seqs = new_remaining
        iteration += 1
    return motifs_found


##############################################
# 4. Main Pipeline
##############################################

def main(bed_file, fasta_file, k=6, chi2_threshold=500, max_motifs=20, flank=7):
    """
    Complete pipeline starting from a BED file and a reference FASTA to extract final methylation motifs.

    1. Load the reference genome.
    2. Parse the BED file to get regions of modified sites.
    3. Extract sequences from these regions.
    4. Use the entire genome as background.
    5. Iteratively extract motifs using the greedy algorithm.

    Prints the final motifs along with their chi2 statistic.
    """
    # Load reference genome (as a dictionary of SeqRecords).
    genome = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))
    # Parse BED file to get regions.
    regions = parse_bed_file(bed_file, flank_size=flank)
    # Extract modified sequences.
    mod_seqs = extract_sequences(regions, genome)
    if not mod_seqs:
        print("No modified sequences extracted from BED.")
        return
    print(f"Extracted {len(mod_seqs)} modified sequences from BED.")
    # Use the entire genome as background.
    ref_seqs = [str(rec.seq).upper() for rec in genome.values()]
    # Run the greedy motif extraction.
    final_motifs = greedy_motif_extraction(mod_seqs, ref_seqs, k, chi2_threshold, max_motifs)
    print("\nFinal extracted motifs:")
    for motif, chi2_val in final_motifs:
        print(f"Motif: {motif}, Chi2: {chi2_val:.2f}")


if __name__ == "__main__":
    # Example usage: adjust the file paths to your BED and FASTA files.
    bed_file = "methylation_sites.bed"  # BED file with modified site coordinates
    fasta_file = "genome.fasta"  # Reference genome FASTA
    main(bed_file, fasta_file, k=6, chi2_threshold=500, max_motifs=20, flank=7)