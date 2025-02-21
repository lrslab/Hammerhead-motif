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
# 4. Motif Merging (Clustering and Consensus)
##############################################

def degenerate_code(bases):
    """
    Map a sorted tuple of bases to an IUPAC degenerate code.
    For example, ('A', 'T') returns 'W'.
    """
    mapping = {
        ('A',): 'A',
        ('C',): 'C',
        ('G',): 'G',
        ('T',): 'T',
        ('A', 'C'): 'M',
        ('A', 'G'): 'R',
        ('A', 'T'): 'W',
        ('C', 'G'): 'S',
        ('C', 'T'): 'Y',
        ('G', 'T'): 'K',
        ('A', 'C', 'G'): 'V',
        ('A', 'C', 'T'): 'H',
        ('A', 'G', 'T'): 'D',
        ('C', 'G', 'T'): 'B',
        ('A', 'C', 'G', 'T'): 'N'
    }
    return mapping.get(tuple(sorted(bases)), 'N')


def consensus_from_cluster(cluster):
    """
    Generate a consensus motif from a list of motifs (all of the same length).
    At each position, use the set of bases from all motifs to produce a degenerate code.
    """
    if not cluster:
        return ""
    length = len(cluster[0])
    consensus = []
    for i in range(length):
        letters = {motif[i] for motif in cluster}
        consensus.append(degenerate_code(letters))
    return ''.join(consensus)


def extract_core(motif, core_length=4):
    """
    Extract the central core of the motif.
    For a motif longer than core_length, return the central core; otherwise, return the motif.
    """
    L = len(motif)
    if L <= core_length:
        return motif
    start = (L - core_length) // 2
    return motif[start:start + core_length]


def merge_by_core(motifs, core_length=4):
    """
    Group motifs by their extracted core (using the central core of length core_length)
    and generate a consensus motif for each group.
    """
    core_groups = {}
    for m in motifs:
        core = extract_core(m, core_length)
        if core not in core_groups:
            core_groups[core] = []
        core_groups[core].append(m)
    merged = [consensus_from_cluster(group) for group in core_groups.values()]
    return merged


def merge_similar_motifs(motifs, max_distance=1):
    """
    Merge similar motifs by clustering those that have Hamming distance <= max_distance,
    and generate a consensus motif for each cluster.
    """
    clusters = []
    used = set()
    motifs = list(motifs)
    for i, m1 in enumerate(motifs):
        if i in used:
            continue
        cluster = [m1]
        used.add(i)
        for j in range(i + 1, len(motifs)):
            m2 = motifs[j]
            if j in used or len(m1) != len(m2):
                continue
            if sum(a != b for a, b in zip(m1, m2)) <= max_distance:
                cluster.append(m2)
                used.add(j)
        clusters.append(cluster)
    merged = [consensus_from_cluster(cluster) for cluster in clusters]
    return merged


def longest_common_substring(strings, min_length=4):
    """
    Given a list of strings, find the longest common substring of at least min_length.
    Returns the longest substring common to all strings, or None if none found.
    """
    if not strings:
        return None
    reference = strings[0]
    common_substrings = set()
    L = len(reference)
    for i in range(L):
        for j in range(i + min_length, L + 1):
            sub = reference[i:j]
            common_substrings.add(sub)
    for s in strings[1:]:
        common_substrings = {sub for sub in common_substrings if sub in s}
        if not common_substrings:
            return None
    # Return the longest common substring
    return max(common_substrings, key=len) if common_substrings else None


def finalize_cluster(cluster, min_common_length=4):
    """
    Given a cluster of motifs (strings), finalize the cluster by extracting the longest
    common substring of length at least min_common_length.
    If found, return it; otherwise, return the consensus from the cluster.
    """
    lcs = longest_common_substring(cluster, min_length=min_common_length)
    if lcs:
        return lcs
    else:
        return consensus_from_cluster(cluster)


def merge_final_motifs(motifs, min_common_length=4, max_distance=1, core_length=4):
    """
    Merge a list of raw motifs using multiple strategies:
      1. Merge by Hamming distance.
      2. Merge by extracted core.
      3. Finalize each cluster by extracting the longest common substring.
    Returns the final list of merged motifs.
    """
    merged_hd = merge_similar_motifs(motifs, max_distance=max_distance)
    merged_core = merge_by_core(motifs, core_length=core_length)
    # Combine both sets and then finalize clusters by common substring.
    combined = list(set(merged_hd + merged_core))
    # If there are multiple motifs, further cluster them by common substring.
    # For simplicity, we treat the entire list as one cluster and finalize.
    final = finalize_cluster(combined, min_common_length=min_common_length)
    return final


##############################################
# 5. Main Pipeline
##############################################

def main(bed_file, fasta_file, k=7, chi2_threshold=500, max_motifs=20, flank=7,
         core_length=4, min_common_length=4):
    """
    Complete pipeline starting from a BED file and a reference FASTA to extract final methylation motifs.
    1. Load the reference genome.
    2. Parse the BED file to get regions of modified sites.
    3. Extract sequences from these regions.
    4. Use the entire genome as background.
    5. Iteratively extract motifs using the greedy algorithm.
    6. Merge similar motif variants using Hamming-distance and core-based methods.
    7. Finalize the merged cluster by extracting the longest common substring.
    Prints the final merged motif.
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
    raw_motifs_with_stats = greedy_motif_extraction(mod_seqs, ref_seqs, k, chi2_threshold, max_motifs)
    print("\nRaw extracted motifs:")
    for motif, chi2_val in raw_motifs_with_stats:
        print(f"Motif: {motif}, Chi2: {chi2_val:.2f}")

    # Extract just the raw motif strings.
    raw_motifs = [m for m, _ in raw_motifs_with_stats]

    # Merge similar motifs using our combined strategy.
    final_motif = merge_final_motifs(raw_motifs, min_common_length=min_common_length,
                                     max_distance=1, core_length=core_length)
    print("\nFinal merged motif:")
    print(final_motif)


if __name__ == "__main__":
    # Example usage: adjust the file paths to your BED and FASTA files.
    bed_file = "methylation_sites.bed"  # BED file with modified site coordinates
    fasta_file = "genome.fasta"  # Reference genome FASTA
    # Adjust k, chi2_threshold, and merging parameters as needed.
    main(bed_file, fasta_file, k=7, chi2_threshold=500, max_motifs=20, flank=7,
         core_length=4, min_common_length=4)