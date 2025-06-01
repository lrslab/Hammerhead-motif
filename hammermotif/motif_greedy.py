#!/usr/bin/env python3
"""
Greedy motif extraction module from hammermotif package.
Optimized for speed with efficient algorithms.
"""
import os
import re
from collections import Counter, defaultdict
import numpy as np
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from scipy.stats import chi2_contingency
import time


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
# 2. Optimized K-mer Counting and Statistics
##############################################

def count_kmers_optimized(seqs, k, seq_indices=None):
    """
    Optimized k-mer counting that tracks which sequences contain each k-mer.
    Returns both counts and sequence indices for each k-mer.
    """
    kmer_data = defaultdict(lambda: {'count': 0, 'seq_indices': set()})
    
    if seq_indices is None:
        seq_indices = list(range(len(seqs)))
    
    # Pre-compile valid bases check
    valid_bases = set("ACGT")
    
    for idx in seq_indices:
        seq = seqs[idx]
        seq_len = len(seq)
        
        # Skip if sequence is too short
        if seq_len < k:
            continue
            
        # Use sliding window
        for i in range(seq_len - k + 1):
            kmer = seq[i:i + k]
            
            # Fast check for valid k-mer
            if all(base in valid_bases for base in kmer):
                kmer_data[kmer]['count'] += 1
                kmer_data[kmer]['seq_indices'].add(idx)
    
    return kmer_data


def count_kmers_background(ref_seqs, k, sample_size=None):
    """
    Count k-mers in background sequences with optional sampling for speed.
    """
    if sample_size and len(ref_seqs) > sample_size:
        # Sample sequences for faster computation
        import random
        sampled = random.sample(ref_seqs, sample_size)
        total_sampled_len = sum(len(s) for s in sampled)
        total_ref_len = sum(len(s) for s in ref_seqs)
        scaling_factor = total_ref_len / total_sampled_len
    else:
        sampled = ref_seqs
        scaling_factor = 1.0
    
    counts = Counter()
    valid_bases = set("ACGT")
    
    for seq in sampled:
        seq_len = len(seq)
        if seq_len < k:
            continue
            
        for i in range(seq_len - k + 1):
            kmer = seq[i:i + k]
            if all(base in valid_bases for base in kmer):
                counts[kmer] += 1
    
    # Scale counts if sampling was used
    if scaling_factor > 1:
        for kmer in counts:
            counts[kmer] = int(counts[kmer] * scaling_factor)
    
    return counts


def total_possible_positions(seqs, k, indices=None):
    """
    Compute the total number of possible positions for a motif of length k.
    """
    if indices is not None:
        return sum(max(0, len(seqs[i]) - k + 1) for i in indices)
    else:
        return sum(max(0, len(s) - k + 1) for s in seqs)


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


##############################################
# 3. Optimized Greedy Motif Extraction
##############################################

def greedy_motif_extraction(mod_seqs, ref_seqs, k, chi2_threshold=500, max_motifs=20):
    """
    Optimized greedy motif extraction with better AT-rich filtering.
    """
    print(f"\nOptimized greedy extraction for k={k}")
    start_time = time.time()
    
    motifs_found = []
    remaining_indices = set(range(len(mod_seqs)))
    
    # Calculate genome GC content for bias correction
    genome_gc = np.mean([calculate_gc_content(seq) for seq in ref_seqs[:5]])
    print(f"  Genome GC content: {genome_gc:.2%}")
    
    # Pre-compute background k-mer counts once (with sampling for speed)
    print("  Computing background model...")
    bg_counts = count_kmers_background(ref_seqs, k, sample_size=min(len(ref_seqs), 5))
    total_bg = sum(max(0, len(s) - k + 1) for s in ref_seqs[:5]) * (len(ref_seqs) / 5)
    
    iteration = 1
    while len(remaining_indices) > 0 and len(motifs_found) < max_motifs:
        if iteration % 5 == 1:
            print(f"  Iteration {iteration}: {len(remaining_indices)} sequences remaining")
        
        # Count k-mers only in remaining sequences
        kmer_data = count_kmers_optimized(mod_seqs, k, remaining_indices)
        
        if not kmer_data:
            break
        
        total_obs = total_possible_positions(mod_seqs, k, remaining_indices)
        
        best_motif = None
        best_chi2 = chi2_threshold
        best_seq_indices = set()
        
        # Evaluate only top k-mers by count for speed
        top_kmers = sorted(kmer_data.items(), 
                          key=lambda x: x[1]['count'], 
                          reverse=True)[:100]  # Only check top 100
        
        for kmer, data in top_kmers:
            count = data['count']
            seq_indices = data['seq_indices']
            
            # Skip if too few occurrences
            if count < 5:
                continue
            
            # Filter out low complexity sequences
            if is_low_complexity(kmer):
                continue
            
            # Calculate GC content of k-mer
            kmer_gc = calculate_gc_content(kmer)
            
            # Skip extremely AT-rich or GC-rich k-mers unless they're very enriched
            if (kmer_gc < 0.2 or kmer_gc > 0.8) and count < 20:
                continue
            
            # Calculate GC bias correction factor
            # If k-mer GC is very different from genome GC, apply correction
            gc_diff = abs(kmer_gc - genome_gc)
            gc_bias_correction = 1.0
            if gc_diff > 0.3:  # More than 30% difference
                gc_bias_correction = 1.5  # Require higher enrichment
            
            bg_count = bg_counts.get(kmer, 1)  # Avoid zero
            chi2, p = compute_chi2_fast(count, total_obs, bg_count, total_bg, gc_bias_correction)
            
            # Additional filter: require higher chi2 for AT-rich motifs
            if kmer_gc < 0.3:  # AT-rich
                chi2 = chi2 * 0.5  # Reduce chi2 score
            
            # Choose the candidate with highest chi2
            if chi2 > best_chi2:
                best_chi2 = chi2
                best_motif = kmer
                best_seq_indices = seq_indices
        
        if best_motif is None:
            break
        
        motifs_found.append((best_motif, best_chi2))
        if iteration <= 10:
            motif_gc = calculate_gc_content(best_motif)
            print(f"    Selected: {best_motif} (χ²={best_chi2:.0f}, GC={motif_gc:.2f}, {len(best_seq_indices)} seqs)")
        
        # Remove sequences containing the motif
        remaining_indices -= best_seq_indices
        iteration += 1
    
    elapsed = time.time() - start_time
    print(f"  Completed in {elapsed:.1f} seconds")
    
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