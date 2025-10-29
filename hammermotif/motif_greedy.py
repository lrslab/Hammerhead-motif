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
from .utils import calculate_gc_content

def greedy_motif_extraction(mod_seqs, ref_seqs, k, chi2_threshold=50, max_motifs=20):
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
    # Use more samples for better background model
    sample_size = min(len(ref_seqs), max(10, len(ref_seqs) // 10))
    bg_counts = count_kmers_background(ref_seqs, k, sample_size=sample_size)
    # Calculate total background positions more accurately
    if sample_size < len(ref_seqs):
        sampled_length = sum(len(s) for s in ref_seqs[:sample_size])
        total_ref_length = sum(len(s) for s in ref_seqs)
        total_bg = sum(max(0, len(s) - k + 1) for s in ref_seqs[:sample_size]) * (total_ref_length / sampled_length)
    else:
        total_bg = sum(max(0, len(s) - k + 1) for s in ref_seqs)
    
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
    
    # Use regex to find valid ACGT stretches
    import re
    valid_pattern = re.compile(r'[ACGT]+')
    
    for idx in seq_indices:
        seq = seqs[idx]
        
        # Find all valid ACGT stretches in the sequence
        for match in valid_pattern.finditer(seq):
            valid_stretch = match.group()
            stretch_len = len(valid_stretch)
            
            # Extract all k-mers from this valid stretch
            if stretch_len >= k:
                for i in range(stretch_len - k + 1):
                    kmer = valid_stretch[i:i + k]
                    kmer_data[kmer]['count'] += 1
                    kmer_data[kmer]['seq_indices'].add(idx)
    
    return kmer_data


def count_kmers_background(ref_seqs, k, sample_size=None):
    """
    Count k-mers in background sequences with optional sampling for speed.
    Optimized version using regex for fast validation.
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
    # Use regex to find all continuous ACGT stretches
    import re
    valid_pattern = re.compile(r'[ACGT]+')
    
    for seq in sampled:
        # Find all valid ACGT stretches in the sequence
        for match in valid_pattern.finditer(seq):
            valid_stretch = match.group()
            stretch_len = len(valid_stretch)
            
            # Extract all k-mers from this valid stretch
            if stretch_len >= k:
                for i in range(stretch_len - k + 1):
                    kmer = valid_stretch[i:i + k]
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

def greedy_motif_extraction(mod_seqs, ref_seqs, k, chi2_threshold=50, max_motifs=20):
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
    # Use more samples for better background model
    sample_size = min(len(ref_seqs), max(10, len(ref_seqs) // 10))
    bg_counts = count_kmers_background(ref_seqs, k, sample_size=sample_size)
    # Calculate total background positions more accurately
    if sample_size < len(ref_seqs):
        sampled_length = sum(len(s) for s in ref_seqs[:sample_size])
        total_ref_length = sum(len(s) for s in ref_seqs)
        total_bg = sum(max(0, len(s) - k + 1) for s in ref_seqs[:sample_size]) * (total_ref_length / sampled_length)
    else:
        total_bg = sum(max(0, len(s) - k + 1) for s in ref_seqs)
    
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
# 4. Main Pipeline (Greedy extraction only)
##############################################

def main(bed_file, fasta_file, k=7, chi2_threshold=50, max_motifs=20, flank=7):
    """
    Greedy motif extraction pipeline. 
    Note: Motif merging is now handled separately in motif_merger.py
    
    1. Load the reference genome.
    2. Parse the BED file to get regions of modified sites.
    3. Extract sequences from these regions.
    4. Use the entire genome as background.
    5. Iteratively extract motifs using the greedy algorithm.
    
    Returns raw motifs with chi2 scores for further processing.
    """
    from .motif import parse_bed_file, extract_sequences
    
    # Load reference genome (as a dictionary of SeqRecords).
    genome = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))
    # Parse BED file to get regions.
    regions = parse_bed_file(bed_file, flank_size=flank)
    # Extract modified sequences.
    mod_seqs = extract_sequences(regions, genome)
    if not mod_seqs:
        print("No modified sequences extracted from BED.")
        return []
    print(f"Extracted {len(mod_seqs)} modified sequences from BED.")
    # Use the entire genome as background.
    ref_seqs = [str(rec.seq).upper() for rec in genome.values()]
    # Run the greedy motif extraction.
    raw_motifs_with_stats = greedy_motif_extraction(mod_seqs, ref_seqs, k, chi2_threshold, max_motifs)
    print("\nRaw extracted motifs:")
    for motif, chi2_val in raw_motifs_with_stats:
        print(f"Motif: {motif}, Chi2: {chi2_val:.2f}")

    return raw_motifs_with_stats


if __name__ == "__main__":
    # Example usage: adjust the file paths to your BED and FASTA files.
    bed_file = "methylation_sites.bed"  # BED file with modified site coordinates
    fasta_file = "genome.fasta"  # Reference genome FASTA
    # Run greedy extraction only
    raw_motifs = main(bed_file, fasta_file, k=7, chi2_threshold=500, max_motifs=20, flank=7)
    print(f"\nExtracted {len(raw_motifs)} raw motifs. Use motif_merger.py for merging.")