import os
import random
import re
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import defaultdict, Counter
from itertools import product
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests
from sklearn.cluster import AffinityPropagation
from tqdm import tqdm  # For progress bars

# IUPAC nucleotide ambiguity codes
IUPAC_CODES = {
    'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T',
    'R': 'AG', 'Y': 'CT', 'S': 'GC', 'W': 'AT',
    'K': 'GT', 'M': 'AC', 'B': 'CGT', 'D': 'AGT',
    'H': 'ACT', 'V': 'ACG', 'N': 'ACGT'
}


def parse_bed_file(bed_file_path, flank_size=7):
    """
    Parse a BED file to extract genomic locations.
    :param bed_file_path: Path to the BED file
    :param flank_size: Number of bases to extend on each side of the region
    :return: List of dictionaries with genomic locations
    """
    locations = []
    with open(bed_file_path) as bed_file:
        for line in bed_file:
            if line.startswith('#') or line.strip() == '':
                continue
            parts = line.strip().split()
            chrom, start, end = parts[0], int(parts[1]), int(parts[2])
            strand = parts[3] if len(parts) > 3 else '+'
            start = max(0, start - flank_size)
            end += flank_size
            locations.append({
                'chrom': chrom,
                'start': start,
                'end': end,
                'strand': strand
            })
    return locations


def extract_sequences(genomic_locations, genome_fasta):
    """
    Extract sequences from the genome based on BED file locations.
    :param genomic_locations: List of locations from parse_bed_file
    :param genome_fasta: Dictionary of chromosome sequences (from SeqIO.to_dict)
    :return: List of SeqRecord objects
    """
    sequences = []
    for loc in genomic_locations:
        chrom = loc['chrom']
        start = loc['start']
        end = loc['end']
        strand = loc['strand']
        if chrom not in genome_fasta:
            continue
        sequence = genome_fasta[chrom].seq[start:end]
        if strand == '-':
            sequence = sequence.reverse_complement()
        sequences.append(SeqRecord(sequence, id=f"{chrom}:{start}-{end}", description=""))
    return sequences


def generate_random_background(genome_fasta, num_sequences, seq_length):
    """
    Generate random background sequences from the longest chromosome.
    :param genome_fasta: Dictionary of chromosome sequences
    :param num_sequences: Number of background sequences to generate
    :param seq_length: Length of each background sequence
    :return: List of SeqRecord objects
    """
    background = []
    # Find the longest chromosome
    longest_chrom = max(genome_fasta.keys(), key=lambda x: len(genome_fasta[x]))
    chrom_seq = genome_fasta[longest_chrom].seq

    for _ in range(num_sequences):
        max_start = len(chrom_seq) - seq_length
        start = random.randint(0, max_start)
        end = start + seq_length
        sequence = chrom_seq[start:end]
        background.append(SeqRecord(sequence, id=f"{longest_chrom}:{start}-{end}", description=""))
    return background


def degenerate_expansion(kmer):
    """
    Expand a degenerate k-mer into all possible non-degenerate combinations.
    :param kmer: A k-mer potentially containing IUPAC codes
    :return: List of all possible non-degenerate k-mers
    """
    return [''.join(p) for p in product(*[IUPAC_CODES[base] for base in kmer])]


def motif_to_regex(motif):
    """
    Convert a motif (potentially with IUPAC codes) to a regular expression.
    :param motif: The motif to convert
    :return: A regex pattern string
    """
    regex = []
    for base in motif.upper():
        bases = IUPAC_CODES.get(base, base)
        if len(bases) == 1:
            regex.append(bases)
        else:
            regex.append(f'[{bases}]')
    return ''.join(regex)


def enhanced_candidate_discovery(reg_seqs, bg_seqs, min_len, max_len):
    """
    Discover candidate motifs from regulatory sequences.
    :param reg_seqs: List of SeqRecord objects (regulatory sequences)
    :param bg_seqs: List of SeqRecord objects (background sequences)
    :param min_len: Minimum motif length to consider
    :param max_len: Maximum motif length to consider
    :return: List of candidate motifs
    """
    candidates = []
    for k in range(min_len, max_len + 1):
        # Scan all possible k-mers (including degenerate ones)
        all_kmers = set()
        for seq in reg_seqs:
            seq_str = str(seq.seq)
            for i in range(len(seq_str) - k + 1):
                kmer = seq_str[i:i + k]
                if all(c in IUPAC_CODES for c in kmer):
                    all_kmers.update(degenerate_expansion(kmer))

        # Calculate weighted enrichment scores
        enrich_scores = {}
        for kmer in all_kmers:
            reg = sum(1 for seq in reg_seqs if kmer in str(seq.seq))
            bg = sum(1 for seq in bg_seqs if kmer in str(seq.seq))
            enrich = ((reg + 0.1) / len(reg_seqs)) / ((bg + 0.1) / len(bg_seqs)) * np.sqrt(k)
            enrich_scores[kmer] = enrich

        # Retain top candidates for this length
        sorted_kmers = sorted(enrich_scores.items(), key=lambda x: -x[1])
        candidates += [kmer for kmer, _ in sorted_kmers[:10]]

    return list(set(candidates))


def precise_enrichment_with_pseudo(reg_seqs, bg_seqs, candidates, fdr_threshold):
    """
    Perform precise enrichment analysis with pseudo-counts.
    :param reg_seqs: List of SeqRecord objects (regulatory sequences)
    :param bg_seqs: List of SeqRecord objects (background sequences)
    :param candidates: List of candidate motifs
    :param fdr_threshold: Significance threshold after FDR correction
    :return: Dictionary of significant motifs and their statistics
    """
    results = []
    for motif in tqdm(candidates, desc="Analyzing motifs"):
        k = len(motif)
        reg_obs = sum(str(seq.seq).count(motif) for seq in reg_seqs)
        bg_obs = sum(str(seq.seq).count(motif) for seq in bg_seqs)

        # Calculate total possible sites
        reg_total = sum(len(seq.seq) - k + 1 for seq in reg_seqs)
        bg_total = sum(len(seq.seq) - k + 1 for seq in bg_seqs)

        # Fisher's exact test with pseudo-counts
        table = [[reg_obs + 0.5, bg_obs + 0.5],
                 [reg_total - reg_obs + 0.5, bg_total - bg_obs + 0.5]]
        _, p = fisher_exact(table, alternative='greater')
        fold = ((reg_obs + 0.5) / reg_total) / ((bg_obs + 0.5) / bg_total)
        results.append((motif, p, fold, reg_obs, bg_obs, reg_total, bg_total))

    # Multiple testing correction
    p_values = [x[1] for x in results]
    _, adj_p_values, _, _ = multipletests(p_values, alpha=fdr_threshold, method='fdr_bh')

    # Filter significant results
    significant = {}
    for (motif, p, fold, reg_obs, bg_obs, reg_total, bg_total), adj_p in zip(results, adj_p_values):
        if adj_p <= fdr_threshold:
            significant[motif] = {
                'p': p,
                'adj_p': adj_p,
                'fold': fold,
                'reg_obs': reg_obs,
                'bg_obs': bg_obs,
                'reg_total': reg_total,
                'bg_total': bg_total
            }
    return significant


def cluster_motifs(motifs):
    """
    Cluster motifs based on sequence similarity.
    :param motifs: List of motifs to cluster
    :return: Dictionary mapping cluster IDs to lists of motifs
    """
    # Build distance matrix using Levenshtein distance
    distance_matrix = np.zeros((len(motifs), len(motifs)))
    for i, m1 in enumerate(motifs):
        for j, m2 in enumerate(motifs):
            distance_matrix[i, j] = levenshtein_distance(m1, m2)

    # Affinity Propagation clustering
    clusterer = AffinityPropagation(affinity='precomputed')
    labels = clusterer.fit_predict(-distance_matrix)

    # Group motifs by cluster
    clustered = defaultdict(list)
    for motif, label in zip(motifs, labels):
        clustered[label].append(motif)
    return clustered


def levenshtein_distance(s1, s2):
    """
    Compute the Levenshtein distance between two strings.
    :param s1: First string
    :param s2: Second string
    :return: Integer distance
    """
    if len(s1) < len(s2):
        return levenshtein_distance(s2, s1)
    if len(s2) == 0:
        return len(s1)
    previous_row = range(len(s2) + 1)
    for i, c1 in enumerate(s1):
        current_row = [i + 1]
        for j, c2 in enumerate(s2):
            insertions = previous_row[j + 1] + 1
            deletions = current_row[j] + 1
            substitutions = previous_row[j] + (c1 != c2)
            current_row.append(min(insertions, deletions, substitutions))
        previous_row = current_row
    return previous_row[-1]


def print_clustered_results(clusters, results_dict):
    """
    Print clustered results with representative motifs.
    :param clusters: Dictionary of motif clusters
    :param results_dict: Dictionary of motif statistics
    """
    print("Significant Motif Clusters:")
    for cluster_id, motifs in clusters.items():
        # Select the motif with the lowest p-value as the representative
        rep_motif = min(motifs, key=lambda x: results_dict[x]['p'])
        stats = results_dict[rep_motif]
        print(f"Cluster {cluster_id} (Representative: {rep_motif})")
        print(f"  Motifs: {', '.join(motifs)}")
        print(f"  p-value: {stats['p']:.2e}, Adjusted p-value: {stats['adj_p']:.2e}")
        print(f"  Fold Enrichment: {stats['fold']:.1f}x")
        print(f"  Observed in Regions: {stats['reg_obs']}/{stats['reg_total']}")
        print(f"  Observed in Background: {stats['bg_obs']}/{stats['bg_total']}\n")


def motif_enrichment_pipeline(bed_file, fasta_file, motifs=None, min_len=4, max_len=16, fdr_threshold=0.05):
    """
    Main pipeline for motif enrichment analysis.
    :param bed_file: Path to the BED file
    :param fasta_file: Path to the genome FASTA file
    :param motifs: Optional list of known motifs to include
    :param min_len: Minimum motif length to consider
    :param max_len: Maximum motif length to consider
    :param fdr_threshold: Significance threshold after FDR correction
    """
    # Load genome and extract sequences
    genome = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))
    reg_seqs = extract_sequences(parse_bed_file(bed_file), genome)
    bg_seqs = generate_random_background(genome, 5000, len(reg_seqs[0].seq))

    # Discover candidates and perform enrichment analysis
    candidate_motifs = enhanced_candidate_discovery(reg_seqs, bg_seqs, min_len, max_len)
    if motifs:
        candidate_motifs += [m.upper() for m in motifs]
    significant = precise_enrichment_with_pseudo(reg_seqs, bg_seqs, candidate_motifs, fdr_threshold)

    # Cluster and print results
    clustered = cluster_motifs(list(significant.keys()))
    print_clustered_results(clustered, significant)


# Example usage
if __name__ == "__main__":
    motif_enrichment_pipeline(
        "methyl_sites.bed",
        "ecoli.fasta",
        motifs=['CCWGG', 'GATC'],  # Known motifs to include
        min_len=4,
        max_len=8,
        fdr_threshold=0.01
    )