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
from tqdm import tqdm

# IUPAC nucleotide ambiguity codes
IUPAC_CODES = {
    'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T',
    'R': 'AG', 'Y': 'CT', 'S': 'GC', 'W': 'AT',
    'K': 'GT', 'M': 'AC', 'B': 'CGT', 'D': 'AGT',
    'H': 'ACT', 'V': 'ACG', 'N': 'ACGT'
}


def parse_bed_file(bed_file_path, flank_size=7):
    """Parse BED file to extract genomic regions with flanking sequences."""
    locations = []
    with open(bed_file_path) as bed_file:
        for line in bed_file:
            if line.startswith('#') or not line.strip():
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
    """Extract sequences from genome based on BED coordinates."""
    sequences = []
    for loc in genomic_locations:
        chrom = loc['chrom']
        if chrom not in genome_fasta:
            continue
        seq = genome_fasta[chrom].seq[loc['start']:loc['end']]
        if loc['strand'] == '-':
            seq = seq.reverse_complement()
        sequences.append(SeqRecord(seq, id=f"{chrom}:{loc['start']}-{loc['end']}"))
    return sequences


def generate_random_background(genome_fasta, num_sequences, seq_length):
    """Generate background sequences from all chromosomes ≥1.5x seq_length."""
    valid_chroms = [c for c in genome_fasta if len(genome_fasta[c]) >= 1.5 * seq_length]
    if not valid_chroms:
        raise ValueError("No suitable chromosomes for background generation")

    # Weight selection by chromosome length
    chrom_weights = [len(genome_fasta[c]) for c in valid_chroms]
    total_weight = sum(chrom_weights)

    background = []
    for _ in range(num_sequences):
        # Select chromosome proportional to its length
        selected_chrom = random.choices(valid_chroms, weights=chrom_weights, k=1)[0]
        max_start = len(genome_fasta[selected_chrom]) - seq_length
        start = random.randint(0, max_start)
        seq = genome_fasta[selected_chrom].seq[start:start + seq_length]
        background.append(SeqRecord(seq, id=f"bg_{selected_chrom}:{start}-{start + seq_length}"))
    return background


def degenerate_expansion(kmer):
    """Expand degenerate motifs to all possible non-degenerate forms."""
    return [''.join(p) for p in product(*[IUPAC_CODES[base] for base in kmer])]


def motif_to_regex(motif):
    """Convert IUPAC motif to regex pattern."""
    return ''.join([f'[{IUPAC_CODES[base]}]' if base in IUPAC_CODES else base
                    for base in motif.upper()])


def enhanced_candidate_discovery(reg_seqs, bg_seqs, min_len, max_len):
    """Discover candidate motifs with length-aware scoring."""
    candidates = []

    for k in range(min_len, max_len + 1):
        # Collect all possible k-mers from regulatory sequences
        reg_kmers = Counter()
        for seq in reg_seqs:
            seq_str = str(seq.seq).upper()
            for i in range(len(seq_str) - k + 1):
                kmer = seq_str[i:i + k]
                if all(c in IUPAC_CODES for c in kmer):
                    reg_kmers.update(degenerate_expansion(kmer))

        # Calculate background frequencies
        bg_kmers = Counter()
        for seq in bg_seqs:
            seq_str = str(seq.seq).upper()
            for i in range(len(seq_str) - k + 1):
                kmer = seq_str[i:i + k]
                if all(c in IUPAC_CODES for c in kmer):
                    bg_kmers.update(degenerate_expansion(kmer))

        # Calculate length-normalized enrichment
        enrich_scores = {}
        total_reg = sum(reg_kmers.values()) + 1e-6
        total_bg = sum(bg_kmers.values()) + 1e-6

        for kmer in reg_kmers:
            reg_count = reg_kmers[kmer]
            bg_count = bg_kmers.get(kmer, 0)
            enrich = (reg_count / total_reg) / ((bg_count / total_bg) * np.log1p(k) +0.000001)
            enrich_scores[kmer] = enrich

        # Retain top candidates per length
        sorted_kmers = sorted(enrich_scores.items(), key=lambda x: -x[1])
        candidates += [kmer for kmer, _ in sorted_kmers[:10]]

    return list(set(candidates))


def precise_enrichment_analysis(reg_seqs, bg_seqs, candidates, alpha=0.01):
    """Perform rigorous enrichment analysis with multiple testing correction."""
    results = []

    for motif in tqdm(candidates, desc="Analyzing motifs"):
        k = len(motif)
        pattern = re.compile(motif_to_regex(motif))

        # Count occurrences
        reg_obs = sum(len(pattern.findall(str(seq.seq))) for seq in reg_seqs)
        bg_obs = sum(len(pattern.findall(str(seq.seq))) for seq in bg_seqs)

        # Calculate possible sites
        reg_total = sum(len(seq.seq) - k + 1 for seq in reg_seqs)
        bg_total = sum(len(seq.seq) - k + 1 for seq in bg_seqs)

        # Fisher's exact test with pseudo-counts
        table = [[reg_obs + 0.5, bg_obs + 0.5],
                 [reg_total - reg_obs + 0.5, bg_total - bg_obs + 0.5]]
        _, p = fisher_exact(table, alternative='greater')
        fold = ((reg_obs + 0.5) / reg_total) / ((bg_obs + 0.5) / bg_total)
        results.append((motif, p, fold, reg_obs, bg_obs, reg_total, bg_total))

    # Multiple testing correction
    pvals = [x[1] for x in results]
    _, adj_pvals, _, _ = multipletests(pvals, alpha=alpha, method='fdr_bh')

    # Filter significant results
    significant = {}
    for (motif, p, fold, r_obs, b_obs, r_tot, b_tot), adj_p in zip(results, adj_pvals):
        if adj_p <= alpha and fold > 2:  # Minimum 2-fold enrichment
            significant[motif] = {
                'p': p,
                'adj_p': adj_p,
                'fold': fold,
                'reg_obs': f"{r_obs}/{r_tot}",
                'bg_obs': f"{b_obs}/{b_tot}"
            }
    return significant


def cluster_motifs(motifs, damping=0.85):
    """Cluster motifs using affinity propagation with adjusted parameters."""
    # Create distance matrix
    dist_matrix = np.zeros((len(motifs), len(motifs)))
    for i, m1 in enumerate(motifs):
        for j, m2 in enumerate(motifs):
            dist_matrix[i, j] = levenshtein_distance(m1, m2)

    # Apply clustering with adjusted parameters
    clusterer = AffinityPropagation(
        affinity='precomputed',
        damping=damping,
        max_iter=500,
        random_state=42
    )
    labels = clusterer.fit_predict(-dist_matrix)

    # Organize clusters
    clusters = defaultdict(list)
    for motif, label in zip(motifs, labels):
        clusters[label].append(motif)

    return clusters


def levenshtein_distance(s1, s2):
    """Calculate Levenshtein edit distance between two strings."""
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


def print_clustered_results(clusters, results):
    """Format output for clustered motifs."""
    print("\nSignificant Motif Clusters:")
    for cluster_id, motifs in clusters.items():
        # Find best representative
        rep = min(motifs, key=lambda x: results[x]['p'])
        stats = results[rep]
        print(f"Cluster {cluster_id} (Representative: {rep})")
        print(f"  Motifs: {', '.join(sorted(motifs, key=len, reverse=True))}")
        print(f"  p-value: {stats['p']:.2e}, Adj.p: {stats['adj_p']:.2e}")
        print(f"  Fold: {stats['fold']:.1f}x")
        print(f"  Regions: {stats['reg_obs']} vs Background: {stats['bg_obs']}\n")


def motif_enrichment_pipeline(bed_file, fasta_file, motifs=None,
                              min_len=3, max_len=5, fdr_threshold=0.01):
    """Main analysis workflow."""
    # Load data
    genome = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))
    reg_seqs = extract_sequences(parse_bed_file(bed_file), genome)
    if not reg_seqs:
        raise ValueError("No regulatory sequences extracted")

    # Generate background
    seq_len = len(reg_seqs[0].seq)
    bg_seqs = generate_random_background(genome, 5000, seq_len)

    # Candidate discovery
    candidates = enhanced_candidate_discovery(reg_seqs, bg_seqs, min_len, max_len)
    if motifs:  # Add known motifs if provided
        candidates += [m.upper() for m in motifs]

    # Enrichment analysis
    significant = precise_enrichment_analysis(reg_seqs, bg_seqs, candidates, fdr_threshold)
    if not significant:
        print("No significant motifs found")
        return

    # Clustering and output
    clusters = cluster_motifs(list(significant.keys()), damping=0.9)
    print_clustered_results(clusters, significant)


if __name__ == "__main__":
    # Example usage matching your parameters
    motif_enrichment_pipeline(
        "methyl_sites.bed",
        "ecoli.fasta",
        motifs=None,
        min_len=3,
        max_len=5,
        fdr_threshold=0.01
    )