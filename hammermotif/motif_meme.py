# !/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 13/2/2025 6:26
# @Author  : Runsheng
# @File    : motif_meme.py

# !/usr/bin/env python3
import os
import re
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


import numpy as np
from collections import defaultdict
from scipy.stats import multinomial
import logging
from typing import List, Dict, Tuple, Optional
import random
from multiprocessing import Pool, cpu_count
import itertools
from dataclasses import dataclass
from Bio.Seq import Seq


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
        chrom_len = len(genome[chrom].seq)
        if end > chrom_len:
            end = chrom_len
        seq = genome[chrom].seq[start:end]
        if strand == '-':
            seq = seq.reverse_complement()
        sequences.append(str(seq).upper())
    return sequences


##############################################
# 2. Simplified MEME-style Motif Discovery (EM algorithm)
##############################################

def meme_motif_discovery(seqs, motif_length, max_iter=100, tol=1e-4):
    """
    A simplified MEME-style motif discovery using the EM algorithm.
    Assumes each sequence may contain at most one motif instance.

    Parameters:
      seqs         - List of input sequences.
      motif_length - Desired motif length.
      max_iter     - Maximum number of EM iterations.
      tol          - Convergence threshold for PWM change.

    Returns:
      pwm          - Position Weight Matrix of shape (motif_length, 4) for A, C, G, T.
      consensus    - Consensus sequence derived from the PWM.
    """
    base2idx = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    num_bases = 4

    # Compute background frequency from all input sequences.
    all_seq = "".join(seqs)
    bg_freq = np.array([all_seq.count(b) for b in "ACGT"], dtype=float)
    bg_freq /= bg_freq.sum()

    # Initialize PWM with random values (and normalize rows)
    pwm = np.random.rand(motif_length, num_bases)
    pwm /= pwm.sum(axis=1, keepdims=True)

    def calc_log_score(seq, pwm):
        """
        Calculate log-probability scores for each possible window in a sequence.
        If any base is invalid, score is set to -inf.
        """
        L = len(seq)
        scores = []
        for i in range(L - motif_length + 1):
            logp = 0.0
            valid = True
            for j in range(motif_length):
                base = seq[i + j]
                if base not in base2idx:
                    valid = False
                    break
                logp += np.log(pwm[j, base2idx[base]] + 1e-10)
            scores.append(logp if valid else -np.inf)
        return np.array(scores)

    # EM iterations
    for it in range(max_iter):
        all_weights = []
        for seq in seqs:
            scores = calc_log_score(seq, pwm)
            if np.all(np.isneginf(scores)):
                # If all scores are -inf, assign uniform weights.
                weights = np.ones_like(scores) / len(scores)
            else:
                max_score = np.max(scores)
                exp_scores = np.exp(scores - max_score)
                weights = exp_scores / (exp_scores.sum() + 1e-10)
            all_weights.append(weights)
        new_pwm = np.zeros_like(pwm)
        for seq, weights in zip(seqs, all_weights):
            L = len(seq)
            for i in range(L - motif_length + 1):
                w = weights[i]
                for j in range(motif_length):
                    base = seq[i + j]
                    if base in base2idx:
                        new_pwm[j, base2idx[base]] += w
        new_pwm += 1e-3  # smoothing
        new_pwm /= new_pwm.sum(axis=1, keepdims=True)
        diff = np.abs(new_pwm - pwm).max()
        pwm = new_pwm
        print(f"EM iteration {it + 1}: max PWM diff = {diff:.2e}")
        if diff < tol:
            print(f"EM converged after {it + 1} iterations (diff={diff:.2e}).")
            break
    else:
        print("EM did not converge within maximum iterations.")
    consensus = "".join(["ACGT"[np.argmax(pwm[i])] for i in range(motif_length)])
    return pwm, consensus


##############################################
# 3. Iterative Masking EM Motif Discovery
##############################################

def mask_motif_in_sequences(seqs, motif):
    """
    For each sequence, mask (replace with 'N') all occurrences of the motif.
    Returns a new list of sequences.
    """
    masked_seqs = []
    pattern = re.compile(motif)
    for s in seqs:
        masked_seqs.append(pattern.sub("N" * len(motif), s))
    return masked_seqs


def iterative_meme_motif_discovery(seqs, motif_length, max_motifs=5, max_iter=100, tol=1e-4):
    """
    Iteratively run the EM algorithm to discover motifs.
    After each motif is discovered, mask its occurrences in the sequences.
    Returns a list of tuples (pwm, consensus) for each discovered motif.
    """
    discovered = []
    current_seqs = seqs.copy()
    iteration = 1
    while current_seqs and len(discovered) < max_motifs:
        print(f"EM Iteration {iteration}: {len(current_seqs)} sequences")
        pwm, consensus = meme_motif_discovery(current_seqs, motif_length, max_iter=max_iter, tol=tol)
        if consensus == "" or set(consensus) == {"N"}:
            print("No valid motif found in this iteration.")
            break
        print(f"Discovered motif: {consensus}")
        discovered.append((pwm, consensus))
        current_seqs = mask_motif_in_sequences(current_seqs, consensus)
        iteration += 1
    return discovered


##############################################
# 4. Main Pipeline: Run EM on different motif lengths
##############################################

def main(bed_file, fasta_file,
         motif_length_small=4, motif_length_large=5,
         max_motifs_small=5, max_motifs_large=5,
         max_iter=100, tol=1e-4, flank=7):
    """
    Complete pipeline:
      1. Load the reference genome.
      2. Parse the BED file and extract modified sequences.
      3. Run iterative EM motif discovery for small motif length (e.g., 4) to capture motifs like GATC.
      4. Run iterative EM motif discovery for large motif length (e.g., 5) to capture motifs like CCWGG.
      5. Output both sets of discovered motifs.
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

    # Run iterative EM motif discovery with small motif length.
    print("\n=== Running EM for small motif (length = {}) ===".format(motif_length_small))
    motifs_small = iterative_meme_motif_discovery(mod_seqs, motif_length_small,
                                                  max_motifs=max_motifs_small,
                                                  max_iter=max_iter, tol=tol)
    print("\nFinal discovered small motifs:")
    for i, (pwm, consensus) in enumerate(motifs_small, start=1):
        print(f"Small Motif {i}: {consensus}")
        print("PWM:")
        print(pwm)

    # Run iterative EM motif discovery with large motif length.
    print("\n=== Running EM for large motif (length = {}) ===".format(motif_length_large))
    motifs_large = iterative_meme_motif_discovery(mod_seqs, motif_length_large,
                                                  max_motifs=max_motifs_large,
                                                  max_iter=max_iter, tol=tol)
    print("\nFinal discovered large motifs:")
    for i, (pwm, consensus) in enumerate(motifs_large, start=1):
        print(f"Large Motif {i}: {consensus}")
        print("PWM:")
        print(pwm)


import numpy as np
from collections import defaultdict
from scipy.stats import multinomial
import logging
from typing import List, Dict, Tuple, Optional, Set
import random
from multiprocessing import Pool, cpu_count
import itertools
from dataclasses import dataclass
from Bio.Seq import Seq
from datetime import datetime

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# IUPAC nucleotide codes
IUPAC_CODES = {
    'A': ['A'],
    'C': ['C'],
    'G': ['G'],
    'T': ['T'],
    'R': ['A', 'G'],
    'Y': ['C', 'T'],
    'M': ['A', 'C'],
    'K': ['G', 'T'],
    'S': ['C', 'G'],
    'W': ['A', 'T'],
    'H': ['A', 'C', 'T'],
    'B': ['C', 'G', 'T'],
    'V': ['A', 'C', 'G'],
    'D': ['A', 'G', 'T'],
    'N': ['A', 'C', 'G', 'T']
}

REVERSE_IUPAC = {
    frozenset(['A']): 'A',
    frozenset(['C']): 'C',
    frozenset(['G']): 'G',
    frozenset(['T']): 'T',
    frozenset(['A', 'G']): 'R',
    frozenset(['C', 'T']): 'Y',
    frozenset(['A', 'C']): 'M',
    frozenset(['G', 'T']): 'K',
    frozenset(['C', 'G']): 'S',
    frozenset(['A', 'T']): 'W',
    frozenset(['A', 'C', 'T']): 'H',
    frozenset(['C', 'G', 'T']): 'B',
    frozenset(['A', 'C', 'G']): 'V',
    frozenset(['A', 'G', 'T']): 'D',
    frozenset(['A', 'C', 'G', 'T']): 'N'
}


@dataclass
class Motif:
    """Class to store motif information"""
    consensus: str
    iupac: str
    pwm: np.ndarray
    width: int
    sites: List[str]
    nsites: int
    evalue: float

    def to_dict(self) -> Dict:
        """Convert motif to dictionary format for JSON serialization"""
        return {
            'consensus': self.consensus,
            'iupac': self.iupac,
            'pwm': self.pwm.tolist(),
            'width': self.width,
            'sites': self.sites,
            'nsites': self.nsites,
            'evalue': self.evalue,
            'timestamp': datetime.utcnow().isoformat()
        }

    def __str__(self) -> str:
        return f"Motif(consensus={self.consensus}, iupac={self.iupac}, nsites={self.nsites}, evalue={self.evalue:.2e})"


import numpy as np
from collections import defaultdict
import logging
from typing import List, Dict, Tuple, Optional
import random
from multiprocessing import Pool, cpu_count
from datetime import datetime

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Default parameters
DEFAULT_MIN_WIDTH = 6
DEFAULT_MAX_WIDTH = 15
DEFAULT_BACKGROUND_FREQ = {'A': 0.25, 'C': 0.25, 'G': 0.25, 'T': 0.25}

class MotifEnrichment:
    def __init__(self, min_width: int = DEFAULT_MIN_WIDTH,
                 max_width: int = DEFAULT_MAX_WIDTH,
                 background_freq: Optional[Dict[str, float]] = None,
                 n_processes: int = None):
        """
        Initialize the MotifEnrichment class.

        Args:
            min_width: Minimum width of motifs to search for
            max_width: Maximum width of motifs to search for
            background_freq: Background nucleotide frequencies
            n_processes: Number of processes for parallel processing
        """
        self.min_width = min_width
        self.max_width = max_width
        self.background_freq = background_freq or DEFAULT_BACKGROUND_FREQ.copy()
        self.n_processes = n_processes or max(1, cpu_count() - 1)
        self.motifs = []

        # Validate input parameters
        if min_width > max_width:
            raise ValueError("min_width must be less than or equal to max_width")
        if background_freq and abs(sum(background_freq.values()) - 1.0) > 1e-6:
            raise ValueError("Background frequencies must sum to 1.0")

        logger.info(f"Initialized MotifEnrichment with min_width={min_width}, "
                    f"max_width={max_width}, n_processes={self.n_processes}")

    def create_position_weight_matrix(self, sequences: List[str],
                                      motif_width: int) -> np.ndarray:
        """
        Create a position weight matrix from a set of aligned sequences.
        """
        counts = np.zeros((4, motif_width))
        nuc_to_idx = {'A': 0, 'C': 1, 'G': 2, 'T': 3}

        for seq in sequences:
            for i, nuc in enumerate(seq[:motif_width]):
                if nuc in nuc_to_idx:
                    counts[nuc_to_idx[nuc], i] += 1

        pseudocount = 0.1
        counts += pseudocount
        pwm = counts / counts.sum(axis=0)[np.newaxis, :]
        return pwm

    def calculate_background_frequencies(self, sequences: List[str]) -> Dict[str, float]:
        """Calculate background nucleotide frequencies from input sequences."""
        total_count = 0
        counts = defaultdict(int)

        for seq in sequences:
            for nuc in seq.upper():
                if nuc in 'ACGT':
                    counts[nuc] += 1
                    total_count += 1

        if total_count == 0:
            return self.background_freq.copy()

        return {nuc: count / total_count for nuc, count in counts.items()}

    def score_sequence(self, sequence: str, pwm: np.ndarray) -> float:
        """Score a sequence against a position weight matrix."""
        score = 0.0
        nuc_to_idx = {'A': 0, 'C': 1, 'G': 2, 'T': 3}

        for i, nuc in enumerate(sequence):
            if nuc in nuc_to_idx and i < pwm.shape[1]:
                if pwm[nuc_to_idx[nuc], i] > 0:
                    score += np.log2(pwm[nuc_to_idx[nuc], i] / self.background_freq[nuc])

        return score

    def initialize_sites(self, sequences: List[str],
                         motif_width: int,
                         min_sites: int) -> List[str]:
        """Initialize motif sites by random selection."""
        sites = []
        for seq in sequences:
            if len(seq) >= motif_width:
                possible_starts = len(seq) - motif_width + 1
                if possible_starts > 0:
                    start = random.randint(0, possible_starts - 1)
                    site = seq[start:start + motif_width].upper()
                    if set(site).issubset(set('ACGT')):
                        sites.append(site)

        return sites if len(sites) >= min_sites else []

    def run_em_algorithm(self, sequences: List[str],
                         initial_pwm: np.ndarray,
                         motif_width: int,
                         max_iterations: int = 100,
                         tolerance: float = 1e-6) -> Tuple[np.ndarray, List[str]]:
        """Run the expectation-maximization algorithm to refine the motif model."""
        current_pwm = initial_pwm
        prev_likelihood = float('-inf')

        for iteration in range(max_iterations):
            sites = []
            likelihood = 0

            for seq in sequences:
                if len(seq) >= motif_width:
                    best_score = float('-inf')
                    best_site = None

                    for i in range(len(seq) - motif_width + 1):
                        subseq = seq[i:i + motif_width].upper()
                        if not set(subseq).issubset(set('ACGT')):
                            continue

                        score = self.score_sequence(subseq, current_pwm)
                        if score > best_score:
                            best_score = score
                            best_site = subseq

                    if best_site:
                        sites.append(best_site)
                        likelihood += best_score

            new_pwm = self.create_position_weight_matrix(sites, motif_width)

            if abs(likelihood - prev_likelihood) < tolerance:
                break

            prev_likelihood = likelihood
            current_pwm = new_pwm

        return current_pwm, sites

    def calculate_evalue(self, sites: List[str], pwm: np.ndarray) -> float:
        """Calculate the E-value for a motif."""
        information_content = 0
        for i in range(pwm.shape[1]):
            for j, nuc in enumerate('ACGT'):
                if pwm[j, i] > 0:
                    information_content += pwm[j, i] * np.log2(pwm[j, i] / self.background_freq[nuc])

        return np.exp(-information_content * len(sites))

    def pwm_to_consensus(self, pwm: np.ndarray, threshold: float = 0.7) -> Tuple[str, str]:
        """Convert PWM to consensus and IUPAC sequences."""
        consensus = ""
        iupac = ""
        nucleotides = ['A', 'C', 'G', 'T']

        for i in range(pwm.shape[1]):
            max_idx = np.argmax(pwm[:, i])
            consensus += nucleotides[max_idx]

            significant_nucs = [nuc for nuc, freq in zip(nucleotides, pwm[:, i])
                                if freq >= threshold]
            if not significant_nucs:
                significant_nucs = [nucleotides[max_idx]]

            iupac += REVERSE_IUPAC[frozenset(significant_nucs)]

        return consensus, iupac

    def process_width(self, args: Tuple[int, List[str], int]) -> Optional[Motif]:
        """Process a single motif width (used for parallel processing)."""
        motif_width, sequences, min_sites = args

        sites = self.initialize_sites(sequences, motif_width, min_sites)
        if not sites:
            return None

        initial_pwm = self.create_position_weight_matrix(sites, motif_width)
        pwm, sites = self.run_em_algorithm(sequences, initial_pwm, motif_width)

        if len(sites) >= min_sites:
            consensus, iupac = self.pwm_to_consensus(pwm)
            return Motif(
                consensus=consensus,
                iupac=iupac,
                pwm=pwm,
                width=motif_width,
                sites=sites,
                nsites=len(sites),
                evalue=self.calculate_evalue(sites, pwm)
            )
        return None

    def find_enriched_motifs(self, sequences: List[str],
                             n_motifs: int = 3,
                             min_sites: int = 10) -> List[Motif]:
        """Find enriched motifs in parallel using multiple processes."""
        logger.info(f"Starting motif enrichment analysis with {len(sequences)} sequences")

        self.background_freq = self.calculate_background_frequencies(sequences)
        width_range = range(self.min_width, self.max_width + 1)
        args = [(w, sequences, min_sites) for w in width_range]

        enriched_motifs = []
        with Pool(processes=self.n_processes) as pool:
            enriched_motifs = pool.map(self.process_width, args)

        enriched_motifs = [m for m in enriched_motifs if m is not None]
        enriched_motifs.sort(key=lambda x: x.evalue)

        self.motifs = enriched_motifs[:n_motifs]
        logger.info(f"Found {len(self.motifs)} enriched motifs")
        return self.motifs