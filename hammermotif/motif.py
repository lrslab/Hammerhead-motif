import os
import math
import random
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# =========================================
# 1) BED/FASTA loading & sequence extraction
# =========================================

def parse_bed_file(bed_file_path, flank_size=7):
    """
    Parse a BED file. For each site, add flank_size upstream and downstream.
    Return a list of dict with keys: chrom, start, end, strand.
    """
    locations = []
    with open(bed_file_path) as bed_file:
        for line in bed_file:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.strip().split()
            chrom = parts[0]
            start = int(parts[1])
            end   = int(parts[2])
            strand = parts[3] if len(parts) > 3 else '+'
            start = max(0, start - flank_size)
            end   = end + flank_size
            locations.append({
                'chrom': chrom,
                'start': start,
                'end': end,
                'strand': strand
            })
    return locations

def extract_sequences(genomic_locations, genome_fasta):
    """
    Extract specified regions from the reference genome (FASTA).
    If the strand is '-', take the reverse complement.
    Return a list of SeqRecord objects.
    """
    sequences = []
    for loc in genomic_locations:
        chrom = loc['chrom']
        start = loc['start']
        end   = loc['end']
        strand= loc['strand']
        if chrom not in genome_fasta:
            continue
        seq_region = genome_fasta[chrom].seq[start:end]
        if strand == '-':
            seq_region = seq_region.reverse_complement()
        rec_id = f"{chrom}:{start}-{end}({strand})"
        sequences.append(SeqRecord(seq_region, id=rec_id, description=""))
    return sequences

# =========================================
# 2) Single-motif EM approach
# =========================================

NUC2IDX = {'A':0, 'C':1, 'G':2, 'T':3}
IDX2NUC = {0:'A', 1:'C', 2:'G', 3:'T'}

def dna_to_numeric(seq):
    """
    Convert an ACGT sequence to a list of numeric indices.
    Unknown bases (N, ...) become -1.
    """
    return [NUC2IDX.get(base, -1) for base in seq.upper()]

def build_pwm(aligned_windows, motif_len, pseudocount=0.5):
    """
    Given aligned motif windows (each a numeric list of length motif_len),
    build a PWM matrix of shape (motif_len x 4).
    """
    counts = [[0.0]*4 for _ in range(motif_len)]
    for window in aligned_windows:
        for i, base_idx in enumerate(window):
            if base_idx >= 0:
                counts[i][base_idx] += 1.0

    pwm = []
    for i in range(motif_len):
        row_sum = sum(counts[i]) + 4*pseudocount
        row = [(counts[i][j] + pseudocount)/row_sum for j in range(4)]
        pwm.append(row)
    return pwm

def score_sequence(seq_numeric, pwm):
    """
    Compute log-prob scores for aligning 'pwm' at each possible start in seq_numeric.
    Return a list of log-prob scores, one for each start.
    """
    motif_len = len(pwm)
    seq_len = len(seq_numeric)
    scores = []
    for start in range(seq_len - motif_len + 1):
        logp = 0.0
        for i in range(motif_len):
            base = seq_numeric[start + i]
            if base == -1:  # unknown base
                logp += math.log(0.25)
            else:
                logp += math.log(pwm[i][base])
        scores.append(logp)
    return scores

def compute_log_likelihood(numeric_seqs, positions, pwm):
    """
    Compute total log-likelihood = sum of log p(motif window) over all sequences,
    based on the assigned positions.
    """
    ll = 0.0
    motif_len = len(pwm)
    for seq_idx, s_num in enumerate(numeric_seqs):
        start = positions[seq_idx]
        if start is None:
            continue
        for i in range(motif_len):
            base = s_num[start + i]
            if base == -1:
                ll += math.log(0.25)
            else:
                ll += math.log(pwm[i][base])
    return ll

def run_em_motif_search(seq_list, motif_len, max_iter=20, pseudocount=0.5, n_starts=3):
    """
    Perform multiple random starts of hard-EM for a single motif of fixed length = motif_len.
    Return (best_pwm, best_loglike, best_positions).
    If no valid alignment is found, best_pwm=None, best_loglike=-inf.
    """
    # Convert sequences to numeric
    numeric_seqs = [dna_to_numeric(s) for s in seq_list]

    # Check if enough sequences can accommodate motif_len
    valid_count = sum([1 for s_num in numeric_seqs if len(s_num) >= motif_len])
    if valid_count == 0:
        return None, float('-inf'), None

    best_global_ll = float('-inf')
    best_global_pwm = None
    best_global_positions = None

    for _ in range(n_starts):
        # Random initialization of positions
        positions = []
        for s_num in numeric_seqs:
            if len(s_num) < motif_len:
                positions.append(None)
            else:
                max_start = len(s_num) - motif_len
                rnd_start = random.randint(0, max_start)
                positions.append(rnd_start)

        current_ll = float('-inf')
        current_pwm = None

        for iteration in range(max_iter):
            # M-step: build PWM from the current positions
            aligned = []
            for seq_idx, s_num in enumerate(numeric_seqs):
                start = positions[seq_idx]
                if start is not None:
                    aligned.append(s_num[start:start+motif_len])

            if not aligned:
                # No sequence was assigned => invalid
                current_ll = float('-inf')
                break

            current_pwm = build_pwm(aligned, motif_len, pseudocount)

            # E-step: reassign positions by maximum log-prob
            changed = False
            for seq_idx, s_num in enumerate(numeric_seqs):
                if len(s_num) < motif_len:
                    positions[seq_idx] = None
                    continue
                s_scores = score_sequence(s_num, current_pwm)
                max_val = float('-inf')
                max_pos = None
                for pos_id, val in enumerate(s_scores):
                    if val > max_val:
                        max_val = val
                        max_pos = pos_id
                old_pos = positions[seq_idx]
                positions[seq_idx] = max_pos
                if old_pos != max_pos:
                    changed = True

            new_ll = compute_log_likelihood(numeric_seqs, positions, current_pwm)
            # If no improvement or no changes, consider it converged
            if new_ll <= current_ll or not changed:
                break
            current_ll = new_ll

        # Compare with global best
        if current_ll > best_global_ll:
            best_global_ll = current_ll
            best_global_pwm = current_pwm
            best_global_positions = positions[:]

    return best_global_pwm, best_global_ll, best_global_positions

def iupac_consensus(pwm, threshold=0.5):
    """
    Convert a PWM to an IUPAC consensus sequence.
    If max base prob >= threshold => single base, else collect bases with prob > 0.01
    """
    if pwm is None:
        return ""
    iupac_map = {
        frozenset(['A']): 'A',
        frozenset(['C']): 'C',
        frozenset(['G']): 'G',
        frozenset(['T']): 'T',
        frozenset(['A','G']): 'R',
        frozenset(['C','T']): 'Y',
        frozenset(['G','C']): 'S',
        frozenset(['A','T']): 'W',
        frozenset(['G','T']): 'K',
        frozenset(['A','C']): 'M',
        frozenset(['C','G','T']): 'B',
        frozenset(['A','G','T']): 'D',
        frozenset(['A','C','T']): 'H',
        frozenset(['A','C','G']): 'V',
        frozenset(['A','C','G','T']): 'N'
    }

    consensus = []
    for row in pwm:
        max_p = max(row)
        max_idx = row.index(max_p)
        if max_p >= threshold:
            consensus.append(IDX2NUC[max_idx])
        else:
            base_set = []
            for i, p in enumerate(row):
                if p > 0.01:
                    base_set.append(IDX2NUC[i])
            if not base_set:
                base_set = ['N']
            consensus.append(iupac_map[frozenset(base_set)])
    return "".join(consensus)

# =========================================
# 3) Multi-motif search pipeline
# =========================================

def mask_motif_sites(seq_list, pwm, positions, motif_len, mask_char='N'):
    """
    Given the best alignment (positions) of the discovered motif in each sequence,
    replace those motif sites with 'N' so that the next search won't keep reinforcing the same motif.

    seq_list: list of DNA strings
    pwm: discovered motif PWM
    positions: the assigned start positions for each sequence
    motif_len: length of the motif
    mask_char: default 'N'

    Return a new list of masked sequences.
    """
    new_seqs = []
    for i, seq in enumerate(seq_list):
        start = positions[i]
        if start is None:
            new_seqs.append(seq)
        else:
            seq_array = list(seq)  # convert to list for easy modification
            # mask that region
            for k in range(motif_len):
                seq_array[start + k] = mask_char
            new_seqs.append("".join(seq_array))
    return new_seqs

def multi_motif_search(seq_list,
                       motif_length_list=[5,6,7,8],
                       max_motifs=2,
                       max_iter=20,
                       n_starts=3,
                       pseudocount=0.5):
    """
    Attempt to discover multiple motifs one by one using the "masking" approach:
    1) For each iteration i in [1..max_motifs]:
       - Run single motif search for lengths in motif_length_list, pick the best loglike
       - If no good motif is found, break
       - Otherwise, record that motif, and mask those sites in the sequences
    2) Return the list of found motifs (PWM, consensus, positions).

    seq_list: list of strings (DNA sequences)
    motif_length_list: which motif lengths to try each time
    max_motifs: how many motifs to find
    """
    results = []         # store tuples: (pwm, loglike, positions, iupac_string)
    masked_seq_list = seq_list[:]

    for m_idx in range(max_motifs):
        best_ll_motif = float('-inf')
        best_len = None
        best_pwm = None
        best_positions = None

        # Try each length in motif_length_list
        for L in motif_length_list:
            pwm, ll, positions = run_em_motif_search(
                seq_list=masked_seq_list,
                motif_len=L,
                max_iter=max_iter,
                pseudocount=pseudocount,
                n_starts=n_starts
            )
            if pwm is not None and ll > best_ll_motif:
                best_ll_motif = ll
                best_len = L
                best_pwm = pwm
                best_positions = positions

        # If we couldn't find any valid motif (best_pwm is None or loglike is -inf)
        if best_pwm is None or best_ll_motif == float('-inf'):
            print(f"No further motif found at iteration {m_idx+1}. Stopping.")
            break

        # Convert best PWM to consensus
        consensus = iupac_consensus(best_pwm, threshold=0.5)
        print(f"Motif #{m_idx+1} - length={best_len}, loglike={best_ll_motif:.2f}, consensus={consensus}")

        # Save the result
        results.append((best_pwm, best_ll_motif, best_positions, consensus))

        # Now mask the found motif from the sequences
        masked_seq_list = mask_motif_sites(masked_seq_list, best_pwm, best_positions, best_len)

    return results

# =========================================
# 4) Complete pipeline: from BED+FASTA to multi-motif
# =========================================

def multi_motif_pipeline(bed_file, fasta_file,
                         flank_size=7,
                         motif_length_list=[5,6,7,8],
                         max_motifs=2,
                         max_iter=20,
                         n_starts=3,
                         pseudocount=0.5):
    """
    1) Read the genome from fasta_file
    2) Parse BED and extract sequences with the given flank_size
    3) Run multi_motif_search on these sequences
    4) Print final results
    """
    # Load genome
    genome = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))

    # Extract target sequences from BED
    locs = parse_bed_file(bed_file, flank_size=flank_size)
    reg_seqs = extract_sequences(locs, genome)
    if not reg_seqs:
        print("No valid sequences from BED.")
        return

    # Convert them to plain strings
    seq_list = [str(rec.seq) for rec in reg_seqs]

    # Run multi motif
    results = multi_motif_search(seq_list,
                                 motif_length_list=motif_length_list,
                                 max_motifs=max_motifs,
                                 max_iter=max_iter,
                                 n_starts=n_starts,
                                 pseudocount=pseudocount)

    # Print final summary
    print("\n=== Final Results ===")
    for i, (pwm, ll, positions, consensus) in enumerate(results, start=1):
        print(f"[Motif {i}] Loglike={ll:.2f}, consensus={consensus}")
        # If you want, print PWM or positions details
        # For example:
        print("PWM:")
        for row_idx, row in enumerate(pwm):
            print(f"  Pos{row_idx+1}  A={row[0]:.3f} C={row[1]:.3f} G={row[2]:.3f} T={row[3]:.3f}")

# =========================================
# Usage example
# =========================================

if __name__ == "__main__":
    # You can replace with your own paths
    bed_file_path = "methyl_sites.bed"
    genome_fasta_path = "ecoli.fasta"

    multi_motif_pipeline(
        bed_file = bed_file_path,
        fasta_file = genome_fasta_path,
        flank_size=7,
        motif_length_list=[5,6,7,8],
        max_motifs=2,      # how many motifs you want to discover
        max_iter=30,
        n_starts=5,
        pseudocount=0.5
    )