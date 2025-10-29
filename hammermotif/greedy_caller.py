#!/usr/bin/env python3
"""
Greedy motif discovery module for bacterial methylation motif calling.

This module implements optimized greedy algorithms for discovering methylation motifs,
including support for gapped and degenerate motifs.
"""

import logging
import re
from collections import Counter, defaultdict
from typing import List, Dict, Tuple, Optional
import numpy as np
from scipy.stats import chi2_contingency, hypergeom
import time
from .motif_greedy import greedy_motif_extraction

logger = logging.getLogger(__name__)


class GreedyMotifCaller:
    """Optimized greedy motif discovery for bacterial methylation sites."""
    
    def __init__(self, genome: Dict[str, str], mod_sites: List[Dict], genome_size: int):
        """
        Initialize the greedy motif caller.
        
        Parameters:
        -----------
        genome : Dict[str, str]
            Dictionary of chromosome sequences
        mod_sites : List[Dict]
            List of modification sites from BED file
        genome_size : int
            Total genome size in base pairs
        """
        self.genome = genome
        self.mod_sites = mod_sites
        self.genome_size = genome_size
        self.greedy_motifs = []
        self.gapped_motifs = []
        self.degenerate_motifs = []
    
    def extract_modified_sequences(self, window_size: int = 20) -> List[str]:
        """Extract sequences around modification sites with optimized performance."""
        sequences = []
        valid_bases = frozenset('ACGT')
        
        for site in self.mod_sites:
            chrom = site['chrom']
            if chrom not in self.genome:
                logger.warning(f"Chromosome {chrom} not found in genome")
                continue
                
            center = (site['start'] + site['end']) // 2
            start = max(0, center - window_size)
            end = min(len(self.genome[chrom]), center + window_size + 1)
            
            sequence = str(self.genome[chrom][start:end]).upper()
            
            # Skip sequences with too many non-ACGT bases
            if sum(1 for b in sequence if b not in valid_bases) > len(sequence) * 0.1:
                continue
            
            # Handle strand orientation
            if site.get('strand') == '-':
                # Simple reverse complement
                complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
                sequence = ''.join(complement.get(b, b) for b in sequence[::-1])
            
            sequences.append(sequence)
        
        logger.info(f"Extracted {len(sequences)} sequences from modification sites")
        return sequences
    
    def run_greedy_discovery(self, k_values: List[int] = None, 
                           chi2_threshold: float = 50,
                           use_parallel: bool = True) -> List[Tuple[str, float]]:
        """
        Run optimized greedy motif discovery algorithm.
        
        Parameters:
        -----------
        k_values : List[int]
            k-mer sizes to search (default: [4, 5, 6, 7, 8])
        chi2_threshold : float
            Chi-square threshold for significance
        use_parallel : bool
            Use parallel processing for multiple k values
        
        Returns:
        --------
        List of (motif, chi2_score) tuples
        """
        if k_values is None:
            k_values = [4, 5, 6, 7, 8]
            
        logger.info("Running optimized greedy motif discovery...")
        
        # Extract sequences once with optimal window size
        mod_seqs = self.extract_modified_sequences(window_size=20)
        
        # Use genome as background
        ref_seqs = [str(seq).upper() for seq in self.genome.values()]
        
        all_motifs = []
        
        if use_parallel and len(k_values) > 1:
            # Parallel processing for multiple k values
            from concurrent.futures import ProcessPoolExecutor, as_completed
            
            with ProcessPoolExecutor(max_workers=min(len(k_values), 4)) as executor:
                futures = {
                    executor.submit(
                        greedy_motif_extraction,
                        mod_seqs, ref_seqs, k=k,
                        chi2_threshold=chi2_threshold,
                        max_motifs=15
                    ): k for k in k_values
                }
                
                for future in as_completed(futures):
                    k = futures[future]
                    try:
                        motifs = future.result()
                        all_motifs.extend(motifs)
                        logger.info(f"Found {len(motifs)} motifs for k={k}")
                    except Exception as e:
                        logger.error(f"Error processing k={k}: {e}")
        else:
            # Sequential processing
            for k in k_values:
                logger.info(f"Searching for {k}-mer motifs...")
                motifs = greedy_motif_extraction(
                    mod_seqs, ref_seqs, k=k,
                    chi2_threshold=chi2_threshold,
                    max_motifs=15
                )
                all_motifs.extend(motifs)
        
        # Remove duplicates and sort by chi2 value
        unique_motifs = {}
        for motif, chi2 in all_motifs:
            if motif not in unique_motifs or chi2 > unique_motifs[motif]:
                unique_motifs[motif] = chi2
        
        all_motifs = [(m, c) for m, c in unique_motifs.items()]
        all_motifs.sort(key=lambda x: x[1], reverse=True)
        
        self.greedy_motifs = all_motifs
        logger.info(f"Found {len(all_motifs)} unique motifs using greedy algorithm")
        
        return all_motifs
    
    def discover_gapped_motifs_optimized(self, max_gap: int = 10,
                                       min_count: int = 10,
                                       use_regex: bool = True) -> List[Tuple[str, float]]:
        """
        Discover motifs with gaps using optimized pattern matching.
        
        Parameters:
        -----------
        max_gap : int
            Maximum gap size (default: 10)
        min_count : int
            Minimum occurrence count (default: 10)
        use_regex : bool
            Use regex for faster pattern matching
        """
        logger.info("Discovering gapped motifs with optimized algorithm...")
        
        mod_seqs = self.extract_modified_sequences(window_size=25)
        ref_seqs = [str(seq).upper() for seq in self.genome.values()]
        
        # Sample reference sequences for speed
        sample_size = min(len(ref_seqs), max(5, len(ref_seqs) // 20))
        ref_sample = ref_seqs[:sample_size]
        scaling_factor = len(ref_seqs) / sample_size
        
        gapped_patterns = []
        
        # Use regex for efficient pattern finding
        if use_regex:
            # Pre-compile regex patterns for common gapped motifs
            pattern_configs = [
                (3, 3, 1, 5),   # XXX_N{1,5}_YYY
                (4, 4, 1, 5),   # XXXX_N{1,5}_YYYY
                (3, 4, 1, 7),   # XXX_N{1,7}_YYYY
                (4, 3, 1, 7),   # XXXX_N{1,7}_YYY
                (5, 3, 1, 10),  # XXXXX_N{1,10}_YYY
            ]
            
            for left_len, right_len, min_gap_len, max_gap_len in pattern_configs:
                # Create regex pattern
                regex_pattern = re.compile(
                    f"([ACGT]{{{left_len}}}).{{{min_gap_len},{max_gap_len}}}([ACGT]{{{right_len}}})"
                )
                
                # Find patterns in modified sequences
                pattern_counts = Counter()
                for seq in mod_seqs:
                    for match in regex_pattern.finditer(seq):
                        left_part = match.group(1)
                        right_part = match.group(2)
                        gap_len = match.end(2) - match.start(1) - left_len - right_len
                        pattern = f"{left_part}{'N' * gap_len}{right_part}"
                        pattern_counts[pattern] += 1
                
                # Filter by minimum count
                for pattern, count in pattern_counts.items():
                    if count >= min_count:
                        # Count in background
                        bg_count = 0
                        pattern_re = re.compile(pattern.replace('N', '.'))
                        for ref_seq in ref_sample:
                            bg_count += len(pattern_re.findall(ref_seq))
                        
                        # Scale background count
                        bg_count = int(bg_count * scaling_factor)
                        
                        # Calculate chi-square
                        if bg_count > 0:
                            chi2 = self._calculate_chi2_optimized(
                                count, len(mod_seqs), bg_count, self.genome_size
                            )
                            if chi2 > 100:  # Higher threshold for gapped motifs
                                gapped_patterns.append((pattern, chi2))
        
        else:
            # Original sliding window approach (slower but more flexible)
            for left_len in [3, 4, 5]:
                for right_len in [3, 4, 5]:
                    for gap_len in range(1, min(max_gap + 1, 11)):
                        pattern_counts = Counter()
                        
                        for seq in mod_seqs:
                            for i in range(len(seq) - left_len - gap_len - right_len + 1):
                                left = seq[i:i+left_len]
                                right = seq[i+left_len+gap_len:i+left_len+gap_len+right_len]
                                
                                if all(b in 'ACGT' for b in left + right):
                                    pattern = f"{left}{'N'*gap_len}{right}"
                                    pattern_counts[pattern] += 1
                        
                        # Process top patterns
                        for pattern, count in pattern_counts.most_common(20):
                            if count < min_count:
                                continue
                            
                            # Count in background (sample)
                            bg_count = 0
                            pattern_re = pattern.replace('N', '.')
                            for ref_seq in ref_sample:
                                bg_count += len(re.findall(pattern_re, ref_seq))
                            
                            # Scale and calculate chi-square
                            bg_count = int(bg_count * scaling_factor)
                            if bg_count > 0:
                                chi2 = self._calculate_chi2_optimized(
                                    count, len(mod_seqs), bg_count, self.genome_size
                                )
                                if chi2 > 100:
                                    gapped_patterns.append((pattern, chi2))
        
        # Sort by chi2 and remove redundant patterns
        gapped_patterns.sort(key=lambda x: x[1], reverse=True)
        
        # Filter out redundant patterns
        filtered_patterns = []
        for pattern, chi2 in gapped_patterns:
            is_redundant = False
            for existing, _ in filtered_patterns:
                if self._is_pattern_redundant(pattern, existing):
                    is_redundant = True
                    break
            if not is_redundant:
                filtered_patterns.append((pattern, chi2))
        
        self.gapped_motifs = filtered_patterns[:10]
        logger.info(f"Found {len(self.gapped_motifs)} gapped motifs")
        
        return self.gapped_motifs
    
    def _calculate_chi2_optimized(self, obs: int, total_obs: int, 
                                 obs_bg: int, total_bg: int) -> float:
        """Optimized chi-square calculation."""
        if total_obs == 0 or total_bg == 0 or obs_bg == 0:
            return 0
        
        # Expected frequency
        expected = (obs + obs_bg) * total_obs / (total_obs + total_bg)
        if expected <= 0:
            return 0
        
        # Chi-square statistic
        chi2 = (obs - expected) ** 2 / expected
        return chi2
    
    def _is_pattern_redundant(self, pattern1: str, pattern2: str) -> bool:
        """Check if two patterns are redundant."""
        # Remove Ns for comparison
        core1 = pattern1.replace('N', '')
        core2 = pattern2.replace('N', '')
        
        # Check if one is substring of another
        if core1 in core2 or core2 in core1:
            return True
        
        # Check if patterns are very similar
        if len(pattern1) == len(pattern2):
            diff = sum(1 for a, b in zip(pattern1, pattern2) 
                      if a != b and a != 'N' and b != 'N')
            if diff <= 1:
                return True
        
        return False
    
    def discover_degenerate_motifs_fast(self, min_length: int = 4, 
                                       min_ic: float = 3.0,
                                       min_sites: int = 10) -> List[Tuple[str, float]]:
        """
        Fast discovery of degenerate motifs using information content filtering.
        """
        logger.info("Fast discovery of degenerate motifs...")
        
        mod_seqs = self.extract_modified_sequences(window_size=30)
        
        # Use position weight matrix approach for speed
        motif_candidates = defaultdict(list)
        
        # Find conserved regions
        for length in range(min_length, min(13, max(len(s) for s in mod_seqs))):
            # Sample sequences for initial pattern discovery
            sample_size = min(len(mod_seqs), 500)
            seq_sample = mod_seqs[:sample_size] if len(mod_seqs) > sample_size else mod_seqs
            
            # Build position frequency matrix
            for i in range(len(seq_sample) - 10):
                for j in range(i + 1, min(i + 50, len(seq_sample))):
                    seq1 = seq_sample[i]
                    seq2 = seq_sample[j]
                    
                    # Find best alignment
                    best_score = 0
                    best_pos1 = 0
                    best_pos2 = 0
                    
                    for pos1 in range(len(seq1) - length + 1):
                        subseq1 = seq1[pos1:pos1 + length]
                        if not all(b in 'ACGT' for b in subseq1):
                            continue
                            
                        for pos2 in range(len(seq2) - length + 1):
                            subseq2 = seq2[pos2:pos2 + length]
                            if not all(b in 'ACGT' for b in subseq2):
                                continue
                            
                            # Calculate similarity
                            score = sum(1 for a, b in zip(subseq1, subseq2) if a == b)
                            if score >= length * 0.7:  # 70% similarity
                                if score > best_score:
                                    best_score = score
                                    best_pos1 = pos1
                                    best_pos2 = pos2
                    
                    # If good alignment found, add to candidates
                    if best_score >= length * 0.7:
                        aligned_seqs = [
                            seq_sample[i][best_pos1:best_pos1 + length],
                            seq_sample[j][best_pos2:best_pos2 + length]
                        ]
                        key = (length, best_score)
                        motif_candidates[key].extend(aligned_seqs)
        
        # Build consensus motifs from candidates
        degenerate_motifs = []
        
        for (length, score), seqs in motif_candidates.items():
            if len(seqs) < min_sites:
                continue
            
            # Build consensus
            consensus = self._build_fast_consensus(seqs)
            
            # Calculate information content
            ic = self._calculate_ic(consensus)
            if ic < min_ic:
                continue
            
            # Count occurrences
            pattern_re = self._pattern_to_regex(consensus)
            count = sum(len(pattern_re.findall(seq)) for seq in mod_seqs)
            
            if count >= min_sites:
                # Calculate chi-square against background
                ref_seqs = [str(seq).upper() for seq in self.genome.values()]
                bg_count = 0
                sample_size = min(3, len(ref_seqs))
                for seq in ref_seqs[:sample_size]:
                    bg_count += len(pattern_re.findall(seq))
                
                bg_count = int(bg_count * len(ref_seqs) / sample_size)
                
                if bg_count > 0:
                    chi2 = self._calculate_chi2_optimized(
                        count, sum(len(s) for s in mod_seqs),
                        bg_count, self.genome_size
                    )
                    
                    if chi2 > 50:
                        degenerate_motifs.append((consensus, chi2))
        
        # Sort and deduplicate
        degenerate_motifs.sort(key=lambda x: x[1], reverse=True)
        self.degenerate_motifs = degenerate_motifs[:20]
        
        logger.info(f"Found {len(self.degenerate_motifs)} degenerate motifs")
        return self.degenerate_motifs
    
    def _build_fast_consensus(self, sequences: List[str]) -> str:
        """Build consensus sequence from aligned sequences."""
        if not sequences:
            return ""
        
        length = len(sequences[0])
        consensus = []
        
        for pos in range(length):
            bases = [seq[pos] for seq in sequences if pos < len(seq)]
            consensus_base = self._get_consensus_base(bases)
            consensus.append(consensus_base)
        
        return ''.join(consensus)
    
    def _get_consensus_base(self, bases: List[str]) -> str:
        """Get IUPAC consensus base."""
        if not bases:
            return 'N'
        
        base_counts = Counter(bases)
        total = len(bases)
        
        # If one base is dominant (>80%), use it
        for base, count in base_counts.items():
            if count / total > 0.8:
                return base
        
        # Use degenerate codes
        unique_bases = set(bases)
        
        if len(unique_bases) == 2:
            bases_sorted = sorted(unique_bases)
            degenerate_map = {
                ('A', 'G'): 'R', ('C', 'T'): 'Y',
                ('A', 'T'): 'W', ('C', 'G'): 'S',
                ('G', 'T'): 'K', ('A', 'C'): 'M'
            }
            return degenerate_map.get(tuple(bases_sorted), 'N')
        
        elif len(unique_bases) == 3:
            missing = {'A', 'C', 'G', 'T'} - unique_bases
            if missing == {'A'}:
                return 'B'
            elif missing == {'C'}:
                return 'D'
            elif missing == {'G'}:
                return 'H'
            elif missing == {'T'}:
                return 'V'
        
        return 'N'
    
    def _calculate_ic(self, motif: str) -> float:
        """Calculate information content of a motif."""
        ic_values = {
            'A': 1.0, 'T': 1.0, 'C': 1.0, 'G': 1.0,
            'R': 0.5, 'Y': 0.5, 'W': 0.5, 'S': 0.5, 'K': 0.5, 'M': 0.5,
            'B': 1/3, 'D': 1/3, 'H': 1/3, 'V': 1/3,
            'N': 0.0
        }
        return sum(ic_values.get(base.upper(), 0) for base in motif)
    
    def _pattern_to_regex(self, pattern: str) -> re.Pattern:
        """Convert IUPAC pattern to regex."""
        iupac_map = {
            'R': '[AG]', 'Y': '[CT]', 'W': '[AT]', 'S': '[GC]',
            'K': '[GT]', 'M': '[AC]', 'B': '[CGT]', 'D': '[AGT]',
            'H': '[ACT]', 'V': '[ACG]', 'N': '.'
        }
        
        regex_pattern = ''
        for base in pattern.upper():
            if base in 'ACGT':
                regex_pattern += base
            else:
                regex_pattern += iupac_map.get(base, '.')
        
        return re.compile(regex_pattern)
    
    def get_all_motifs(self) -> Dict[str, List[Tuple[str, float]]]:
        """Get all discovered motifs organized by type."""
        return {
            'greedy': self.greedy_motifs,
            'gapped': self.gapped_motifs,
            'degenerate': self.degenerate_motifs
        } 