#!/usr/bin/env python3
"""
Motif merger module for bacterial methylation motif calling.

This module handles merging motifs from different discovery methods and
creating consensus sequences.
"""

import logging
import re
from collections import Counter, defaultdict
from typing import List, Dict, Tuple, Optional, Set
import numpy as np

logger = logging.getLogger(__name__)


class MotifMerger:
    """Merge and consolidate motifs from different discovery methods."""
    
    def __init__(self):
        """Initialize the motif merger."""
        self.iupac_compatible = self._get_iupac_compatible()
        self.final_motifs = []
    
    def _get_iupac_compatible(self) -> Dict[str, Set[str]]:
        """Get IUPAC base compatibility mapping."""
        return {
            'A': set('ARWMDHVN'),
            'C': set('CYSMHVBN'),
            'G': set('GRSKDVBN'),
            'T': set('YTWKDHBN'),
            'R': set('AGRN'),
            'Y': set('CTYN'),
            'W': set('ATWN'),
            'S': set('CGSN'),
            'K': set('GTKN'),
            'M': set('ACMN'),
            'B': set('CGTYSKBN'),
            'D': set('AGTRWKDN'),
            'H': set('ACTYWMHN'),
            'V': set('ACGRSMVN'),
            'N': set('ACGTRYWSKMDHVBN')
        }
    
    def merge_motifs(self, 
                    greedy_motifs: List[Tuple[str, float]],
                    meme_motifs: List[Dict],
                    chi2_threshold: float = 30,
                    evalue_threshold: float = 0.01) -> List[str]:
        """
        Score-aware motif merging that preserves high-scoring motifs.
        Now includes substring-based merging with score priority.
        """
        logger.info("Merging motifs from different methods...")
        
        # Collect motifs with their scores
        scored_motifs = []
        
        # Add significant greedy motifs with scores
        for motif, chi2 in greedy_motifs:
            if chi2 >= chi2_threshold:
                scored_motifs.append((motif, chi2))
        
        # Add significant MEME motifs
        for motif_dict in meme_motifs:
            if motif_dict.get('evalue', 1.0) <= evalue_threshold:
                score = -np.log10(motif_dict.get('evalue', 1e-10))
                
                if 'degenerate_consensus' in motif_dict:
                    scored_motifs.append((motif_dict['degenerate_consensus'], score))
                elif 'consensus' in motif_dict:
                    scored_motifs.append((motif_dict['consensus'], score))
                elif 'regex' in motif_dict:
                    iupac = self._regex_to_iupac(motif_dict['regex'])
                    if iupac:
                        scored_motifs.append((iupac, score))
        
        # New: Apply substring-based merging with score priority
        merged_motifs = self._substring_aware_merging(scored_motifs)
        
        self.final_motifs = merged_motifs
        logger.info(f"Final merged motifs: {merged_motifs}")
        
        return merged_motifs
    
    def _group_non_redundant_motifs(self, motifs: List[str], 
                                   motifs_with_scores: List[Tuple[str, float, str]]) -> List[str]:
        """Group similar motifs but preserve high-scoring longer motifs."""
        # Create score lookup
        score_map = {motif: score for motif, score, _ in motifs_with_scores}
        
        # Preserve motifs with high scores and length >= 6
        protected_motifs = []
        remaining_motifs = []
        
        for motif in motifs:
            score = score_map.get(motif, 0)
            if len(motif) >= 6 and score > 1000:  # High-scoring longer motifs
                protected_motifs.append(motif)
            else:
                remaining_motifs.append(motif)
        
        # Group only the remaining shorter/lower-scoring motifs
        if remaining_motifs:
            grouped_remaining = self._group_and_merge_similar(remaining_motifs)
        else:
            grouped_remaining = []
        
        # Combine protected and grouped motifs
        final_result = protected_motifs + grouped_remaining
        
        # Remove any duplicates while preserving order
        seen = set()
        unique_result = []
        for motif in final_result:
            if motif not in seen:
                unique_result.append(motif)
                seen.add(motif)
        
        return unique_result
    
    def _group_and_merge_similar(self, motifs: List[str]) -> List[str]:
        """Group similar motifs and create consensus."""
        if not motifs:
            return []
        
        # Group similar motifs
        groups = []
        used = set()
        
        for i, motif1 in enumerate(motifs):
            if i in used:
                continue
                
            group = [motif1]
            used.add(i)
            
            for j, motif2 in enumerate(motifs[i+1:], i+1):
                if j in used:
                    continue
                
                if self._are_motifs_similar(motif1, motif2):
                    group.append(motif2)
                    used.add(j)
            
            groups.append(group)
        
        # Create consensus for each group
        merged = []
        for group in groups:
            if len(group) == 1:
                merged.append(group[0])
            else:
                consensus = self._create_group_consensus(group)
                merged.append(consensus)
        
        return merged
    
    def _are_motifs_similar(self, motif1: str, motif2: str) -> bool:
        """Check if two motifs are similar."""
        motif1 = motif1.upper()
        motif2 = motif2.upper()
        
        # Quick substring check
        if motif1 in motif2 or motif2 in motif1:
            return True
        
        # Check for common core motifs (length >= 4)
        if self._have_common_core(motif1, motif2):
            return True
        
        # Check circular permutations for short motifs
        if len(motif1) == len(motif2) and len(motif1) <= 10:
            for i in range(len(motif1)):
                if motif1[i:] + motif1[:i] == motif2:
                    return True
        
        # Handle gapped motifs
        if 'N' in motif1 or 'N' in motif2:
            return self._compare_gapped_motifs(motif1, motif2)
        
        # Check similarity for similar length motifs (lower threshold)
        len_diff = abs(len(motif1) - len(motif2))
        if len_diff <= 2:
            score = self._calculate_similarity_score(motif1, motif2)
            if score >= 0.6:  # Lowered from 0.75
                return True
        
        # Check sliding window for different lengths (lower threshold)
        if len_diff > 2:
            shorter, longer = (motif1, motif2) if len(motif1) < len(motif2) else (motif2, motif1)
            if len(shorter) >= 4:
                for i in range(len(longer) - len(shorter) + 1):
                    window = longer[i:i + len(shorter)]
                    score = self._calculate_similarity_score(shorter, window)
                    if score >= 0.75:  # Lowered from 0.85
                        return True
        
        return False
    
    def _have_common_core(self, motif1: str, motif2: str) -> bool:
        """Check if two motifs share a common core of at least 4 bases."""
        # Find all substrings of length 4+ in both motifs
        min_core_length = 4
        max_len = min(len(motif1), len(motif2))
        
        for core_len in range(min_core_length, max_len + 1):
            # Get all substrings of this length from motif1
            for i in range(len(motif1) - core_len + 1):
                substr1 = motif1[i:i + core_len]
                # Skip if contains too many degenerate bases
                if substr1.count('N') + sum(substr1.count(x) for x in 'RYWSKMDHVB') > core_len // 2:
                    continue
                
                # Check if this substring appears in motif2
                if substr1 in motif2:
                    return True
                
                # Check with one mismatch allowed for longer cores
                if core_len >= 6:
                    for j in range(len(motif2) - core_len + 1):
                        substr2 = motif2[j:j + core_len]
                        mismatches = sum(1 for a, b in zip(substr1, substr2) if a != b)
                        if mismatches <= 1:
                            return True
        
        return False
    
    def _calculate_similarity_score(self, motif1: str, motif2: str) -> float:
        """Calculate similarity score between two motifs."""
        if len(motif1) != len(motif2):
            # Pad shorter motif with N's
            if len(motif1) < len(motif2):
                motif1 = motif1 + 'N' * (len(motif2) - len(motif1))
            else:
                motif2 = motif2 + 'N' * (len(motif1) - len(motif2))
        
        matches = 0
        for base1, base2 in zip(motif1, motif2):
            if base1 == base2:
                matches += 1.0
            elif base2 in self.iupac_compatible.get(base1, set()) or \
                 base1 in self.iupac_compatible.get(base2, set()):
                matches += 0.8
            elif base1 == 'N' or base2 == 'N':
                matches += 0.5
        
        return matches / len(motif1)
    
    def _compare_gapped_motifs(self, motif1: str, motif2: str) -> bool:
        """Compare motifs containing gaps."""
        # Extract non-N positions
        core1 = [(i, b) for i, b in enumerate(motif1) if b != 'N']
        core2 = [(i, b) for i, b in enumerate(motif2) if b != 'N']
        
        if abs(len(core1) - len(core2)) > 2:
            return False
        
        # Compare cores
        core_seq1 = ''.join(b for _, b in core1)
        core_seq2 = ''.join(b for _, b in core2)
        
        if core_seq1 in core_seq2 or core_seq2 in core_seq1:
            return True
        
        # Check similarity
        score = self._calculate_similarity_score(core_seq1, core_seq2)
        return score >= 0.75
    
    def _create_group_consensus(self, group: List[str]) -> str:
        """Create consensus from a group of similar motifs."""
        if len(group) == 1:
            return group[0]
        
        # Convert all to uppercase
        group = [m.upper() for m in group]
        
        # Handle gapped motifs specially
        if any('N' in m for m in group):
            return self._create_gapped_consensus(group)
        
        # First, try to find common core and build around it
        core = self._find_common_core(group)
        if core and len(core) >= 4:
            # Extend the core with flanking regions
            extended = self._extend_core_with_flanks(core, group)
            if extended:
                return extended
        
        # Fallback to alignment-based consensus
        aligned = self._align_motifs(group)
        return self._build_consensus(aligned)
    
    def _find_common_core(self, motifs: List[str]) -> Optional[str]:
        """Find the longest common substring among all motifs."""
        if not motifs:
            return None
        
        reference = motifs[0]
        common_substrings = set()
        
        # Find all substrings of length 4+ in reference
        for length in range(4, len(reference) + 1):
            for start in range(len(reference) - length + 1):
                substring = reference[start:start + length]
                # Check if this substring appears in all other motifs
                if all(substring in motif for motif in motifs[1:]):
                    common_substrings.add(substring)
        
        # Return the longest common substring
        if common_substrings:
            return max(common_substrings, key=len)
        
        return None
    
    def _extend_core_with_flanks(self, core: str, motifs: List[str]) -> Optional[str]:
        """Extend a core motif with common flanking regions."""
        # Find positions of core in each motif
        core_positions = []
        for motif in motifs:
            pos = motif.find(core)
            if pos != -1:
                core_positions.append((motif, pos))
            else:
                return None  # Core not found in all motifs
        
        # Determine common left and right flanks
        max_left = min(pos for _, pos in core_positions)
        max_right = min(len(motif) - pos - len(core) for motif, pos in core_positions)
        
        # Build consensus with flanks
        extended = ""
        
        # Left flank
        for i in range(max_left):
            left_pos = max_left - 1 - i
            bases = [motif[pos - left_pos - 1] for motif, pos in core_positions if pos > left_pos]
            if bases:
                consensus_base = self._get_consensus_base(bases)
                extended = consensus_base + extended
        
        # Core
        extended += core
        
        # Right flank
        for i in range(max_right):
            bases = [motif[pos + len(core) + i] for motif, pos in core_positions 
                    if pos + len(core) + i < len(motif)]
            if bases:
                consensus_base = self._get_consensus_base(bases)
                extended += consensus_base
        
        return extended if len(extended) >= len(core) else core
    
    def _align_motifs(self, motifs: List[str]) -> List[str]:
        """Align motifs for consensus building."""
        if not motifs:
            return []
        
        # Find the best reference (most similar to all others)
        best_ref = motifs[0]
        best_score = 0
        
        for candidate in motifs:
            total_score = sum(
                self._calculate_similarity_score(candidate, other)
                for other in motifs if other != candidate
            )
            if total_score > best_score:
                best_score = total_score
                best_ref = candidate
        
        # Align all to reference
        aligned = []
        for motif in motifs:
            if motif == best_ref:
                aligned.append(motif)
            else:
                aligned_motif = self._align_to_reference(best_ref, motif)
                aligned.append(aligned_motif)
        
        return aligned
    
    def _align_to_reference(self, reference: str, motif: str) -> str:
        """Align a motif to reference."""
        if len(motif) == len(reference):
            return motif
        
        if len(motif) < len(reference):
            # Find best position to align
            best_pos = 0
            best_score = 0
            
            for pos in range(len(reference) - len(motif) + 1):
                padded = 'N' * pos + motif + 'N' * (len(reference) - len(motif) - pos)
                score = self._calculate_similarity_score(reference, padded)
                if score > best_score:
                    best_score = score
                    best_pos = pos
            
            return 'N' * best_pos + motif + 'N' * (len(reference) - len(motif) - best_pos)
        
        else:
            # Motif is longer - find best substring
            best_start = 0
            best_score = 0
            
            for start in range(len(motif) - len(reference) + 1):
                substring = motif[start:start + len(reference)]
                score = self._calculate_similarity_score(reference, substring)
                if score > best_score:
                    best_score = score
                    best_start = start
            
            return motif[best_start:best_start + len(reference)]
    
    def _build_consensus(self, aligned_motifs: List[str]) -> str:
        """Build consensus from aligned motifs."""
        if not aligned_motifs:
            return ""
        
        consensus = []
        length = max(len(m) for m in aligned_motifs)
        
        for pos in range(length):
            bases = []
            for motif in aligned_motifs:
                if pos < len(motif) and motif[pos] != 'N':
                    bases.append(motif[pos])
            
            if bases:
                consensus_base = self._get_consensus_base(bases)
                consensus.append(consensus_base)
            else:
                consensus.append('N')
        
        # Clean up - remove leading/trailing N's
        consensus_str = ''.join(consensus).strip('N')
        
        return consensus_str if consensus_str else aligned_motifs[0]
    
    def _create_gapped_consensus(self, group: List[str]) -> str:
        """Create consensus for gapped motifs."""
        # Find most common gap pattern
        gap_patterns = Counter()
        
        for motif in group:
            pattern = tuple(i for i, b in enumerate(motif) if b == 'N')
            gap_patterns[pattern] += 1
        
        if gap_patterns:
            # Use most common pattern
            common_pattern = gap_patterns.most_common(1)[0][0]
            
            # Filter motifs with this pattern
            matching = [m for m in group if 
                       tuple(i for i, b in enumerate(m) if b == 'N') == common_pattern]
            
            if matching:
                return self._build_consensus(matching)
        
        # Fallback to motif with highest IC
        return max(group, key=self._calculate_ic)
    
    def _get_consensus_base(self, bases: List[str]) -> str:
        """Get IUPAC consensus base from list of bases."""
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
            bases_sorted = tuple(sorted(unique_bases))
            degenerate_map = {
                ('A', 'G'): 'R', ('C', 'T'): 'Y',
                ('A', 'T'): 'W', ('C', 'G'): 'S',
                ('G', 'T'): 'K', ('A', 'C'): 'M'
            }
            return degenerate_map.get(bases_sorted, 'N')
        
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
    
    def _regex_to_iupac(self, regex: str) -> Optional[str]:
        """Convert regex pattern to IUPAC notation."""
        conversions = {
            '[AG]': 'R', '[GA]': 'R',
            '[CT]': 'Y', '[TC]': 'Y',
            '[AT]': 'W', '[TA]': 'W',
            '[GC]': 'S', '[CG]': 'S',
            '[GT]': 'K', '[TG]': 'K',
            '[AC]': 'M', '[CA]': 'M',
            '.': 'N'
        }
        
        result = regex
        for pattern, iupac in conversions.items():
            result = result.replace(pattern, iupac)
        
        # Check if fully converted
        if re.match(r'^[ACGTRYWSKMDHVBN]+$', result):
            return result
        
        return None
    
    def _simple_core_clustering(self, motifs: List[str]) -> List[str]:
        """
        Simple clustering based on identifying core motifs.
        
        1. Identify high-quality core motifs (length 6+, good IC)
        2. Group other motifs that contain these cores
        3. Handle remaining motifs with similarity
        """
        if not motifs:
            return []
        
        # Sort by length and information content
        motifs.sort(key=lambda m: (len(m), self._calculate_ic(m)), reverse=True)
        
        # Step 1: Identify core motifs
        cores = []
        remaining = []
        
        for motif in motifs:
            # Good cores: length 6+, high information content, no degenerates
            if (len(motif) >= 6 and 
                self._calculate_ic(motif) >= 4.0 and
                all(base in 'ACGT' for base in motif)):
                cores.append(motif)
            else:
                remaining.append(motif)
        
        # Step 2: Group remaining motifs with cores
        used_motifs = set()
        final_motifs = []
        
        for core in cores:
            if core in used_motifs:
                continue
            
            # Check if any other core contains this core or vice versa
            best_core = core
            for other_core in cores:
                if other_core != core and other_core not in used_motifs:
                    # If this core is contained in another core, use the larger one
                    if core in other_core:
                        best_core = other_core
                        break
                    # If another core is contained in this core, this is better
                    elif other_core in core:
                        used_motifs.add(other_core)
            
            # Group remaining motifs that contain this core
            for motif in remaining:
                if motif not in used_motifs and best_core in motif:
                    used_motifs.add(motif)
            
            final_motifs.append(best_core)
            used_motifs.add(best_core)
        
        # Step 3: Handle ungrouped motifs with similarity clustering
        ungrouped = [m for m in remaining if m not in used_motifs]
        
        if ungrouped:
            # Simple similarity grouping for remaining motifs
            similarity_groups = self._group_and_merge_similar(ungrouped)
            final_motifs.extend(similarity_groups)
        
        return final_motifs
    
    def _smart_motif_clustering(self, motifs: List[str]) -> List[str]:
        """
        Smart motif clustering that identifies core motifs and merges variants.
        
        Strategy:
        1. Find the best representative motif for each cluster
        2. Remove redundant motifs that are substrings of better motifs
        3. Group similar remaining motifs
        """
        if not motifs:
            return []
        
        # Step 1: Remove substring redundancy - prefer longer, more informative motifs
        filtered_motifs = self._remove_substring_redundancy(motifs)
        
        # Step 2: Group remaining similar motifs
        final_motifs = self._group_and_merge_similar(filtered_motifs)
        
        return final_motifs
    
    def _remove_substring_redundancy(self, motifs: List[str]) -> List[str]:
        """
        Remove motifs that are substrings of longer, more significant motifs.
        Keep the most informative representatives.
        """
        # Sort by length descending, then by information content
        motifs_sorted = sorted(motifs, key=lambda m: (len(m), self._calculate_ic(m)), reverse=True)
        
        filtered = []
        
        for motif in motifs_sorted:
            # Check if this motif is a substring of any already selected motif
            is_redundant = False
            
            for selected in filtered:
                # If current motif is a substring of a selected longer motif, skip it
                if motif in selected:
                    is_redundant = True
                    break
                # If selected motif is a substring of current motif, replace it
                elif selected in motif:
                    filtered.remove(selected)
                    break
            
            if not is_redundant:
                filtered.append(motif)
        
        return filtered
    
    def _score_based_merging(self, motif_scores: List[Tuple[str, float]]) -> List[str]:
        """
        Merge motifs based on scores and similarity.
        
        Strategy:
        1. Identify high-scoring motifs as anchors
        2. Group related motifs around these anchors
        3. Create consensus for similar motifs without clear anchors
        """
        if not motif_scores:
            return []
        
        # Sort by score descending
        motif_scores.sort(key=lambda x: x[1], reverse=True)
        
        # Step 1: Identify anchor motifs (high score and good length)
        anchors = []
        remaining = []
        
        for motif, score in motif_scores:
            # High-scoring motifs with length 6+ become anchors
            if score > 1000 and len(motif) >= 6:
                anchors.append((motif, score))
            # Also consider very high-scoring shorter motifs
            elif score > 1500 and len(motif) >= 4:
                anchors.append((motif, score))
            else:
                remaining.append((motif, score))
        
        # Step 2: Group remaining motifs with anchors
        final_motifs = []
        used_motifs = set()
        
        for anchor_motif, anchor_score in anchors:
            if anchor_motif in used_motifs:
                continue
            
            # Find motifs that should be merged with this anchor
            group = [anchor_motif]
            
            for motif, score in remaining:
                if motif in used_motifs:
                    continue
                
                # Check if motif contains the anchor or vice versa
                if anchor_motif in motif or motif in anchor_motif:
                    group.append(motif)
                    used_motifs.add(motif)
                # Check if motifs are very similar
                elif self._are_motifs_similar(anchor_motif, motif):
                    group.append(motif)
                    used_motifs.add(motif)
            
            # Use the anchor as representative (it has the highest score)
            final_motifs.append(anchor_motif)
            used_motifs.add(anchor_motif)
        
        # Step 3: Handle remaining ungrouped motifs
        ungrouped = [(m, s) for m, s in remaining if m not in used_motifs]
        
        if ungrouped:
            # Group similar ungrouped motifs
            ungrouped_motifs = [m for m, s in ungrouped]
            grouped = self._group_and_merge_similar(ungrouped_motifs)
            final_motifs.extend(grouped)
        
        return final_motifs
    

    
    def _is_very_similar_or_contained(self, motif1: str, motif2: str) -> bool:
        """Check if motifs are very similar or one contains the other completely."""
        # Exact match
        if motif1 == motif2:
            return True
        
        # Check for core substring relationships
        # If one motif is a core substring of another (length >= 6)
        if len(motif1) >= 6 and motif1 in motif2:
            return True
        elif len(motif2) >= 6 and motif2 in motif1:
            return True
        
        # For shorter motifs, require more difference in length
        if motif1 in motif2 and len(motif1) < len(motif2) - 2:
            return True
        elif motif2 in motif1 and len(motif2) < len(motif1) - 2:
            return True
        
        # Check for shared core (e.g., CAGAAG in CAGAAGTA and CAGAAGG)
        if len(motif1) >= 6 and len(motif2) >= 6:
            # Find longest common substring
            for length in range(min(len(motif1), len(motif2)), 5, -1):
                for i in range(len(motif1) - length + 1):
                    substr = motif1[i:i + length]
                    if substr in motif2 and length >= 6:
                        return True
        
        # Very high similarity for same-length motifs
        if len(motif1) == len(motif2):
            diff_count = sum(1 for a, b in zip(motif1, motif2) if a != b)
            if diff_count <= 1:  # At most 1 difference
                return True
        
        return False 

    def _substring_aware_merging(self, scored_motifs: List[Tuple[str, float]]) -> List[str]:
        """
        Merge motifs where shorter, higher-scoring motifs take priority over
        longer motifs that contain them.
        
        Example: GATC (chi2=4587) vs CGATC (chi2=2916) -> keep GATC
        """
        if not scored_motifs:
            return []
        
        # Remove duplicates, keep highest score for each motif
        motif_scores = {}
        for motif, score in scored_motifs:
            if motif not in motif_scores or score > motif_scores[motif]:
                motif_scores[motif] = score
        
        # Sort by score descending
        sorted_motifs = sorted(motif_scores.items(), key=lambda x: x[1], reverse=True)
        
        # Separate standard and gapped motifs for different treatment
        standard_motifs = [(m, s) for m, s in sorted_motifs if 'N' not in m]
        gapped_motifs = [(m, s) for m, s in sorted_motifs if 'N' in m]
        
        logger.info("Performing substring-aware merging...")
        logger.info(f"Input: {len(standard_motifs)} standard, {len(gapped_motifs)} gapped motifs")
        
        # Process standard motifs with existing logic
        final_standard = self._merge_standard_motifs(standard_motifs, motif_scores)
        
        # Process gapped motifs with more aggressive merging
        final_gapped = self._merge_gapped_motifs_aggressively(gapped_motifs, motif_scores)
        
        # Combine and limit total number more conservatively
        final_motifs = final_standard + final_gapped
        
        # Only limit if we have too many motifs (more flexible threshold)
        if len(final_motifs) > 8:
            # Sort by score and keep top motifs, but be more generous
            motif_scores_final = {m: motif_scores[m] for m in final_motifs}
            final_sorted = sorted(motif_scores_final.items(), key=lambda x: x[1], reverse=True)
            final_motifs = [m for m, s in final_sorted[:8]]
            logger.info(f"Limited final motifs to top 8 by score")
        
        logger.info(f"Final: {len(final_motifs)} motifs")
        return final_motifs
    
    def _merge_standard_motifs(self, standard_motifs: List[Tuple[str, float]], 
                              motif_scores: Dict[str, float]) -> List[str]:
        """Merge standard motifs using conservative substring logic."""
        final_motifs = []
        used_motifs = set()
        
        for motif, score in standard_motifs:
            if motif in used_motifs:
                continue
            
            # Check if this motif should be the representative for a group
            motifs_to_merge = [motif]
            should_skip_current = False
            
            # Find all other motifs that contain this motif or are contained by this motif
            for other_motif, other_score in standard_motifs:
                if other_motif == motif or other_motif in used_motifs:
                    continue
                
                # Case 1: Other motif contains this motif (e.g., CGGATCC contains GATC)
                if motif in other_motif:
                    # Be much more conservative - only merge if score difference is very significant
                    score_ratio = other_score / score if score > 0 else 1
                    length_diff = len(other_motif) - len(motif)
                    
                    # Only merge if longer motif has significantly higher score relative to length increase
                    if score_ratio > (1.5 + 0.1 * length_diff):  # More stringent threshold
                        # Longer motif has much higher score, skip current motif
                        used_motifs.add(motif)
                        should_skip_current = True
                        logger.info(f"  Skipping {motif}, keeping {other_motif} (score ratio {score_ratio:.2f})")
                        break
                    elif score >= other_score:
                        # Current motif (shorter) has equal or higher score, merge longer into shorter
                        motifs_to_merge.append(other_motif)
                        used_motifs.add(other_motif)
                        logger.info(f"  Merging {other_motif} -> {motif} (preserving shorter: {score:.0f} >= {other_score:.0f})")
                
                # Case 2: This motif contains other motif (e.g., CGGATCC contains GATC)
                elif other_motif in motif:
                    # Be conservative about merging shorter motifs into longer ones
                    score_ratio = score / other_score if other_score > 0 else 1
                    length_diff = len(motif) - len(other_motif)
                    
                    # Only merge if current (longer) motif has significantly higher score
                    if score_ratio > (1.5 + 0.1 * length_diff):
                        # Current motif (longer) has much higher score, merge shorter
                        motifs_to_merge.append(other_motif)
                        used_motifs.add(other_motif)
                        logger.info(f"  Merging {other_motif} -> {motif} (score ratio {score_ratio:.2f})")
                    else:
                        # Shorter motif should be preserved separately, skip current
                        used_motifs.add(motif)
                        should_skip_current = True
                        logger.info(f"  Preserving shorter motif {other_motif}, skipping {motif}")
                        break
            
            if not should_skip_current and motif not in used_motifs:
                final_motifs.append(motif)
                used_motifs.add(motif)
                if len(motifs_to_merge) > 1:
                    logger.info(f"  Selected representative: {motif} (merged {len(motifs_to_merge)} motifs)")
        
        return final_motifs
    
    def _merge_gapped_motifs_aggressively(self, gapped_motifs: List[Tuple[str, float]], 
                                         motif_scores: Dict[str, float]) -> List[str]:
        """
        Merge gapped motifs more aggressively to reduce redundancy.
        Only keep the top 2-3 most significant gapped motifs.
        """
        if not gapped_motifs:
            return []
        
        logger.info("Applying aggressive gapped motif merging...")
        
        # Apply higher threshold for gapped motifs
        high_score_gapped = [(m, s) for m, s in gapped_motifs if s > 1000]
        
        if not high_score_gapped:
            # If no high-scoring gapped motifs, take top 2
            high_score_gapped = gapped_motifs[:2]
        
        # Group similar gapped motifs by pattern similarity
        grouped_gapped = self._group_gapped_by_pattern(high_score_gapped)
        
        # Take best representative from each group, max 3 groups
        final_gapped = []
        for group in grouped_gapped[:3]:  # Max 3 gapped motif groups
            # Choose highest scoring motif from group
            best_motif = max(group, key=lambda x: motif_scores[x[0]])
            final_gapped.append(best_motif[0])
            logger.info(f"  Selected gapped representative: {best_motif[0]} (score {best_motif[1]:.0f})")
        
        return final_gapped
    
    def _group_gapped_by_pattern(self, gapped_motifs: List[Tuple[str, float]]) -> List[List[Tuple[str, float]]]:
        """Group gapped motifs by similar patterns."""
        groups = []
        used = set()
        
        for i, (motif1, score1) in enumerate(gapped_motifs):
            if i in used:
                continue
            
            group = [(motif1, score1)]
            used.add(i)
            
            # Find similar gapped motifs
            for j, (motif2, score2) in enumerate(gapped_motifs[i+1:], i+1):
                if j in used:
                    continue
                
                if self._are_gapped_motifs_similar(motif1, motif2):
                    group.append((motif2, score2))
                    used.add(j)
                    logger.info(f"  Grouping gapped motifs: {motif1} ~ {motif2}")
            
            groups.append(group)
        
        # Sort groups by best score in each group
        groups.sort(key=lambda g: max(score for _, score in g), reverse=True)
        return groups
    
    def _are_gapped_motifs_similar(self, motif1: str, motif2: str) -> bool:
        """Check if two gapped motifs are similar enough to merge."""
        # Extract non-N parts (cores)
        core1 = motif1.replace('N', '')
        core2 = motif2.replace('N', '')
        
        # If cores are very similar or one contains the other
        if core1 in core2 or core2 in core1:
            return True
        
        # If cores have high similarity
        if len(core1) == len(core2):
            mismatches = sum(1 for a, b in zip(core1, core2) if a != b)
            if mismatches <= 1:  # At most 1 mismatch
                return True
        
        # Check if gap patterns are similar
        gap_pattern1 = self._extract_gap_pattern(motif1)
        gap_pattern2 = self._extract_gap_pattern(motif2)
        
        if gap_pattern1 == gap_pattern2:
            return True
        
        return False
    
    def _extract_gap_pattern(self, motif: str) -> tuple:
        """Extract the gap pattern from a gapped motif."""
        # Return positions and lengths of N regions
        gaps = []
        in_gap = False
        gap_start = 0
        gap_length = 0
        
        for i, char in enumerate(motif):
            if char == 'N':
                if not in_gap:
                    in_gap = True
                    gap_start = i
                    gap_length = 1
                else:
                    gap_length += 1
            else:
                if in_gap:
                    gaps.append((gap_start, gap_length))
                    in_gap = False
        
        if in_gap:  # Gap at the end
            gaps.append((gap_start, gap_length))
        
        return tuple(gaps)
    
    def degenerate_code(self, bases):
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
    
    def consensus_from_cluster(self, cluster):
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
            consensus.append(self.degenerate_code(letters))
        return ''.join(consensus)
    
    def extract_core(self, motif, core_length=4):
        """
        Extract the central core of the motif.
        For a motif longer than core_length, return the central core; otherwise, return the motif.
        """
        L = len(motif)
        if L <= core_length:
            return motif
        start = (L - core_length) // 2
        return motif[start:start + core_length]
    
    def merge_by_core(self, motifs, core_length=4):
        """
        Group motifs by their extracted core (using the central core of length core_length)
        and generate a consensus motif for each group.
        """
        core_groups = {}
        for m in motifs:
            core = self.extract_core(m, core_length)
            if core not in core_groups:
                core_groups[core] = []
            core_groups[core].append(m)
        merged = [self.consensus_from_cluster(group) for group in core_groups.values()]
        return merged
    
    def merge_similar_motifs(self, motifs, max_distance=1):
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
        merged = [self.consensus_from_cluster(cluster) for cluster in clusters]
        return merged
    
    def longest_common_substring(self, strings, min_length=4):
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
    
    def finalize_cluster(self, cluster, min_common_length=4):
        """
        Given a cluster of motifs (strings), finalize the cluster by extracting the longest
        common substring of length at least min_common_length.
        If found, return it; otherwise, return the consensus from the cluster.
        """
        lcs = self.longest_common_substring(cluster, min_length=min_common_length)
        if lcs:
            return lcs
        else:
            return self.consensus_from_cluster(cluster)
    
    def merge_final_motifs(self, motifs, min_common_length=4, max_distance=1, core_length=4):
        """
        Merge a list of raw motifs using multiple strategies:
          1. Merge by Hamming distance.
          2. Merge by extracted core.
          3. Finalize each cluster by extracting the longest common substring.
        Returns the final list of merged motifs.
        """
        merged_hd = self.merge_similar_motifs(motifs, max_distance=max_distance)
        merged_core = self.merge_by_core(motifs, core_length=core_length)
        # Combine both sets and then finalize clusters by common substring.
        combined = list(set(merged_hd + merged_core))
        # If there are multiple motifs, further cluster them by common substring.
        # For simplicity, we treat the entire list as one cluster and finalize.
        final = self.finalize_cluster(combined, min_common_length=min_common_length)
        return final 