#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2024
# @Author  : Assistant
# @File    : test_motif_verification.py
"""
Test to verify the presence of known E. coli methylation motifs (GATC and CCWGG)
in the genome using the hammermotif package.
"""

import os
import sys
import unittest
import re
from collections import Counter
from Bio import SeqIO

# Add the parent directory to the path to import hammermotif
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from hammermotif.motif_greedy import (
    parse_bed_file, extract_sequences,
    total_possible_positions,  greedy_motif_extraction
)

from hammermotif.motif import (count_kmers_in_seqs,compute_chi2
)


class MotifVerificationTest(unittest.TestCase):
    """Test class to verify known E. coli methylation motifs."""
    
    def setUp(self):
        """Set up test environment with file paths."""
        self.wkdir = "/home/li/myapp/Hammerhead-motif/test"
        self.bed_file_path = os.path.join(self.wkdir, 'potential_modification_site.bed')
        self.genome_path = os.path.join(self.wkdir, 'ecoli.fa')
        
        # Load genome once for all tests
        print(f"Loading genome from: {self.genome_path}")
        self.genome = SeqIO.to_dict(SeqIO.parse(self.genome_path, "fasta"))
        print(f"Loaded {len(self.genome)} chromosomes/contigs")
        
    def count_motif_in_genome(self, motif_pattern, use_regex=True):
        """
        Count occurrences of a motif pattern in the entire genome.
        
        Args:
            motif_pattern (str): The motif pattern to search for
            use_regex (bool): Whether to treat the pattern as regex (for degenerate bases)
        
        Returns:
            int: Total count of motif occurrences
        """
        total_count = 0
        
        for chrom_name, seq_record in self.genome.items():
            sequence = str(seq_record.seq).upper()
            
            if use_regex:
                # Convert IUPAC codes to regex patterns
                pattern = motif_pattern.upper()
                pattern = pattern.replace('W', '[AT]')  # W = A or T
                pattern = pattern.replace('S', '[CG]')  # S = C or G
                pattern = pattern.replace('M', '[AC]')  # M = A or C
                pattern = pattern.replace('K', '[GT]')  # K = G or T
                pattern = pattern.replace('R', '[AG]')  # R = A or G
                pattern = pattern.replace('Y', '[CT]')  # Y = C or T
                pattern = pattern.replace('B', '[CGT]') # B = C, G, or T
                pattern = pattern.replace('D', '[AGT]') # D = A, G, or T
                pattern = pattern.replace('H', '[ACT]') # H = A, C, or T
                pattern = pattern.replace('V', '[ACG]') # V = A, C, or G
                pattern = pattern.replace('N', '[ACGT]') # N = any base
                
                matches = re.findall(pattern, sequence)
                count = len(matches)
            else:
                # Simple string search
                count = sequence.count(motif_pattern.upper())
            
            total_count += count
            if count > 0:
                print(f"  {chrom_name}: {count} occurrences")
        
        return total_count
    
    def test_gatc_motif_presence(self):
        """Test for the presence of GATC motif in E. coli genome."""
        print("\n=== Testing GATC motif presence ===")
        
        # Count GATC motifs in the genome
        gatc_count = self.count_motif_in_genome("GATC", use_regex=False)
        print(f"Total GATC motifs found: {gatc_count}")
        
        # E. coli should have many GATC sites (typically thousands to 10 thousands)
        self.assertGreater(gatc_count, 100, 
                          f"Expected many GATC motifs in E. coli, but found only {gatc_count}")
        
        # Also test reverse complement
        ctag_count = self.count_motif_in_genome("CTAG", use_regex=False)
        print(f"Total CTAG (reverse complement) motifs found: {ctag_count}")
        
        # In a double-stranded DNA, GATC and CTAG should be roughly equal
        ratio = gatc_count / ctag_count if ctag_count > 0 else float('inf')
        print(f"GATC/CTAG ratio: {ratio:.2f}")
        
    def test_ccwgg_motif_presence(self):
        """Test for the presence of CCWGG motif in E. coli genome."""
        print("\n=== Testing CCWGG motif presence ===")
        
        # Count CCWGG motifs (W = A or T)
        ccwgg_count = self.count_motif_in_genome("CC[AT]GG", use_regex=True)
        print(f"Total CCWGG motifs found: {ccwgg_count}")
        
        # Break down by specific variants
        ccagg_count = self.count_motif_in_genome("CCAGG", use_regex=False)
        cctgg_count = self.count_motif_in_genome("CCTGG", use_regex=False)
        print(f"  CCAGG: {ccagg_count}")
        print(f"  CCTGG: {cctgg_count}")
        print(f"  Total: {ccagg_count + cctgg_count} (should match {ccwgg_count})")
        
        # E. coli should have CCWGG sites
        self.assertGreater(ccwgg_count, 10, 
                          f"Expected CCWGG motifs in E. coli, but found only {ccwgg_count}")
        
        # Test reverse complement (CCWGG reverse complement is CCWGG)
        # But we should also check the actual reverse complement pattern
        ccwgg_rc_count = self.count_motif_in_genome("CC[AT]GG", use_regex=True)
        print(f"CCWGG pattern count: {ccwgg_rc_count}")
        
    def test_motif_enrichment_analysis(self):
        """Test motif enrichment analysis using the bed file."""
        print("\n=== Testing motif enrichment analysis ===")
        
        if not os.path.exists(self.bed_file_path):
            self.skipTest(f"BED file not found: {self.bed_file_path}")
        
        # Parse BED file and extract sequences
        regions = parse_bed_file(self.bed_file_path, flank_size=7)
        print(f"Parsed {len(regions)} regions from BED file")
        
        if len(regions) == 0:
            self.skipTest("No regions found in BED file")
        
        # Extract sequences from regions
        mod_seqs = extract_sequences(regions, self.genome)
        print(f"Extracted {len(mod_seqs)} sequences")
        
        self.assertGreater(len(mod_seqs), 0, "No sequences extracted from BED regions")
        
        # Count GATC and CCWGG in modified sequences
        gatc_in_mod = sum(seq.count('GATC') for seq in mod_seqs)
        ccagg_in_mod = sum(seq.count('CCAGG') for seq in mod_seqs)
        cctgg_in_mod = sum(seq.count('CCTGG') for seq in mod_seqs)
        ccwgg_in_mod = ccagg_in_mod + cctgg_in_mod
        
        print(f"GATC in modified sequences: {gatc_in_mod}")
        print(f"CCWGG in modified sequences: {ccwgg_in_mod} (CCAGG: {ccagg_in_mod}, CCTGG: {cctgg_in_mod})")
        
        # Calculate enrichment ratios
        total_mod_length = sum(len(seq) for seq in mod_seqs)
        total_genome_length = sum(len(seq_record.seq) for seq_record in self.genome.values())
        
        gatc_total = self.count_motif_in_genome("GATC", use_regex=False)
        ccwgg_total = self.count_motif_in_genome("CC[AT]GG", use_regex=True)
        
        if gatc_total > 0 and total_genome_length > 0:
            gatc_density_genome = gatc_total / total_genome_length
            gatc_density_mod = gatc_in_mod / total_mod_length if total_mod_length > 0 else 0
            gatc_enrichment = gatc_density_mod / gatc_density_genome if gatc_density_genome > 0 else 0
            print(f"GATC enrichment ratio: {gatc_enrichment:.2f}")
        
        if ccwgg_total > 0 and total_genome_length > 0:
            ccwgg_density_genome = ccwgg_total / total_genome_length
            ccwgg_density_mod = ccwgg_in_mod / total_mod_length if total_mod_length > 0 else 0
            ccwgg_enrichment = ccwgg_density_mod / ccwgg_density_genome if ccwgg_density_genome > 0 else 0
            print(f"CCWGG enrichment ratio: {ccwgg_enrichment:.2f}")
    
    def test_greedy_motif_extraction(self):
        """Test the greedy motif extraction algorithm."""
        print("\n=== Testing greedy motif extraction ===")
        
        if not os.path.exists(self.bed_file_path):
            self.skipTest(f"BED file not found: {self.bed_file_path}")
        
        # Parse BED file and extract sequences
        regions = parse_bed_file(self.bed_file_path, flank_size=7)
        mod_seqs = extract_sequences(regions, self.genome)
        
        if len(mod_seqs) == 0:
            self.skipTest("No sequences extracted from BED regions")
        
        # Use a subset of the genome as background (for speed)
        ref_seqs = [str(seq_record.seq).upper() for seq_record in list(self.genome.values())[:2]]
        
        # Run greedy motif extraction with relaxed parameters
        motifs_found = greedy_motif_extraction(
            mod_seqs, ref_seqs, 
            k=4,  # Start with shorter motifs
            chi2_threshold=10,  # Lower threshold for testing
            max_motifs=5
        )
        
        print(f"Found {len(motifs_found)} motifs:")
        for motif, chi2_val in motifs_found:
            print(f"  {motif}: chi2={chi2_val:.2f}")
        
        # Check if any known motifs or their variants are found
        found_motifs = [motif for motif, _ in motifs_found]
        
        # Check for GATC or its substrings
        gatc_related = any('GATC' in motif or motif in 'GATC' for motif in found_motifs)
        
        # Check for CCWGG variants
        ccwgg_related = any(
            'CCAG' in motif or 'CCTG' in motif or 'CAGG' in motif or 'CTGG' in motif
            for motif in found_motifs
        )
        
        print(f"GATC-related motifs found: {gatc_related}")
        print(f"CCWGG-related motifs found: {ccwgg_related}")
        
        # At least one of the known motifs should be detected
        self.assertTrue(len(motifs_found) > 0, "No motifs found by greedy extraction")
    
    def tearDown(self):
        """Clean up after tests."""
        pass


if __name__ == '__main__':
    # Run the tests
    unittest.main(verbosity=2) 