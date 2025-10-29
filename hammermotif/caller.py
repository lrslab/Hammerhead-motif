#!/usr/bin/env python3
"""
Comprehensive Bacterial DNA Methylation Motif Caller

This is the main integration script that coordinates all motif discovery methods
and generates comprehensive reports.

Features:
- Integrates greedy, MEME-based motif discovery
- Analyzes known motifs for enrichment
- Generates comprehensive reports
- Supports batch processing
"""

import os
import sys
import json
import logging
import argparse
import re
import numpy as np
from pathlib import Path
from datetime import datetime
from typing import List, Dict, Tuple, Optional

# Import split modules
from .utils import read_fasta_file, read_bed_file
from .greedy_caller import GreedyMotifCaller
from .meme_caller import MEMEMotifCaller
from .meme_parser import MEMEParser
from .motif_merger import MotifMerger

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


class MethylationMotifCaller:
    """Main integration class for bacterial methylation motif discovery."""
    
    def __init__(self, genome_file: str, bed_file: str, output_dir: str = None):
        """
        Initialize the motif caller.
        
        Parameters:
        -----------
        genome_file : str
            Path to genome FASTA file
        bed_file : str
            Path to BED file with modification sites
        output_dir : str
            Output directory for results
        """
        self.genome_file = genome_file
        self.bed_file = bed_file
        self.genome_name = Path(genome_file).stem
        
        if output_dir is None:
            output_dir = f"results_{self.genome_name}_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True, parents=True)
        
        # Load genome and modification sites
        logger.info(f"Loading genome from {genome_file}")
        self.genome = read_fasta_file(genome_file)
        self.genome_size = sum(len(seq) for seq in self.genome.values())
        
        logger.info(f"Loading modification sites from {bed_file}")
        self.mod_sites = read_bed_file(bed_file)
        
        logger.info(f"Genome: {len(self.genome)} chromosomes, {self.genome_size:,} bp")
        logger.info(f"Modification sites: {len(self.mod_sites):,}")
        
        # Initialize components
        self.greedy_caller = GreedyMotifCaller(self.genome, self.mod_sites, self.genome_size)
        self.meme_caller = MEMEMotifCaller(self.genome, self.mod_sites)
        self.meme_parser = MEMEParser()
        self.motif_merger = MotifMerger()
        
        # Storage for results
        self.greedy_results = {}
        self.meme_results = []
        self.final_motifs = []
        self.known_motif_analysis = {}
    
    def run_greedy_discovery(self, k_values: List[int] = None, 
                           chi2_threshold: float = 50,
                           discover_gapped: bool = True,
                           discover_degenerate: bool = True,
                           use_parallel: bool = True) -> Dict:
        """
        Run all greedy-based discovery methods.
        
        Returns dictionary with all greedy results.
        """
        logger.info("="*60)
        logger.info("Running Greedy Motif Discovery")
        logger.info("="*60)
        
        results = {}
        
        # Standard greedy discovery
        logger.info("\n1. Standard k-mer motifs:")
        greedy_motifs = self.greedy_caller.run_greedy_discovery(
            k_values=k_values,
            chi2_threshold=chi2_threshold,
            use_parallel=use_parallel
        )
        results['standard'] = greedy_motifs
        
        # Gapped motifs
        if discover_gapped:
            logger.info("\n2. Gapped motifs:")
            gapped_motifs = self.greedy_caller.discover_gapped_motifs_optimized(
                max_gap=10,
                min_count=10,
                use_regex=True
            )
            results['gapped'] = gapped_motifs
        
        # Degenerate motifs
        if discover_degenerate:
            logger.info("\n3. Degenerate motifs:")
            degenerate_motifs = self.greedy_caller.discover_degenerate_motifs_fast(
                min_length=4,
                min_ic=3.0,
                min_sites=10
            )
            results['degenerate'] = degenerate_motifs
        
        self.greedy_results = results
        return results
    
    def run_meme_discovery(self, nmotifs: int = 10, 
                          minw: int = 4, 
                          maxw: int = 12,
                          skip_if_no_meme: bool = True) -> List[Dict]:
        """
        Run MEME-based motif discovery.
        
        Returns list of MEME motifs.
        """
        logger.info("\n" + "="*60)
        logger.info("Running MEME Motif Discovery")
        logger.info("="*60)
        
        try:
            # Check if MEME is available
            import subprocess
            result = subprocess.run(['which', 'meme'], capture_output=True)
            if result.returncode != 0:
                if skip_if_no_meme:
                    logger.warning("MEME not found in PATH, skipping MEME analysis")
                    return []
                else:
                    raise RuntimeError("MEME not found in PATH")
            
            # Run MEME
            output_dir = self.meme_caller.discover_motifs(
                nmotifs=nmotifs,
                minw=minw,
                maxw=maxw,
                use_parallel=True
            )
            
            if output_dir:
                # Parse results
                meme_txt = os.path.join(output_dir, "meme.txt")
                meme_xml = os.path.join(output_dir, "meme.xml")
                
                if os.path.exists(meme_txt):
                    self.meme_results = self.meme_parser.parse_meme_output(meme_txt)
                elif os.path.exists(meme_xml):
                    self.meme_results = self.meme_parser.parse_meme_xml(meme_xml)
                
                logger.info(f"Found {len(self.meme_results)} MEME motifs")
                
                # Cleanup
                self.meme_caller.cleanup()
            
        except Exception as e:
            logger.error(f"MEME analysis failed: {e}")
            if not skip_if_no_meme:
                raise
        
        return self.meme_results
    
    def merge_and_finalize_motifs(self, chi2_threshold: float = 30,
                                 evalue_threshold: float = 0.01) -> List[str]:
        """
        Merge motifs from all discovery methods.
        """
        logger.info("\n" + "="*60)
        logger.info("Merging and Finalizing Motifs")
        logger.info("="*60)
        
        # Collect all greedy motifs
        all_greedy = []
        for motif_type, motifs in self.greedy_results.items():
            all_greedy.extend(motifs)
        
        # Merge using the merger module
        self.final_motifs = self.motif_merger.merge_motifs(
            all_greedy,
            self.meme_results,
            chi2_threshold=chi2_threshold,
            evalue_threshold=evalue_threshold
        )
        
        logger.info(f"\nFinal motifs ({len(self.final_motifs)}):")
        for i, motif in enumerate(self.final_motifs, 1):
            logger.info(f"  {i}. {motif}")
        
        return self.final_motifs
    
    def analyze_known_motifs(self, known_motifs: List[str]) -> Dict[str, Dict]:
        """
        Analyze enrichment of known motifs in the dataset.
        """
        logger.info("\n" + "="*60)
        logger.info("Analyzing Known Motifs")
        logger.info("="*60)
        
        from scipy.stats import hypergeom
        import re
        
        results = {}
        
        # Extract sequences around modification sites
        mod_seqs = self.greedy_caller.extract_modified_sequences(window_size=25)
        total_mod_length = sum(len(seq) for seq in mod_seqs)
        
        for motif in known_motifs:
            logger.info(f"\nAnalyzing: {motif}")
            
            # Convert to regex pattern
            pattern = self._motif_to_regex(motif)
            
            # Count in modified sequences
            mod_count = sum(len(pattern.findall(seq)) for seq in mod_seqs)
            
            # Count in genome
            genome_count = 0
            for chrom_seq in self.genome.values():
                genome_count += len(pattern.findall(str(chrom_seq).upper()))
            
            # Calculate enrichment
            if genome_count > 0:
                expected = (genome_count / self.genome_size) * total_mod_length
                enrichment = mod_count / expected if expected > 0 else 0
                
                # Calculate p-value
                p_value = hypergeom.sf(
                    mod_count - 1,
                    self.genome_size,
                    genome_count,
                    total_mod_length
                )
            else:
                enrichment = 0
                p_value = 1.0
            
            results[motif] = {
                'mod_count': int(mod_count),
                'genome_count': int(genome_count),
                'enrichment': float(enrichment),
                'p_value': float(p_value),
                'significant': bool(p_value < 0.01 and enrichment > 2)
            }
            
            logger.info(f"  Modified sites: {mod_count}")
            logger.info(f"  Genome total: {genome_count}")
            logger.info(f"  Enrichment: {enrichment:.2f}x")
            logger.info(f"  P-value: {p_value:.2e}")
            logger.info(f"  Significant: {'Yes' if results[motif]['significant'] else 'No'}")
        
        self.known_motif_analysis = results
        return results
    
    def _motif_to_regex(self, motif: str) -> re.Pattern:
        """Convert IUPAC motif to regex pattern."""
        iupac_map = {
            'R': '[AG]', 'Y': '[CT]', 'W': '[AT]', 'S': '[GC]',
            'K': '[GT]', 'M': '[AC]', 'B': '[CGT]', 'D': '[AGT]',
            'H': '[ACT]', 'V': '[ACG]', 'N': '.'
        }
        
        pattern = motif.upper()
        for base, regex in iupac_map.items():
            pattern = pattern.replace(base, regex)
        
        return re.compile(pattern)
    
    def generate_report(self, include_sites: bool = False):
        """Generate comprehensive analysis report."""
        logger.info("\n" + "="*60)
        logger.info("Generating Report")
        logger.info("="*60)
        
        # Text report
        report_file = self.output_dir / f"{self.genome_name}_motif_report.txt"
        
        with open(report_file, 'w') as f:
            f.write("Bacterial Methylation Motif Analysis Report\n")
            f.write("=" * 60 + "\n\n")
            
            f.write(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"Genome: {self.genome_name}\n")
            f.write(f"Genome file: {self.genome_file}\n")
            f.write(f"BED file: {self.bed_file}\n")
            f.write(f"Genome size: {self.genome_size:,} bp\n")
            f.write(f"Modification sites: {len(self.mod_sites):,}\n\n")
            
            # Greedy results
            f.write("GREEDY MOTIF DISCOVERY\n")
            f.write("-" * 60 + "\n\n")
            
            if 'standard' in self.greedy_results:
                f.write("Standard motifs:\n")
                for motif, chi2 in self.greedy_results['standard'][:20]:
                    f.write(f"  {motif:<20} χ² = {chi2:>10.2f}\n")
                f.write("\n")
            
            if 'gapped' in self.greedy_results:
                f.write("Gapped motifs:\n")
                for motif, chi2 in self.greedy_results['gapped'][:10]:
                    f.write(f"  {motif:<20} χ² = {chi2:>10.2f}\n")
                f.write("\n")
            
            if 'degenerate' in self.greedy_results:
                f.write("Degenerate motifs:\n")
                for motif, chi2 in self.greedy_results['degenerate'][:10]:
                    f.write(f"  {motif:<20} χ² = {chi2:>10.2f}\n")
                f.write("\n")
            
            # MEME results
            if self.meme_results:
                f.write("MEME MOTIF DISCOVERY\n")
                f.write("-" * 60 + "\n\n")
                for motif in self.meme_results[:10]:
                    consensus = motif.get('degenerate_consensus', motif.get('consensus', 'N/A'))
                    f.write(f"  {consensus:<20} E-value = {motif.get('evalue', 0):.2e}, "
                           f"Sites = {motif.get('nsites', 0)}\n")
                f.write("\n")
            
            # Final motifs
            f.write("FINAL MERGED MOTIFS\n")
            f.write("-" * 60 + "\n")
            for i, motif in enumerate(self.final_motifs, 1):
                f.write(f"  {i}. {motif}\n")
            f.write("\n")
            
            # Known motif analysis
            if self.known_motif_analysis:
                f.write("KNOWN MOTIF ANALYSIS\n")
                f.write("-" * 60 + "\n")
                f.write(f"{'Motif':<15} {'Enrichment':>10} {'P-value':>12} {'Significant':>12}\n")
                f.write("-" * 50 + "\n")
                for motif, stats in self.known_motif_analysis.items():
                    f.write(f"{motif:<15} {stats['enrichment']:>10.2f}x "
                           f"{stats['p_value']:>12.2e} "
                           f"{'Yes' if stats['significant'] else 'No':>12}\n")
        
        # JSON report for programmatic access
        json_report = self.output_dir / f"{self.genome_name}_motif_report.json"
        
        report_data = {
            'metadata': {
                'date': datetime.now().isoformat(),
                'genome': self.genome_name,
                'genome_file': str(self.genome_file),
                'bed_file': str(self.bed_file),
                'genome_size': self.genome_size,
                'modification_sites': len(self.mod_sites)
            },
            'greedy_motifs': {
                motif_type: [(m, float(s)) for m, s in motifs[:20]]
                for motif_type, motifs in self.greedy_results.items()
            },
            'meme_motifs': [
                {
                    'consensus': m.get('consensus', ''),
                    'degenerate_consensus': m.get('degenerate_consensus', ''),
                    'evalue': float(m.get('evalue', 0)),
                    'nsites': m.get('nsites', 0),
                    'width': m.get('width', 0)
                }
                for m in self.meme_results
            ],
            'final_motifs': self.final_motifs,
            'known_motif_analysis': {
                motif: {
                    'mod_count': int(stats['mod_count']),
                    'genome_count': int(stats['genome_count']),
                    'enrichment': float(stats['enrichment']),
                    'p_value': float(stats['p_value']),
                    'significant': bool(stats['significant'])  # Convert numpy bool to Python bool
                }
                for motif, stats in self.known_motif_analysis.items()
            } if self.known_motif_analysis else {}
        }
        
        with open(json_report, 'w') as f:
            json.dump(report_data, f, indent=2)
        
        logger.info(f"Reports saved to:")
        logger.info(f"  Text: {report_file}")
        logger.info(f"  JSON: {json_report}")
    
    def run_complete_analysis(self, 
                            known_motifs: List[str] = None,
                            use_meme: bool = True,
                            use_parallel: bool = True):
        """
        Run complete motif analysis pipeline.
        
        Parameters:
        -----------
        known_motifs : List[str]
            List of known motifs to analyze
        use_meme : bool
            Whether to run MEME analysis
        use_parallel : bool
            Whether to use parallel processing
        """
        logger.info(f"\nStarting complete analysis for {self.genome_name}")
        logger.info("=" * 60)
        
        # 1. Run greedy discovery
        self.run_greedy_discovery(use_parallel=use_parallel)
        
        # 2. Run MEME if requested
        if use_meme:
            self.run_meme_discovery()
        
        # 3. Merge and finalize motifs
        self.merge_and_finalize_motifs()
        
        # 4. Analyze known motifs if provided
        if known_motifs:
            self.analyze_known_motifs(known_motifs)
        
        # 5. Generate report
        self.generate_report()
        
        logger.info("\nAnalysis complete!")
        
        return self.final_motifs


def main():
    """Command-line interface for the motif caller."""
    parser = argparse.ArgumentParser(
        description='Bacterial DNA Methylation Motif Discovery',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage
  python -m hammermotif.caller genome.fa modifications.bed
  
  # With known motifs
  python -m hammermotif.caller genome.fa modifications.bed --known-motifs GATC CCWGG
  
  # Without MEME
  python -m hammermotif.caller genome.fa modifications.bed --no-meme
  
  # Custom output directory
  python -m hammermotif.caller genome.fa modifications.bed -o results/
        """
    )
    
    parser.add_argument('genome', help='Genome FASTA file')
    parser.add_argument('bed', help='BED file with modification sites')
    parser.add_argument('-o', '--output', help='Output directory', default=None)
    parser.add_argument('--known-motifs', nargs='+', help='Known motifs to analyze')
    parser.add_argument('--no-meme', action='store_true', help='Skip MEME analysis')
    parser.add_argument('--no-parallel', action='store_true', help='Disable parallel processing')
    parser.add_argument('--chi2-threshold', type=float, default=50, 
                       help='Chi-square threshold for motif significance')
    parser.add_argument('--k-values', nargs='+', type=int, default=None,
                       help='k-mer sizes to test (default: 4 5 6 7 8)')
    
    args = parser.parse_args()
    
    # Create caller instance
    caller = MethylationMotifCaller(args.genome, args.bed, args.output)
    
    # Run analysis
    caller.run_complete_analysis(
        known_motifs=args.known_motifs,
        use_meme=not args.no_meme,
        use_parallel=not args.no_parallel
    )


if __name__ == '__main__':
    main()
