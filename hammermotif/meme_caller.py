#!/usr/bin/env python3
"""
MEME motif discovery module for bacterial methylation motif calling.

This module handles running MEME and managing the MEME process.
"""

import os
import logging
import subprocess
import tempfile
import random
from typing import List, Dict, Optional
from pathlib import Path

logger = logging.getLogger(__name__)


class MEMEMotifCaller:
    """MEME-based motif discovery for bacterial methylation sites."""
    
    def __init__(self, genome: Dict[str, str], mod_sites: List[Dict]):
        """
        Initialize the MEME motif caller.
        
        Parameters:
        -----------
        genome : Dict[str, str]
            Dictionary of chromosome sequences
        mod_sites : List[Dict]
            List of modification sites from BED file
        """
        self.genome = genome
        self.mod_sites = mod_sites
        self.meme_motifs = []
        self.temp_files = []
    
    def extract_sequences_for_meme(self, window_size: int = 25) -> List[str]:
        """
        Extract sequences around modification sites for MEME analysis.
        
        Parameters:
        -----------
        window_size : int
            Window size around modification sites (default: 25)
        """
        sequences = []
        
        for site in self.mod_sites:
            chrom = site['chrom']
            if chrom not in self.genome:
                logger.warning(f"Chromosome {chrom} not found in genome")
                continue
                
            center = (site['start'] + site['end']) // 2
            start = max(0, center - window_size)
            end = min(len(self.genome[chrom]), center + window_size + 1)
            
            sequence = self.genome[chrom][start:end].upper()
            
            # Skip sequences with too many N's
            if sequence.count('N') > len(sequence) * 0.1:
                continue
            
            # Handle strand orientation
            if site.get('strand') == '-':
                # Simple reverse complement
                complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
                sequence = ''.join(complement.get(b, 'N') for b in sequence[::-1])
            
            sequences.append(sequence)
        
        return sequences
    
    def prepare_meme_input(self, sequences: List[str], 
                          max_sequences: int = 20000,
                          min_length: int = 20) -> str:
        """
        Prepare input FASTA file for MEME.
        
        Parameters:
        -----------
        sequences : List[str]
            Sequences to analyze
        max_sequences : int
            Maximum number of sequences to use
        min_length : int
            Minimum sequence length
            
        Returns:
        --------
        Path to temporary FASTA file
        """
        # Filter sequences by length
        sequences = [seq for seq in sequences if len(seq) >= min_length]
        
        # Sample if too many sequences
        if len(sequences) > max_sequences:
            random.seed(42)  # For reproducibility
            sequences = random.sample(sequences, max_sequences)
        
        # Create temporary FASTA file
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', 
                                       delete=False) as f:
            for i, seq in enumerate(sequences):
                f.write(f">seq_{i}\n{seq}\n")
            temp_fasta = f.name
        
        self.temp_files.append(temp_fasta)
        logger.info(f"Created MEME input file with {len(sequences)} sequences")
        
        return temp_fasta
    
    def run_meme(self, nmotifs: int = 10, minw: int = 4, maxw: int = 12,
                 mod: str = "zoops", objfun: str = "classic",
                 markov_order: int = 0, time: int = 7200,
                 nostatus: bool = True, use_parallel: bool = True) -> Optional[str]:
        """
        Run MEME for motif discovery.
        
        Parameters:
        -----------
        nmotifs : int
            Number of motifs to find
        minw : int
            Minimum motif width
        maxw : int
            Maximum motif width
        mod : str
            Motif distribution model (oops, zoops, anr)
        objfun : str
            Objective function (classic, cd)
        markov_order : int
            Order of Markov background model
        time : int
            Maximum runtime in seconds
        nostatus : bool
            Suppress status messages
        use_parallel : bool
            Use parallel processing if available
            
        Returns:
        --------
        Path to MEME output directory or None if failed
        """
        logger.info(f"Running MEME with parameters: nmotifs={nmotifs}, "
                   f"width={minw}-{maxw}, mod={mod}")
        
        # Extract sequences
        sequences = self.extract_sequences_for_meme(window_size=25)
        if not sequences:
            logger.error("No sequences extracted for MEME")
            return None
        
        # Prepare input
        input_fasta = self.prepare_meme_input(sequences)
        
        # Create output directory
        output_dir = tempfile.mkdtemp(prefix="meme_output_")
        self.temp_files.append(output_dir)
        
        # Build MEME command
        cmd = [
            'meme', input_fasta,
            '-o', output_dir,
            '-dna',
            '-revcomp',
            '-nmotifs', str(nmotifs),
            '-minw', str(minw),
            '-maxw', str(maxw),
            '-mod', mod,
            '-objfun', objfun,
            '-markov_order', str(markov_order),
            '-time', str(time)
        ]
        
        if nostatus:
            cmd.append('-nostatus')
        
        if use_parallel:
            # Check if MEME supports parallel processing
            try:
                result = subprocess.run(['meme', '-h'], 
                                      capture_output=True, text=True)
                if '-p' in result.stdout:
                    import multiprocessing
                    ncpus = min(4, multiprocessing.cpu_count())
                    cmd.extend(['-p', str(ncpus)])
            except:
                pass
        
        # Run MEME
        try:
            logger.info(f"Running command: {' '.join(cmd)}")
            result = subprocess.run(cmd, capture_output=True, text=True,
                                  timeout=time + 300)  # Add buffer to timeout
            
            if result.returncode == 0:
                logger.info("MEME completed successfully")
                return output_dir
            else:
                logger.error(f"MEME failed with return code {result.returncode}")
                logger.error(f"STDERR: {result.stderr}")
                return None
                
        except subprocess.TimeoutExpired:
            logger.error(f"MEME timed out after {time} seconds")
            return None
        except Exception as e:
            logger.error(f"Error running MEME: {e}")
            return None
    
    def cleanup(self):
        """Clean up temporary files."""
        for temp_file in self.temp_files:
            try:
                if os.path.isdir(temp_file):
                    subprocess.run(['rm', '-rf', temp_file], check=False)
                elif os.path.isfile(temp_file):
                    os.unlink(temp_file)
            except Exception as e:
                logger.warning(f"Failed to clean up {temp_file}: {e}")
        
        self.temp_files = []
    
    def discover_motifs(self, **kwargs) -> Optional[str]:
        """
        Main method to discover motifs using MEME.
        
        Returns path to MEME output directory or None if failed.
        """
        try:
            output_dir = self.run_meme(**kwargs)
            return output_dir
        finally:
            # Always cleanup temporary input files
            # Keep output directory for parsing
            input_files = [f for f in self.temp_files if f.endswith('.fa')]
            for f in input_files:
                try:
                    os.unlink(f)
                    self.temp_files.remove(f)
                except:
                    pass 