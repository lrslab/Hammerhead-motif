#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 14/3/2025 9:19
# @Author  : Runsheng
# @File    : meme_runner.py

import subprocess
import logging
import re
import numpy as np
from pathlib import Path
from Bio import motifs
from Bio.motifs import meme
import os


def run_meme(input_fasta, output_dir, meme_path="meme", **kwargs):
    """
    Run MEME with optimized parameters for methylation motifs.
    
    Parameters match the successful configuration:
    meme inputFile.fa -dna -oc . -nostatus -time 14400 -mod zoops 
    -nmotifs 10 -minw 4 -maxw 16 -objfun classic -revcomp -markov_order 0
    """
    cmd = [meme_path, input_fasta, "-oc", output_dir]

    # Add parameters from kwargs with proper formatting
    for key, value in kwargs.items():
        if key in ['dna', 'revcomp', 'nostatus']:
            # Boolean flags
            if value:
                cmd.append(f"-{key}")
        elif key == 'markov_order':
            # Special case for markov_order parameter
            cmd.extend(["-markov_order", str(value)])
        elif key == 'objfun':
            # Special case for objfun parameter
            cmd.extend(["-objfun", str(value)])
        elif isinstance(value, bool):
            # Other boolean parameters
            if value:
                cmd.append(f"-{key}")
        else:
            # Regular parameters with values
            cmd.extend([f"-{key}", str(value)])

    try:
        # Log the command being run
        logging.info(f"Running MEME command: {' '.join(cmd)}")
        
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        
        # Log successful completion
        logging.info("MEME completed successfully")
        
        # If there's any stderr output that's not an error, log it as info
        if result.stderr:
            logging.info(f"MEME stderr output: {result.stderr}")
            
    except subprocess.CalledProcessError as e:
        logging.error(f"MEME failed with return code {e.returncode}")
        logging.error(f"MEME stderr: {e.stderr}")
        logging.error(f"MEME stdout: {e.stdout}")
        raise


def parse_meme_output(meme_output_file):
    """
    Parse MEME output file using BioPython and extract motif information.
    Tries XML first (if available), then falls back to text parsing.
    """
    motifs_list = []
    
    # Check if XML file exists alongside the text file
    if meme_output_file.endswith('.txt'):
        xml_file = meme_output_file.replace('.txt', '.xml')
    else:
        xml_file = meme_output_file
    
    # Try XML parsing first if available
    if os.path.exists(xml_file) and xml_file.endswith('.xml'):
        try:
            logging.info(f"Attempting to parse MEME XML file: {xml_file}")
            with open(xml_file) as f:
                record = meme.read(f)
            
            # Extract motifs from the record
            for i, motif in enumerate(record):
                # Get consensus sequence
                consensus = str(motif.consensus)
                
                # Get degenerate consensus
                degenerate_consensus = str(motif.degenerate_consensus) if hasattr(motif, 'degenerate_consensus') else consensus
                
                # Get E-value
                evalue = motif.evalue if hasattr(motif, 'evalue') else 0.0
                
                # Get number of sites
                num_sites = len(motif.instances) if hasattr(motif, 'instances') else 0
                
                # Get motif length
                length = motif.length if hasattr(motif, 'length') else len(consensus)
                
                # Get information content
                ic = motif.ic if hasattr(motif, 'ic') else 0.0
                
                # Get position weight matrix
                pwm = []
                if hasattr(motif, 'counts'):
                    # Convert counts to PWM
                    for pos in range(length):
                        counts = {}
                        total = 0
                        for base in 'ACGT':
                            count = motif.counts[base][pos] if base in motif.counts else 0
                            counts[base] = count
                            total += count
                        
                        if total > 0:
                            probs = [counts[base]/total for base in 'ACGT']
                        else:
                            probs = [0.25, 0.25, 0.25, 0.25]
                        pwm.append(probs)
                
                motif_info = {
                    'consensus': consensus,
                    'degenerate_consensus': degenerate_consensus,
                    'evalue': evalue,
                    'nsites': num_sites,
                    'ic': ic,
                    'length': length,
                    'pwm': np.array(pwm) if pwm else np.array([])
                }
                
                motifs_list.append(motif_info)
                
                logging.info(f"Parsed motif {i+1}: {consensus} (E-value: {evalue:.2e}, Sites: {num_sites})")
            
            return motifs_list
            
        except Exception as e:
            logging.warning(f"Error parsing MEME XML with BioPython: {e}")
            logging.info("Falling back to text file parsing...")
    
    # Fall back to text parsing
    try:
        # Use BioPython's MEME text parser
        text_file = meme_output_file if meme_output_file.endswith('.txt') else meme_output_file.replace('.xml', '.txt')
        
        if os.path.exists(text_file):
            logging.info(f"Attempting to parse MEME text file: {text_file}")
            with open(text_file) as f:
                record = meme.read(f)
            
            # Extract motifs from the record
            for i, motif in enumerate(record):
                # Get consensus sequence
                consensus = str(motif.consensus)
                
                # Get E-value
                evalue = motif.evalue if hasattr(motif, 'evalue') else 0.0
                
                # Get number of sites
                num_sites = len(motif.instances) if hasattr(motif, 'instances') else 0
                
                # Get information content
                ic = motif.ic if hasattr(motif, 'ic') else 0.0
                
                # Get position weight matrix
                pwm = []
                if hasattr(motif, 'counts'):
                    # Convert counts to PWM
                    for pos in range(len(motif)):
                        counts = [motif.counts[base][pos] for base in 'ACGT']
                        total = sum(counts)
                        if total > 0:
                            probs = [count/total for count in counts]
                        else:
                            probs = [0.25, 0.25, 0.25, 0.25]
                        pwm.append(probs)
                
                motif_info = {
                    'consensus': consensus,
                    'evalue': evalue,
                    'nsites': num_sites,
                    'ic': ic,
                    'pwm': np.array(pwm) if pwm else np.array([])
                }
                
                motifs_list.append(motif_info)
                
                logging.info(f"Parsed motif {i+1}: {consensus} (E-value: {evalue:.2e}, Sites: {num_sites})")
    
    except Exception as e:
        logging.error(f"Error parsing MEME output with BioPython: {e}")
        
        # Fallback to manual parsing if BioPython fails
        logging.info("Attempting manual parsing of MEME output...")
        motifs_list = parse_meme_output_manual(meme_output_file)
    
    return motifs_list


def parse_meme_output_manual(meme_output_file):
    """
    Manual parsing of MEME output as fallback.
    """
    motifs_list = []
    current_motif = None
    in_sites_section = False
    in_pwm_section = False
    
    with open(meme_output_file) as f:
        lines = f.readlines()
        
        i = 0
        while i < len(lines):
            line = lines[i].strip()
            
            # Start of a new motif
            if line.startswith("MOTIF") and "MEME-" in line:
                if current_motif and current_motif.get('consensus'):
                    motifs_list.append(current_motif)
                
                parts = line.split()
                motif_id = parts[1]
                width = int(parts[-1]) if parts[-1].isdigit() else 0
                
                current_motif = {
                    'id': motif_id,
                    'width': width,
                    'consensus': '',
                    'evalue': 0.0,
                    'nsites': 0,
                    'ic': 0.0,
                    'pwm': []
                }
                in_sites_section = False
                in_pwm_section = False
            
            # E-value line
            elif current_motif and "E-value" in line:
                match = re.search(r'E-value[=:]\s*([0-9.e+-]+)', line)
                if match:
                    try:
                        current_motif['evalue'] = float(match.group(1))
                    except:
                        current_motif['evalue'] = 0.0
            
            # Sites line
            elif current_motif and "sites sorted by position p-value" in line:
                in_sites_section = True
                current_motif['sites'] = []
            
            # Count sites
            elif current_motif and in_sites_section and line and not line.startswith('-'):
                parts = line.split()
                if len(parts) >= 5 and parts[0].isalnum():
                    current_motif['sites'].append(parts[-1])  # The actual sequence
                elif not line:
                    in_sites_section = False
                    current_motif['nsites'] = len(current_motif.get('sites', []))
            
            # Letter-probability matrix
            elif current_motif and "letter-probability matrix" in line:
                in_pwm_section = True
                current_motif['pwm'] = []
                # Extract width and nsites from the line if available
                match = re.search(r'w=\s*(\d+)', line)
                if match:
                    current_motif['width'] = int(match.group(1))
                match = re.search(r'nsites=\s*(\d+)', line)
                if match:
                    current_motif['nsites'] = int(match.group(1))
            
            # PWM values
            elif current_motif and in_pwm_section:
                parts = line.split()
                if len(parts) == 4:
                    try:
                        probs = [float(p) for p in parts]
                        current_motif['pwm'].append(probs)
                    except:
                        in_pwm_section = False
                elif line.startswith('-') or not line:
                    in_pwm_section = False
                    # Generate consensus from PWM
                    if current_motif['pwm']:
                        consensus = generate_consensus_from_pwm(current_motif['pwm'])
                        current_motif['consensus'] = consensus
            
            i += 1
    
    # Don't forget the last motif
    if current_motif and current_motif.get('consensus'):
        motifs_list.append(current_motif)
    
    # If no consensus was generated from PWM, try to extract from sites
    for motif in motifs_list:
        if not motif.get('consensus') and motif.get('sites'):
            # Use the most common pattern from sites
            motif['consensus'] = get_consensus_from_sites(motif['sites'])
    
    return motifs_list


def generate_consensus_from_pwm(pwm):
    """Generate consensus sequence from position weight matrix."""
    bases = ['A', 'C', 'G', 'T']
    consensus = ''
    
    for position in pwm:
        max_idx = position.index(max(position))
        consensus += bases[max_idx]
    
    return consensus


def get_consensus_from_sites(sites):
    """Generate consensus from aligned sites."""
    if not sites:
        return ''
    
    # Find the most common length
    lengths = [len(s) for s in sites]
    common_length = max(set(lengths), key=lengths.count)
    
    # Filter sites to common length
    filtered_sites = [s for s in sites if len(s) == common_length]
    
    if not filtered_sites:
        return ''
    
    # Build consensus
    consensus = ''
    for i in range(common_length):
        bases = [s[i] for s in filtered_sites if i < len(s)]
        if bases:
            # Most common base at this position
            most_common = max(set(bases), key=bases.count)
            consensus += most_common
    
    return consensus