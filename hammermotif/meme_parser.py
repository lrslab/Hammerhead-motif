#!/usr/bin/env python3
"""
MEME output parser module for bacterial methylation motif calling.

This module handles parsing MEME output files and extracting motif information.
"""

import os
import re
import logging
from typing import List, Dict, Optional, Tuple
from pathlib import Path

logger = logging.getLogger(__name__)


class MEMEParser:
    """Parser for MEME output files."""
    
    def __init__(self):
        self.motifs = []
    
    def parse_meme_output(self, meme_file: str) -> List[Dict]:
        """
        Parse MEME output text file.
        
        Parameters:
        -----------
        meme_file : str
            Path to meme.txt output file
            
        Returns:
        --------
        List of motif dictionaries
        """
        if not os.path.exists(meme_file):
            logger.error(f"MEME output file not found: {meme_file}")
            return []
        
        motifs = []
        
        with open(meme_file, 'r') as f:
            content = f.read()
        
        # Find all motif blocks
        motif_pattern = re.compile(
            r'MOTIF\s+(\S+)\s+MEME-(\d+).*?'
            r'letter-probability matrix:.*?'
            r'log-odds matrix:.*?'
            r'((?:[\s\d.-]+\n)+)',
            re.DOTALL
        )
        
        for match in motif_pattern.finditer(content):
            motif_name = match.group(1)
            motif_num = int(match.group(2))
            
            # Extract motif information
            motif_info = self._extract_motif_info(content, motif_name)
            
            if motif_info:
                motif_info['name'] = motif_name
                motif_info['number'] = motif_num
                motifs.append(motif_info)
        
        # Extract additional information
        for motif in motifs:
            # Get consensus sequence
            consensus = self._extract_consensus(content, motif['name'])
            if consensus:
                motif['consensus'] = consensus
            
            # Get degenerate consensus
            degenerate = self._build_degenerate_consensus(motif)
            if degenerate:
                motif['degenerate_consensus'] = degenerate
        
        self.motifs = motifs
        logger.info(f"Parsed {len(motifs)} motifs from MEME output")
        
        return motifs
    
    def _extract_motif_info(self, content: str, motif_name: str) -> Dict:
        """Extract detailed information for a specific motif."""
        info = {}
        
        # Find the motif block
        motif_block_pattern = re.compile(
            rf'MOTIF\s+{re.escape(motif_name)}\s+MEME.*?'
            rf'width\s*=\s*(\d+).*?'
            rf'sites\s*=\s*(\d+).*?'
            rf'E-value\s*=\s*([\d.e+-]+)',
            re.DOTALL | re.IGNORECASE
        )
        
        match = motif_block_pattern.search(content)
        if match:
            info['width'] = int(match.group(1))
            info['nsites'] = int(match.group(2))
            info['evalue'] = float(match.group(3))
        
        # Extract position weight matrix
        pwm = self._extract_pwm(content, motif_name)
        if pwm:
            info['pwm'] = pwm
        
        # Extract sites
        sites = self._extract_sites(content, motif_name)
        if sites:
            info['sites'] = sites
        
        return info
    
    def _extract_pwm(self, content: str, motif_name: str) -> Optional[Dict[str, List[float]]]:
        """Extract position weight matrix for a motif."""
        # Find the letter-probability matrix section
        pattern = re.compile(
            rf'MOTIF\s+{re.escape(motif_name)}.*?'
            rf'letter-probability matrix:.*?\n'
            rf'((?:\s*[\d.]+\s+[\d.]+\s+[\d.]+\s+[\d.]+\s*\n)+)',
            re.DOTALL
        )
        
        match = pattern.search(content)
        if not match:
            return None
        
        matrix_text = match.group(1).strip()
        
        # Parse the matrix
        pwm = {'A': [], 'C': [], 'G': [], 'T': []}
        
        for line in matrix_text.split('\n'):
            if line.strip():
                values = line.strip().split()
                if len(values) == 4:
                    try:
                        pwm['A'].append(float(values[0]))
                        pwm['C'].append(float(values[1]))
                        pwm['G'].append(float(values[2]))
                        pwm['T'].append(float(values[3]))
                    except ValueError:
                        continue
        
        return pwm if pwm['A'] else None
    
    def _extract_consensus(self, content: str, motif_name: str) -> Optional[str]:
        """Extract consensus sequence for a motif."""
        # Look for consensus in the summary section
        pattern = re.compile(
            rf'{re.escape(motif_name)}\s+(\d+)\s+([ACGT]+)',
            re.IGNORECASE
        )
        
        match = pattern.search(content)
        if match:
            return match.group(2)
        
        # Try to find in multilevel consensus
        pattern2 = re.compile(
            rf'Multilevel\s+{re.escape(motif_name)}\s+consensus\s+([ACGTRYWSKMDHVBN]+)',
            re.IGNORECASE
        )
        
        match = pattern2.search(content)
        if match:
            return match.group(1)
        
        return None
    
    def _extract_sites(self, content: str, motif_name: str) -> List[Dict]:
        """Extract binding sites for a motif."""
        sites = []
        
        # Find the sites section for this motif
        pattern = re.compile(
            rf'MOTIF\s+{re.escape(motif_name)}.*?'
            rf'sites sorted by position p-value.*?\n'
            rf'((?:.*?\n)*?)(?:--------|\n\n)',
            re.DOTALL
        )
        
        match = pattern.search(content)
        if not match:
            return sites
        
        sites_text = match.group(1)
        
        # Parse individual sites
        site_pattern = re.compile(
            r'(\S+)\s+([+-])\s+(\d+)\s+[\d.e+-]+\s+([ACGT]+)'
        )
        
        for site_match in site_pattern.finditer(sites_text):
            site = {
                'sequence_name': site_match.group(1),
                'strand': site_match.group(2),
                'start': int(site_match.group(3)),
                'site_sequence': site_match.group(4)
            }
            sites.append(site)
        
        return sites
    
    def _build_degenerate_consensus(self, motif: Dict) -> Optional[str]:
        """Build degenerate consensus from PWM."""
        if 'pwm' not in motif:
            return None
        
        pwm = motif['pwm']
        if not pwm or not pwm['A']:
            return None
        
        consensus = []
        length = len(pwm['A'])
        
        for pos in range(length):
            # Get probabilities at this position
            probs = {
                'A': pwm['A'][pos],
                'C': pwm['C'][pos],
                'G': pwm['G'][pos],
                'T': pwm['T'][pos]
            }
            
            # Determine consensus base
            consensus_base = self._get_consensus_from_probs(probs)
            consensus.append(consensus_base)
        
        return ''.join(consensus)
    
    def _get_consensus_from_probs(self, probs: Dict[str, float], 
                                 threshold: float = 0.25) -> str:
        """Get IUPAC consensus base from probabilities."""
        # Get bases above threshold
        significant_bases = []
        for base, prob in probs.items():
            if prob >= threshold:
                significant_bases.append(base)
        
        if not significant_bases:
            return 'N'
        
        # Sort for consistent output
        significant_bases.sort()
        
        # Map to IUPAC code
        if len(significant_bases) == 1:
            return significant_bases[0]
        elif len(significant_bases) == 2:
            pair = tuple(significant_bases)
            iupac_map = {
                ('A', 'G'): 'R', ('C', 'T'): 'Y',
                ('A', 'T'): 'W', ('C', 'G'): 'S',
                ('G', 'T'): 'K', ('A', 'C'): 'M'
            }
            return iupac_map.get(pair, 'N')
        elif len(significant_bases) == 3:
            missing = {'A', 'C', 'G', 'T'} - set(significant_bases)
            if missing == {'A'}:
                return 'B'
            elif missing == {'C'}:
                return 'D'
            elif missing == {'G'}:
                return 'H'
            elif missing == {'T'}:
                return 'V'
        
        return 'N'
    
    def parse_meme_xml(self, xml_file: str) -> List[Dict]:
        """
        Parse MEME XML output for more structured data.
        
        Parameters:
        -----------
        xml_file : str
            Path to meme.xml file
            
        Returns:
        --------
        List of motif dictionaries
        """
        try:
            import xml.etree.ElementTree as ET
            
            tree = ET.parse(xml_file)
            root = tree.getroot()
            
            motifs = []
            
            # Find all motif elements
            for motif_elem in root.findall('.//motif'):
                motif = {
                    'id': motif_elem.get('id'),
                    'name': motif_elem.get('name'),
                    'width': int(motif_elem.get('width', 0)),
                    'sites': int(motif_elem.get('sites', 0)),
                    'evalue': float(motif_elem.get('e_value', 0))
                }
                
                # Extract PWM from scores
                scores_elem = motif_elem.find('scores')
                if scores_elem is not None:
                    pwm = self._parse_xml_pwm(scores_elem)
                    if pwm:
                        motif['pwm'] = pwm
                        motif['degenerate_consensus'] = self._build_degenerate_consensus(motif)
                
                # Extract regular expression
                re_elem = motif_elem.find('regular_expression')
                if re_elem is not None and re_elem.text:
                    motif['regex'] = re_elem.text.strip()
                
                motifs.append(motif)
            
            return motifs
            
        except Exception as e:
            logger.error(f"Error parsing MEME XML: {e}")
            return []
    
    def _parse_xml_pwm(self, scores_elem) -> Optional[Dict[str, List[float]]]:
        """Parse PWM from XML scores element."""
        pwm = {'A': [], 'C': [], 'G': [], 'T': []}
        
        alphabet_matrix = scores_elem.find('alphabet_matrix')
        if alphabet_matrix is None:
            return None
        
        for column in alphabet_matrix.findall('alphabet_array'):
            values = {}
            for letter in column.findall('value'):
                letter_id = letter.get('letter_id')
                if letter_id in pwm:
                    values[letter_id] = float(letter.text)
            
            # Add values in order
            for base in ['A', 'C', 'G', 'T']:
                if base in values:
                    pwm[base].append(values[base])
        
        return pwm if pwm['A'] else None
    
    def get_consensus_motifs(self) -> List[str]:
        """Get list of consensus sequences from parsed motifs."""
        consensus_list = []
        
        for motif in self.motifs:
            # Prefer degenerate consensus
            if 'degenerate_consensus' in motif:
                consensus_list.append(motif['degenerate_consensus'])
            elif 'consensus' in motif:
                consensus_list.append(motif['consensus'])
            elif 'regex' in motif:
                # Convert regex to IUPAC if possible
                iupac = self._regex_to_iupac(motif['regex'])
                if iupac:
                    consensus_list.append(iupac)
        
        return consensus_list
    
    def _regex_to_iupac(self, regex: str) -> Optional[str]:
        """Convert simple regex pattern to IUPAC notation."""
        # Handle simple character classes
        conversions = {
            '[AG]': 'R', '[GA]': 'R',
            '[CT]': 'Y', '[TC]': 'Y',
            '[AT]': 'W', '[TA]': 'W',
            '[GC]': 'S', '[CG]': 'S',
            '[GT]': 'K', '[TG]': 'K',
            '[AC]': 'M', '[CA]': 'M',
            '[CGT]': 'B', '[CTG]': 'B', '[GCT]': 'B', '[GTC]': 'B', '[TCG]': 'B', '[TGC]': 'B',
            '[AGT]': 'D', '[ATG]': 'D', '[GAT]': 'D', '[GTA]': 'D', '[TAG]': 'D', '[TGA]': 'D',
            '[ACT]': 'H', '[ATC]': 'H', '[CAT]': 'H', '[CTA]': 'H', '[TAC]': 'H', '[TCA]': 'H',
            '[ACG]': 'V', '[AGC]': 'V', '[CAG]': 'V', '[CGA]': 'V', '[GAC]': 'V', '[GCA]': 'V',
            '.': 'N'
        }
        
        # Simple conversion
        result = regex
        for pattern, iupac in conversions.items():
            result = result.replace(pattern, iupac)
        
        # Check if fully converted
        if re.match(r'^[ACGTRYWSKMDHVBN]+$', result):
            return result
        
        return None 