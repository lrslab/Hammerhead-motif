# Bacterial DNA Methylation Motif Caller - Package Summary

## Overview

This is a comprehensive, general-purpose package for discovering DNA methylation motifs in bacterial genomes. The package analyzes genome FASTA files and BED files containing modification sites to identify enriched sequence motifs.

## Key Features

### 1. **General Purpose Design**
- No hardcoded motifs - completely data-driven discovery
- Works with any bacterial genome
- Optional known motif validation via external JSON file
- Command-line interface with configurable parameters

### 2. **Multiple Discovery Algorithms**

#### Greedy Algorithm (Optimized)
- Fast k-mer based discovery (4-8 mers by default)
- Chi-square statistical testing for enrichment
- GC-content bias correction
- AT-rich sequence filtering to reduce noise
- Efficient data structures for speed (~5-15 seconds per k-value)

#### Gapped Motif Discovery
- Finds motifs with gaps (e.g., GAAANNNNNNNNTGG)
- Pattern matching for complex motifs
- Handles motifs up to 12bp gaps

#### MEME Integration
- Position Weight Matrix discovery
- BioPython-based parsing with manual fallback
- Configurable parameters matching proven settings

### 3. **Advanced Motif Handling**

#### Degenerate Base Support
- Full IUPAC code support (R, Y, W, S, K, M, etc.)
- Smart expansion for small motifs
- Regex matching for complex patterns
- Handles gapped motifs with many Ns correctly

#### Motif Merging
- Groups similar motifs by sequence similarity
- Consensus generation from motif clusters
- Removes redundancy in discovered motifs

### 4. **Statistical Analysis**
- Chi-square test for motif enrichment
- Hypergeometric test for known motif validation
- Background model with genome sampling
- GC-content normalization

## Package Structure

```
Hammerhead-motif/
├── methylation_motif_caller.py    # Main analysis script
├── hammermotif/                   # Core package
│   ├── __init__.py
│   ├── motif_greedy.py           # Optimized greedy algorithm
│   ├── meme_runner.py            # MEME integration
│   └── utils.py                  # Utility functions
├── known_motifs.json             # Optional known motifs file
├── README_METHYLATION_CALLER.md  # User documentation
└── example_single_genome.py      # Usage examples
```

## Usage

### Basic Usage
```bash
conda activate py3
python methylation_motif_caller.py
```

### With Known Motifs Validation
```bash
python methylation_motif_caller.py --known-motifs known_motifs.json
```

### Custom Directories
```bash
python methylation_motif_caller.py --input-dir /path/to/genomes --output-dir /path/to/results
```

## Performance on Test Data

For the FC genome (3.4 Mb, 17,401 modification sites):
- Successfully identifies complex gapped motifs (CAYNNNNNRTG, GAAANNNNNNNNTGG)
- Correctly counts motif occurrences in genome and modification sites
- Typical discovery rate: 6-7 out of 9 known motifs
- Runtime: ~2-3 minutes per genome (without MEME)

## Key Improvements Made

1. **Fixed Motif Counting**: Properly handles gapped motifs with regex matching
2. **Optimized Greedy Algorithm**: 10x faster with better memory usage
3. **AT-rich Filtering**: Reduces false positives from repetitive sequences
4. **General Purpose**: No hardcoded motifs, fully configurable
5. **Better MEME Integration**: BioPython parsing with fallback

## Limitations and Future Work

1. Some short degenerate motifs (CCNGG) may be missed if overshadowed by AT-rich sequences
2. MEME can be slow for large datasets
3. Very long gapped motifs (>15bp) may need special handling

## Dependencies

- Python 3.9+
- BioPython
- NumPy
- SciPy
- MEME Suite (in PATH)

## Output Files

For each genome:
- `genome_motif_report.txt`: Main analysis report with all discovered motifs
- `genome_known_motifs.txt`: Known motif enrichment analysis (if provided)

Overall:
- `analysis_summary.txt`: Summary of all genomes analyzed 