# Bacterial DNA Methylation Motif Caller

A comprehensive Python package for discovering DNA methylation motifs in bacterial genomes using multiple algorithms and statistical approaches.

## Features

- **Multiple Discovery Algorithms**:
  - Greedy k-mer extraction with chi-square statistics
  - MEME integration for position weight matrix discovery
  - Statistical enrichment analysis

- **Comprehensive Analysis**:
  - Batch processing of multiple bacterial genomes
  - Known motif verification and enrichment testing
  - Motif merging and consensus generation
  - Detailed reporting with statistics

- **Flexible Input/Output**:
  - Supports standard BED format for modification sites
  - Works with FASTA genome files
  - Generates organized output with reports for each genome

## Installation

### Prerequisites

1. Activate the conda environment:
```bash
conda activate py3
```

2. Ensure MEME is installed and in PATH (already configured in the environment)

3. Required Python packages (should be already installed):
- biopython
- numpy
- scipy
- pandas

## Usage

### Quick Start

To analyze all bacterial genomes in the test directory:

```bash
python methylation_motif_caller.py
```

This will:
1. Process all genome/BED file pairs in `./test/all/`
2. Run multiple motif discovery algorithms
3. Generate results in the `results/` directory

### Example Usage

For a complete example of analyzing a single genome (FC genome with known complex motifs):

```bash
cd test
python example_single_genome.py
```

For running the complete analysis on FC genome:

```bash
cd test
python test_fc_complete_analysis.py
```

### Input Files

The package expects pairs of files for each genome:
- Genome FASTA file: `genome_name.fa` or `genome_name.fasta`
- BED file with modification sites: `genome_name.bed`

BED file format:
```
chrom   start   end     score   strand
chr1    1000    1001    0.95    +
chr1    2000    2001    0.87    -
```

### Output Structure

```
results/
├── genome1/
│   ├── genome1_motif_report.txt       # Main analysis report
│   └── genome1_known_motifs.txt       # Known motif enrichment analysis
├── genome2/
│   └── ...
└── analysis_summary.txt               # Overall summary of all genomes
```

## Algorithm Details

### 1. Greedy Motif Discovery

The greedy algorithm iteratively:
- Extracts k-mers from sequences around modification sites
- Calculates chi-square statistics vs. background
- Selects the most enriched motif
- Removes sequences containing that motif
- Repeats until no significant motifs remain

### 2. MEME Integration

Uses MEME (Multiple EM for Motif Elicitation) to:
- Find position weight matrices
- Discover motifs of variable length
- Provide E-values for statistical significance

MEME is run with optimized parameters for bacterial methylation motifs:
- `-mod zoops`: Zero or one occurrence per sequence
- `-objfun classic`: Classic objective function for better sensitivity
- `-markov_order 0`: 0-order Markov background model
- `-minw 4 -maxw 16`: Search for motifs between 4-16 bp
- `-revcomp`: Search both DNA strands
- `-time 14400`: 4-hour time limit for thorough search

### 3. Motif Merging

Combines similar motifs using:
- Hamming distance clustering
- Core sequence extraction
- Longest common substring identification
- IUPAC degenerate code consensus

## Example Results

For the provided bacterial genomes, the package successfully identifies known methylation motifs:

- **E. coli**: GATC (Dam), CCWGG (Dcm)
- **S. aureus**: GATC, GGANNNNNNTGG
- **Salmonella**: GATC, CAGAG
- **K. pneumoniae**: GATC, CCWGG, GTGANNNNNNTGG

## Advanced Usage

### Custom Analysis

```python
from methylation_motif_caller import MethylationMotifCaller

# Initialize caller
caller = MethylationMotifCaller(
    genome_file="path/to/genome.fa",
    bed_file="path/to/modifications.bed",
    output_dir="custom_output"
)

# Run complete analysis
motifs = caller.run_complete_analysis()

# Or run specific steps
caller.run_greedy_discovery(k_values=[4, 5, 6, 7])
caller.run_meme_discovery(nmotifs=20)
caller.merge_and_finalize_motifs()
```

### Parameters

Key parameters can be adjusted:

- `k_values`: List of k-mer sizes to test (default: [4, 5, 6, 7])
- `chi2_threshold`: Minimum chi-square value for significance (default: 100)
- `window_size`: Size of sequence window around modifications (default: 20)
- `nmotifs`: Number of motifs for MEME to find (default: 10)

## Interpreting Results

### Motif Report

The main report includes:
- Genome statistics
- Greedy algorithm results with chi-square values
- MEME results with E-values
- Final merged motifs

### Enrichment Analysis

For known motifs, the package calculates:
- Occurrence counts in modification sites vs. genome
- Enrichment ratio
- P-value using hypergeometric test
- Significance determination

### Quality Indicators

Good motifs typically have:
- Chi-square value > 100
- E-value < 0.01
- Enrichment ratio > 2
- Consistent appearance across algorithms

## Troubleshooting

### Common Issues

1. **No motifs found**:
   - Check BED file has sufficient modification sites
   - Lower chi2_threshold parameter
   - Try different k-mer sizes

2. **MEME fails**:
   - Ensure MEME is properly installed
   - Check sequence quality and length
   - Adjust minw/maxw parameters

3. **Memory issues**:
   - Process genomes individually
   - Reduce window_size parameter
   - Use subset of modification sites

## Citation

If you use this package, please cite:
- The MEME Suite: Bailey et al. (2015) Nucleic Acids Research
- BioPython: Cock et al. (2009) Bioinformatics

## Contact

For questions or issues, please refer to the package documentation or create an issue on the repository. 