# Hammerhead Motif Discovery Package

A comprehensive Python package for discovering DNA methylation motifs in bacterial genomes from Hammerhead output.

## Overview

This package provides multiple algorithms for discovering DNA methylation motifs from bacterial genome sequences and modification site data:

- **Greedy k-mer extraction** with chi-square statistics
- **MEME integration** for position weight matrix discovery  
- **Gapped motif discovery** for complex motifs like CAYNNNNNRTG
- **Statistical enrichment analysis** with hypergeometric tests
- **Batch processing** for multiple genomes

## Installation

```bash
conda activate py3  # Or your preferred environment
# Ensure MEME suite is installed and in PATH
```

## Quick Start

```bash
# Analyze all genomes in test/all directory
python methylation_motif_caller.py

# Or analyze a single genome
python methylation_motif_caller.py genome.fa modifications.bed output_dir
```

## Documentation

- **Detailed usage and algorithms**: See [README_METHYLATION_CALLER.md](README_METHYLATION_CALLER.md)
- **Package components**: See [PACKAGE_SUMMARY.md](PACKAGE_SUMMARY.md)
- **Example scripts**: See the `test/` directory

## Package Structure

```
.
├── methylation_motif_caller.py   # Main general-purpose script
├── hammermotif/                  # Core package modules
│   ├── meme_runner.py           # MEME integration
│   ├── motif_greedy.py          # Greedy algorithm
│   └── ...                      # Other modules
└── test/                        # Test data and example scripts
    ├── all/                     # Test genomes and BED files
    ├── example_single_genome.py # Example usage script
    └── test_*.py               # Various test scripts
```

## Key Features

- Discovers both simple (GATC) and complex gapped motifs (CAYNNNNNRTG)
- Handles IUPAC degenerate codes
- Optimized for speed (5-15 seconds per genome)
- Comprehensive statistical analysis
- Detailed reports with enrichment statistics


### get example data and run hammerhead to get the bed file indicating the modifided sites
```bash 
wget https://figshare.com/ndownloader/files/46437190 -O ecoli.fa

# small test data, only have 642 sites
# The data will give 642 modified sequences from BED, which only contains the GATC sites, not the CCWGG sites.

wget https://figshare.com/ndownloader/files/46437193 -O test.fastq.gz

# full test data, have 30000 sites for GATC and 1600 sites for CCWGG
wget https://figshare.com/ndownloader/files/42651082 -O ecoli_re.fastq.gz

# run the full set
#hammerhead --ref ecoli.fa --read test.fastq.gz --min_depth 5 --min_depth_strand 3
hammerhead --ref ecoli.fa --read  ecoli_re.fastq.gz --min_depth 5 --min_depth_strand 3


```
The code will generate the files blow:
```bash
ls
ecoli.fa      enrichment.bed  mapping.mpileup.txt  potential_modification_site.bed  test.fastq.gz
ecoli.fa.fai  enrichment.fa   mapping.sort.bam     potential_modification_site.txt
```


### start the motif enrichment analysis using two files:
-  potential_modification_site.bed 
-  ecoli.fa
```bash


### the ways we tried for the motif enrichment methods
1. Greedy method: give a K, caculate the enrichment fold and P value (or chi-sequre value) for each K-mer, get\
the top K-mer with the highest enrichment fold and lowest P value one by one. Use the merge kmer function to \
merge the similar one to IUPAC code to get the final motif.

2. Include MEME wrapper for motif enrichment 

3. Gibbs sampling method: Gibbs sampling algorithm

4. Homer method: mask and run