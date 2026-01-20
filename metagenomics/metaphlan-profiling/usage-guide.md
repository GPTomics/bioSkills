# MetaPhlAn Profiling - Usage Guide

## Overview
MetaPhlAn uses clade-specific marker genes to provide accurate taxonomic profiling of metagenomic samples, outputting relative abundances that sum to 100%.

## Prerequisites
```bash
conda install -c bioconda metaphlan

# Database downloads automatically on first run (~1GB)
# Or download manually:
metaphlan --install
```

## Quick Start
Tell your AI agent what you want to do:
- "Profile the taxonomic composition of my metagenome"
- "Get species-level abundances for my microbial community"
- "Merge MetaPhlAn profiles from multiple samples"

## Example Prompts
### Basic Profiling
> "Run MetaPhlAn on sample.fastq.gz and output species abundances"

> "Profile my metagenome and save the intermediate mapping file"

### Multi-sample Analysis
> "Process all fastq files through MetaPhlAn and merge into one table"

> "Create a merged abundance table for downstream visualization"

### Specific Taxonomic Levels
> "Extract only species-level results from my MetaPhlAn output"

> "Get genus-level abundances from the merged table"

### Comparison with Kraken2
> "I have both MetaPhlAn and Kraken2 results - help me compare them"

## What the Agent Will Do
1. Verify MetaPhlAn installation and database availability
2. Run profiling with appropriate parameters for input format
3. Generate taxonomic profile with relative abundances
4. Merge multiple samples if requested
5. Filter to specific taxonomic levels as needed

## Tips
- MetaPhlAn outputs relative abundance (all values sum to 100% at each level)
- Low mapping rate is normal - only marker genes are targeted
- Save the mapping file (`--mapout`) for faster re-analysis
- UNCLASSIFIED means reads didn't match any marker gene
- Use `merge_metaphlan_tables.py` to combine multiple profiles

## MetaPhlAn vs Kraken2

| Feature | MetaPhlAn | Kraken2 |
|---------|-----------|---------|
| Method | Marker genes | K-mers |
| Output | Relative abundance | Read counts |
| Accuracy | Higher at species | Good overall |
| Database | Smaller (~1GB) | Larger (8-50GB) |
| Speed | Slower | Very fast |

## Common Issues

### No Database Found
```bash
metaphlan --install
```

### Low Mapping Rate
Normal for some samples. MetaPhlAn only considers marker genes. Check for host contamination if rate is very low.

### Output All Zeros
- Check input file is not empty
- Verify `--input_type` matches file format
- Sample may have very low microbial content

## Resources
- [MetaPhlAn GitHub](https://github.com/biobakery/MetaPhlAn)
- [MetaPhlAn Wiki](https://github.com/biobakery/MetaPhlAn/wiki)
