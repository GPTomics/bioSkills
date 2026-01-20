# Kraken2 Classification - Usage Guide

## Overview
Kraken2 is a fast taxonomic classifier that uses exact k-mer matches to assign reads to taxonomic nodes. It's highly accurate for well-represented taxa and ideal for screening large datasets.

## Prerequisites
```bash
conda install -c bioconda kraken2

# Download pre-built database
wget https://genome-idx.s3.amazonaws.com/kraken/k2_standard_08gb_20230605.tar.gz
mkdir kraken2_db && tar -xzf k2_standard_08gb_20230605.tar.gz -C kraken2_db

# Or build standard database (~100GB disk, ~50GB final)
kraken2-build --standard --db kraken2_standard_db --threads 16
```

## Quick Start
Tell your AI agent what you want to do:
- "Classify my metagenomic reads to identify organisms"
- "Run Kraken2 on paired-end reads and generate a report"
- "Filter out human reads from my metagenome"

## Example Prompts
### Basic Classification
> "Run Kraken2 on reads_R1.fastq.gz and reads_R2.fastq.gz using the standard database"

> "Classify my shotgun metagenomics data and show me the top 20 species"

### Database Selection
> "Use the viral database to screen for viral sequences in my sample"

> "Run Kraken2 with the standard-8 database since I have limited RAM"

### Post-classification Analysis
> "Extract only bacterial classifications from my Kraken2 output"

> "Get species-level results and sort by abundance"

## What the Agent Will Do
1. Verify Kraken2 installation and database availability
2. Run classification with appropriate parameters for paired/single reads
3. Generate per-read classifications and taxonomic report
4. Parse and summarize results by taxonomic level

## Tips
- Use `--memory-mapping` for low-memory systems (slower but works)
- Classification rate of 30-70% is normal for environmental samples
- Store databases on SSD for better performance
- Run Bracken after Kraken2 for improved abundance estimates
- Use `--threads` to parallelize (8-16 typically optimal)

## Output Files

| File | Description |
|------|-------------|
| output.kraken | Per-read classifications |
| report.txt | Taxonomic summary |

## Resources
- [Kraken2 GitHub](https://github.com/DerrickWood/kraken2)
- [Pre-built Databases](https://benlangmead.github.io/aws-indexes/k2)
