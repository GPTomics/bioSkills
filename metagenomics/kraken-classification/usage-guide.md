# Kraken2 Classification - Usage Guide

## Overview

Kraken2 is a fast taxonomic classifier that uses exact k-mer matches to assign reads to taxonomic nodes. It's highly accurate for well-represented taxa and very fast, but requires substantial memory for the database.

## When to Use This Skill

- You have metagenomic reads and want to know what organisms are present
- You need fast classification for screening
- You want to classify reads for downstream assembly filtering
- You're looking for viral or rare taxa (better database coverage)

## Installation

```bash
conda install -c bioconda kraken2
```

## Database Setup

### Option 1: Download Pre-built Database

```bash
# Download from https://benlangmead.github.io/aws-indexes/k2
wget https://genome-idx.s3.amazonaws.com/kraken/k2_standard_08gb_20230605.tar.gz
mkdir kraken2_db && tar -xzf k2_standard_08gb_20230605.tar.gz -C kraken2_db
```

### Option 2: Build Standard Database

```bash
kraken2-build --standard --db kraken2_standard_db --threads 16
```

This downloads RefSeq bacterial, archaeal, viral, and human genomes. Requires ~100GB disk space during build, ~50GB final.

## Basic Workflow

```bash
# 1. Classify reads
kraken2 --db kraken2_db \
    --threads 8 \
    --paired \
    --output classifications.kraken \
    --report kraken_report.txt \
    reads_R1.fastq.gz reads_R2.fastq.gz

# 2. View top taxa
head kraken_report.txt

# 3. Get species only
awk '$4 == "S"' kraken_report.txt | sort -k1 -nr | head -20
```

## Common Issues

### Database Loading Fails

- Check available RAM (`free -h`)
- Use `--memory-mapping` for low-memory systems
- Use smaller database (standard-8)

### Low Classification Rate

- Normal for environmental samples (30-70% is common)
- May indicate novel organisms
- Try different database with broader coverage

### Slow Performance

- Use more threads (`--threads`)
- Ensure database is on fast storage (SSD)
- Pre-load database into memory cache

## Output Files

| File | Description |
|------|-------------|
| output.kraken | Per-read classifications |
| report.txt | Taxonomic summary |

## Next Steps

After Kraken2 classification, typically run Bracken for species abundance estimation.

## Resources

- [Kraken2 GitHub](https://github.com/DerrickWood/kraken2)
- [Pre-built Databases](https://benlangmead.github.io/aws-indexes/k2)
