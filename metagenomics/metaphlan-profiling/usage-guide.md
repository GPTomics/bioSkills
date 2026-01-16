# MetaPhlAn Profiling - Usage Guide

## Overview

MetaPhlAn (Metagenomic Phylogenetic Analysis) uses clade-specific marker genes to provide accurate taxonomic profiling of metagenomic samples. It outputs relative abundances that sum to 100% and is particularly accurate at species level.

## When to Use This Skill

- You want accurate species-level abundances
- You need relative abundances (percentages)
- You're comparing microbial communities across samples
- You want a standard profiling method for publication

## MetaPhlAn vs Kraken2

| Feature | MetaPhlAn | Kraken2 |
|---------|-----------|---------|
| Method | Marker genes | K-mers |
| Output | Relative abundance | Read counts |
| Accuracy | Higher at species | Good overall |
| Database | Smaller (~1GB) | Larger (8-50GB) |
| Speed | Slower | Very fast |

## Installation

```bash
conda install -c bioconda metaphlan
```

Database downloads automatically on first run (~1GB).

## Basic Workflow

```bash
# 1. Profile sample
metaphlan sample.fastq.gz \
    --input_type fastq \
    --nproc 8 \
    --output_file profile.txt \
    --mapout sample.map.bz2

# 2. View results
head profile.txt

# 3. Extract species
grep "s__" profile.txt | grep -v "t__"
```

## Processing Multiple Samples

```bash
# Process all samples
for fq in *.fastq.gz; do
    sample=$(basename $fq .fastq.gz)
    metaphlan $fq --input_type fastq --nproc 4 -o profiles/${sample}.txt
done

# Merge into one table
merge_metaphlan_tables.py profiles/*.txt > merged_profiles.txt
```

## Common Issues

### No Database Found

```bash
# Download database manually
metaphlan --install
```

### Low Mapping Rate

Normal for some samples. MetaPhlAn only considers marker genes, so many reads won't map. Check for host contamination.

### Output All Zeros

- Check input file is not empty
- Verify input_type matches file format
- Sample may have very low microbial content

## Output Interpretation

- Values are percentages (relative abundance)
- All values sum to 100% at each taxonomic level
- UNCLASSIFIED means reads didn't match any marker

## Resources

- [MetaPhlAn GitHub](https://github.com/biobakery/MetaPhlAn)
- [MetaPhlAn Wiki](https://github.com/biobakery/MetaPhlAn/wiki)
