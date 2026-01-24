# Microbiome Pipeline Usage Guide

## Overview

Complete 16S rRNA amplicon sequencing workflow from raw FASTQ reads to differential abundance testing using compositionally-aware methods.

## Prerequisites

```r
BiocManager::install(c('dada2', 'phyloseq', 'ALDEx2'))
install.packages(c('vegan', 'ggplot2'))
```

## Quick Start

Tell your AI agent what you want to do:
- "Run the microbiome pipeline on my 16S FASTQ files"
- "Process my amplicon data with DADA2 and run diversity analysis"
- "Find differentially abundant taxa between my groups"

## Example Prompts

### Basic Analysis
> "I have 16S rRNA sequencing data, run the full pipeline"

> "Process my paired-end amplicon reads and assign taxonomy"

### Diversity Analysis
> "Calculate alpha and beta diversity for my microbiome samples"

> "Run PERMANOVA to test if microbial communities differ between groups"

### Differential Abundance
> "Find taxa that differ between treatment and control using ALDEx2"

> "Run compositional differential abundance analysis on my microbiome data"

## Pipeline Stages

### 1. Quality Filtering (DADA2)
- Inspect quality profiles
- Trim to quality thresholds
- Remove PhiX contamination

### 2. ASV Inference
- Learn error rates from data
- Denoise to exact sequences
- Merge forward/reverse pairs

### 3. Chimera Removal
- Consensus-based detection
- Remove bimeric sequences

### 4. Taxonomy Assignment
- SILVA 138 classifier
- Species-level matching
- Confidence filtering

### 5. Diversity Analysis
- Alpha: Richness, Shannon, Simpson
- Beta: Bray-Curtis, UniFrac
- PERMANOVA for group testing

### 6. Differential Abundance
- ALDEx2 for compositionality
- Effect size + FDR filtering
- Visualize significant taxa

## Input Requirements

### FASTQ Files
```
raw_reads/
├── Sample1_R1_001.fastq.gz
├── Sample1_R2_001.fastq.gz
├── Sample2_R1_001.fastq.gz
└── ...
```

### Metadata
```csv
SampleID,Group,Subject,Timepoint
Sample1,Control,S1,T0
Sample2,Treatment,S2,T0
```

### Reference Databases
- SILVA 138.1 training set
- SILVA species assignment
- Download from: https://zenodo.org/record/4587955

## Parameter Selection

### truncLen
Based on quality profiles:
- V4 (515F-806R): c(240, 160) typical
- V3-V4: c(280, 200) typical
- Ensure 20bp overlap for merging

### maxEE
Maximum expected errors:
- Strict: c(2, 2)
- Permissive: c(5, 5)

## Expected Results

| Metric | Typical Range |
|--------|---------------|
| ASVs | 500-5000 |
| Reads/sample | 10k-100k |
| Genus assignment | 70-90% |
| Differential taxa | 5-20% |

## Common Issues

### Low merge rate
- Amplicon too long for read length
- Poor reverse read quality
- Solution: Adjust truncLen or use longer reads

### Many unassigned
- Novel taxa not in database
- Poor amplification
- Solution: Try GTDB database

### No differential taxa
- Insufficient replication
- High inter-individual variation
- Solution: Increase n, reduce FDR stringency

## Tips

- **Primer trimming**: Remove primer sequences before DADA2 processing
- **truncLen**: Set based on quality profiles, ensure 20bp overlap for merging
- **Replicates**: Microbiome studies typically need more replicates (n>=5) due to high variability
- **Compositionality**: Use ALDEx2 or ANCOM-BC, not standard differential tests
- **Reference database**: SILVA 138 or GTDB for taxonomy assignment
