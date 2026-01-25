---
name: bio-ribo-seq-riboseq-preprocessing
description: Preprocess ribosome profiling data including adapter trimming, size selection, rRNA removal, and alignment. Use when preparing Ribo-seq reads for downstream analysis of translation.
tool_type: cli
primary_tool: bowtie2
---

# Ribo-seq Preprocessing

## Workflow Overview

```
Raw Ribo-seq FASTQ
    |
    v
Adapter trimming (cutadapt)
    |
    v
Size selection (25-35 nt typical)
    |
    v
rRNA removal (SortMeRNA/bowtie2)
    |
    v
Alignment to transcriptome
    |
    v
Quality filtered BAM
```

## Adapter Trimming

```bash
# Trim 3' adapter
cutadapt \
    -a CTGTAGGCACCATCAAT \
    -m 20 \
    -M 40 \
    -o trimmed.fastq.gz \
    input.fastq.gz
```

## Size Selection

```bash
# Select ribosome footprint size range
# Typical: 28-32 nt (protected by ribosome)
cutadapt \
    -m 28 \
    -M 32 \
    -o size_selected.fastq.gz \
    trimmed.fastq.gz
```

## rRNA Removal

```bash
# Option 1: SortMeRNA (comprehensive)
sortmerna \
    --ref rRNA_databases/silva-bac-16s-id90.fasta \
    --ref rRNA_databases/silva-euk-18s-id95.fasta \
    --ref rRNA_databases/silva-euk-28s-id98.fasta \
    --reads size_selected.fastq.gz \
    --aligned rRNA_reads \
    --other non_rRNA_reads \
    --fastx \
    --threads 8

# Option 2: Bowtie2 to rRNA index
bowtie2 -x rRNA_index \
    -U size_selected.fastq.gz \
    --un non_rRNA.fastq.gz \
    -S /dev/null \
    -p 8
```

## Alignment to Transcriptome

```bash
# STAR alignment (spliced)
STAR --runMode alignReads \
    --genomeDir STAR_index \
    --readFilesIn non_rRNA.fastq.gz \
    --readFilesCommand zcat \
    --outFilterMultimapNmax 1 \
    --outFilterMismatchNmax 2 \
    --alignIntronMax 1 \
    --outSAMtype BAM SortedByCoordinate \
    --outFileNamePrefix riboseq_

# Or bowtie2 to transcriptome
bowtie2 -x transcriptome_index \
    -U non_rRNA.fastq.gz \
    -S aligned.sam \
    --no-unal \
    -p 8
```

## Quality Metrics

```bash
# Check read length distribution
samtools view aligned.bam | \
    awk '{print length($10)}' | \
    sort | uniq -c | sort -k2n

# Expected: Peak at 28-30 nt

# Check mapping rate
samtools flagstat aligned.bam
```

## Python Preprocessing

```python
import pysam
import numpy as np
from collections import Counter

def get_length_distribution(bam_path):
    '''Get read length distribution from BAM'''
    lengths = Counter()
    with pysam.AlignmentFile(bam_path, 'rb') as bam:
        for read in bam:
            if not read.is_unmapped:
                lengths[read.query_length] += 1
    return lengths

def filter_by_length(bam_in, bam_out, min_len=28, max_len=32):
    '''Filter BAM by read length'''
    with pysam.AlignmentFile(bam_in, 'rb') as infile:
        with pysam.AlignmentFile(bam_out, 'wb', template=infile) as outfile:
            for read in infile:
                if min_len <= read.query_length <= max_len:
                    outfile.write(read)
```

## Related Skills

- ribosome-periodicity - Validate preprocessing quality
- read-qc - General quality control
- read-alignment - Alignment concepts
