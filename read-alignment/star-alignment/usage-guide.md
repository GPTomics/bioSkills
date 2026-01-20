# STAR RNA-seq Alignment - Usage Guide

## Overview
STAR (Spliced Transcripts Alignment to a Reference) is the most widely used RNA-seq aligner. It's extremely fast, supports splice-aware alignment, and can detect novel junctions with two-pass mode.

## Prerequisites
```bash
conda install -c bioconda star samtools
```

## Quick Start
Tell your AI agent what you want to do:
- "Align my RNA-seq reads using STAR"
- "Generate a STAR index for my reference genome and GTF annotation"
- "Run two-pass STAR alignment with gene counting"

## Example Prompts

### Index Generation
> "Create a STAR index for genome.fa with genes.gtf annotation for 150bp reads"

> "Generate STAR genome index with appropriate sjdbOverhang for my 100bp reads"

### Basic Alignment
> "Align my paired-end RNA-seq reads to the STAR index"

> "Run STAR alignment on sample_R1.fq.gz and sample_R2.fq.gz"

### Two-Pass Mode
> "Align RNA-seq reads using STAR two-pass mode for novel junction discovery"

> "Run STAR with twopassMode Basic and generate gene counts"

### Multi-Sample Studies
> "Collect splice junctions from all samples and realign with combined junctions"

> "Run two-pass alignment across multiple samples for my population study"

## What the Agent Will Do
1. Generate STAR index if not present (with correct sjdbOverhang)
2. Align reads with appropriate parameters
3. Run two-pass mode if requested for novel junction discovery
4. Generate gene counts if quantMode is enabled
5. Sort and index output BAM file

## Tips
- Set sjdbOverhang to (read length - 1) for optimal junction detection
- Use --twopassMode Basic for most RNA-seq analyses
- STAR requires ~30GB RAM for human genome; use HISAT2 for low-memory systems
- Use --quantMode GeneCounts for built-in gene counting
- Check Log.final.out for alignment quality metrics

## Index Parameters

| Read Length | sjdbOverhang |
|-------------|--------------|
| 50bp | 49 |
| 100bp | 99 |
| 150bp | 149 |

## Complete RNA-seq Pipeline

```bash
GENOME=genome.fa
GTF=genes.gtf
INDEX_DIR=star_index
R1=sample_R1.fq.gz
R2=sample_R2.fq.gz
PREFIX=sample_
THREADS=16

if [ ! -d "$INDEX_DIR" ]; then
    mkdir -p $INDEX_DIR
    STAR --runMode genomeGenerate \
        --runThreadN $THREADS \
        --genomeDir $INDEX_DIR \
        --genomeFastaFiles $GENOME \
        --sjdbGTFfile $GTF \
        --sjdbOverhang 149
fi

STAR --runThreadN $THREADS \
    --genomeDir $INDEX_DIR \
    --readFilesIn $R1 $R2 \
    --readFilesCommand zcat \
    --outFileNamePrefix $PREFIX \
    --outSAMtype BAM SortedByCoordinate \
    --twopassMode Basic \
    --quantMode GeneCounts

samtools index ${PREFIX}Aligned.sortedByCoord.out.bam
```

## Stranded Libraries

| Library Type | STAR quantMode Column |
|--------------|----------------------|
| Unstranded | Column 2 |
| Forward (Ligation) | Column 3 |
| Reverse (dUTP, TruSeq) | Column 4 |

## Troubleshooting

### Low Unique Mapping
- Check read quality
- Verify correct reference genome
- Check for contamination

### Memory Errors
```bash
--limitBAMsortRAM 10000000000
--outSAMtype BAM Unsorted
```

### STAR vs HISAT2

| Feature | STAR | HISAT2 |
|---------|------|--------|
| Memory | High (~30GB) | Low (~8GB) |
| Speed | Very fast | Fast |
| Two-pass | Built-in | Manual |
| Quantification | Built-in | No |
