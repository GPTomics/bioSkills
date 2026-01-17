# STAR RNA-seq Alignment Usage Guide

## Overview

STAR (Spliced Transcripts Alignment to a Reference) is the most widely used RNA-seq aligner. It's extremely fast, supports splice-aware alignment, and can detect novel junctions with two-pass mode.

## Installation

```bash
conda install -c bioconda star samtools
```

## Quick Start

```bash
# 1. Generate index (once per reference/annotation)
STAR --runMode genomeGenerate \
    --runThreadN 8 \
    --genomeDir star_index/ \
    --genomeFastaFiles genome.fa \
    --sjdbGTFfile genes.gtf \
    --sjdbOverhang 99  # Read length - 1

# 2. Align reads
STAR --runThreadN 8 \
    --genomeDir star_index/ \
    --readFilesIn sample_R1.fq.gz sample_R2.fq.gz \
    --readFilesCommand zcat \
    --outFileNamePrefix sample_ \
    --outSAMtype BAM SortedByCoordinate

# 3. Index BAM
samtools index sample_Aligned.sortedByCoord.out.bam
```

## Index Considerations

### sjdbOverhang
Set to (read length - 1). For mixed read lengths, use the longest:
- 50bp reads: `--sjdbOverhang 49`
- 100bp reads: `--sjdbOverhang 99`
- 150bp reads: `--sjdbOverhang 149`

### Memory Requirements
- Human genome: ~30GB RAM for index generation, ~30GB for alignment
- Smaller genomes: proportionally less

### Pre-built Indices
ENCODE and other projects provide pre-built indices:
- https://www.encodeproject.org/references/

## Complete RNA-seq Pipeline

```bash
# Variables
GENOME=genome.fa
GTF=genes.gtf
INDEX_DIR=star_index
R1=sample_R1.fq.gz
R2=sample_R2.fq.gz
PREFIX=sample_
THREADS=16

# Generate index if needed
if [ ! -d "$INDEX_DIR" ]; then
    mkdir -p $INDEX_DIR
    STAR --runMode genomeGenerate \
        --runThreadN $THREADS \
        --genomeDir $INDEX_DIR \
        --genomeFastaFiles $GENOME \
        --sjdbGTFfile $GTF \
        --sjdbOverhang 149
fi

# Two-pass alignment with counting
STAR --runThreadN $THREADS \
    --genomeDir $INDEX_DIR \
    --readFilesIn $R1 $R2 \
    --readFilesCommand zcat \
    --outFileNamePrefix $PREFIX \
    --outSAMtype BAM SortedByCoordinate \
    --twopassMode Basic \
    --quantMode GeneCounts \
    --outSAMunmapped Within \
    --outSAMattributes NH HI AS NM MD

# Index BAM
samtools index ${PREFIX}Aligned.sortedByCoord.out.bam
```

## Understanding Two-Pass Mode

### Single Pass (Default)
- Uses only annotated junctions from GTF
- Faster
- Good when annotation is comprehensive

### Two-Pass Basic
- First pass: discovers novel junctions
- Second pass: realigns with all junctions
- Better for novel transcript discovery
- Recommended for most RNA-seq

```bash
--twopassMode Basic
```

### Two-Pass Per Sample
For population studies, first collect junctions from all samples:
```bash
# Pass 1: Collect junctions from all samples
for sample in *.fq.gz; do
    STAR --genomeDir index --readFilesIn $sample \
        --outFileNamePrefix ${sample}_ --outSAMtype None
done

# Combine junctions
cat *SJ.out.tab | awk '$7 > 2' > filtered_junctions.tab

# Pass 2: Realign with combined junctions
STAR --genomeDir index --sjdbFileChrStartEnd filtered_junctions.tab ...
```

## Stranded Libraries

| Library Type | STAR quantMode Column |
|--------------|----------------------|
| Unstranded | Column 2 |
| Forward (Ligation) | Column 3 |
| Reverse (dUTP, TruSeq) | Column 4 |

Most Illumina libraries are reverse stranded (column 4).

## Interpreting Log.final.out

Key metrics to check:
```
                     Number of input reads: 50000000
          Uniquely mapped reads number: 45000000
               Uniquely mapped reads %: 90%
        % of reads mapped to multiple loci: 5%
        % of reads unmapped: too short: 3%
```

Good alignment:
- Uniquely mapped > 70-80%
- Multi-mapped < 10%
- Unmapped too short < 10%

## Troubleshooting

### Low Unique Mapping
- Check read quality
- Verify correct reference genome
- Check for contamination

### High Multi-mapping
- Normal for genes with paralogs
- Consider using featureCounts with `-M --fraction`

### Memory Errors
```bash
# Reduce sorting memory
--limitBAMsortRAM 10000000000

# Or output unsorted BAM
--outSAMtype BAM Unsorted
# Then sort with samtools
samtools sort -m 4G -@ 8 Aligned.out.bam -o Aligned.sorted.bam
```

### Slow Performance
```bash
# Use shared memory for multiple samples
STAR --genomeLoad LoadAndKeep ...
```

## STAR vs HISAT2

| Feature | STAR | HISAT2 |
|---------|------|--------|
| Speed | Very fast | Fast |
| Memory | High (~30GB) | Low (~8GB) |
| Two-pass | Built-in | Manual |
| Quantification | Built-in | No |
| Fusion detection | Yes | No |

STAR is preferred when memory is available; HISAT2 for memory-constrained systems.
