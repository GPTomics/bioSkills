# HISAT2 RNA-seq Alignment Usage Guide

## Overview

HISAT2 is a memory-efficient splice-aware aligner for RNA-seq data. It uses a graph-based index that enables fast alignment with low memory usage (~8GB for human genome vs ~30GB for STAR).

## Installation

```bash
conda install -c bioconda hisat2 samtools
```

## Quick Start

```bash
# 1. Build index with annotation
hisat2_extract_splice_sites.py genes.gtf > splicesites.txt
hisat2_extract_exons.py genes.gtf > exons.txt
hisat2-build -p 8 --ss splicesites.txt --exon exons.txt genome.fa hs2_index

# 2. Align reads
hisat2 -p 8 -x hs2_index -1 sample_R1.fq.gz -2 sample_R2.fq.gz | \
    samtools sort -o sample.bam -
samtools index sample.bam
```

## Pre-built Indices

HISAT2 provides pre-built indices with annotation:
- http://daehwankimlab.github.io/hisat2/download/

These include splice sites and exons, saving index build time.

## Complete Pipeline

```bash
# Variables
GENOME=genome.fa
GTF=genes.gtf
INDEX=hisat2_index
R1=sample_R1.fq.gz
R2=sample_R2.fq.gz
OUT=sample
THREADS=8

# Build index if needed
if [ ! -f "${INDEX}.1.ht2" ]; then
    hisat2_extract_splice_sites.py $GTF > splicesites.txt
    hisat2_extract_exons.py $GTF > exons.txt
    hisat2-build -p $THREADS --ss splicesites.txt --exon exons.txt $GENOME $INDEX
fi

# Align with strandedness (most Illumina is RF)
hisat2 -p $THREADS -x $INDEX \
    --rna-strandness RF \
    -1 $R1 -2 $R2 2> ${OUT}_hisat2.log | \
    samtools sort -@ 4 -o ${OUT}.bam -

samtools index ${OUT}.bam
```

## Strandedness Detection

Determine strandedness from your library prep:

| Kit/Method | Strandedness | HISAT2 Flag |
|------------|--------------|-------------|
| Unstranded | None | (default) |
| TruSeq Stranded | Reverse | --rna-strandness RF |
| dUTP method | Reverse | --rna-strandness RF |
| Ligation method | Forward | --rna-strandness FR |

Use RSeQC's `infer_experiment.py` if unknown:
```bash
infer_experiment.py -r genes.bed -i aligned.bam
```

## Two-Pass Mode for Novel Junctions

Unlike STAR, HISAT2 doesn't have built-in two-pass mode, but you can implement it:

```bash
# First pass: collect novel junctions
for sample in sample1 sample2 sample3; do
    hisat2 -p 8 -x index \
        --novel-splicesite-outfile ${sample}_junctions.txt \
        -1 ${sample}_R1.fq.gz -2 ${sample}_R2.fq.gz \
        -S /dev/null 2> ${sample}_pass1.log
done

# Combine junctions (filter by read support)
cat *_junctions.txt | \
    awk '{key=$1"\t"$2"\t"$3"\t"$4; count[key]+=$5} END {for(k in count) if(count[k]>=5) print k"\t"count[k]}' \
    > combined_junctions.txt

# Second pass: align with combined junctions
for sample in sample1 sample2 sample3; do
    hisat2 -p 8 -x index \
        --novel-splicesite-infile combined_junctions.txt \
        -1 ${sample}_R1.fq.gz -2 ${sample}_R2.fq.gz | \
        samtools sort -@ 4 -o ${sample}.bam -
    samtools index ${sample}.bam
done
```

## Downstream Analysis

### For featureCounts
```bash
# Coordinate-sorted BAM (default)
hisat2 -p 8 -x index -1 r1.fq.gz -2 r2.fq.gz | samtools sort -o aligned.bam -
featureCounts -p -a genes.gtf -o counts.txt aligned.bam
```

### For StringTie
```bash
# Use --dta flag
hisat2 -p 8 -x index --dta -1 r1.fq.gz -2 r2.fq.gz | samtools sort -o aligned.bam -
stringtie aligned.bam -G genes.gtf -o transcripts.gtf
```

## HISAT2 vs STAR Decision Guide

| Factor | Choose HISAT2 | Choose STAR |
|--------|---------------|-------------|
| Memory | Limited (<16GB) | Available (>32GB) |
| Built-in counting | Not needed | Needed |
| Fusion detection | Not needed | Needed |
| Two-pass mode | Manual OK | Want built-in |
| Accuracy | Similar | Similar |
| Speed | Fast | Faster |

## Troubleshooting

### Low Alignment Rate
1. Check read quality
2. Verify reference genome version
3. Check strandedness setting
4. Look for contamination

### Index Building Fails
```bash
# Check available memory
free -h

# Try with fewer threads
hisat2-build -p 4 genome.fa index
```

### Interpreting Alignment Summary
```
50000000 reads; of these:
  50000000 (100.00%) were paired; of these:
    5000000 (10.00%) aligned concordantly 0 times  # Check if high
    43000000 (86.00%) aligned concordantly exactly 1 time  # Good
    2000000 (4.00%) aligned concordantly >1 times  # Normal
90.00% overall alignment rate  # Target >80%
```

If concordant 0 times is high:
- Wrong reference
- Poor quality reads
- Wrong strandedness setting
