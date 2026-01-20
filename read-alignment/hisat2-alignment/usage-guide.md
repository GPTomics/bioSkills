# HISAT2 RNA-seq Alignment - Usage Guide

## Overview
HISAT2 is a memory-efficient splice-aware aligner for RNA-seq data. It uses a graph-based index that enables fast alignment with low memory usage (~8GB for human genome vs ~30GB for STAR).

## Prerequisites
```bash
conda install -c bioconda hisat2 samtools
```

## Quick Start
Tell your AI agent what you want to do:
- "Align my RNA-seq reads using HISAT2"
- "Build a HISAT2 index with splice site annotation"
- "Run HISAT2 with the correct strandedness for my TruSeq library"

## Example Prompts

### Index Building
> "Build a HISAT2 index for genome.fa with splice sites from genes.gtf"

> "Create HISAT2 index with exon and splice site annotations"

### Basic Alignment
> "Align my paired-end RNA-seq reads to the HISAT2 index"

> "Run HISAT2 alignment on sample_R1.fq.gz and sample_R2.fq.gz"

### Stranded Libraries
> "Align RNA-seq reads with reverse strandedness (TruSeq/dUTP library)"

> "Run HISAT2 with --rna-strandness RF for my stranded library"

### Two-Pass Mode
> "Collect novel junctions from all samples and realign with combined junctions"

> "Run two-pass HISAT2 alignment for better novel junction detection"

### Downstream Tools
> "Align reads for StringTie transcript assembly using the --dta flag"

> "Run HISAT2 alignment optimized for featureCounts"

## What the Agent Will Do
1. Extract splice sites and exons from GTF if building index
2. Build or locate HISAT2 index
3. Align reads with appropriate strandedness setting
4. Apply downstream-tool-specific flags (--dta for StringTie)
5. Sort and index output BAM file

## Tips
- Use pre-built indices from HISAT2 website to save time
- Set --rna-strandness correctly (RF for most Illumina stranded libraries)
- Use --dta flag if planning to use StringTie downstream
- HISAT2 is ideal for memory-constrained systems (<16GB RAM)
- Use RSeQC's infer_experiment.py if strandedness is unknown

## Strandedness Settings

| Kit/Method | HISAT2 Flag |
|------------|-------------|
| Unstranded | (default) |
| TruSeq Stranded | --rna-strandness RF |
| dUTP method | --rna-strandness RF |
| Ligation method | --rna-strandness FR |

## Complete Pipeline

```bash
GENOME=genome.fa
GTF=genes.gtf
INDEX=hisat2_index
R1=sample_R1.fq.gz
R2=sample_R2.fq.gz
OUT=sample
THREADS=8

if [ ! -f "${INDEX}.1.ht2" ]; then
    hisat2_extract_splice_sites.py $GTF > splicesites.txt
    hisat2_extract_exons.py $GTF > exons.txt
    hisat2-build -p $THREADS --ss splicesites.txt --exon exons.txt $GENOME $INDEX
fi

hisat2 -p $THREADS -x $INDEX \
    --rna-strandness RF \
    -1 $R1 -2 $R2 2> ${OUT}_hisat2.log | \
    samtools sort -@ 4 -o ${OUT}.bam -

samtools index ${OUT}.bam
```

## HISAT2 vs STAR Decision Guide

| Factor | Choose HISAT2 | Choose STAR |
|--------|---------------|-------------|
| Memory | Limited (<16GB) | Available (>32GB) |
| Built-in counting | Not needed | Needed |
| Fusion detection | Not needed | Needed |
| Two-pass mode | Manual OK | Want built-in |

## Troubleshooting

### Low Alignment Rate
- Check read quality
- Verify reference genome version
- Check strandedness setting
- Look for contamination

### Interpreting Alignment Summary
```
50000000 reads; of these:
  43000000 (86.00%) aligned concordantly exactly 1 time  # Good
  5000000 (10.00%) aligned concordantly 0 times  # Check if high
90.00% overall alignment rate  # Target >80%
```
