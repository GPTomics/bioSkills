# Adapter Trimming - Usage Guide

## Overview

Sequencing adapters must be removed before alignment to prevent misalignment and artifacts. Cutadapt offers precise control with error-tolerant matching, while Trimmomatic provides efficient paired-end handling with palindrome mode.

## When to Use This Skill

- FastQC shows adapter contamination
- Reads contain sequencing artifacts
- Preparing reads for alignment
- Trimming custom primers or barcodes

## Installation

```bash
# Cutadapt
pip install cutadapt
# or
conda install -c bioconda cutadapt

# Trimmomatic
conda install -c bioconda trimmomatic
```

## Choosing a Tool

| Feature | Cutadapt | Trimmomatic |
|---------|----------|-------------|
| Flexibility | High | Medium |
| Paired-end | Good | Excellent (palindrome) |
| Speed | Fast | Fast |
| Error tolerance | Configurable | Fixed |
| Documentation | Excellent | Good |

**Recommendation**: Use Cutadapt for most cases. Use Trimmomatic when you need palindrome mode for paired-end data with adapter read-through.

## Common Workflows

### Illumina TruSeq (Cutadapt)

```bash
cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
         -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
         -m 20 -j 4 \
         -o trimmed_R1.fastq.gz -p trimmed_R2.fastq.gz \
         raw_R1.fastq.gz raw_R2.fastq.gz
```

### Nextera (Cutadapt)

```bash
cutadapt -a CTGTCTCTTATACACATCT -A CTGTCTCTTATACACATCT \
         -m 20 -j 4 \
         -o trimmed_R1.fastq.gz -p trimmed_R2.fastq.gz \
         raw_R1.fastq.gz raw_R2.fastq.gz
```

### TruSeq with Trimmomatic (Palindrome Mode)

```bash
trimmomatic PE -phred33 -threads 4 \
    raw_R1.fastq.gz raw_R2.fastq.gz \
    trimmed_R1.fastq.gz unpaired_R1.fastq.gz \
    trimmed_R2.fastq.gz unpaired_R2.fastq.gz \
    ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10:2:keepBothReads \
    MINLEN:20
```

## Identifying Your Adapter

1. **Check FastQC report** - "Overrepresented sequences" and "Adapter Content"
2. **Check sequencing facility** - They should provide adapter info
3. **Common patterns**:
   - `AGATCGGAAGAGC` - Illumina universal
   - `CTGTCTCTTATACACATCT` - Nextera
   - `TGGAATTCTCGGGTGCCAAGG` - Small RNA

## Troubleshooting

### Adapters Still Present

Increase error rate or reduce overlap requirement:
```bash
cutadapt -a ADAPTER -e 0.15 -O 3 -o out.fq in.fq
```

### Too Many Reads Discarded

Reduce minimum length:
```bash
cutadapt -a ADAPTER -m 15 -o out.fq in.fq
```

### Unknown Adapter

Use FastQC "Overrepresented sequences" to identify, or try common adapters:
```bash
# Try all common Illumina adapters
cutadapt -a AGATCGGAAGAGC \
         -a CTGTCTCTTATACACATCT \
         -a TGGAATTCTCGGGTGCCAAGG \
         -o out.fq in.fq
```

## Resources

- [Cutadapt Documentation](https://cutadapt.readthedocs.io/)
- [Trimmomatic Manual](http://www.usadellab.org/cms/?page=trimmomatic)
- [Illumina Adapter Sequences](https://support.illumina.com/downloads/illumina-adapter-sequences-document-1000000002694.html)
