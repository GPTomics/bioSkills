# read-alignment

## Overview

Align short reads to reference genomes using standard aligners. Covers DNA alignment with bwa-mem2/bowtie2 and RNA-seq spliced alignment with STAR/HISAT2.

**Tool type:** cli | **Primary tools:** bwa-mem2, bowtie2, STAR, HISAT2

## Skills

| Skill | Description |
|-------|-------------|
| bwa-alignment | Align DNA reads with bwa-mem2 (BWA successor) |
| bowtie2-alignment | Align DNA/RNA reads with Bowtie2 |
| star-alignment | Spliced alignment for RNA-seq with STAR |
| hisat2-alignment | Spliced alignment for RNA-seq with HISAT2 |

## Example Prompts

- "Align my DNA reads to the human genome"
- "Build a bwa-mem2 index for my reference"
- "Run bwa-mem2 on paired-end reads"
- "Align RNA-seq reads with STAR"
- "Generate a STAR genome index with my GTF"
- "Run STAR with two-pass mode"
- "Align reads with bowtie2 in local mode"
- "Run HISAT2 on my RNA-seq data"
- "Create a HISAT2 index with splice sites"
- "Align single-end reads and output sorted BAM"

## Requirements

```bash
# bwa-mem2 (faster successor to bwa)
conda install -c bioconda bwa-mem2

# bowtie2
conda install -c bioconda bowtie2

# STAR
conda install -c bioconda star

# HISAT2
conda install -c bioconda hisat2

# samtools (for BAM processing)
conda install -c bioconda samtools
```

## Related Skills

- **read-qc** - Upstream quality control and trimming
- **alignment-files** - Post-alignment BAM processing
- **rna-quantification** - Downstream quantification of aligned reads
- **variant-calling** - Variant calling from aligned reads
