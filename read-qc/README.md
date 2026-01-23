# read-qc

## Overview

Read quality control and preprocessing - the first step in any NGS workflow. Covers quality assessment with FastQC/MultiQC, adapter trimming, quality filtering, and contamination screening.

**Tool type:** cli | **Primary tools:** FastQC, MultiQC, fastp, Trimmomatic, Cutadapt

## Skills

| Skill | Description |
|-------|-------------|
| quality-reports | Generate and interpret QC reports with FastQC/MultiQC |
| adapter-trimming | Remove sequencing adapters with Cutadapt/Trimmomatic |
| quality-filtering | Quality/length trimming, N removal with fastp/Trimmomatic |
| fastp-workflow | All-in-one modern preprocessing with fastp |
| contamination-screening | Detect sample contamination with FastQ Screen |
| umi-processing | Extract and deduplicate reads using UMIs with umi_tools |
| rnaseq-qc | RNA-seq specific QC: rRNA, strandedness, gene body coverage |

## Example Prompts

- "How do I check the quality of my FASTQ files?"
- "Run FastQC on all my samples"
- "Generate a combined QC report with MultiQC"
- "Remove Illumina adapters from paired-end reads"
- "Trim adapters from my single-end reads with Cutadapt"
- "Trim low quality bases from the 3' end"
- "Filter reads shorter than 50bp"
- "Run fastp on my RNA-seq data"
- "Use fastp for adapter trimming and quality filtering"
- "Check if my samples have contamination"
- "Screen my reads against multiple genomes"
- "What is the GC content distribution in my reads?"
- "Extract UMIs from my reads and deduplicate"
- "Process single-cell data with UMIs"
- "Check rRNA contamination in my RNA-seq"
- "Verify strandedness of my library"
- "Check gene body coverage"
- "Calculate transcript integrity (TIN)"

## Requirements

```bash
# FastQC and MultiQC
conda install -c bioconda fastqc multiqc

# fastp (recommended all-in-one)
conda install -c bioconda fastp

# Cutadapt
conda install -c bioconda cutadapt

# Trimmomatic
conda install -c bioconda trimmomatic

# FastQ Screen
conda install -c bioconda fastq-screen

# umi_tools
conda install -c bioconda umi_tools

# RNA-seq QC
conda install -c bioconda sortmerna rseqc
```

## Related Skills

- **sequence-io** - FASTQ file reading and writing
- **alignment-files** - Downstream BAM processing
- **metagenomics** - QC before taxonomic classification
