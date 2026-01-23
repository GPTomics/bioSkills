# small-rna-seq

## Overview

Analyze small RNA sequencing data including miRNA, piRNA, and snoRNA for discovery, quantification, and differential expression analysis.

**Tool type:** mixed | **Primary tools:** miRDeep2, miRge3, cutadapt, DESeq2

## Skills

| Skill | Description |
|-------|-------------|
| smrna-preprocessing | Adapter trimming and size selection for small RNA reads |
| mirdeep2-analysis | De novo miRNA discovery and quantification |
| mirge3-analysis | Fast miRNA quantification with isomiR detection |
| differential-mirna | Differential expression analysis of small RNAs |
| target-prediction | miRNA target gene prediction |

## Example Prompts

- "Trim adapters from my small RNA-seq FASTQ files"
- "Quantify miRNAs in my samples using miRge3"
- "Discover novel miRNAs with miRDeep2"
- "Find differentially expressed miRNAs between conditions"
- "Predict target genes for my DE miRNAs"
- "Detect isomiRs and A-to-I editing events"
- "Analyze piRNA expression in germline samples"

## Requirements

```bash
# Preprocessing
pip install cutadapt

# miRDeep2 (Perl-based)
conda install -c bioconda mirdeep2

# miRge3
pip install mirge3

# Differential expression (R)
BiocManager::install(c('DESeq2', 'edgeR'))

# Target prediction
conda install -c bioconda miranda
```

## Related Skills

- **read-qc** - General read quality control
- **differential-expression** - Bulk RNA-seq DE analysis framework
- **rna-quantification** - RNA quantification concepts
