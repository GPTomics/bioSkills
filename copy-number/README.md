# copy-number

## Overview

Detect and analyze copy number variants (CNVs) from whole genome, exome, or targeted sequencing data. Covers read-depth based CNV detection with CNVkit and GATK, visualization, and functional annotation.

**Tool type:** mixed | **Primary tools:** CNVkit, GATK, Python

## Skills

| Skill | Description |
|-------|-------------|
| cnvkit-analysis | Detect CNVs from targeted/exome sequencing with CNVkit |
| gatk-cnv | CNV calling using GATK best practices workflow |
| cnv-visualization | Plot copy number profiles, segments, and genome-wide views |
| cnv-annotation | Annotate CNVs with genes, pathways, and clinical significance |

## Example Prompts

- "Detect CNVs from my exome BAM files"
- "Run CNVkit on tumor-normal pairs"
- "Build a CNVkit reference from normal samples"
- "Call CNVs with GATK"
- "Plot copy number profile for a sample"
- "Visualize CNVs as a heatmap across samples"
- "Annotate CNVs with affected genes"
- "Find CNVs overlapping cancer driver genes"
- "Compare CNV calls between callers"
- "Filter CNVs by size and quality"

## Requirements

```bash
# CNVkit
conda install -c bioconda cnvkit

# GATK (for GATK CNV workflow)
conda install -c bioconda gatk4

# Python libraries
pip install cnvlib matplotlib pandas seaborn

# R (for some visualizations)
# install.packages(c('ggplot2', 'DNAcopy'))
```

## Related Skills

- **alignment-files** - Input BAM processing
- **variant-calling** - SNV/indel calling (complementary)
- **long-read-sequencing** - SV detection from long reads
- **population-genetics** - Population-level CNV analysis
