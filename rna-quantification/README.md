# rna-quantification

## Overview

Quantify gene and transcript expression from RNA-seq data. Covers BAM-based counting with featureCounts and alignment-free quantification with Salmon/kallisto, plus import to R for differential expression.

**Tool type:** mixed | **Primary tools:** featureCounts, Salmon, kallisto, tximport

## Skills

| Skill | Description |
|-------|-------------|
| featurecounts-counting | Count reads per gene from BAM files |
| alignment-free-quant | Pseudo-alignment quantification with Salmon/kallisto |
| tximport-workflow | Import transcript estimates to R for DESeq2/edgeR |
| count-matrix-qc | QC and combine sample counts, detect outliers |

## Example Prompts

- "Count reads per gene from my BAM files"
- "Run featureCounts with paired-end data"
- "Count multi-mapping reads fractionally"
- "Quantify transcripts with Salmon"
- "Build a Salmon index from my transcriptome"
- "Run kallisto on paired-end RNA-seq"
- "Import Salmon results into R for DESeq2"
- "Use tximport to prepare data for edgeR"
- "Combine multiple sample counts into a matrix"
- "Check for sample outliers before DE analysis"
- "Generate PCA of my count matrix"
- "Normalize counts for visualization"

## Requirements

```bash
# featureCounts (part of Subread)
conda install -c bioconda subread

# Salmon
conda install -c bioconda salmon

# kallisto
conda install -c bioconda kallisto
```

```r
BiocManager::install(c('tximport', 'tximeta'))
```

## Related Skills

- **read-qc** - Upstream quality control
- **alignment-files** - BAM file processing
- **differential-expression** - Downstream DE analysis
- **genome-intervals** - GTF/GFF annotation handling
