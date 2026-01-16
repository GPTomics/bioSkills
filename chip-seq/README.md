# chip-seq

## Overview

ChIP-seq analysis using MACS3 for peak calling, ChIPseeker for annotation, and DiffBind for differential binding. Used to study protein-DNA interactions including transcription factor binding and histone modifications.

**Tool type:** mixed | **Primary tools:** MACS3 (CLI), ChIPseeker (R), DiffBind (R)

## Skills

| Skill | Description |
|-------|-------------|
| peak-calling | Call peaks with MACS3 from ChIP-seq BAM files |
| peak-annotation | Annotate peaks to genes with ChIPseeker |
| differential-binding | Differential binding analysis with DiffBind |
| chipseq-visualization | Visualize ChIP-seq data and peaks |

## Example Prompts

- "Call peaks from my ChIP-seq BAM with MACS3"
- "Run MACS3 with input control"
- "Call broad peaks for H3K27me3"
- "Annotate my peaks with nearby genes"
- "Find peaks in promoters vs enhancers"
- "What genes have peaks in their promoter?"
- "Find differential peaks between treatment and control"
- "Compare binding between two conditions with DiffBind"
- "Identify lost and gained peaks"
- "Create a heatmap of peak signal"
- "Plot peak distribution around TSS"
- "Visualize peaks in a genomic region"

## Requirements

```bash
# MACS3
conda install -c bioconda macs3
```

```r
BiocManager::install(c('ChIPseeker', 'DiffBind', 'TxDb.Hsapiens.UCSC.hg38.knownGene'))
```

## Related Skills

- **alignment-files** - BAM preparation before peak calling
- **pathway-analysis** - Functional enrichment of peak-associated genes
- **genome-intervals** - BED file operations for peak regions
- **methylation-analysis** - Combine with DNA methylation data
