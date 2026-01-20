# ATAC-seq QC - Usage Guide

## Overview
Assess ATAC-seq library quality using TSS enrichment, FRiP scores, fragment size distributions, mitochondrial contamination, and library complexity metrics.

## Prerequisites
```bash
pip install pysam pyBigWig numpy pandas matplotlib
conda install -c bioconda samtools picard bedtools
```

## Quick Start
Tell your AI agent what you want to do:
- "Run QC on my ATAC-seq BAM file"
- "Calculate TSS enrichment score for my sample"

## Example Prompts
### Comprehensive QC
> "Generate a full QC report for my ATAC-seq sample including TSS enrichment, FRiP, and fragment distribution"

### TSS Enrichment
> "Calculate the TSS enrichment score for my ATAC-seq data using ENCODE standards"

### Fragment Analysis
> "Plot the fragment size distribution from my ATAC-seq BAM to check nucleosome periodicity"

### Library Complexity
> "Calculate library complexity metrics (NRF, PBC1, PBC2) for my ATAC-seq sample"

### Batch QC
> "Run MultiQC to aggregate QC metrics across all my ATAC-seq samples"

## What the Agent Will Do
1. Calculate alignment statistics and mitochondrial contamination fraction
2. Compute FRiP (Fraction of Reads in Peaks) using peak file
3. Generate fragment size distribution plot to verify nucleosome periodicity
4. Calculate TSS enrichment score using gene annotations
5. Assess library complexity (NRF, PBC1, PBC2)
6. Compile metrics into a QC report

## Key QC Metrics

| Metric | What It Measures | Target |
|--------|------------------|--------|
| TSS Enrichment | Signal at promoters | >7 (ENCODE standard) |
| FRiP | Signal in peaks vs background | >0.2 |
| Fragment Size | Nucleosome periodicity | Clear pattern |
| MT Fraction | Contamination | <20% |
| NRF | Library complexity | >0.9 |

## Fragment Size Interpretation

```
Fragment Size    Peak Type
<100 bp         Nucleosome-free regions (NFR)
180-247 bp      Mono-nucleosome
315-473 bp      Di-nucleosome
558-615 bp      Tri-nucleosome
```

Good QC: Clear peaks at ~50bp (NFR), ~200bp (mono), ~400bp (di)
Poor QC: Single peak or no periodicity

## Tips
- TSS enrichment >7 is the ENCODE standard for high-quality ATAC-seq
- High MT fraction (>20%) indicates incomplete nuclear isolation
- Poor fragment distribution may indicate over-transposition or degraded chromatin
- Use MultiQC to compare QC metrics across multiple samples
- Low correlation between replicates may indicate technical issues or mislabeling
