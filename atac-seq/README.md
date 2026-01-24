# atac-seq

## Overview

Analyze ATAC-seq data for chromatin accessibility profiling. Covers peak calling with MACS3, quality control metrics, differential accessibility analysis, and transcription factor footprinting.

**Tool type:** mixed | **Primary tools:** MACS3, DiffBind, chromVAR, TOBIAS

## Skills

| Skill | Description |
|-------|-------------|
| atac-peak-calling | Call accessible regions with MACS3 using ATAC-specific parameters |
| atac-qc | Assess data quality: fragment sizes, TSS enrichment, FRiP |
| differential-accessibility | Find differentially accessible regions between conditions |
| footprinting | Detect transcription factor binding sites within accessible regions |
| nucleosome-positioning | Extract nucleosome positions with NucleoATAC, ATACseqQC |
| motif-deviation | TF motif accessibility variability with chromVAR |

## Example Prompts

- "Call peaks from my ATAC-seq BAM files"
- "Run MACS3 on ATAC-seq data"
- "Check the fragment size distribution"
- "Calculate TSS enrichment score"
- "Find peaks that differ between conditions"
- "Run differential accessibility with DiffBind"
- "Perform TF footprinting analysis"
- "Run TOBIAS for footprinting"
- "Generate ATAC-seq QC report"
- "Separate nucleosome-free and mono-nucleosomal reads"
- "Run chromVAR to find variable TF motifs"
- "Identify differential motif accessibility between conditions"

## Requirements

```bash
# Peak calling
conda install -c bioconda macs3

# QC tools
pip install atacseqqc
conda install -c bioconda deeptools

# Differential accessibility
# R packages: DiffBind, DESeq2

# Footprinting
conda install -c bioconda tobias
```

```r
BiocManager::install(c('DiffBind', 'ATACseqQC', 'ChIPseeker', 'chromVAR', 'motifmatchr'))
```

## Related Skills

- **read-alignment** - Align ATAC-seq reads (bowtie2)
- **chip-seq** - Similar peak-based analysis
- **genome-intervals** - BED file operations
- **alignment-files** - BAM preprocessing
