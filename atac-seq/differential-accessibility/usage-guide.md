# Differential Accessibility - Usage Guide

## Overview
Identify chromatin regions that change accessibility between conditions using DiffBind or DESeq2, revealing condition-specific regulatory elements.

## Prerequisites
```r
BiocManager::install(c('DiffBind', 'DESeq2', 'ChIPseeker', 'clusterProfiler'))
BiocManager::install('TxDb.Hsapiens.UCSC.hg38.knownGene')
BiocManager::install('org.Hs.eg.db')
```

```bash
conda install -c bioconda bedtools homer
```

## Quick Start
Tell your AI agent what you want to do:
- "Run differential accessibility analysis between treatment and control"
- "Find regions that gain or lose accessibility after drug treatment"

## Example Prompts
### Basic Differential Analysis
> "Run DiffBind analysis to compare chromatin accessibility between treatment and control ATAC-seq samples"

### Multi-Condition Comparison
> "Compare accessibility across WT, KO, and rescue conditions with multiple contrasts"

### Batch-Corrected Analysis
> "Run differential accessibility analysis accounting for batch effects in my experimental design"

### Downstream Annotation
> "Annotate differential peaks with nearby genes and run GO enrichment analysis"

### Motif Enrichment
> "Find enriched TF motifs in regions that gain accessibility after treatment"

## What the Agent Will Do
1. Create consensus peak set across all samples
2. Count reads in peaks for each sample
3. Normalize counts using DiffBind or DESeq2
4. Generate QC visualizations (PCA, correlation heatmap)
5. Run statistical testing for differential accessibility
6. Annotate significant peaks with nearby genes
7. Optionally run motif enrichment analysis

## Sample Sheet Format
```csv
SampleID,Condition,Replicate,bamReads,Peaks,PeakCaller
ctrl_rep1,control,1,ctrl_rep1.bam,ctrl_rep1_peaks.narrowPeak,macs
ctrl_rep2,control,2,ctrl_rep2.bam,ctrl_rep2_peaks.narrowPeak,macs
treat_rep1,treated,1,treat_rep1.bam,treat_rep1_peaks.narrowPeak,macs
treat_rep2,treated,2,treat_rep2.bam,treat_rep2_peaks.narrowPeak,macs
```

## Interpreting Results

### Key Columns
- **Conc**: Mean concentration (log2)
- **Fold**: Log2 fold change
- **p.value**: Raw p-value
- **FDR**: Adjusted p-value

### Standard Thresholds
- Significant: FDR < 0.05, |Fold| > 1 (2-fold change)
- Stringent: FDR < 0.01, |Fold| > 2 (4-fold change)

## Tips
- Use at least 2 biological replicates per condition
- Normalize to reads in peaks (not total library size) for ATAC-seq
- Check PCA and correlation plots to identify batch effects or outliers
- Use design formula (~Batch + Condition) to correct for known batch effects
- Opened regions (positive fold change) often indicate activated regulatory elements
- Follow up with motif analysis to identify TFs driving accessibility changes
