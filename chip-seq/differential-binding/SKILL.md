---
name: bio-chipseq-differential-binding
description: Differential binding analysis using DiffBind. Compare ChIP-seq peaks between conditions with statistical rigor. Requires replicate samples. Outputs differentially bound regions with fold changes and p-values.
tool_type: r
primary_tool: DiffBind
---

# Differential Binding with DiffBind

## Create Sample Sheet

```r
# Create sample sheet as data frame or CSV
samples <- data.frame(
    SampleID = c('ctrl_1', 'ctrl_2', 'treat_1', 'treat_2'),
    Tissue = c('cell', 'cell', 'cell', 'cell'),
    Factor = c('H3K4me3', 'H3K4me3', 'H3K4me3', 'H3K4me3'),
    Condition = c('control', 'control', 'treatment', 'treatment'),
    Replicate = c(1, 2, 1, 2),
    bamReads = c('ctrl1.bam', 'ctrl2.bam', 'treat1.bam', 'treat2.bam'),
    Peaks = c('ctrl1_peaks.narrowPeak', 'ctrl2_peaks.narrowPeak',
              'treat1_peaks.narrowPeak', 'treat2_peaks.narrowPeak'),
    PeakCaller = c('macs', 'macs', 'macs', 'macs')
)

write.csv(samples, 'samples.csv', row.names = FALSE)
```

## Load Data

```r
library(DiffBind)

# From sample sheet
dba_obj <- dba(sampleSheet = 'samples.csv')

# View summary
dba_obj
```

## Count Reads in Peaks

```r
# Count reads in consensus peaks (DiffBind 3.0+ defaults)
# summits=250 and bUseSummarizeOverlaps=TRUE are now defaults
dba_obj <- dba.count(dba_obj)

# With specific parameters
dba_obj <- dba.count(
    dba_obj,
    summits = 250,         # Re-center peaks around summits (default in 3.0)
    minOverlap = 2         # Peak must be in at least 2 samples
)
```

## Normalize Data

```r
# Normalize (required before analysis)
dba_obj <- dba.normalize(dba_obj)

# Check normalization
dba.normalize(dba_obj, bRetrieve = TRUE)
```

## Set Up Contrast (DiffBind 3.0+)

```r
# Recommended: design formula approach (DiffBind 3.0+)
dba_obj <- dba.contrast(dba_obj, design = '~ Condition')

# Or use categories for automatic contrast
dba_obj <- dba.contrast(dba_obj, categories = DBA_CONDITION)

# Legacy approach (retained for backward compatibility, not recommended)
# dba_obj <- dba.contrast(dba_obj, group1 = dba_obj$masks$control,
#                         group2 = dba_obj$masks$treatment)
```

## Run Differential Analysis

```r
# Analyze with DESeq2 (default)
dba_obj <- dba.analyze(dba_obj, method = DBA_DESEQ2)

# Or with edgeR
dba_obj <- dba.analyze(dba_obj, method = DBA_EDGER)
```

## View Results

```r
# Summary of differential peaks
dba.show(dba_obj, bContrasts = TRUE)

# Retrieve differential binding results
db_results <- dba.report(dba_obj)
db_results
```

## Filter Results

```r
# Get significant peaks (FDR < 0.05, |FC| > 2)
db_sig <- dba.report(dba_obj, th = 0.05, fold = 2)

# Get all results for custom filtering
db_all <- dba.report(dba_obj, th = 1)
```

## Export Results

```r
# To data frame
results_df <- as.data.frame(dba.report(dba_obj, th = 1))

# Export to CSV
write.csv(results_df, 'differential_binding.csv', row.names = FALSE)

# Export to BED
library(rtracklayer)
export(db_sig, 'diff_peaks.bed', format = 'BED')
```

## Visualization

```r
# PCA plot
dba.plotPCA(dba_obj, DBA_CONDITION, label = DBA_ID)

# Correlation heatmap
dba.plotHeatmap(dba_obj)

# MA plot
dba.plotMA(dba_obj)

# Volcano plot
dba.plotVolcano(dba_obj)

# Heatmap of differential peaks
dba.plotHeatmap(dba_obj, contrast = 1, correlations = FALSE)
```

## Venn Diagram of Peaks

```r
# Overlap between conditions
dba.plotVenn(dba_obj, dba_obj$masks$control)
dba.plotVenn(dba_obj, dba_obj$masks$treatment)
```

## Profile Plots

```r
# Average signal profile
profiles <- dba.plotProfile(dba_obj)
```

## Get Consensus Peaks

```r
# Export consensus peakset
consensus <- dba.peakset(dba_obj, bRetrieve = TRUE)
export(consensus, 'consensus_peaks.bed', format = 'BED')
```

## Multi-Factor Design

```r
# With blocking factor (e.g., batch correction)
dba_obj <- dba.contrast(dba_obj, design = '~ Batch + Condition')
dba_obj <- dba.analyze(dba_obj)
```

## DiffBind 3.0 Notes

DiffBind 3.0+ introduced significant changes:
- `dba.normalize()` is now required before analysis
- Default `summits=250` recenters peaks (was FALSE in older versions)
- Use design formulas instead of group1/group2 for contrasts
- Blacklist filtering is applied by default

## Sample Sheet Columns

| Column | Required | Description |
|--------|----------|-------------|
| SampleID | Yes | Unique identifier |
| Tissue | No | Tissue/cell type |
| Factor | No | ChIP target |
| Condition | Yes | Experimental condition |
| Treatment | No | Additional grouping |
| Replicate | Yes | Replicate number |
| bamReads | Yes | Path to BAM file |
| Peaks | Yes | Path to peak file |
| PeakCaller | Yes | macs, bed, narrow |
| bamControl | No | Path to input BAM |

## Key Functions

| Function | Purpose |
|----------|---------|
| dba | Create DBA object |
| dba.count | Count reads in peaks |
| dba.normalize | Normalize counts |
| dba.contrast | Set up comparison |
| dba.analyze | Run differential analysis |
| dba.report | Get results |
| dba.plotPCA | PCA visualization |
| dba.plotMA | MA plot |
| dba.plotHeatmap | Heatmap |

## Related Skills

- peak-calling - Generate input peak files
- peak-annotation - Annotate differential peaks
- differential-expression - Compare with RNA-seq
- pathway-analysis - Functional enrichment
