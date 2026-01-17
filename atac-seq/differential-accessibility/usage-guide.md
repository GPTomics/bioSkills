# Differential Accessibility Usage Guide

## Overview

Differential accessibility analysis identifies chromatin regions that change between conditions. This reveals condition-specific regulatory elements and helps understand how chromatin state relates to gene expression changes.

## Workflow Overview

```
1. Call peaks per sample
2. Create consensus peak set
3. Count reads in peaks
4. Normalize counts
5. Statistical testing
6. Annotate and interpret
```

## Complete DiffBind Pipeline

```r
library(DiffBind)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)

# 1. Create sample sheet
# Required columns: SampleID, Condition, Replicate, bamReads, Peaks
samples <- read.csv('samples.csv')

# 2. Load samples
dba <- dba(sampleSheet=samples)
print(dba)  # Summary of peaks

# 3. Establish consensus peak set and count reads
dba <- dba.count(dba,
    summits=250,       # Re-center on summit, extend 250bp each side
    minOverlap=2,      # Peak must be in >=2 samples
    bParallel=TRUE)    # Use parallel processing

# 4. Normalization (critical step)
dba <- dba.normalize(dba,
    normalize=DBA_NORM_NATIVE,     # Use native DESeq2/edgeR normalization
    library=DBA_LIBSIZE_PEAKREADS) # Normalize to reads in peaks

# 5. QC visualizations
pdf('qc_plots.pdf')
dba.plotPCA(dba, attributes=DBA_CONDITION, label=DBA_ID)
plot(dba)  # Correlation heatmap
dev.off()

# 6. Set up contrast
dba <- dba.contrast(dba,
    contrast=c('Condition', 'treated', 'control'))

# 7. Differential analysis
dba <- dba.analyze(dba, method=DBA_DESEQ2)

# 8. Get results
results <- dba.report(dba, th=1)  # Get all peaks
sig_results <- dba.report(dba, th=0.05, fold=1)  # Significant only

print(paste('Total peaks:', length(results)))
print(paste('Significant (FDR<0.05, |FC|>2):', length(sig_results)))

# 9. Visualizations
pdf('differential_plots.pdf')
dba.plotMA(dba)
dba.plotVolcano(dba)
dba.plotHeatmap(dba, contrast=1, correlations=FALSE)
dev.off()

# 10. Annotate peaks
peakAnno <- annotatePeak(sig_results,
    TxDb=TxDb.Hsapiens.UCSC.hg38.knownGene,
    annoDb='org.Hs.eg.db')

# 11. Save results
results_df <- as.data.frame(results)
results_df$annotation <- NA
results_df$gene <- NA

anno_df <- as.data.frame(peakAnno)
write.csv(anno_df, 'differential_accessibility.csv', row.names=FALSE)
```

## Sample Sheet Format

```csv
SampleID,Condition,Replicate,bamReads,Peaks,PeakCaller
ctrl_rep1,control,1,ctrl_rep1.bam,ctrl_rep1_peaks.narrowPeak,macs
ctrl_rep2,control,2,ctrl_rep2.bam,ctrl_rep2_peaks.narrowPeak,macs
treat_rep1,treated,1,treat_rep1.bam,treat_rep1_peaks.narrowPeak,macs
treat_rep2,treated,2,treat_rep2.bam,treat_rep2_peaks.narrowPeak,macs
```

## Handling Multiple Conditions

```r
# Three-way comparison
samples <- data.frame(
    SampleID = c('wt1', 'wt2', 'ko1', 'ko2', 'rescue1', 'rescue2'),
    Condition = c('WT', 'WT', 'KO', 'KO', 'Rescue', 'Rescue'),
    ...
)

dba <- dba(sampleSheet=samples)
dba <- dba.count(dba)
dba <- dba.normalize(dba)

# Multiple contrasts
dba <- dba.contrast(dba, contrast=c('Condition', 'KO', 'WT'))
dba <- dba.contrast(dba, contrast=c('Condition', 'Rescue', 'KO'))
dba <- dba.analyze(dba)

# Get each contrast
ko_vs_wt <- dba.report(dba, contrast=1)
rescue_vs_ko <- dba.report(dba, contrast=2)
```

## Batch Correction

```r
# Include batch in design
samples$Batch <- factor(c('A', 'B', 'A', 'B', ...))

dba <- dba(sampleSheet=samples)
dba <- dba.count(dba)
dba <- dba.normalize(dba)

# Use design formula
dba <- dba.contrast(dba, design='~Batch + Condition')
dba <- dba.analyze(dba)
```

## Alternative: DESeq2 Directly

For more control over the analysis:

```r
library(DESeq2)

# Count reads in consensus peaks
# Use featureCounts or bedtools multicov

counts <- read.delim('peak_counts.txt', row.names=1)
coldata <- read.csv('sample_info.csv', row.names=1)

# Ensure matching order
counts <- counts[, rownames(coldata)]

dds <- DESeqDataSetFromMatrix(
    countData=round(counts),
    colData=coldata,
    design=~condition)

# Filter low counts
keep <- rowSums(counts(dds) >= 10) >= min(table(coldata$condition))
dds <- dds[keep,]

# Run analysis
dds <- DESeq(dds)
res <- results(dds, contrast=c('condition', 'treated', 'control'))

# Shrink fold changes for visualization
res_shrunk <- lfcShrink(dds, coef='condition_treated_vs_control', type='apeglm')
```

## Interpreting Results

### Key Columns in DiffBind Output
- **Conc**: Mean concentration (log2)
- **Conc_control/treated**: Mean in each condition
- **Fold**: Log2 fold change
- **p.value**: Raw p-value
- **FDR**: Adjusted p-value

### Thresholds
```r
# Standard thresholds
sig <- results[results$FDR < 0.05 & abs(results$Fold) > 1, ]

# More stringent
sig_strict <- results[results$FDR < 0.01 & abs(results$Fold) > 2, ]

# More accessible in treatment
opened <- sig[sig$Fold > 0, ]

# Less accessible in treatment
closed <- sig[sig$Fold < 0, ]
```

## Downstream Analysis

### GO Enrichment
```r
library(clusterProfiler)

# Get genes near differential peaks
genes <- unique(anno_df$SYMBOL[!is.na(anno_df$SYMBOL)])

# Enrichment
go <- enrichGO(genes, OrgDb=org.Hs.eg.db, keyType='SYMBOL', ont='BP')
dotplot(go)
```

### Motif Analysis
```bash
# Extract sequences from differential peaks
bedtools getfasta -fi genome.fa -bed opened_peaks.bed > opened.fa
bedtools getfasta -fi genome.fa -bed closed_peaks.bed > closed.fa

# Run HOMER for motif enrichment
findMotifsGenome.pl opened_peaks.bed hg38 motif_results/ -size 200
```

## Troubleshooting

### No Significant Peaks
- Check QC metrics first
- Verify biological replicates are consistent
- Consider relaxing thresholds for exploratory analysis

### Too Many Significant Peaks
- Check for batch effects
- Verify normalization is appropriate
- Consider more stringent thresholds

### Low Correlation Between Replicates
- May indicate technical issues
- Consider removing problematic samples
- Check for mislabeling
