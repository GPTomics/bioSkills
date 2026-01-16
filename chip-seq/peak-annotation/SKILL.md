---
name: bio-chipseq-peak-annotation
description: Annotate ChIP-seq peaks to genomic features and genes using ChIPseeker. Assign peaks to promoters, exons, introns, and intergenic regions. Find nearest genes and calculate distance to TSS. Generate annotation plots and statistics.
tool_type: r
primary_tool: ChIPseeker
---

# Peak Annotation with ChIPseeker

## Load Peaks and Annotations

```r
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

# Read peaks from MACS2
peaks <- readPeakFile('sample_peaks.narrowPeak')
```

## Annotate Peaks

```r
# Annotate with default settings
peak_anno <- annotatePeak(
    peaks,
    TxDb = txdb,
    annoDb = 'org.Hs.eg.db'
)

# View annotation summary
peak_anno
```

## Custom Promoter Definition

```r
# Define promoter region (-3kb to +3kb from TSS)
peak_anno <- annotatePeak(
    peaks,
    TxDb = txdb,
    tssRegion = c(-3000, 3000),  # Promoter definition
    annoDb = 'org.Hs.eg.db'
)
```

## Extract Annotated Data Frame

```r
# Convert to data frame
anno_df <- as.data.frame(peak_anno)

# Key columns: seqnames, start, end, annotation, distanceToTSS, SYMBOL, GENENAME
head(anno_df)

# Export to CSV
write.csv(anno_df, 'annotated_peaks.csv', row.names = FALSE)
```

## Get Genes with Peaks in Promoter

```r
# Filter for promoter peaks
promoter_peaks <- anno_df[grep('Promoter', anno_df$annotation), ]

# Get unique genes
promoter_genes <- unique(promoter_peaks$SYMBOL)
```

## Annotation Pie Chart

```r
# Pie chart of genomic feature distribution
plotAnnoPie(peak_anno)

# Bar plot alternative
plotAnnoBar(peak_anno)
```

## Distance to TSS Plot

```r
# Distribution of peaks relative to TSS
plotDistToTSS(peak_anno, title = 'Distribution of peaks relative to TSS')
```

## Compare Multiple Peak Sets

```r
# Read multiple peak files
peak_files <- list(
    H3K4me3 = 'H3K4me3_peaks.narrowPeak',
    H3K27ac = 'H3K27ac_peaks.narrowPeak',
    H3K27me3 = 'H3K27me3_peaks.broadPeak'
)

peak_list <- lapply(peak_files, readPeakFile)

# Annotate all
anno_list <- lapply(peak_list, annotatePeak, TxDb = txdb, annoDb = 'org.Hs.eg.db')

# Compare annotations
plotAnnoBar(anno_list)
plotDistToTSS(anno_list)
```

## Venn Diagram of Peak Overlap

```r
# Find overlapping peaks
genes_list <- lapply(anno_list, function(x) as.data.frame(x)$SYMBOL)
vennplot(genes_list)
```

## Coverage Plot

```r
# Plot peak coverage around TSS
covplot(peaks, weightCol = 'V5')  # V5 is score column in narrowPeak
```

## Profile Heatmap Around TSS

```r
# Get promoter coordinates
promoter <- getPromoters(TxDb = txdb, upstream = 3000, downstream = 3000)

# Get tag matrix
tagMatrix <- getTagMatrix(peaks, windows = promoter)

# Plot heatmap
tagHeatmap(tagMatrix, xlim = c(-3000, 3000), color = 'red')

# Average profile
plotAvgProf(tagMatrix, xlim = c(-3000, 3000), xlab = 'Distance from TSS')
```

## Functional Enrichment of Peak Genes

```r
library(clusterProfiler)

# Get genes from peaks
genes <- unique(anno_df$ENTREZID)

# GO enrichment
ego <- enrichGO(
    gene = genes,
    OrgDb = org.Hs.eg.db,
    ont = 'BP',
    pAdjustMethod = 'BH',
    pvalueCutoff = 0.05
)
```

## Seq2Gene - All Genes in Peak Regions

```r
# Find all genes overlapping peak regions (not just nearest)
genes_in_peaks <- seq2gene(peaks, tssRegion = c(-1000, 1000), flankDistance = 3000, TxDb = txdb)
```

## Different Organisms

```r
# Mouse
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
peak_anno_mm <- annotatePeak(peaks, TxDb = TxDb.Mmusculus.UCSC.mm10.knownGene, annoDb = 'org.Mm.eg.db')

# Zebrafish
library(TxDb.Drerio.UCSC.danRer11.refGene)
library(org.Dr.eg.db)
```

## Key Functions

| Function | Purpose |
|----------|---------|
| readPeakFile | Read peak file (BED, narrowPeak) |
| annotatePeak | Annotate peaks to genes |
| plotAnnoPie | Pie chart of annotations |
| plotAnnoBar | Bar plot of annotations |
| plotDistToTSS | Distance to TSS distribution |
| getPromoters | Get promoter regions |
| getTagMatrix | Coverage matrix around regions |
| tagHeatmap | Heatmap of signal |
| plotAvgProf | Average profile plot |
| seq2gene | Map peaks to all overlapping genes |

## Annotation Categories

| Category | Description |
|----------|-------------|
| Promoter | Within tssRegion of TSS |
| 5' UTR | 5' untranslated region |
| 3' UTR | 3' untranslated region |
| Exon | Coding exon |
| Intron | Intronic region |
| Downstream | Within 3kb downstream |
| Distal Intergenic | Beyond gene regions |

## Related Skills

- peak-calling - Generate peak files with MACS2
- differential-binding - Find differential peaks
- pathway-analysis - Functional enrichment
- chipseq-visualization - Additional visualizations
