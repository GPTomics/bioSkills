---
name: bio-epitranscriptomics-modification-visualization
description: Create metagene plots and browser tracks for RNA modification data. Use when visualizing m6A distribution patterns around genomic features like stop codons.
tool_type: r
primary_tool: Guitar
---

## Version Compatibility

Reference examples tested with: deepTools 3.5+

Before using code patterns, verify installed versions match. If versions differ:
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Modification Visualization

**"Visualize m6A distribution around stop codons"** → Create metagene plots and genome browser tracks showing RNA modification patterns relative to transcript landmarks (5'UTR, CDS, 3'UTR, stop codon).
- R: `Guitar::GuitarPlot()` for metagene distribution plots
- CLI: `deeptools computeMatrix` → `plotHeatmap` for modification heatmaps

## Metagene Plots with Guitar

```r
library(Guitar)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

# Load m6A peaks
peaks <- import('m6a_peaks.bed')

# Create metagene plot
# Shows distribution relative to transcript features
GuitarPlot(
    peaks,
    txdb = TxDb.Hsapiens.UCSC.hg38.knownGene,
    saveToPDFprefix = 'm6a_metagene'
)
```

## Custom Metagene with deepTools

**Goal:** Create a metagene profile showing m6A enrichment distribution relative to gene body landmarks (TSS, TES).

**Approach:** Compute the log2 IP/input ratio as a bigWig track with bamCompare, then build a signal matrix over scaled gene regions with computeMatrix and render as a profile plot.

```bash
# Create bigWig from IP/Input ratio
bamCompare -b1 IP.bam -b2 Input.bam \
    --scaleFactors 1:1 \
    --ratio log2 \
    -o IP_over_Input.bw

# Metagene around stop codons
computeMatrix scale-regions \
    -S IP_over_Input.bw \
    -R genes.bed \
    --regionBodyLength 2000 \
    -a 500 -b 500 \
    -o matrix.gz

plotProfile -m matrix.gz -o metagene.pdf
```

## Browser Tracks

```bash
# Create normalized bigWig for genome browser
bamCoverage -b IP.bam \
    --normalizeUsing CPM \
    -o IP_normalized.bw

# Peak BED to bigBed
bedToBigBed m6a_peaks.bed chrom.sizes m6a_peaks.bb
```

## Heatmaps

```r
library(ComplexHeatmap)

# m6A signal around peaks
Heatmap(
    signal_matrix,
    name = 'm6A signal',
    cluster_rows = TRUE,
    show_row_names = FALSE
)
```

## Related Skills

- epitranscriptomics/m6a-peak-calling - Generate peaks for visualization
- data-visualization/genome-tracks - IGV, UCSC integration
- chip-seq/chipseq-visualization - Similar techniques
