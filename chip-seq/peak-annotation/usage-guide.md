# Peak Annotation - Usage Guide

## Overview

Annotate ChIP-seq peaks to genomic features (promoters, exons, introns, intergenic) and nearest genes. Three approaches available: ChIPseeker (R), HOMER annotatePeaks.pl (CLI), and Python (pandas/pyranges). All support custom GTF annotation files alongside pre-built annotation databases.

## Prerequisites

```r
# R (ChIPseeker with custom GTF support)
BiocManager::install(c('ChIPseeker', 'GenomicFeatures', 'rtracklayer'))

# For standard genomes only (not needed with custom GTF)
BiocManager::install(c('TxDb.Hsapiens.UCSC.hg38.knownGene', 'org.Hs.eg.db'))

# For functional enrichment
BiocManager::install('clusterProfiler')
```

```bash
# CLI (HOMER)
conda install -c bioconda homer
```

```bash
# Python
pip install pandas pyranges
```

## Quick Start

Tell your AI agent what you want to do:
- "Annotate my peaks to the nearest genes using the provided GTF"
- "Classify each peak as promoter, exon, intron, or intergenic"
- "Find the nearest gene and distance to TSS for each peak"
- "Show the distribution of peaks across genomic features"

## Example Prompts

### Annotation with Custom GTF
> "Annotate my ChIP-seq peaks using the GENCODE GTF file I provided. Define promoters as 2kb from TSS. Report nearest gene, distance to TSS, and feature category for each peak."

> "Use the GTF annotations to classify each peak in peaks.bed as promoter, exon, intron, or intergenic. Export as a TSV with columns: chr, start, end, nearest_gene, distance_to_tss, feature."

### Standard Genome Annotation
> "Annotate my narrowPeak file to hg38 genes using ChIPseeker"

> "Run HOMER annotatePeaks.pl on my peaks against hg38"

### Feature Analysis
> "What fraction of my H3K4me3 peaks are at promoters?"

> "Find which genes have peaks in their promoters and export as a gene list"

> "Compare peak annotation distributions between my two ChIP-seq samples"

### Visualization
> "Create a pie chart showing peak distribution across genomic features"

> "Plot the distance of peaks to the nearest TSS"

## What the Agent Will Do

1. Determine whether to use a pre-built annotation database or a provided GTF file
2. Parse the annotation and extract gene models, exon coordinates, and TSS positions
3. Handle coordinate system differences (GTF 1-based vs BED 0-based)
4. Choose a gene assignment convention (host-gene or nearest-TSS) based on the biological context
5. Classify peaks as promoter, exon, intron, or intergenic with proper priority
6. Collapse detailed subcategories (5' UTR, 3' UTR -> exon; downstream -> intergenic)
7. Compute signed distance to TSS (negative = upstream, positive = downstream)
8. Export annotated peaks with gene symbols and feature assignments

## Tips

- When a specific GTF is provided, always use it directly instead of a pre-built TxDb package to avoid gene name mismatches between annotation sources
- Promoter window size significantly affects results; match the definition to the analysis requirements (commonly 1-3kb)
- Ensure chromosome naming is consistent between peaks and annotations (chr1 vs 1); use `seqlevelsStyle()` in R to convert
- ChIPseeker's `annoDb` parameter does not work with custom TxDb from `makeTxDbFromGFF()`. Map gene symbols from the GTF manually using rtracklayer
- HOMER's built-in promoter definition (-1kb to +100bp) is not configurable; reclassify by Distance to TSS column for custom windows
- Gene assignment and feature classification can be decoupled: by default, HOMER and ChIPseeker assign the nearest-TSS gene independently of which gene's feature the peak overlaps. A peak inside gene A's intron near gene B's TSS will be labeled gene B / intron. The intron belongs to gene A, not B. Use ChIPseeker's `overlap='all'` or the host-gene Python approach to couple gene and feature. This is important when gene-feature consistency matters (e.g., functional enrichment of annotated genes)
- Feature priority matters when categories overlap: promoter > exon > intron > intergenic
- Peak center = (start + end) / 2 for BED coordinates; this is the reference point for all distance calculations
- For GENCODE GTFs, gene IDs include version suffixes (e.g., ENSG00000142192.25); strip versions with `sub('\\..*', '', id)` when matching

## Related Skills

- peak-calling - Generate peak files with MACS3 or HOMER
- motif-analysis - De novo and known motif enrichment
- differential-binding - Compare peaks between conditions
- chipseq-visualization - Signal tracks, heatmaps, profile plots
- genome-intervals/gtf-gff-handling - Parse and convert GTF/GFF files
- genome-intervals/proximity-operations - bedtools closest and window operations
- pathway-analysis/go-enrichment - Functional enrichment of peak-associated genes
