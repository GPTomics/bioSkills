# RNA Modification Visualisation - Usage Guide

## Overview

Visualise RNA-modification data with the canonical 5'UTR / CDS / 3'UTR transcript-feature metagene plot (Guitar), the genome-coordinate metagene (deepTools), peak-centred heatmaps clustered by condition (ComplexHeatmap / deepTools), IP-vs-input paired browser tracks at specific loci (pyGenomeTracks / ggcoverage / Gviz), DRACH sequence-logo plots (ggseqlogo), and the 5'UTR / CDS / 3'UTR stacked-bar feature-distribution summary (ChIPseeker + ggplot2). Establishes stop-codon enrichment in the metagene as the biological QC anchor: a MeRIP dataset without it is suspect.

## Prerequisites

```r
BiocManager::install(c('Guitar', 'TxDb.Hsapiens.UCSC.hg38.knownGene', 'ComplexHeatmap',
                       'rtracklayer', 'GenomicFeatures', 'BSgenome.Hsapiens.UCSC.hg38',
                       'ChIPseeker', 'ggseqlogo', 'Gviz'))
install.packages(c('ggcoverage', 'ggplot2', 'dplyr', 'circlize'))
```

```bash
conda install -c bioconda deeptools pygenometracks ucsc-bedtobigbed homer
```

Reference inputs:

- Called peak BED from m6a-peak-calling
- IP-over-Input log2 bigWig tracks from merip-preprocessing
- Matching GENCODE / Ensembl GTF + TxDb
- BSgenome for sequence extraction (motif logos)

## Quick Start

- "Render the canonical Guitar metagene plot with stop-codon enrichment"
- "Build a peak-centred heatmap clustered by condition"
- "Make a pyGenomeTracks figure of paired IP / Input / log2 tracks at a specific locus"
- "Render a DRACH sequence logo from peak-centre 5-mers as a sanity check"
- "Build the 5'UTR / CDS / 3'UTR stacked-bar peak distribution"

## Example Prompts

### Metagene (THE Canonical Plot)

> "Render the Guitar GuitarPlot of my exomePeak2 m6A peaks using TxDb.Hsapiens.UCSC.hg38.knownGene; expect stop-codon-proximal enrichment as the biological QC anchor."

> "Confirm my MeRIP shows the canonical Dominissini 2012 / Meyer 2012 stop-codon enrichment pattern before downstream analysis."

> "Build a per-condition Guitar metagene with WT vs KO peak BEDs side-by-side."

### Browser Tracks

> "Build a pyGenomeTracks INI for paired IP and Input bigWig + log2 IP/Input + m6A peaks BED + GENCODE GTF; render at chr19:54,792,000-54,799,000."

> "Render an ggcoverage browser plot of IP-over-Input log2 with peaks highlighted at the METTL3 locus."

### Heatmaps

> "Build a peak-centred heatmap with deepTools computeMatrix reference-point + plotHeatmap; +/-500 bp window; k-means k=3."

> "Render a ComplexHeatmap of multi-condition IP/Input signal at peaks with row annotation showing per-condition mean."

### DRACH Logo

> "Render a ggseqlogo of peak-centre 5-mers from my exomePeak2 BED; confirm DRACH consensus visually."

> "Run HOMER findMotifsGenome.pl -rna on my peaks and report DRACH enrichment E-value."

### Stacked Bar

> "Render the Figure 1 5'UTR / CDS / 3'UTR / intron / intergenic stacked bar via ChIPseeker annotatePeak + ggplot2."

> "Compare peak feature distribution between WT and KO conditions as side-by-side stacked bars."

### Genome-Coordinate Metagene

> "Build a deepTools scale-regions metagene over protein-coding genes; +/-500 bp flanks; render as plotProfile + plotHeatmap."

## What the Agent Will Do

1. Build the Guitar metagene FIRST as a smoke test for stop-codon enrichment (the biological QC anchor)
2. If metagene does NOT show stop-codon enrichment, halt and recommend re-inspecting IP enrichment in merip-preprocessing
3. Render the 5'UTR / CDS / 3'UTR stacked bar paired with the metagene
4. Render the DRACH sequence logo as a sanity check on antibody specificity
5. Build pyGenomeTracks INI for paired IP / Input / log2 browser tracks at specific loci
6. Render peak-centred heatmap with deepTools computeMatrix + plotHeatmap (single-condition) or ComplexHeatmap (multi-condition)
7. For differential m6A analyses, render volcano + MA plots (defers to m6a-differential)
8. For per-condition or per-cluster comparison, build faceted plots
9. Export publication-ready PDF and PNG outputs

## Tips

- The Guitar metagene with stop-codon enrichment is the biological QC anchor. Build this FIRST after peak calling. If it does NOT show the canonical pattern, do not proceed.
- Pair the metagene with the 5'UTR / CDS / 3'UTR stacked bar; they are complementary biology-QC anchors, not redundant.
- DRACH sequence logo is a sanity check on antibody specificity, NOT a per-peak filter.
- Guitar's older API uses `txdb=`; newer uses `txTxdb=`. Verify with `?GuitarPlot`.
- pyGenomeTracks is the reproducible-figure choice for browser tracks; IGV is for interactive inspection only.
- ggcoverage is the ggplot2-native alternative to pyGenomeTracks; combinable with custom annotation layers.
- For peak-centred heatmaps with multi-condition clustering, use ComplexHeatmap; deepTools plotHeatmap is single-condition view.
- deepTools `--operation log2` is the modern syntax; `--ratio log2` is being phased out.
- The Volcano / MA plot for differential m6A lives in m6a-differential, not here; don't duplicate.
- bigWig signal magnitudes should be CPM-normalised in merip-preprocessing for cross-sample comparability.

## Related Skills

- merip-preprocessing - Generates the bigWig tracks used here
- m6a-peak-calling - Generates the peak BED for metagene / stacked bar / heatmap
- m6a-differential - Volcano / MA plots live there, not here
- m6anet-analysis - Per-site DRS modification calls; metagene visualisation is analogous
- data-visualization/ggplot2-fundamentals - General ggplot2 grammar
- data-visualization/multipanel-figures - Combining metagene + heatmap + volcano
- data-visualization/heatmaps-clustering - General heatmap clustering
- data-visualization/volcano-and-ma-plots - General volcano / MA recipes
- data-visualization/genome-tracks - General browser-track rendering
- data-visualization/sequence-logos - General sequence-logo plotting
- chip-seq/chipseq-visualization - Sibling browser-track + peak-centred heatmap patterns
