# Genome Browser Tracks - Usage Guide

## Overview

Genome-track plots show stacked aligned data (coverage from BigWig, peaks from BED, genes from GTF, Hi-C from cool, loops from BedPE) at a genomic locus. pyGenomeTracks is config-driven and reproducible; Gviz is the R Bioconductor option; IGV batch produces screenshots from the interactive tool. The dominant correctness trap is spike-in normalization - deepTools `bamCoverage --normalizeUsing` cannot implement ChIP-Rx; the scale-factor must be computed externally and applied with `--normalizeUsing None`.

## Prerequisites

```bash
pip install pyGenomeTracks deepTools
# IGV: https://software.broadinstitute.org/software/igv/
```

```r
BiocManager::install(c('Gviz', 'GenomicRanges', 'TxDb.Hsapiens.UCSC.hg38.knownGene'))
```

## Quick Start

Tell your AI agent what you want to do:
- "Write a pyGenomeTracks INI for ChIP-seq + peaks + genes at a locus"
- "Generate BigWigs from BAMs with ChIP-Rx spike-in normalization"
- "Render Gviz multi-track plot with ideogram + coverage + peaks + gene model"
- "IGV batch script to screenshot 20 loci"
- "Hi-C triangle track with loops overlay"

## Example Prompts

### ChIP-seq locus

> "pyGenomeTracks INI: x-axis at top, H3K27ac BigWig (orange, height 3), narrowPeak track, BedPE loops as arcs, GENCODE GTF as UCSC-style merged. Region chr1:1000000-2000000."

### Spike-in BigWig

> "Generate BigWigs for ChIP-Rx samples. Compute spike-in scale factor from samtools view of the spike-aligned BAM. Apply via --scaleFactor with --normalizeUsing None."

### Cross-sample comparison

> "Build pyGenomeTracks figure comparing control vs treatment ChIP-seq. Shared y-axis (min=0, max=50) across samples for visual comparability."

### Hi-C locus

> "Hi-C cool matrix track + ChIP-seq coverage + gene model. depth=1500000 for a 2 Mb region."

### IGV screenshots

> "IGV batch script: snapshot the top 20 differentially accessible peaks at 500 bp each side. Sort by base; max panel height 400."

## What the Agent Will Do

1. Confirm genome build and effective genome size from data.
2. For ChIP-Rx: compute spike-in scale factor externally; pass via `--scaleFactor` + `--normalizeUsing None`.
3. Build pyGenomeTracks INI in top-to-bottom order with explicit per-track styling.
4. Set shared y-axis (min, max) across samples when comparing.
5. For gene tracks: UCSC style + merge_transcripts=true for human/mouse density.
6. For Hi-C: depth ≥ half region width.
7. Invoke `pyGenomeTracks --tracks tracks.ini --region ... --outFileName ... --width 18.3 --dpi 300` (width in CM).
8. Verify output; spot-check spike-in normalization against known reference loci.

## Tips

- **`bamCoverage --normalizeUsing` does NOT implement ChIP-Rx spike-in.** All options (RPKM/CPM/BPM/RPGC) divide by sample-internal reads and undo spike-in. Use `--normalizeUsing None` + `--scaleFactor`.

- **Shared y-axis across samples** for cross-condition comparison. Set explicit `min_value` and `max_value` in pyGenomeTracks INI.

- **pyGenomeTracks `--width` is in CENTIMETERS**, not inches. Nature single col = 8.9 cm; double = 18.3 cm.

- **Gene-track style**: UCSC + `merge_transcripts = true` for dense human/mouse loci; flybase only for sparse organisms.

- **Hi-C `depth` must be ≥ half region width**; smaller collapses the triangle to a thin band.

- **`--decreasingXAxis` flips orientation** for minus-strand-of-interest loci.

- **Track order in pyGenomeTracks INI = render order top-to-bottom.** Genes typically last (bottom).

- **IGV batch is fail-tolerant.** Verify each snapshot was produced; missing snapshots are silent.

- **For published figures, choose pyGenomeTracks or Gviz** over IGV batch - config / code is reproducible; GUI screenshots are not.

- **deepTools `bamCoverage --binSize 10` for resolution** vs `--binSize 50` for file size. Defaults vary by tool version.

## Related Skills

- alignment-files/bam-statistics - BAM-level QC before bigwig
- chip-seq/peak-calling - Peak files for tracks
- chip-seq/chipseq-visualization - ChIP-specific visualization
- hi-c-analysis/hic-visualization - Hi-C contact maps
- alternative-splicing/sashimi-plots - Splice-junction tracks
- data-visualization/multipanel-figures - Combining tracks with other plots
- genome-intervals/bigwig-tracks - BigWig file handling
