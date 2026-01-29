# Genome Browser Tracks - Usage Guide

## Overview

Create publication-quality genome browser visualizations showing coverage, peaks, genes, and Hi-C data using pyGenomeTracks (Python) or IGV batch scripting.

## Prerequisites

```bash
# pyGenomeTracks
pip install pyGenomeTracks

# IGV (download from Broad Institute)
# https://software.broadinstitute.org/software/igv/download
```

```r
# Gviz (R alternative)
BiocManager::install('Gviz')
```

## Quick Start

Tell your AI agent:
- "Create a genome browser figure for chr1:1000000-2000000"
- "Show my ChIP-seq coverage and peaks for this region"
- "Make track plots comparing treatment vs control"
- "Generate IGV snapshots for my list of regions"

## Example Prompts

### Basic Track Plots

> "Plot my bigWig coverage for the MYC locus"

> "Show peaks and gene annotations for chr17:7500000-7700000"

### Multi-Sample

> "Compare ChIP-seq signals between treatment and control"

> "Overlay two samples with shared y-axis"

### Batch Generation

> "Generate snapshots for all regions in my BED file"

> "Create track plots for my top 20 DE genes"

## What the Agent Will Do

1. Identify input files (bigWig, BED, BAM, GTF)
2. Create INI configuration file for pyGenomeTracks
3. Set appropriate track heights and colors
4. Configure axis and spacers
5. Run pyGenomeTracks for specified region
6. Export at publication resolution

## Tool Selection

| Tool | Best For |
|------|----------|
| pyGenomeTracks | Reproducible figures, scripting, Hi-C |
| IGV batch | Quick snapshots, interactive exploration |
| Gviz | R workflows, Bioconductor integration |

## Tips

- **pyGenomeTracks** produces the most publication-ready output
- **overlay_previous = share-y** for comparing samples on same scale
- **min_value/max_value** should match across samples for fair comparison
- **width parameter** controls plot width in cm (default 40)
- **dpi 300** minimum for publication
- **Use BED12** for transcript structures with exon blocks
