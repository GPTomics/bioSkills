# CNV Visualization Usage Guide

## Overview

A copy number figure is an argument: the plot type, the y-axis quantity, and where the diploid baseline sits decide what a reader can conclude. This skill covers genome-wide and per-chromosome log2 scatter plots, B-allele-frequency tracks, ideograms, cohort heatmaps, circos views, and caller-native diagnostic plots, and, critically, when each view is the wrong one. The central rule: a depth-only log2 plot cannot show loss of heterozygosity or purity, and misleads outright if the baseline is centered on a non-diploid mode.

## Prerequisites

```bash
pip install matplotlib pandas numpy seaborn
conda install -c bioconda cnvkit          # for built-in scatter/diagram/heatmap
# R route: install.packages('ggplot2'); BiocManager::install('karyoploteR')
```

Inputs: caller output such as CNVkit `.cnr`/`.cns`, GATK `denoisedCR.tsv`/`modelFinal.seg`, or ASCAT/FACETS/Sequenza segment tables. For allele-specific tracks, a VCF of germline heterozygous SNP B-allele frequencies.

## Quick Start

Tell the AI agent what to do:
- "Create a genome-wide CNV profile with segments for this sample"
- "Make a two-panel log2 + BAF plot so I can see loss of heterozygosity"
- "Build a recurrent-CNV heatmap across my 40-sample cohort"
- "Diagnose whether my whole-genome-loss plot is a baseline centering artifact"
- "Plot chromosome 17 at high resolution to inspect an ERBB2 amplification boundary"

## Example Prompts

### Single-sample figures

> "Create a publication genome-wide log2 profile with colored segments for this CNVkit output, rasterizing the per-bin points so the PDF stays small."

> "Plot a combined log2 and B-allele-frequency figure for this tumor and point out any copy-neutral LOH regions."

### Cohort and context

> "Build a cohort CNV heatmap, and choose a bin size that keeps focal drivers like MYC and ERBB2 visible rather than averaging them away."

> "Make an ideogram showing which cytobands and genes this sample's CNVs overlap."

### Diagnosis

> "My genome-wide plot shows almost everything as deleted. Decide whether this is real or a diploid-baseline centering error, using the BAF track to check."

> "Inspect the ASCAT sunrise plot and FACETS diagnostic plot and tell me whether the purity/ploidy fit is trustworthy."

## What the Agent Will Do

1. Identify the biological question and pick the plot type that answers it
2. Confirm whether LOH or allele-specific state matters; if so, require a BAF track
3. Anchor the y-axis baseline to the caller's ploidy estimate, not the data mode
4. Render the figure (CNVkit built-ins, matplotlib, ggplot2, or caller-native plots)
5. Choose bin sizes and rasterization appropriate to focal vs arm-level questions
6. Flag when a log2-only view is insufficient or potentially misleading

## Tips

- Pair every log2 plot with a BAF track when LOH, allele-specific gain, or whole-genome doubling is in question; log2 alone cannot show them.
- Anchor the diploid baseline to the caller's ploidy; centering on the data mode inverts gain/loss calls in whole-genome-doubled tumors.
- Label log2 axes "log2 copy ratio", not "copy number"; log2 amplitude shrinks with decreasing purity and is not comparable across samples.
- For cohort heatmaps, match the bin size to the question: Mb bins for arm-level events, gene-level or GISTIC-peak bins for focal drivers.
- Rasterize per-bin scatter points (keep segments vector) so publication PDFs stay small.
- Always inspect caller-native diagnostic plots (ASCAT sunrise, Sequenza contour, FACETS plotSample) before trusting downstream calls.

## Related Skills

- copy-number/cnvkit-analysis - Generates .cnr/.cns inputs and built-in plots
- copy-number/gatk-cnv - GATK denoised-ratio and modeled-segment plots
- copy-number/allele-specific-copy-number - Source of BAF tracks and ploidy estimates
- copy-number/recurrent-cnv - Cohort recurrence underlying heatmaps
- data-visualization/ggplot2-fundamentals - Publication-figure grammar
- data-visualization/circos-plots - Circular genome layouts for CNV + SV
