# BigWig Tracks - Usage Guide

## Overview
bigWig is an indexed, compressed, multi-resolution binary container for continuous genomic signal -- coverage, fold-change, conservation (phyloP/phastCons), methylation rate, or signal p-value. It is fast because it answers a wide query from precomputed zoom-level summaries rather than from the base data, which makes the choice of summary statistic the defining expert nuance: a wide `mean` annihilates narrow tall peaks while staying faithful to broad domains, so the same region can read flat zoomed-out and obviously peaked zoomed-in. This skill covers querying a bigWig correctly with pyBigWig and the UCSC Kent tools, picking the right statistic and exactness for the biological question, handling no-data (NaN, not zero), and building a valid track from a sorted bedGraph plus a chrom.sizes file.

## Prerequisites
```bash
# pyBigWig (Python read/write; numpy support for array values)
pip install pyBigWig numpy

# UCSC Kent tools (CLI build/query)
conda install -c bioconda ucsc-bedgraphtobigwig ucsc-bigwigtobedgraph ucsc-bigwiginfo ucsc-bigwigsummary ucsc-bigwigaverageoverbed

# deepTools (compare, correlate, metaprofiles)
conda install -c bioconda deeptools
```

Building any bigWig also needs a chrom.sizes file (`cut -f1,2 reference.fa.fai > chrom.sizes`).

## Quick Start
Tell your AI agent what you want to do:
- "Get the mean signal from my bigWig for each peak in peaks.bed"
- "What is the peak height in chr1:1,000,000-2,000,000 -- the mean looks flat"
- "Build a browser-ready bigWig from my coverage bedGraph"
- "Compare my treatment and control bigWigs as a log2 ratio track"
- "Make a TSS metaprofile heatmap from my signal bigWig"

## Example Prompts

### Extracting Signal Correctly
> "Compute mean signal per gene from coverage.bw over genes.bed. This is a read-depth track, so treat uncovered bases as zero, and use exact values, not zoom approximations, since the numbers go in a table."
> "Get the peak height (not the mean) in chr3:5M-6M from my ChIP bigWig -- a single mean over that window dilutes the summit away."
> "Extract per-base values for chr1:1,000,000-1,001,000 from my methylation bigWig and average only the covered positions (no-data is undefined, not zero)."

### Building and Converting
> "Sort coverage.bedGraph and convert it to bigWig using hg38 chrom.sizes; verify with bigWigInfo."
> "Write a bigWig from these (chrom, start, end, value) intervals with pyBigWig, adding the header before the entries."
> "Convert my bigWig back to bedGraph for just chr1:1M-2M."

### Comparing and Profiling
> "Make a log2(IP/input) track from chip.bw and input.bw with a pseudocount."
> "Compute a correlation matrix across my four replicate bigWigs over a regions BED and plot a PCA."
> "Build a reference-point matrix of signal 2 kb around every TSS and plot a heatmap and profile."

## What the Agent Will Do
1. Name the biological question first -- peak height, total amount, level, or assayed fraction -- to fix the summary statistic.
2. Run `bigWigInfo` to check zoom levels, coverage, and chrom naming before trusting any query.
3. Choose `max`/`sum`/`mean`/`coverage` accordingly and pass `exact=True` (or use `values()`/`bigWigAverageOverBed`) when a number enters a result.
4. Decide gap handling biologically: coverage/depth tracks count gaps as zero (`mean0`); rate/ratio tracks treat gaps as undefined (`mean`/`np.nanmean`).
5. For builds, sort the bedGraph, supply a matching chrom.sizes, and add the header before entries.
6. For comparison/profiling, route to deepTools (`bigwigCompare`, `multiBigwigSummary`, `computeMatrix`).

## Tips
- A wide query is a summary, not the underlying data: `mean` (the default) hides narrow peaks; use `max` for peaks, `sum` for totals, `coverage` for assayed-fraction.
- `exact=False` is the pyBigWig default and reads zoom levels -- pass `exact=True` whenever a reviewer might recompute the number.
- No-data is NaN, not 0: `np.mean` poisons to NaN; pick `np.nanmean` (covered-only) or `np.nan_to_num().mean()` (gaps as zero) by the track's biology.
- `bigWigAverageOverBed` is the purpose-built per-region tool; its `mean` vs `mean0` columns are the same NaN-vs-zero decision.
- "I scanned the track and saw nothing" is not evidence of absence unless the scan was finer than the feature width -- zoom-out erases sharp signal.
- bedGraph must be sorted (`sort -k1,1 -k2,2n`) and non-overlapping before `bedGraphToBigWig`, and chrom.sizes must match the reference naming exactly.
- Use `bamCoverage` (deepTools) to generate a normalized track from a BAM -- that lives upstream in chip-seq/atac-seq.
- Signal -> bigWig; discrete features (peaks, genes) -> bigBed.

## Related Skills

- bedgraph-handling - The text bedGraph this skill converts to/from, and exact-arithmetic alternative
- coverage-analysis - Generates the per-base depth/bedGraph that becomes a bigWig
- bed-file-basics - The region BED files passed to bigWigAverageOverBed/computeMatrix
- chip-seq/chipseq-visualization - Generates normalized tracks (bamCoverage) and renders computeMatrix metaprofiles
- atac-seq/footprinting - Consumes bigWig signal over motif sites
- data-visualization/genome-tracks - Renders the bigWig in a browser figure
