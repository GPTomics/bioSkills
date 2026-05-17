# ChIP-seq Visualization - Usage Guide

## Overview

Generate ChIP-seq heatmaps, profile plots, genome-browser tracks, and sample-correlation visualizations. Covers deepTools (production), pyGenomeTracks (modern config-driven tracks), Gviz / EnrichedHeatmap / ChIPseeker (R), and IGV batch scripts. Embeds bigWig normalization decisions (CPM, BPM, RPGC, spike-in scaled), bamCompare operation choice (log2/subtract/SES), and k-means clustering of heatmaps for biological subgrouping. The most consequential choice is bigWig normalization; it determines whether visual comparison reflects biology.

## Prerequisites

```bash
conda install -c bioconda deeptools samtools bedtools pygenometracks
```

```r
BiocManager::install(c('Gviz', 'EnrichedHeatmap', 'ChIPseeker',
                       'TxDb.Hsapiens.UCSC.hg38.knownGene', 'rtracklayer'))
```

```bash
# IGV (for batch screenshots)
# Download from https://software.broadinstitute.org/software/igv/download
```

## Quick Start

Tell the agent what to do:
- "Generate spike-in-scaled bigWigs from a BAM list and Drosophila spike-in counts"
- "Create a TSS-centered heatmap with k-means clustering across 4 ChIP samples"
- "Compare ChIP signal between treatment and control with log2 bamCompare bigWig"
- "Build a pyGenomeTracks INI for chr1:1M-2M with H3K27ac, H3K4me3, MACS peaks, gene annotation"
- "Make an IGV batch script for 10 published loci to generate consistent screenshots"
- "Compute Spearman correlation between 6 replicates with multiBamSummary and plotCorrelation"
- "Why does my heatmap look identical across samples when I expect treatment-specific signal?"

## Example Prompts

### Standard heatmap + profile (TSS centered)
> "Convert chip.bam and input.bam to bigWigs with CPM normalization. Compute matrix centered on TSS (±3 kb). Generate heatmap with `--zMin -3 --zMax 3 --colorMap RdBu_r` and profile plot."

### k-means clustering for biology
> "I have H3K27ac signal at 5000 promoters. Use plotHeatmap with `--kmeans 3` to identify clusters of promoters by signal pattern. Output the cluster assignment BED."

### log2 ratio of ChIP vs Input
> "Generate a log2 ratio bigWig: bamCompare ChIP vs input with `--operation log2 --pseudocount 1 --skipZeroOverZero`. Use bin size 50."

### Spike-in scaled tracks
> "I have 120k Drosophila reads in ctrl_1 and 85k in treat_1. Compute scale factors (smallest / each) and generate spike-in-scaled bigWigs. Do NOT also pass `--normalizeUsing`."

### Cross-condition comparison
> "Build heatmaps showing H3K27ac signal at union SE regions for DMSO vs BETi. Use RPGC normalization for cross-sample comparability OR spike-in scaled if I have spike-in data."

### Genome browser views
> "Generate pyGenomeTracks figure for chr1:1500000-1700000 showing H3K27ac, H3K4me3, MACS narrowPeak, SE BED, and gene annotation."

### IGV batch
> "Create an IGV batch script to screenshot 5 specific loci with my BAM, BW, and peak files. End with `exit`."

### Diagnostic
> "My heatmap looks pancake-flat. Diagnose: is it normalization (`--zMin --zMax` not set), data quality, or wrong matrix mode?"

## What the Agent Will Do

1. **Choose normalization** based on biology:
   - Within-sample: CPM (standard) or BPM (length-aware)
   - Cross-sample comparable: RPGC with read-length-matched effective genome size
   - Cross-condition with global shifts: spike-in scaled (mutually exclusive with `--normalizeUsing`)
   - ChIP vs input: bamCompare log2 (most common) / subtract / SES
2. **Generate bigWigs**: `bamCoverage` or `bamCompare`; extend single-end reads to fragment length; bin size 10-50 bp
3. **Compute signal matrix**: reference-point (TSS / peak summit) with `-b 3000 -a 3000`, or scale-regions for gene bodies
4. **Generate plots:**
   - Heatmaps with `--zMin/--zMax` set; `--kmeans N` for biological subgrouping
   - Profile plots with `--perGroup` for sample comparison
   - Browser tracks via pyGenomeTracks INI (recommended) or Gviz (R)
5. **Sample correlation**: `multiBamSummary bins` → `plotCorrelation --corMethod spearman`
6. **Validate**: check that scale, normalization, and figure annotations match the biological claim
7. **Export**: PDF (publication) or PNG (presentation); pyGenomeTracks supports both
8. **Document**: normalization method + parameters, bin size, color map, z-scale, region count, kmeans settings if used

## Tips

- **Normalization choice is the single biggest decision.** Get it right before generating any heatmap.
- **`--normalizeUsing` and `--scaleFactor` are mutually exclusive in bamCoverage.** Use one OR the other.
- **For log2 bamCompare, always add `--pseudocount 1 --skipZeroOverZero`** to avoid -Inf in plots.
- **Set `--zMin/--zMax` explicitly.** Default ranges are outlier-driven and produce uninformative heatmaps.
- **k-means clusters by signal in the FIRST `-S` sample only.** Order samples so the most-discriminating one is first.
- **Pyramids of error: read length must match the bigWig effective genome size.** Use deepTools `effectiveGenomeSize` table for accurate RPGC.
- **For broad histone marks, use larger bin sizes (50-100 bp) and broader matrix windows (5-10 kb).**
- **pyGenomeTracks is easier to script than Gviz** for pipeline figure generation across many regions; Gviz is more flexible for one-off publication figures.
- **For spike-in scaled tracks, internal-control sanity check**: blacklist regions should show no signal change post-scaling. See chip-seq/spike-in-normalization.
- **IGV batch scripts MUST end with `exit`** or IGV hangs.

## Troubleshooting

### Heatmap looks pancake-flat

1. `--zMin/--zMax` not set; outliers dominate the color scale → set explicit range
2. Bin size too large for signal pattern → use 10 bp for sharp marks
3. Window too narrow for biology → broaden `-b -a` parameters
4. Signal too low → check FRiP, depth (see chipseq-qc)

### Spike-in scaled tracks identical to CPM tracks

`--normalizeUsing` and `--scaleFactor` both passed; deepTools applies scale first then normalizes (undoes the scale). Use one only.

### Profile plot shows identical curves for treatment and control

1. Wrong bigWig used (sample swap)
2. Composition-bias-driven normalization erased the difference
3. Replicate signal averaged across one good + one bad replicate

### Tracks display "no data" in pyGenomeTracks

1. Region outside chromosome bounds → check `samtools view -H bam | head`
2. INI file path errors → use absolute paths in `file =`
3. Chromosome naming mismatch (chr vs no chr)

### IGV batch script hangs

Missing `exit` at end; IGV waits indefinitely for next command.

### bamCompare log2 ratios show -Inf or NaN

Missing `--pseudocount 1` or `--skipZeroOverZero`. Either fix produces plottable values.

### Heatmap k-means clusters look random

k-means clusters by FIRST sample only. If the first sample has weak signal, clustering becomes noise-driven. Reorder samples so the most-discriminating one is first; or compute matrix per-sample and combine externally for joint clustering.

## Related Skills

- chip-seq/peak-calling - Peak files for reference regions
- chip-seq/chipseq-qc - Replicate correlation, fingerprint plots
- chip-seq/spike-in-normalization - Spike-in scaled bigWig generation
- chip-seq/differential-binding - Visualize differential peaks
- chip-seq/super-enhancers - Heatmaps of SE constituents
- chip-seq/peak-annotation - Annotation-aware visualization (promoter vs enhancer)
- data-visualization/genome-tracks - General genome track patterns
- data-visualization/heatmaps-clustering - Heatmap clustering conventions
