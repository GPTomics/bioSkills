---
name: bio-data-visualization-genome-tracks
description: Build genome-browser-style multi-track figures with pyGenomeTracks (config-driven), Gviz (R), and IGV batch screenshotting. Covers BigWig coverage tracks, BED/peak overlays, gene-model rendering, Hi-C matrix tracks, BedPE link arcs, spike-in-aware normalization, and the bamCoverage --normalizeUsing trap. Use when producing publication figures of genomic loci with stacked aligned tracks (coverage, peaks, genes, interactions) for ChIP-seq, ATAC-seq, RNA-seq, Hi-C, or generic locus visualization.
tool_type: mixed
primary_tool: pyGenomeTracks
---

## Version Compatibility

Reference examples tested with: pyGenomeTracks 3.9+, Gviz 1.46+ (Bioconductor), deepTools 3.5+, GenomicRanges 1.54+, IGV 2.18+ (batch mode).

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)`
- R: `packageVersion('<pkg>')` then `?function_name`
- CLI: `<tool> --version` then `<tool> --help`

If code throws ImportError, AttributeError, or TypeError, introspect the installed package and adapt the example to match the actual API rather than retrying.

# Genome Browser Tracks

**"Plot a genomic locus with multiple tracks"** -> Build a stacked figure where each track (coverage from BigWig, peaks from BED, genes from GTF, Hi-C from cool, loops from BedPE) is aligned to genome coordinates. The decisions that matter: track normalization (especially for ChIP-Rx spike-in), gene-model rendering style (UCSC vs FlyBase), y-axis sharing across samples, and which tool fits the workflow — pyGenomeTracks (config-driven, reproducible, headless), Gviz (R Bioconductor), IGV batch (interactive-tool screenshots).

- Python / CLI: `pyGenomeTracks` (Lopez-Delisle 2021 *Bioinformatics* 37:422)
- R: `Gviz::plotTracks` (Hahne-Ivanek 2016)
- Interactive: IGV (Robinson 2011 *Nat Biotechnol* 29:24) with batch scripting

## The Single Most Important Modern Insight -- Spike-In Normalization Cannot Be Done With --normalizeUsing

deepTools `bamCoverage` is the canonical BigWig generator. Its `--normalizeUsing` flag accepts {RPKM, CPM, BPM, RPGC, None} — none of which implement ChIP-Rx spike-in normalization. All four divide by *sample-internal* mapped read counts and will UNDO any spike-in correction.

For ChIP-Rx (Orlando 2014 *Cell Rep* 9:1163):
1. Compute spike-in scale factor externally: `scale = 1 / (spike_reads_per_million)` OR per Orlando method
2. Pass via `--scaleFactor <value>` with `--normalizeUsing None`
3. **Do NOT combine `--scaleFactor` with `--normalizeUsing CPM/RPGC`** — re-normalizes the signal and undoes spike-in

This is the most common silent error in ChIP-seq visualization. The BigWig looks fine; the cross-sample comparison is wrong by the spike-in factor.

## pyGenomeTracks — Config-Driven, Reproducible

**Goal:** Render a multi-track locus figure from a config file specifying each track's source file, style, height, and color.

**Approach:** Write an `.ini` file with one section per track; invoke `pyGenomeTracks --tracks tracks.ini --region chr1:1000000-2000000 --outFileName out.pdf`.

```ini
# tracks.ini
[x-axis]
where = top
fontsize = 8

[h3k27ac]
file = h3k27ac.bw
title = H3K27ac
height = 3
color = #D55E00
min_value = 0
max_value = 50
number_of_bins = 700
summary_method = mean
nans_to_zeros = true

[spacer]
height = 0.3

[peaks]
file = h3k27ac_peaks.narrowPeak
title = Peaks
height = 0.8
color = #888888
display = collapsed
labels = false
file_type = narrowPeak

[loops]
file = loops.bedpe
title = Loops
height = 2
file_type = links
links_type = arcs
color = '#0072B2'
line_width = 0.5

[hic]
file = matrix.cool
title = Hi-C (KR-normalized)
height = 8
depth = 1000000
min_value = 0
max_value = auto
transform = log1p
colormap = RdYlBu_r

[genes]
file = gencode.v44.gtf
title = Genes
height = 5
fontsize = 8
style = UCSC                            # or 'flybase'; UCSC merges transcripts, flybase shows all
prefered_name = gene_name
merge_transcripts = true
color = '#3C5488'
border_color = black
```

```bash
pyGenomeTracks --tracks tracks.ini \
    --region chr1:1000000-2000000 \
    --outFileName locus.pdf \
    --width 18 \                          # CENTIMETERS not inches; default 40 cm
    --dpi 300

# For multiple regions from a BED:
pyGenomeTracks --tracks tracks.ini --BED regions.bed \
    --outFileName multi.pdf
```

**`--width` is in centimeters**, not inches. Default 40 cm; Nature double-column = 18.3 cm. `--decreasingXAxis` flips orientation for minus-strand loci.

## Gviz (R Bioconductor)

```r
library(Gviz)
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

# Tracks
axTrack <- GenomeAxisTrack()
itrack <- IdeogramTrack(genome = 'hg38', chromosome = 'chr1')

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
grTrack <- GeneRegionTrack(txdb, genome = 'hg38', chromosome = 'chr1',
                            name = 'Genes', transcriptAnnotation = 'symbol',
                            collapseTranscripts = 'meta')

dTrack <- DataTrack(range = 'h3k27ac.bw', type = 'h',
                     chromosome = 'chr1', name = 'H3K27ac',
                     col.histogram = '#D55E00', fill.histogram = '#D55E00')

aTrack <- AnnotationTrack(range = 'peaks.bed', name = 'Peaks',
                           chromosome = 'chr1', fill = '#888888',
                           stacking = 'dense')

# Render
plotTracks(list(itrack, axTrack, dTrack, aTrack, grTrack),
           from = 1000000, to = 2000000,
           sizes = c(1, 1, 3, 1, 4),
           background.title = 'transparent',
           cex.title = 0.7,
           cex.axis = 0.6)
```

## IGV Batch Scripting

For interactive-tool screenshots without launching the GUI:

```bash
# batch.txt
new
genome hg38
load sample.bam
load peaks.bed
snapshotDirectory ./screenshots
goto chr1:1000000-2000000
sort base
maxPanelHeight 500
snapshot region1.png
goto chr2:5000000-6000000
snapshot region2.png
exit
```

```bash
igv -b batch.txt
```

IGV batch is suitable when the workflow requires IGV's specific rendering style (allele frequencies, split-read pairs, soft-clipped sequences) — features pyGenomeTracks and Gviz don't replicate.

## BigWig Generation — The Spike-In Trap

```bash
# WITHOUT spike-in (e.g., RNA-seq, ATAC-seq):
bamCoverage -b sample.bam -o sample.bw \
    --binSize 10 \
    --normalizeUsing BPM \
    --effectiveGenomeSize 2913022398        # hg38 effective; check for build

# CORRECT ChIP-Rx spike-in:
# 1. Compute scale factor externally
SPIKE_RPM=$(samtools view -c sample.spike.bam) 
SCALE_FACTOR=$(echo "scale=10; 1000000 / $SPIKE_RPM" | bc)

# 2. Apply --scaleFactor with --normalizeUsing None
bamCoverage -b sample.bam -o sample.bw \
    --binSize 10 \
    --normalizeUsing None \                  # CRITICAL: None
    --scaleFactor $SCALE_FACTOR

# INCORRECT (silent error):
bamCoverage -b sample.bam -o sample.bw \
    --normalizeUsing CPM \                   # WRONG: undoes spike-in
    --scaleFactor $SCALE_FACTOR
```

## Track Comparison Across Samples

For multi-sample tracks (control vs treatment), set shared y-axis explicitly:

```ini
[sample1_bw]
file = sample1.bw
title = Control
height = 3
color = '#0072B2'
min_value = 0
max_value = 100                              # SHARED max across samples

[sample2_bw]
file = sample2.bw
title = Treatment
height = 3
color = '#D55E00'
min_value = 0
max_value = 100                              # SAME max for visual comparability
overlay_previous = share-y                   # for overlay; omit for stack
```

Without shared y-axis, the "taller" sample is the one with stronger absolute signal — but the figure visually conflates signal magnitude with rendering scale.

## Per-Method Failure Modes

### bamCoverage --normalizeUsing undoes spike-in

**Trigger:** ChIP-Rx workflow using `--normalizeUsing CPM` AND `--scaleFactor`.

**Mechanism:** CPM normalization divides by sample-internal reads; cancels the spike-in factor.

**Symptom:** Spike-in-normalized tracks look the same as un-normalized; cross-condition comparison wrong.

**Fix:** `--normalizeUsing None` with `--scaleFactor`. Validate by examining tracks at known reference loci where signal should match between samples.

### Different y-axis across samples

**Trigger:** Auto-scaled `max_value = auto` per-sample.

**Mechanism:** Each track scales independently to its own max.

**Symptom:** Visual "looks same" across samples that actually differ in magnitude.

**Fix:** Set explicit `min_value` and `max_value` to the same value across samples.

### Wrong gene-model style

**Trigger:** `style = flybase` for human data (or vice versa).

**Mechanism:** UCSC merges overlapping transcripts; flybase shows all isoforms; pile-up of isoforms unreadable for transcript-dense human loci.

**Symptom:** Gene track is a forest of overlapping arrows.

**Fix:** `style = UCSC` for human/mouse; `merge_transcripts = true` to collapse to canonical isoform.

### pyGenomeTracks --width interpreted as inches

**Trigger:** `--width 7` thinking inches.

**Mechanism:** Default unit is centimeters; `--width 7` is 7 cm = 2.75 inches.

**Symptom:** Tiny figure that doesn't match journal column width.

**Fix:** `--width 18.3` for Nature double column (18.3 cm = 183 mm). `--width 8.9` for single column.

### Track order top-down vs bottom-up confusion

**Trigger:** Expecting tracks in config-file order; pyGenomeTracks renders top-to-bottom (config[0] = top).

**Mechanism:** Convention differs across tools (Gviz top-to-bottom; some browsers bottom-to-top).

**Symptom:** Gene model at top instead of bottom.

**Fix:** Verify against config file order; for "genes at bottom" put `[genes]` section last.

### Hi-C matrix track depth too low

**Trigger:** `depth = 100000` for a 2 Mb region.

**Mechanism:** Hi-C matrix track shows interactions up to `depth` distance; smaller than region collapses the triangle.

**Symptom:** Hi-C track shows only a thin band.

**Fix:** `depth` should be ≥ half the region width; for 2 Mb region, `depth = 1000000` minimum.

### IGV batch script silent failures

**Trigger:** Typo in batch command; IGV continues to next command.

**Mechanism:** IGV batch mode doesn't fail-fast.

**Symptom:** Subset of snapshots missing; no error.

**Fix:** Verify each snapshot was produced; small batches and `set echo TRUE` for debugging.

## Reconciliation: When Tracks Disagree

| Pattern | Cause | Action |
|---------|-------|--------|
| Tracks look identical pre/post spike-in | --normalizeUsing canceled spike-in | Switch to None + --scaleFactor |
| Coverage differs between bamCoverage and IGV | Different binning; smoothing default | Specify --binSize explicitly; verify with raw BAM |
| Peaks in different positions across tools | Different peak-caller output (MACS narrowPeak vs broadPeak) | Document caller; cross-reference upstream chip-seq/peak-calling |
| Hi-C matrix orientation flipped | Pre-rotation vs post-rotation convention | Most tools assume upper-triangle; check vendor |

## Quantitative Thresholds

| Threshold | Value | Source |
|-----------|-------|--------|
| pyGenomeTracks --width default | 40 cm | Tool default; Nature ~18.3 cm |
| pyGenomeTracks --dpi recommended | 300 for publication | Standard |
| bamCoverage --binSize typical | 10-50 bp | Resolution vs file size trade-off |
| Hi-C track depth | >= half region width | Tool convention |
| Effective genome size hg38 | 2913022398 | UCSC |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| Spike-in normalized tracks look unnormalized | --normalizeUsing canceled spike-in | --normalizeUsing None + --scaleFactor |
| Y-axis differs across samples | Auto-scaling per-track | Explicit min/max in config |
| Gene track unreadable | flybase style on dense human locus | UCSC + merge_transcripts = true |
| Figure tiny | --width interpreted as inches | --width in CM |
| Hi-C band thin | depth too small | depth >= 0.5 × region width |
| IGV screenshots missing | Batch error silent | Verify per-snapshot; small batches |
| Coverage off by 2x | Strand-specific issue | Use --filterRNAstrand or split strands |

## References

- Hahne F, Ivanek R. 2016. Visualizing genomic data using Gviz and Bioconductor. *Methods Mol Biol* 1418:335-351.
- Lopez-Delisle L, Rabbani L, Wolff J, et al. 2021. pyGenomeTracks: reproducible plots for multivariate genomic datasets. *Bioinformatics* 37(3):422-423.
- Orlando DA, Chen MW, Brown VE, et al. 2014. Quantitative ChIP-seq normalization reveals global modulation of the epigenome. *Cell Rep* 9(3):1163-1170.
- Ramírez F, Ryan DP, Grüning B, et al. 2016. deepTools2: a next generation web server for deep-sequencing data analysis. *Nucleic Acids Res* 44(W1):W160-W165.
- Robinson JT, Thorvaldsdóttir H, Winckler W, et al. 2011. Integrative Genomics Viewer. *Nat Biotechnol* 29(1):24-26.

## Related Skills

- alignment-files/bam-statistics - BAM-level QC before bigwig
- chip-seq/peak-calling - Peak files for tracks
- chip-seq/chipseq-visualization - ChIP-seq-specific tracks
- hi-c-analysis/hic-visualization - Hi-C-specific contact maps
- alternative-splicing/sashimi-plots - Splice-junction tracks
- data-visualization/multipanel-figures - Combining track figures
- genome-intervals/bigwig-tracks - BigWig file handling
