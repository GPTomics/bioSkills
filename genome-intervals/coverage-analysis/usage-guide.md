# Coverage Analysis - Usage Guide

## Overview
Coverage analysis measures and interprets how deeply sequencing reads cover a genome or set of target regions. The central idea is that coverage is a distribution over positions, not a single number: the mean is inflated by repeat/rDNA pileups and blind to GC-extreme or poorly-mappable holes, so the decision-grade summary is the median plus a breadth curve (% of target reaching each depth) plus an evenness number (CV, Fano, Picard fold-80, or Gini). This skill covers mosdepth (the modern fast default for breadth curves and callable BEDs), bedtools genomecov and coverage (bedGraph tracks and per-target stats), samtools depth and coverage (per-base depth and per-contig depth+breadth), and the silent confounders -- duplicates, MAPQ filtering, GC bias, mate-overlap double-counting, and capture unevenness -- that bias depth before any plot is drawn.

## Prerequisites
```bash
# CLI tools
conda install -c bioconda bedtools mosdepth samtools

# Python (optional, for pybedtools and parsing dist files)
pip install pybedtools numpy
# pybedtools requires a bedtools binary on PATH
```

## Quick Start
Tell your AI agent what you want to do:
- "What is the median depth and breadth at 20x for my BAM?"
- "Build a callable-region BED from my alignments"
- "Generate a coverage bedGraph track for the genome browser"
- "Compute per-target coverage stats across my capture panel"
- "Check whether my high-depth amplicon panel has a mate-overlap problem"

## Example Prompts

### Adequacy and breadth
> "Run mosdepth on sample.bam, give me the median depth and the breadth at 1x, 10x, 20x, and 30x, and tell me whether the library is even or spiky."
> "My run reports mean 30x WGS -- is that actually enough? Show me the cumulative coverage curve and the fraction of the genome callable at 20x."
> "Compare per-contig mean depth and breadth across chromosomes and flag any contig with high mean depth but low breadth."

### Tracks and callable regions
> "Create a bedGraph coverage track from alignments.bam that includes zero-coverage regions, then prepare it for conversion to bigWig."
> "Make a callable-region BED that classifies each base as NO_COVERAGE, LOW, CALLABLE, or HIGH using mosdepth quantize."

### Per-target / capture QC
> "Compute coverage statistics for each exon in targets.bed from my exome BAM and list every target below 20x mean depth."
> "Run capture-uniformity QC and report the fold-80 base penalty and on-target percentage."

### Short-insert and spliced data
> "My cfDNA panel shows inflated variant allele fractions -- check for mate-overlap double-counting and recompute depth correcting for it."
> "Calculate exon coverage from my RNA-seq BAM, making sure introns are not counted as covered."

## What the Agent Will Do
1. Confirm the BAM/CRAM is coordinate-sorted and indexed, and that duplicates were marked (depth on an un-deduped BAM is a vanity number).
2. Decide what to count -- MAPQ threshold, duplicate handling, read span vs fragment, mate-overlap correction -- and state it explicitly.
3. Choose a tool: mosdepth for breadth curves and callable BEDs, samtools coverage for a per-contig glance, bedtools for tracks and per-target stats, Picard CollectHsMetrics for capture QC.
4. Report median (not mean), a breadth/cumulative curve at the caller's threshold, and an evenness number; flag mean/median skew, fold-80, or blacklist contamination.
5. Output in the requested form (summary table, bedGraph track, callable BED, per-target stats) and route normalized cross-sample tracks to the deepTools-based visualization skills.

## Tips
- Mark duplicates before measuring; otherwise the DUP filter has nothing to drop and depth is inflated.
- Never report mean depth alone -- report median and breadth at a depth threshold; mean/median > 1.1-1.2 signals a skewed distribution.
- mosdepth `--fast-mode`/`-x` silently disables mate-overlap correction; do not use it for VAF-sensitive short-insert data.
- Pre-1.13 samtools depth caps at 8000 and truncates silently; check `samtools --version` and add `-d 0` on old builds.
- In `samtools coverage`, the `coverage` column is breadth (% bases >=1x); `meandepth` is depth -- never report the bare word "coverage" without a qualifier.
- bedtools coverage reports stats for the `-a` file (A = targets, B = reads) since v2.24.0; older habits get this backwards silently.
- Bare `bedtools genomecov` outputs a histogram, not a track; add `-bg`/`-bga`. Always add `-split` for spliced RNA-seq.
- CRAM input to samtools depth requires `--reference`.
- Mask ENCODE-blacklist / low-mappability regions before summarizing, or the tail and mean are dominated by artifacts.
- A spiky library cannot be rescued by sequencing deeper -- fix the library (PCR-free, better capture, UMIs).

## Related Skills
- bedgraph-handling - bedGraph tracks this skill emits, and their normalization
- bigwig-tracks - Convert the coverage bedGraph to an indexed bigWig for browsers
- interval-arithmetic - Intersect coverage/callable BEDs with target regions
- alignment-files/pileup-generation - Per-base pileup upstream of depth and per-call DP
- alignment-files/bam-statistics - flagstat/idxstats and dup rate that explain coverage confounders
- chip-seq/chipseq-visualization - deepTools normalized coverage tracks for cross-sample comparison
- data-visualization/genome-tracks - Render the coverage tracks built here
