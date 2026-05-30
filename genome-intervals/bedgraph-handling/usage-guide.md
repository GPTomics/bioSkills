# bedGraph Handling - Usage Guide

## Overview
bedGraph is a 4-column (chrom/start/end/value) 0-based half-open text format for continuous genomic signal - coverage, scores, log2 ratios. This skill covers generating signal tracks from a BAM, the normalization step that makes tracks comparable across samples (and the global-perturbation case where library-size normalization is the wrong tool), multi-sample track arithmetic, and converting the text bedGraph into a browser-ready bigWig under the strict sort/overlap/chrom.sizes contract. The central judgment is that a raw bedGraph is a library-size artifact, the wrong normalization is worse than none, and bigWig - not bedGraph - is the deliverable.

## Prerequisites
```bash
# deepTools (bamCoverage/bamCompare/bigwigCompare/multiBigwigSummary)
conda install -c bioconda deeptools

# bedtools (genomecov/unionbedg/merge/map)
conda install -c bioconda bedtools

# UCSC Kent utilities (bedGraphToBigWig/bigWigToBedGraph/fetchChromSizes)
conda install -c bioconda ucsc-bedgraphtobigwig ucsc-bigwigtobedgraph

# pyBigWig (read/extract bigWig in Python)
pip install pyBigWig
```

## Quick Start
Tell your AI agent what you want to do:
- "Generate a normalized coverage bigWig from my BAM file"
- "Make a log2 treatment-vs-input track from two BAMs"
- "Convert my bedGraph to bigWig for the genome browser"
- "Check whether library-size normalization is safe for my knockdown experiment"
- "Build a sample correlation heatmap from my bigWig tracks"

## Example Prompts

### Generating a normalized track
> "Make an RPGC-normalized bigWig from chip.bam for GRCh38, 25 bp bins, extending reads to fragments and excluding chrX and chrM from the scale factor."
> "Generate a CPM-normalized bedGraph from my ATAC BAM so I can inspect the values, then convert it to bigWig."

### Deciding the normalization
> "I knocked down a histone acetyltransferase and want to compare H3K27ac between control and KD. Is CPM/RPGC normalization safe here, or do I need a spike-in?"
> "Which normalization should I use to compare ChIP signal across samples when a few regions have very high coverage?"

### Multi-sample arithmetic
> "Make a log2 ratio track of treatment over input from the two BAMs, normalizing for depth first."
> "Stack my three replicate bedGraphs into a value matrix over a common interval partition."
> "Build a Spearman correlation heatmap across my six bigWig tracks to check sample relatedness."

### Conversion and troubleshooting
> "My bedGraphToBigWig run says the input is not case-sensitive sorted - fix the sort."
> "My bigWig loads in the browser but the heights look wrong - help me find what corrupted it."
> "I have overlapping intervals in my bedGraph; collapse them so conversion succeeds."

## What the Agent Will Do
1. Establish what is being measured and whether a global-level change is plausible - if so, flag that library-size normalization will erase it and route to a spike-in.
2. Generate the track, preferring deepTools bamCoverage (BAM -> normalized bigWig in one step) over a hand-built bedGraph.
3. Pick the normalization from the taxonomy (None for inspection only; BPM for composition-aware cross-sample; RPGC for interpretable browser viewing), supplying the correct effective-genome-size for the build, read length, and filtering.
4. For ChIP/ATAC, extend reads to fragments; for RNA-seq, use split/stranded handling.
5. For comparisons, use bamCompare (normalizes depth then computes) on raw BAMs, or bigwigCompare only on already-normalized tracks.
6. When converting a text bedGraph, C-locale-sort, ensure non-overlapping intervals, derive chrom.sizes from the exact aligned-to FASTA, inspect the text, then run bedGraphToBigWig.

## Tips
- Ship bigWig, never bedGraph - bedGraph is unindexed scratch; bigWig is the random-access artifact.
- Inspect the text bedGraph (`head`, `awk '$4 > 1000'`, `sort -k4 -n`) before converting - it is the last human-readable checkpoint before opaque binary.
- Always `LC_COLLATE=C sort -k1,1 -k2,2n` before bedGraphToBigWig; a locale-aware sort that looks fine interactively fails in the scheduler.
- Derive chrom.sizes from the exact FASTA the reads aligned to (`samtools faidx` + `cut -f1,2`), not a generic download - a mismatched assembly/patch errors or pads phantom genome.
- Choose binSize to match feature width (sharp TF/ATAC 10-25 bp, broad marks 50-200 bp); tracks compared bin-for-bin must share the same binSize.
- Get the effective-genome-size table right: non-N length if multimappers were kept, the read-length unique-k-mer value if you filtered to unique alignments (~7% difference).
- A log2 browser track's apparent fold-changes are partly a pseudocount and bin-size artifact - do not read them as measured.

## Related Skills

- coverage-analysis - Per-base depth generation and distribution-vs-mean diagnostics feeding bedGraph tracks
- bigwig-tracks - Reading, extracting, and writing the bigWig deliverable this skill produces
- chip-seq/spike-in-normalization - The bench-decided external-reference scaling when a global change makes library-size normalization wrong
- chip-seq/chipseq-visualization - Render the normalized signal tracks built here
- atac-seq/footprinting - Consumes high-resolution coverage/bigWig signal over motif sites
- data-visualization/genome-tracks - Render the bedGraph/bigWig tracks for figures
