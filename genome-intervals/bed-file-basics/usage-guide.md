# BED File Basics - Usage Guide

## Overview
BED (Browser Extensible Data) is the standard format for genomic intervals and the coordinate substrate the whole interval category rests on. This skill covers reading, creating, validating, sorting, and converting BED files (BED3 through BED12, narrowPeak/broadPeak) with bedtools (CLI) and pybedtools/pyranges/pandas (Python). Its central concern is the coordinate convention: BED is 0-based half-open while GTF/GFF, SAM, VCF, and wiggle are 1-based closed, and a botched conversion shifts every answer by one base with no error. It also covers the silent traps that make an analysis "run fine but wrong" -- chromosome-name mismatch, CRLF corruption, and the lexicographic-vs-version sort mismatch under bedtools `-sorted`.

## Prerequisites
```bash
conda install -c bioconda bedtools samtools ucsc-bedtobigbed ucsc-liftover crossmap
pip install pybedtools pyranges pandas
```
pybedtools requires a `bedtools` binary on PATH. pyranges is pure-Python (no binary). Generate a genome/chrom.sizes file from the same reference FASTA the rest of the pipeline used: `samtools faidx ref.fa && cut -f1,2 ref.fa.fai > genome.txt`.

## Quick Start
Tell your AI agent what you want to do:
- "Create a BED file from my list of peak coordinates"
- "Convert my VCF variants to BED, handling the coordinate systems"
- "Why does my intersect return zero overlaps?"
- "Sort my BED file the way bedtools expects"
- "Lift my hg19 BED coordinates over to hg38 and tell me what failed to map"
- "Validate this BED file and check the BED12 blocks"

## Example Prompts

### Creating and reading
> "Create a BED6 file from this DataFrame with columns chrom, start, end, name, score, strand, then read it back and report the interval lengths."
> "Build a minimal BED3 from these peak coordinates and save it sorted."

### Coordinate conversion
> "Convert my VCF variant positions to a BED file, subtracting 1 from the start so the coordinates line up, and verify on a known single-base variant."
> "I have a 1-based GTF feature -- give me the equivalent BED interval and confirm the length is unchanged."

### Debugging silent failures
> "My peaks.bed intersects zero genes -- check whether it is a chromosome-naming mismatch (chr1 vs 1) before assuming there is no biological overlap."
> "This BED works in one tool and errors in another on the last column -- check for CRLF line endings."
> "I ran bedtools intersect -sorted and the output only has low-numbered chromosomes -- diagnose the sort-order mismatch."

### Cross-assembly liftover
> "Lift this hg19 peak BED over to hg38 with liftOver and report how many regions failed to map and why I should not ignore them."
> "Use CrossMap to convert my GRCh37 intervals to GRCh38 and record the chain-file provenance with the output."

### Validation and windows
> "Validate this BED file: field-count consistency, negative/inverted intervals, chromosome naming, and whether it is sorted."
> "Reconstruct the absolute exon coordinates from this BED12 and confirm the block invariants hold."
> "Make 10 kb sliding windows with a 5 kb step across hg38 and number them."

## What the Agent Will Do
1. Establish the coordinate convention from the file format (BED 0-based half-open vs GTF/VCF/SAM 1-based closed).
2. Read or create the intervals with bedtools, pybedtools, pyranges, or pandas as appropriate.
3. Check the silent traps: chromosome-name harmonization across files, CRLF, and sort order before any `-sorted` operation.
4. Apply the requested operation (sort, validate, filter, convert), using `start - 1, end unchanged` at any 1-based boundary and a genome.txt derived from the aligned FASTA where lengths are needed.
5. Verify conversions on a known single-base landmark and save the output with consistent field counts.

## Tips
- The first question on any interval file is "what convention is it in?" -- answered from the format, not by looking at the file.
- Harmonize chromosome naming (chr1 vs 1, chrM vs MT) before any cross-file operation; a mismatch gives a valid empty result, not an error.
- Sort with `bedtools sort` or `sort -k1,1 -k2,2n` (lexicographic) before tabix/bigBed; for `-sorted` ops, sort both files identically or pass `-g genome.txt`.
- Never open a BED file in Excel: it date-mangles gene names (SEPT9 -> 9-Sep) and float-truncates large coordinates.
- Generate genome.txt from the exact reference FASTA the BAM was aligned to; a generic chrom.sizes silently rots slop/complement/windows.
- Use `pybedtools.cleanup()` at the end of a pybedtools script to remove temp files.

## Related Skills

- interval-arithmetic - Set operations on the intervals defined here
- gtf-gff-handling - 1-based annotation parsing and gene-model hierarchy
- coverage-analysis - Per-base depth feeding bedGraph intervals
- alignment-files/sam-bam-basics - BAM-to-BED conversion and SAM-vs-BAM coordinates
- variant-calling/vcf-basics - VCF POS-to-BED conversion and indel anchoring
