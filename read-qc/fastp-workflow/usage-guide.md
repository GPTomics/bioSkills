# fastp Workflow - Usage Guide

## Overview
fastp is a modern, all-in-one FASTQ preprocessor that does adapter trimming, quality/length filtering, 2-color poly-G removal, base correction, and QC reporting in a single fast C++ pass, replacing separate Cutadapt, Trimmomatic, and FastQC steps for bulk Illumina data. Its defining feature is paired-end adapter detection by overlap analysis (no adapter sequence needed). It is the default for bulk preprocessing, but not a substitute for cutadapt's precision on small-RNA/amplicon adapters or for proper UMI-based molecule counting.

## Prerequisites
```bash
conda install -c bioconda fastp
```

## Quick Start
Tell your AI agent what you want to do:
- "Process my paired-end FASTQ files with fastp"
- "Clean up my NovaSeq data with poly-G trimming"
- "Run fastp with Q20 filtering and a minimum length of 36 bp"
- "Aggregate all my fastp reports with MultiQC"

## Example Prompts

### Basic processing
> "Run fastp on my paired-end FASTQ files and generate HTML and JSON reports"

> "Process all samples in my data directory with fastp"

### Custom parameters
> "Use fastp with Q20 filtering, a minimum length of 36 bp, and 3' sliding-window trimming"

> "Enable overlap base correction for my paired-end reads"

### Platform-specific
> "Process my NovaSeq data with poly-G trimming"

> "Merge overlapping paired-end reads from my short-insert library"

### Batch processing
> "Run fastp on all my samples and aggregate the JSON reports with MultiQC"

## What the Agent Will Do
1. Run fastp with overlap-based PE adapter trimming (no adapter sequence needed)
2. Apply per-read quality filtering and optional light window trimming with a length gate
3. Handle platform-specific issues (auto poly-G for 2-color instruments)
4. Generate HTML and JSON reports for MultiQC aggregation
5. Optionally merge overlapping reads or extract UMIs (dedup happens later, post-alignment)

## Comparison with the traditional pipeline

| Task | Traditional | fastp |
|------|-------------|-------|
| QC report | FastQC | Built-in HTML/JSON |
| Adapter trim | Cutadapt | Overlap analysis (no sequence needed) |
| Quality trim/filter | Trimmomatic | Built-in |
| Poly-G | Manual | Auto for 2-color |
| Speed | 3 passes | 1 pass |

## Tips
- For paired-end data fastp needs no adapter sequence; it finds the read overlap and trims read-through automatically. Provide `--adapter_sequence` for single-end.
- `--correction` uses the overlap to fix a low-quality base from its high-quality mate (paired-end only).
- Poly-G trimming auto-enables for NextSeq/NovaSeq from the instrument ID; leave it on.
- Do NOT use `--dedup` for RNA-seq, amplicon, or any assay where identical reads are real signal; it is sequence-identity dedup and removes biological duplicates. Use UMIs or coordinate dedup instead.
- `--umi` only extracts the UMI; molecule-accurate dedup/consensus happens after alignment (umi_tools/fgbio).
- For light RNA-seq trimming, skip aggressive quality cutting; the aligner soft-clips tails.
- JSON reports feed MultiQC directly for cohort-level review.

## Resources
- [fastp GitHub](https://github.com/OpenGene/fastp)
- [fastp Publication](https://doi.org/10.1093/bioinformatics/bty560)

## Related Skills
read-qc/adapter-trimming - Precise adapter/primer control for small-RNA and amplicon
read-qc/quality-filtering - Detailed quality/length filtering and the trim-light evidence base
read-qc/quality-reports - Aggregate fastp JSON across samples with MultiQC
read-qc/umi-processing - Molecule-accurate UMI dedup and consensus after alignment
alignment-files/duplicate-handling - Coordinate-based duplicate marking for DNA variant calling
