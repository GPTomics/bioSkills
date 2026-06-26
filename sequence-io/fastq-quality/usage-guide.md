# FASTQ Quality - Usage Guide

## Overview

This skill enables AI agents to help you work with FASTQ quality scores for NGS data: accessing Phred scores, filtering and trimming by quality, profiling quality by position, and converting between the Sanger/Phred+33, Solexa, and Illumina/Phred+64 encodings using Biopython.

## Prerequisites

```bash
pip install biopython
```

## Quick Start

Tell your AI agent what you want to do:

- "Calculate average quality for each read"
- "Filter reads with mean quality below 25"
- "Trim low-quality bases from the 3' end with a sliding window"
- "Show the per-position quality profile for the first 50 bases"
- "Convert this old Illumina 1.3 FASTQ to standard Phred+33"
- "Which quality encoding does this FASTQ use?"

## Example Prompts

### Quality Analysis
> "Show me the quality distribution of reads.fastq and the fraction of reads above Q30"

### Filtering
> "Keep only reads with average quality >= 30 and write them to filtered.fastq"

### Trimming
> "Trim the 3' end where a 5-base sliding window drops below Q20"

### Per-Position Profiling
> "Generate a per-position mean quality profile and tell me where quality starts dropping"

### Encoding and Conversion
> "This is legacy Solexa data - convert it to standard Phred+33 without corrupting the scores"
> "Confirm whether reads.fastq is Phred+33 or Phred+64 before I filter it"

## What the Agent Will Do
1. Determine the FASTQ quality encoding from instrument/run metadata before parsing - never guessing the offset
2. Parse records and read scores from `record.letter_annotations['phred_quality']`
3. Compute per-read or per-position quality metrics
4. Filter, trim, or convert encodings as requested
5. Report summary statistics and flag binning or B-tail artifacts

## Quality Score Reference

- Q20 = 99% accuracy (1 error per 100 bases)
- Q30 = 99.9% accuracy (1 error per 1000 bases)
- Q40 = 99.99% accuracy (1 error per 10000 bases)

## The Four Encodings

| Variant | Format string | Offset | Score type |
|---------|---------------|--------|------------|
| Sanger / Illumina 1.8+ | `'fastq'` / `'fastq-sanger'` | 33 | Phred |
| Solexa / Illumina 1.0 | `'fastq-solexa'` | 64 | Solexa odds (can go negative) |
| Illumina 1.3-1.7 | `'fastq-illumina'` | 64 | Phred |

## Tips

- Never guess the offset. A wrong format string either raises loudly or, worse, silently shifts every score by exactly 31 in the ASCII overlap region. Confirm the encoding from the instrument, not from the data.
- Auto-detection is provably ambiguous: a high-quality Sanger file and a low-quality Illumina-1.3 file can be byte-identical. Scanning the minimum byte can rule out encodings but cannot confirm one.
- Solexa uses an odds score that goes negative; it is stored in `record.letter_annotations['solexa_quality']`, and converting it to Phred and back is lossy at low quality.
- The attribute is `letter_annotations`, not `per_letter_annotations` (the latter does not exist).
- In Illumina 1.5-1.7 data, a trailing run of `B` (Q2) is a quality-control flag, not real Q2 - strip it rather than trimming it as a score.
- NovaSeq bins quality to four levels (Q2, Q12, Q23, Q37), so quality histograms show spikes; this is expected and interacts with GATK BQSR (coarser corrections).
- Process large files as iterators (`SeqIO.parse`), not lists, to keep memory flat.

## Related Skills

- read-sequences - Parse FASTQ records and choose parse vs index for large files
- filter-sequences - Filter reads by length and content alongside quality
- paired-end-fastq - Keep R1/R2 synchronized when filtering paired reads
- sequence-statistics - Summary statistics across read sets
- read-qc/quality-reports - FastQC-style aggregate quality reports
