# Paired-End FASTQ - Usage Guide

## Overview

This skill enables AI agents to help work with paired-end sequencing data (R1/R2 files) using Biopython, keeping mates synchronized so downstream aligners pair the right reads.

## Prerequisites

```bash
pip install biopython
```

## Quick Start

Tell your AI agent what you want to do:

- "Filter paired reads keeping only pairs where both mates pass quality, and put orphans in separate files"
- "Interleave R1 and R2 into a single file for bwa mem -p"
- "Deinterleave this combined FASTQ into separate R1/R2 files"
- "Check whether my R1 and R2 files are still synchronized (matching counts and IDs)"
- "Find the R2 file that pairs with this R1"

## Example Prompts

### Synchronized Filtering
> "Filter my paired FASTQ files keeping pairs where both reads have mean quality >= 30, and route reads whose mate was discarded into orphan files"

### Interleaving
> "Combine reads_R1.fastq and reads_R2.fastq into an interleaved file, R1 then R2 alternating"

### Deinterleaving
> "Split this interleaved FASTQ back into separate R1 and R2 files and confirm the alternation was intact"

### Validation
> "Check if my R1 and R2 files have matching read counts and whether the read names still match as pairs"

### Mate Matching
> "My reads use the CASAVA 1.8 header format; match R1 and R2 by the ID before the space"

## Common Naming Patterns

| Pattern | R1 | R2 |
|---------|-----|-----|
| Illumina | `sample_R1_001.fastq` | `sample_R2_001.fastq` |
| Simple | `sample_1.fastq` | `sample_2.fastq` |
| Underscore | `sample_R1.fastq` | `sample_R2.fastq` |
| Dotted | `sample.R1.fastq.gz` | `sample.R2.fastq.gz` |

## Read-Name Conventions

| Era | Mate marker | Example |
|-----|-------------|---------|
| Pre-CASAVA 1.8 | `/1`, `/2` suffix on the name | `@HWUSI-EAS100R:6:73:941:1973#0/1` |
| CASAVA 1.8+ | second field after a space | `@EAS139:136:FC706VJ:2:2104:15343:197393 1:Y:18:ATCACG` |

In CASAVA 1.8+, the ID before the space is identical for both mates; tools match on that shared ID.

## What the Agent Will Do
1. Open both R1 and R2 FASTQ files (handling gzip if needed)
2. Iterate paired records in lockstep with `zip`
3. Verify pairing by matching the read ID up to the first space (and stripping `/1`/`/2`)
4. Apply any filter to BOTH mates, routing lone survivors to orphan files
5. Interleave, deinterleave, or report paired statistics as requested

## Tips

- Never filter, trim, or subsample one mate without the other; desync causes silent mismapping
- Keep pairs together: a pair survives only if both mates pass, otherwise route the survivor to an orphan file
- `zip` stops at the shorter file, so check that R1 and R2 counts match before processing files of unknown provenance
- CASAVA 1.8+ IDs match only up to the first space; do not compare full header descriptions
- Remove orphans before interleaving, or the strict R1,R2 alternation breaks and desyncs the file
- Use streaming generators for large files; reserve `SeqIO.index` for random access
- R2 mean quality is typically a bit lower than R1; a small gap is expected, not a defect

## Related Skills

- read-sequences - Parse individual FASTQ files and choose parse vs index
- fastq-quality - Phred encoding and quality interpretation before paired filtering
- filter-sequences - Single-file filtering criteria (apply to both mates here)
- compressed-files - gzip vs BGZF handling for paired files
- read-qc/quality-reports - FastQC/MultiQC per-mate quality assessment
- alignment-files/sam-bam-basics - Align paired reads with bwa mem after filtering
