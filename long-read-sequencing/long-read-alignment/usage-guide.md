# Long-Read Alignment - Usage Guide

## Overview
minimap2 is the standard aligner for long reads. It's fast, accurate, and handles the higher error rates of long-read data. Different presets optimize for ONT, PacBio, or RNA sequencing.

## Prerequisites
```bash
conda install -c bioconda minimap2 samtools
```

## Quick Start
Tell your AI agent what you want to do:
- "Align my Nanopore reads to the reference genome"
- "Map PacBio HiFi reads and generate a sorted BAM"

## Example Prompts

### Basic Alignment
> "Align my ONT reads in reads.fastq.gz to reference.fa and output a sorted BAM"

> "Map my PacBio HiFi reads to the human genome using minimap2"

### With Read Groups
> "Align reads with sample ID 'patient1' and add read group information for downstream GATK"

### RNA/cDNA
> "Align my Nanopore direct RNA reads to the transcriptome using the splice preset"

### Performance
> "Build a minimap2 index for the reference genome so future alignments are faster"

## What the Agent Will Do
1. Select the appropriate minimap2 preset based on your data type (map-ont, map-hifi, map-pb, splice)
2. Run alignment with proper threading
3. Pipe output through samtools for sorting and indexing
4. Check alignment statistics with flagstat

## Choosing the Right Preset

| Data Type | Preset |
|-----------|--------|
| Nanopore (any) | map-ont |
| PacBio HiFi (Q20+) | map-hifi |
| PacBio CLR (older) | map-pb |
| cDNA/direct RNA | splice |

## Tips
- Always use the correct preset for your data type - wrong preset will give poor alignments
- Add read groups with `-R '@RG\tID:run1\tSM:sample1'` for downstream tools like GATK
- Pre-build index with `minimap2 -d ref.mmi ref.fa` for multiple alignment runs
- Check alignment rate with `samtools flagstat` - low rates indicate wrong preset or reference mismatch
- Use PAF format (`-x` without `-a`) when you only need alignment coordinates, not full BAM

## Resources
- [minimap2 Manual](https://lh3.github.io/minimap2/minimap2.html)
- [minimap2 Paper](https://doi.org/10.1093/bioinformatics/bty191)
