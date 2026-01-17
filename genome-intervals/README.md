# genome-intervals

## Overview

Genomic interval operations using BEDTools, pybedtools, and pyBigWig. Covers BED file manipulation, interval arithmetic, GTF/GFF parsing, coverage analysis, and bigWig track generation.

**Tool type:** mixed | **Primary tools:** BEDTools, pybedtools, pyBigWig

## Skills

| Skill | Description |
|-------|-------------|
| bed-file-basics | BED format, file creation, validation |
| interval-arithmetic | intersect, subtract, merge, complement operations |
| gtf-gff-handling | Parse and convert GTF/GFF annotation files |
| proximity-operations | closest, window, flank, slop for proximity queries |
| coverage-analysis | genomecov, coverage calculations, bedGraph generation |
| bedgraph-handling | Create, convert, and manipulate bedGraph files |
| bigwig-tracks | Create and read bigWig browser tracks |

## Example Prompts

- "Find peaks that overlap with promoters"
- "Get intervals unique to sample A but not sample B"
- "Merge overlapping intervals in my BED file"
- "Extract gene coordinates from GTF"
- "Find the nearest gene to each peak"
- "Add 500bp flanks to my intervals"
- "Extend intervals by 1kb on each side"
- "Calculate coverage across my BED regions"
- "Create a normalized bedGraph from my BAM"
- "Convert bedGraph to bigWig"
- "Merge bedGraph files from multiple samples"
- "Get the complement of my intervals"
- "Window my genome into 1kb bins"
- "Parse a GTF file and extract exon coordinates"

## Requirements

```bash
# BEDTools (CLI)
conda install -c bioconda bedtools

# pybedtools
pip install pybedtools

# pyBigWig
pip install pyBigWig

# GTF/GFF parsing
pip install gtfparse gffutils
```

## Related Skills

- **alignment-files** - BAM file processing, depth calculation
- **variant-calling** - VCF region filtering
- **chip-seq** - Peak file manipulation
- **rna-quantification** - GTF annotation for counting
