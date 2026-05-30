# genome-intervals

## Overview

Genomic interval algebra and signal tracks: coordinate systems, set operations, overlap-significance testing, GTF/GFF gene-model parsing, proximity queries, coverage, and bedGraph/bigWig track generation with bedtools/pybedtools/pyranges, gffutils, deepTools, and pyBigWig.

**Tool type:** mixed | **Primary tools:** bedtools, pybedtools, pyranges, gffutils, gffread, deepTools, pyBigWig, mosdepth

## Skills

| Skill | Description |
|-------|-------------|
| bed-file-basics | BED format, 0-based vs 1-based coordinate systems, sorting, validation, conversion |
| interval-arithmetic | intersect, subtract, merge, complement, map, groupby; the sorted-input and -split preconditions |
| overlap-significance | Test whether two interval sets overlap more than chance (permutation nulls, GAT, regioneR, LOLA, GREAT) |
| gtf-gff-handling | Parse, traverse, and convert GTF/GFF3 gene-model hierarchies; derive introns/UTRs/TSS |
| proximity-operations | closest, window, flank, slop; strand-aware promoters and nearest-feature assignment |
| coverage-analysis | genomecov, coverage, samtools depth, mosdepth; depth-vs-breadth and the distribution |
| bedgraph-handling | Create, normalize, and convert bedGraph signal tracks |
| bigwig-tracks | Create and read bigWig browser tracks; the zoom-level summary-statistic trap |

## Example Prompts

- "Find peaks that overlap with promoters"
- "Is the overlap between my peaks and enhancers more than expected by chance?"
- "Get intervals unique to sample A but not sample B"
- "Merge overlapping intervals in my BED file"
- "Convert these GTF coordinates to a BED file"
- "Extract gene, exon, and CDS coordinates from a GTF"
- "Find the nearest gene to each peak"
- "Build strand-aware promoter windows from a TSS list"
- "Extend intervals by 1kb on each side"
- "What fraction of my target is covered at 30x?"
- "Create a normalized bigWig from my BAM"
- "Convert bedGraph to bigWig"
- "Get the mean signal per gene from a bigWig"
- "Window my genome into 1kb bins"

## Requirements

```bash
# CLI tools
conda install -c bioconda bedtools gffread ucsc-bedgraphtobigwig ucsc-bigwigaverageoverbed ucsc-liftover crossmap samtools mosdepth deeptools gat

# Python
pip install pybedtools pyranges pyBigWig gffutils gtfparse bioframe

# R (overlap-significance)
# BiocManager::install(c('regioneR', 'LOLA', 'rGREAT'))
```

pybedtools shells out to the bedtools binary, so bedtools must be on PATH.

## Related Skills

- **alignment-files** - BAM pileup/depth feeding coverage-analysis and track generation
- **read-alignment** - Produces the BAMs that coverage and bedGraph operations consume
- **chip-seq** - Peak BED files this category operates on; spike-in normalization and track visualization
- **atac-seq** - ATAC peak sets for arithmetic/proximity; footprinting over coverage signal
- **variant-calling** - VCF positions converted to BED and region-filtered by interval
- **rna-quantification** - Consumes GTF/GFF features parsed here for read counting
- **data-visualization** - Renders the bedGraph/bigWig tracks built here
