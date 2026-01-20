# Nucleosome Positioning Usage Guide

Extract nucleosome positions from ATAC-seq fragment size patterns.

## Prerequisites

```r
# R packages
BiocManager::install('ATACseqQC')
```

```bash
# NucleoATAC
pip install nucleoatac
```

## Background

ATAC-seq fragments reflect chromatin structure:
- **< 100 bp**: Nucleosome-free regions (NFR)
- **180-247 bp**: Mono-nucleosome
- **315-473 bp**: Di-nucleosome
- **558-615 bp**: Tri-nucleosome

## Quick Start

### ATACseqQC (R)

```r
library(ATACseqQC)

# Read BAM
bam <- 'atac.bam'

# Fragment size distribution
fragSize <- fragSizeDist(bam, 'sample')
```

### Separate by Fragment Size

```r
library(ATACseqQC)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

# Split reads by nucleosome occupancy
objs <- splitGAlignmentsByCut(gal, txs=txs, genome=genome)

# NFR: nucleosome-free
# mono: mononucleosome
# di: dinucleosome
```

### NucleoATAC

```bash
# Call nucleosomes
nucleoatac run \
    --bed peaks.bed \
    --bam atac.bam \
    --fasta genome.fa \
    --out nucleoatac_output
```

## Fragment Size Analysis

```r
library(ATACseqQC)
library(ggplot2)

# Get fragment sizes
sizes <- fragSizeDist(bamfile, 'sample')

# Plot
ggplot(sizes, aes(x=fragLen, y=proportion)) +
    geom_line() +
    geom_vline(xintercept=c(100, 180, 247), linetype='dashed') +
    labs(x='Fragment Length', y='Proportion')
```

## Nucleosome Occupancy

```bash
# Generate nucleosome signal track
nucleoatac nuc \
    --bed peaks.bed \
    --bam atac.bam \
    --fasta genome.fa \
    --out nuc_signal
```

## TSS Nucleosome Patterns

```r
# V-plot around TSS
library(ATACseqQC)

tsse <- TSSEscore(gal, txs)
# NFR at TSS with flanking +1/-1 nucleosomes
```

## Common Analyses

### NFR at Promoters

```r
# Extract NFR reads
nfr_reads <- objs$NFR

# Coverage at promoters
promoter_cov <- coverage(nfr_reads)
```

### Nucleosome Phasing

```bash
# NucleoATAC outputs phased nucleosome positions
# Look for regular ~180bp spacing
```

## Interpretation

| Pattern | Meaning |
|---------|---------|
| Strong NFR peak | Active promoter/enhancer |
| Regular spacing | Well-positioned nucleosomes |
| Fuzzy positioning | Dynamic chromatin |

## Tips

- Use properly deduplicated BAM files
- Paired-end data required for accurate sizing
- Filter for high-quality fragments (MAPQ > 30)
- Consider mitochondrial contamination

## See Also

- [ATACseqQC vignette](https://bioconductor.org/packages/ATACseqQC/)
- [NucleoATAC paper](https://genome.cshlp.org/content/25/11/1757)
