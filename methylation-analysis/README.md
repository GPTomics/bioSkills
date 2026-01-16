# Methylation Analysis

DNA methylation analysis using Bismark for bisulfite sequencing alignment and methylKit/bsseq for downstream analysis.

## Overview

This category covers bisulfite sequencing data analysis: Bismark handles alignment and methylation calling from raw reads, while methylKit and bsseq provide differential methylation analysis. DNA methylation is a key epigenetic modification involved in gene regulation, development, and disease.

**Tool type:** `mixed`
**Primary tools:** Bismark (CLI), methylKit (R), bsseq (R)

## Skills

| Skill | Description |
|-------|-------------|
| [bismark-alignment](bismark-alignment/) | Bisulfite read alignment with Bismark |
| [methylation-calling](methylation-calling/) | Extract methylation calls from Bismark output |
| [methylkit-analysis](methylkit-analysis/) | Methylation analysis with methylKit in R |
| [dmr-detection](dmr-detection/) | Differentially methylated region detection |

## Workflow

```
Raw Bisulfite Reads (FASTQ)
    |
    v
[bismark-alignment] -----> Align to bisulfite-converted genome
    |
    v
[methylation-calling] ---> Extract CpG methylation levels
    |
    v
[methylkit-analysis] ----> Import, filter, normalize
    |
    v
[dmr-detection] ---------> Find differentially methylated regions
    |
    v
[pathway-analysis] ------> Annotate and interpret DMRs
```

## Bisulfite Sequencing Types

| Method | Coverage | Cost | Use Case |
|--------|----------|------|----------|
| WGBS | Whole genome | High | Comprehensive profiling |
| RRBS | CpG-rich regions | Medium | Cost-effective, CpG islands |
| Targeted | Specific regions | Low | Validation, specific loci |

## Methylation Context

| Context | Pattern | Description |
|---------|---------|-------------|
| CpG | 5'-CG-3' | Most common, often in islands |
| CHG | 5'-CHG-3' | Plant-specific |
| CHH | 5'-CHH-3' | Plant-specific |

(H = A, C, or T)

## Example Prompts

### Bismark Alignment
- "Align my bisulfite sequencing reads with Bismark"
- "How do I prepare a genome for bisulfite alignment?"
- "Run Bismark on paired-end RRBS data"

### Methylation Calling
- "Extract CpG methylation levels from my BAM file"
- "Get methylation calls from Bismark output"
- "Create a coverage file for my methylation data"

### methylKit Analysis
- "Load my methylation data into methylKit"
- "Normalize my bisulfite sequencing samples"
- "Compare methylation between treatment groups"

### DMR Detection
- "Find differentially methylated regions between conditions"
- "Identify DMRs with at least 25% methylation difference"
- "Annotate my DMRs with gene information"

## Requirements

### Bismark (CLI)

```bash
# Conda installation (recommended)
conda install -c bioconda bismark

# Also requires bowtie2 or hisat2
conda install -c bioconda bowtie2

# Verify installation
bismark --version
```

### R Packages

```r
if (!require('BiocManager', quietly = TRUE))
    install.packages('BiocManager')

BiocManager::install(c('methylKit', 'bsseq', 'GenomicRanges'))

# For annotation
BiocManager::install(c('annotatr', 'TxDb.Hsapiens.UCSC.hg38.knownGene'))
```

## Key Functions

| Function | Tool | Purpose |
|----------|------|---------|
| bismark_genome_preparation | Bismark | Prepare bisulfite genome index |
| bismark | Bismark | Align bisulfite reads |
| bismark_methylation_extractor | Bismark | Extract methylation calls |
| methRead | methylKit | Read methylation data |
| filterByCoverage | methylKit | Filter low coverage CpGs |
| normalizeCoverage | methylKit | Normalize between samples |
| calculateDiffMeth | methylKit | Differential methylation |
| getMethylDiff | methylKit | Extract significant DMCs |
| tileMethylCounts | methylKit | Tile-based DMR analysis |
| read.bismark | bsseq | Read Bismark output |
| BSmooth | bsseq | Smooth methylation data |

## Notes

- **Genome preparation required** - Must run bismark_genome_preparation before alignment
- **Directional libraries** - Most protocols are directional; use `--non_directional` if not
- **Coverage filtering** - Remove CpGs with <10x coverage for reliability
- **Multiple testing** - Use q-value < 0.05 for DMC/DMR calls
- **Strand collapse** - CpG methylation is usually symmetric; collapse strands

## Related Skills

- **alignment-files** - BAM file manipulation after alignment
- **sequence-io** - FASTQ handling before alignment
- **differential-expression** - Similar statistical approaches for DMRs
- **pathway-analysis** - Functional annotation of genes near DMRs

## References

- [Bismark User Guide](https://rawgit.com/FelixKrueger/Bismark/master/Docs/Bismark_User_Guide.html)
- [methylKit Bioconductor](https://bioconductor.org/packages/methylKit/)
- [bsseq Bioconductor](https://bioconductor.org/packages/bsseq/)
