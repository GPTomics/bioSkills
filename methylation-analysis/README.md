# methylation-analysis

## Overview

DNA methylation analysis using Bismark for bisulfite sequencing alignment and methylKit/bsseq for downstream analysis. Covers alignment, methylation calling, and detection of differentially methylated regions (DMRs).

**Tool type:** mixed | **Primary tools:** Bismark (CLI), methylKit (R), bsseq (R)

## Skills

| Skill | Description |
|-------|-------------|
| bismark-alignment | Bisulfite read alignment with Bismark |
| methylation-calling | Extract methylation calls from Bismark output |
| methylkit-analysis | Methylation analysis with methylKit in R |
| dmr-detection | Differentially methylated region detection |

## Example Prompts

- "Align my bisulfite sequencing reads with Bismark"
- "How do I prepare a genome for bisulfite alignment?"
- "Run Bismark on paired-end RRBS data"
- "Extract CpG methylation levels from my BAM file"
- "Get methylation calls from Bismark output"
- "Create a coverage file for my methylation data"
- "Load my methylation data into methylKit"
- "Normalize my bisulfite sequencing samples"
- "Compare methylation between treatment groups"
- "Find differentially methylated regions between conditions"
- "Identify DMRs with at least 25% methylation difference"
- "Annotate my DMRs with gene information"

## Requirements

```bash
# Bismark
conda install -c bioconda bismark bowtie2
```

```r
BiocManager::install(c('methylKit', 'bsseq', 'GenomicRanges'))
```

## Related Skills

- **alignment-files** - BAM file manipulation after alignment
- **sequence-io** - FASTQ handling before alignment
- **pathway-analysis** - Functional annotation of genes near DMRs
