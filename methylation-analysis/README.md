# methylation-analysis

## Overview

DNA methylation analysis using Bismark for bisulfite sequencing alignment, methylKit/bsseq for downstream analysis, and scipy/limma for per-CpG differential testing. Covers alignment, methylation calling, per-CpG statistical testing, and detection of differentially methylated regions (DMRs).

**Tool type:** mixed | **Primary tools:** Bismark (CLI), methylKit (R), bsseq (R), scipy (Python)

## Skills

| Skill | Description |
|-------|-------------|
| bismark-alignment | Bisulfite read alignment with Bismark |
| methylation-calling | Extract methylation calls from Bismark output |
| methylkit-analysis | Methylation analysis with methylKit in R |
| differential-cpg-testing | Per-CpG differential methylation testing |
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
- "Run per-CpG differential methylation testing with Welch's t-test"
- "Test individual CpG sites for methylation differences using limma on M-values"
- "Compute beta values and delta-beta effect sizes from my BS-seq counts"
- "Find differentially methylated regions between conditions"
- "Identify DMRs with at least 25% methylation difference"
- "Annotate my DMRs with gene information"

## Requirements

```bash
# Bismark
conda install -c bioconda bismark bowtie2

# Python (for per-CpG testing)
pip install numpy pandas scipy statsmodels
```

```r
BiocManager::install(c('methylKit', 'bsseq', 'GenomicRanges', 'limma', 'DSS'))
```

## Related Skills

- **alignment-files** - BAM file manipulation after alignment
- **sequence-io** - FASTQ handling before alignment
- **pathway-analysis** - Functional annotation of genes near DMRs
