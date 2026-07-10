# Methylation Pipeline - Usage Guide

## Overview

This workflow processes bisulfite sequencing data from FASTQ to differentially methylated regions (DMRs), covering alignment, methylation calling, and statistical analysis.

## Prerequisites

```bash
# CLI tools
conda install -c bioconda bismark bowtie2 trim-galore samtools

# R packages
BiocManager::install(c('methylKit', 'genomation'))
```

## Quick Start

Tell your AI agent what you want to do:
- "Run the methylation pipeline on my bisulfite-seq data"
- "Find differentially methylated regions between treatment and control"
- "Align my WGBS data with Bismark and call methylation"

## Example Prompts

### Starting from FASTQ
> "Process my RRBS data through methylation calling"

> "Align bisulfite sequencing reads and extract methylation"

> "Run Bismark on my whole-genome bisulfite data"

### Analysis
> "Find DMRs with methylKit"

> "Compare methylation between tumor and normal"

> "Run per-CpG differential methylation with a beta-binomial count model (DSS) instead of methylKit"

> "Annotate my DMRs with gene features"

### Pipeline-level decisions
> "How do I check bisulfite conversion, and why does it matter before I compute methylation levels?"

> "My coverage looks doubled in mate-overlap regions - did I call methylation correctly?"

> "Should I use methylKit tiles or dmrseq/DSS for region-level DMRs?"

> "Is a t-test on beta values good enough, or do I need a count model?"

## Input Requirements

| Input | Format | Description |
|-------|--------|-------------|
| FASTQ files | .fastq.gz | Paired-end bisulfite-treated reads |
| Reference | FASTA | Genome (Bismark will prepare) |

## What the Workflow Does

1. **Quality Control** - Trim adapters and low-quality bases
2. **Alignment** - Map bisulfite-converted reads with Bismark
3. **Deduplication** - Remove PCR duplicates
4. **Methylation Calling** - Extract methylation status per CpG
5. **Per-CpG Analysis** - Statistical testing with methylKit (R) or scipy (Python)
6. **DMR Detection** - Find differentially methylated regions

## Tips

- **Coverage**: WGBS needs 10-30x coverage; RRBS can work with less
- **RRBS: skip deduplication** - reads legitimately stack at MspI cut sites, so positional dedup destroys real coverage (dedup WGBS/EM-seq only)
- **Conversion rate**: Should be >99%; check with spike-in controls
- **M-bias**: Check for position bias and trim if needed
- **Replicates**: Minimum 2-3 per condition for reliable DMR calling
- **Per-CpG test choice**: Sequencing counts carry coverage (precision), so route them to a beta-binomial / overdispersion count model (DSS, or methylKit overdispersion='MN'); a bare-beta t-test discards coverage and is only a quick look. Array/continuous data uses limma on M-values. See methylation-analysis/differential-cpg-testing.

## Related Skills

- methylation-analysis/bismark-alignment - Bisulfite/EM-seq alignment, library/strand model, conversion QC
- methylation-analysis/methylation-calling - Per-CpG calling from BAM (Bismark/MethylDackel), contexts, variant-aware
- methylation-analysis/methylkit-analysis - methylKit object model and overdispersion gotchas
- methylation-analysis/differential-cpg-testing - Per-CpG testing (count-vs-continuous fork)
- methylation-analysis/dmr-detection - Selection-aware region callers (dmrseq/DSS) and PMD segmentation
- methylation-analysis/array-preprocessing - Alternate entry: Infinium IDAT to beta/M matrix
- methylation-analysis/cell-type-deconvolution - Cell-fraction covariates for bulk-tissue EWAS
- methylation-analysis/epigenetic-clocks - DNAm age and age acceleration
- methylation-analysis/ewas-design - EWAS confounding, batch, inflation, and replication
