# Small RNA-seq Pipeline - Usage Guide

## Overview

Complete workflow from small RNA-seq FASTQ to differential miRNA expression and target prediction, with the small-RNA-specific decisions surfaced at each step. The pipeline trims kit-specific adapters, quantifies known miRNAs and isomiRs (miRge3) or discovers novel miRNAs (miRDeep2), runs compositionally-aware differential expression (DESeq2 on raw counts), and predicts targets filtered by expression evidence. When abundant non-miRNA classes appear, it branches to tRF/piRNA profiling.

## Prerequisites

```bash
conda install -c bioconda cutadapt mirdeep2 bowtie viennarna miranda
pip install mirge3 umi_tools
# BiocManager::install(c('DESeq2', 'apeglm'))
```

## Quick Start

- "Analyze my small RNA-seq from FASTQ to differential miRNAs and targets"
- "Quantify known miRNAs with miRge3, then run DESeq2 on the raw counts"
- "Discover novel miRNAs with miRDeep2 and validate the candidates"
- "My library is mostly tRFs, not miRNAs - what should I run?"

## Example Prompts

### Full Pipeline

> "Run the complete small RNA-seq pipeline from trimming to target prediction"

> "Find DE miRNAs between my conditions and predict their targets"

### Specific Steps

> "Just quantify known miRNAs and isomiRs with miRge3"

> "Predict targets for my DE miRNAs and intersect with my mRNA-seq"

> "Profile the tRFs and piRNAs instead of miRNAs"

### Pipeline-level decisions

> "Which normalization should I use - my library is dominated by a few miRNAs?"

> "My samples are plasma/serum - how do I normalize when there's no stable endogenous miRNA?"

> "Should I quantify known miRNAs (miRge3) or discover novel ones (miRDeep2)?"

> "My size factors look skewed and lots of miRNAs are 'down' - is that real?"

## What the Agent Will Do

1. Trim the kit adapter (discarding untrimmed reads) and handle any 4N spacer or UMI
2. Quantify known miRNAs/isomiRs with miRge3, or discover novel miRNAs with miRDeep2 when needed
3. Check the RNA-class composition and branch to tRF/piRNA profiling if non-miRNA classes dominate
4. Run DESeq2 on raw counts, inspecting size factors for compositional distortion
5. Predict targets and filter them by anti-correlated mRNA DE and validated databases

## Tips

- The 3' adapter is on every real read, so `--discard-untrimmed` is correct; strip NEXTflex 4N after trimming, extract QIAseq UMIs before/after alignment
- The length histogram is QC: a 21-23 nt peak is a good miRNA library; a 26-32 nt peak is piRNA; a 30+ nt smear is degradation or contamination
- Ligation bias makes absolute cross-miRNA abundance untrustworthy; compare the same miRNA across samples, never across kits
- Feed RAW counts to DESeq2 (RPM is display only) and inspect size factors - a few dominant miRNAs distort normalization
- A miRDeep2 score is a hypothesis; pick a cutoff from survey.pl signal-to-noise and validate novel calls
- A predicted target is a hypothesis; filter by anti-correlated mRNA DE before trusting it
- If non-miRNA classes dominate, this is a tRF/piRNA analysis - switch tools accordingly

## Related Skills

- small-rna-seq/smrna-preprocessing - Adapter, UMI, and 4N handling
- small-rna-seq/mirdeep2-analysis - Novel miRNA discovery
- small-rna-seq/mirge3-analysis - Known-miRNA and isomiR quantification
- small-rna-seq/differential-mirna - Compositionally-aware DE
- small-rna-seq/target-prediction - Seed prediction filtered by expression
- small-rna-seq/trf-pirna-profiling - tRF and piRNA profiling
- workflow-management/nf-core-pipelines - Run nf-core/smrnaseq as a curated, reproducible pipeline instead of chaining the steps by hand
