# Differential miRNA Expression - Usage Guide

## Overview

Identify differentially expressed miRNAs between conditions with DESeq2 or edgeR, accounting for the fact that miRNA data is compositionally harder than mRNA. A handful of hyper-abundant miRNAs can be more than half of all reads, so global-scaling normalization is fragile and the normalization choice - not the DE engine - is the dominant decision. The workflow filters near-zero noise first, inspects whether a few miRNAs dominate, tests RAW counts (RPM is display only), and shrinks fold-changes with apeglm for the many low-count miRNAs. Biofluid data has no trustworthy endogenous normalizer and needs spike-ins plus hemolysis and batch control, and a miRNA dropping can reflect target-directed degradation (TDMD) rather than transcriptional repression.

## Prerequisites

```r
BiocManager::install(c('DESeq2', 'edgeR', 'apeglm', 'EnhancedVolcano', 'pheatmap'))
# Optional compositional sensitivity analysis: BiocManager::install('ALDEx2')
```

## Quick Start

Tell your AI agent:
- "Run DESeq2 on my raw miRge3 count matrix between treatment and control"
- "Check whether a few dominant miRNAs are distorting my normalization"
- "Cross-check the result with edgeR"
- "My data is plasma - which normalizer should I use?"
- "Shrink the fold-changes with apeglm and report base mean for each hit"

## Example Prompts

### Basic DE Analysis

> "Find DE miRNAs between tumor and normal from my raw count matrix"

> "My size factors are far from 1 and one miRNA is 60% of reads - what normalization should I use?"

> "Run edgeR quasi-likelihood as a cross-check of the DESeq2 calls"

### Compositional and Biofluid Caveats

> "Run an ALDEx2 compositional sensitivity analysis alongside DESeq2"

> "Normalize my plasma miRNA-seq with cel-miR-39 spike-ins and flag hemolysis"

> "A perturbation may have changed the whole miRNA pool - how do I even detect that?"

### Visualization and Export

> "Make a volcano plot and a heatmap of significant miRNAs"

> "Export significant miRNAs with their base mean and shrunk log2FC"

## What the Agent Will Do

1. Load the RAW count matrix (not RPM) and build sample metadata
2. Prefilter at a lower threshold than mRNA and inspect for compositional dominance (size factors, top-miRNA fraction)
3. Run DESeq2 (or edgeR TMM) and apply apeglm LFC shrinkage for low-count miRNAs
4. Filter on FDR and effect size while keeping expression level visible, to avoid calling tiny miRNAs DE
5. Visualize and export, and flag biofluid/TDMD caveats where relevant

## Tips

- The normalization choice changes the result more than the DE engine does - filter first, inspect dominance, and report the normalizer
- RPM is for display only; feed RAW counts to DESeq2/edgeR
- Use a lower prefilter than mRNA (miRNAs have fewer total counts and most miRBase entries are noise), but justify the threshold
- apeglm shrinkage matters more here because few features weaken the dispersion prior
- edgeR names its FDR column `FDR`, not `padj`; pass `group` to `filterByExpr`
- Biofluids have no endogenous normalizer - use cel-miR-39 spike-ins, flag hemolysis (miR-451a:miR-23a-3p), and model batch
- A global pool shift (e.g. Dicer/Drosha loss) is invisible to all internal normalizers; only spike-ins or cell-number normalization detect it
- A miRNA going down can be TDMD (target-driven decay), not transcriptional repression - check pri/pre-miRNA before claiming repression
- Decide mature-level vs isomiR-level DE: collapse to mature for a well-powered "which miRNAs changed" question; isomiR-level is sparser and 5' isomiRs can move opposite the canonical form, so never silently sum them
- apeglm shrinks a named coef; use ashr for arbitrary or multi-level contrasts

## Related Skills

- mirge3-analysis - Produces the raw count matrix
- mirdeep2-analysis - Alternative quantification
- target-prediction - Targets of DE miRNAs
- differential-expression/deseq2-basics - General DESeq2 mechanics
- differential-expression/edger-basics - General edgeR mechanics
