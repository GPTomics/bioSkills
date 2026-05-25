# differential-expression

## Overview

Postdoc-grade differential expression analysis of bulk RNA-seq count data with DESeq2, edgeR, and limma-voom. Encodes the decision space (Wald vs LRT, apeglm vs ashr, TREAT vs post-hoc filtering, design-inclusion vs ComBat-then-test, factor vs spline for time), the canonical failure modes (Nygaard 2016 batch-correction cardinal sin, padj=NA three-meanings, alphabetical reference level inverting fold changes, edgeR v4 legacy=FALSE flip), and the modern alternatives (pseudobulk DE, transcript-DE via catchSalmon, DREAM mixed models for repeated measures).

**Tool type:** r | **Primary tools:** DESeq2, edgeR, limma, apeglm, sva, IHW

## Skills

| Skill | Description |
|-------|-------------|
| deseq2-basics | DESeq2 negative-binomial GLM, Wald and LRT testing, apeglm/ashr/normal shrinkage, Cook's outliers, independent filtering, vst/rlog, tximport |
| edger-basics | edgeR QL F-test (modern default), TMM/TMMwsp, TREAT, voomWithQualityWeights, edgeR v4 legacy flip, transcript-DE via catchSalmon |
| de-results | padj=NA three meanings, BH/Storey/IHW/lfsr, TREAT vs post-hoc, p-value histogram diagnostics, gene annotation, GSEA ranked input, ORA proper background, SABV |
| de-visualization | DESeq2/edgeR built-in plots; MA and volcano with shrunken-LFC compression; PCA on VST/rlog never raw; heatmap row-scaling trap; blind=TRUE vs FALSE |
| batch-correction | Design-inclusion vs Nygaard 2016 cardinal sin; ComBat vs ComBat-seq; SVA-captures-biology; RUVg/s/r; removeBatchEffect visualization-only; single-cell boundary |
| timeseries-de | DESeq2 LRT, splines, maSigPro (Nueda 2014), ImpulseDE2 impulse-assumption failure modes, DREAM for repeated measures, pseudoreplication, trajectory clustering (DPGP, Mfuzz) |

## Example Prompts

- "Run DESeq2 on treated vs control and give me the top 20 genes at padj<0.05 with apeglm-shrunken fold changes"
- "Apply log fold change shrinkage with apeglm for my GSEA ranked list"
- "I have paired tumor/normal from 8 patients -- set up the design correctly"
- "Why is padj NA for my master regulator? Help me diagnose"
- "Test whether the drug effect differs between WT and KO genotypes"
- "Run edgeR QL F-test with robust dispersion and TREAT for a 1.5-fold biological threshold"
- "Include batch as a covariate in DESeq2 -- do NOT subtract before DE"
- "Build a volcano plot with apeglm-shrunken LFC on x and unshrunken Wald p on y"
- "Diagnose this p-value histogram and tell me what's wrong with the model"
- "Detect any sample swaps in my cohort before I run DE"
- "Fit a time-course DE with splines and test the treatment-time interaction"
- "Cluster my time-course DE genes into trajectory shapes with Mfuzz"
- "Reproduce a pre-2025 edgeR result -- set legacy=TRUE"
- "Prepare a ranked gene list for GSEA from DESeq2 results"
- "Run DE on prokaryotic RNA-seq under stress -- TMM/RLE assumption may fail"
- "Run DE using PyDESeq2 in Python"

## Requirements

```r
BiocManager::install(c('DESeq2', 'edgeR', 'limma', 'apeglm', 'ashr',
                        'IHW', 'tximport', 'tximeta', 'sva', 'RUVSeq',
                        'variancePartition', 'ImpulseDE2', 'maSigPro', 'Mfuzz',
                        'AnnotationDbi', 'org.Hs.eg.db', 'clusterProfiler',
                        'EnhancedVolcano', 'qvalue'))
install.packages(c('ggplot2', 'pheatmap', 'ggrepel', 'dplyr', 'splines',
                    'matrixStats', 'UpSetR', 'openxlsx'))
```

```bash
pip install pydeseq2 pandas numpy
```

## Related Skills

- **expression-matrix** - Count matrix loading, normalization choice, gene ID mapping, metadata joins, sparse handling
- **rna-quantification** - featureCounts, salmon, kallisto, tximport upstream of DE
- **pathway-analysis** - ORA (with proper background) and GSEA (with shrunken LFC or Wald stat) downstream of DE
- **data-visualization** - Custom volcano/MA with ggrepel, PCA/UMAP, heatmap customization, UpSet plots
- **read-qc** - RIN/DV200 as covariates; degradation effects on DE
- **single-cell** - Pseudobulk DE; cell-level DE FDR-inflation problem
- **temporal-genomics** - Circadian (JTK_CYCLE, MetaCycle), trajectory modeling beyond DE
- **alignment-files** - Process BAM files for counting
- **workflows** - End-to-end pipelines wrapping these skills
