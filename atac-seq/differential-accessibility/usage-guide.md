# Differential Accessibility - Usage Guide

## Overview

Identify chromatin regions with significantly different accessibility between conditions. Choose between consensus-peak workflows (DiffBind, DESeq2 directly) and peak-free sliding-window workflows (csaw). Pick the correct normalization (reads-in-peaks vs full-library vs spike-in) for the biological setting, correct batch effects with SVA/RUVseq, and report effect sizes with calibrated FDR.

## Prerequisites

```r
BiocManager::install(c('DiffBind', 'DESeq2', 'edgeR', 'csaw', 'limma',
                       'ChIPseeker', 'sva', 'RUVSeq', 'GenomicRanges',
                       'TxDb.Hsapiens.UCSC.hg38.knownGene', 'org.Hs.eg.db'))
```

```bash
conda install -c bioconda subread bedtools homer
```

Inputs: per-sample BAM files (deduplicated, chrM-stripped), per-sample peak files (narrowPeak), and a sample sheet defining conditions and any covariates.

## Quick Start

Tell your AI agent what you want to do:
- "Run DiffBind with default DESeq2 backend on 3 control vs 3 treated ATAC samples"
- "Use edgeR QL F-test because I only have 2 replicates per condition"
- "Switch DiffBind normalization from reads-in-peaks to full-library because the treatment causes global compaction"
- "Run csaw sliding-window analysis since differentiation collapses peak calls"
- "Add SVA surrogate variables for hidden batch and re-fit"
- "Filter to FDR < 0.05 AND abs(log2FC) >= 1, then annotate to nearest gene with ChIPseeker"

## Example Prompts

### Standard 3v3 Comparison
> "Run DiffBind on three replicates per group with DESeq2 backend, `summits=250` for fixed-width counting, native normalization (`DBA_NORM_NATIVE`, i.e. DESeq2 RLE), and report peaks with FDR < 0.05 and abs(log2FC) >= 1."

### Low Replicate Power
> "I have 2 reps per condition. Use DiffBind with edgeR backend (better at low n) or skip DiffBind and run edgeR QL directly. Do not use DESeq2 + apeglm shrinkage because shrinkage at n=2 is over-aggressive."

### Global Accessibility Shift
> "Treatment is an HDAC inhibitor that globally opens chromatin. Run DiffBind with `normalize=DBA_NORM_LIB` (full library; also DiffBind's default) instead of reads-in-peaks (`library=DBA_LIBSIZE_PEAKREADS`), since RiP normalization will erase the global biology."

### Window-based Analysis
> "Conditions are pre- and post-differentiation. Peak structure differs dramatically. Use csaw with `width=150` windows, `filterWindowsGlobal` at 3-fold over background (`$filter > log2(3)`), edgeR QL F-test, then merge significant adjacent windows."

### Hidden Batch Correction
> "Three batches of samples processed weeks apart. Run svaseq with n.sv=2 to estimate hidden surrogates, then re-fit DESeq2 with `~SV1 + SV2 + condition` design."

### Multi-factor Design
> "Adjust for batch and donor as covariates. Use DiffBind `dba.contrast(..., design='~Batch + Donor + Condition')` with explicit factor levels."

### Annotate and Enrich
> "After DA analysis, annotate peaks to genes using ChIPseeker with promoter region `-2000 to +500`, then run GO enrichment on opened-vs-closed gene lists separately."

## What the Agent Will Do

1. Verify input BAMs are deduplicated, MAPQ-filtered, chrM-stripped
2. Build consensus peakset (DiffBind built-in OR custom iterative overlap)
3. Re-center peaks on summits +/- 250 bp for fixed-width counting
4. Choose normalization based on whether treatment causes global accessibility shift
5. Pick statistical method by replicate count and peak-structure stability
6. Add batch / donor / time as covariates if applicable
7. Run SVA or RUVseq for hidden batch when needed
8. Apply FDR < 0.05 and abs(log2FC) >= 1 thresholds (or stricter for high-confidence)
9. Annotate differential peaks with ChIPseeker; split opened vs closed lists
10. Pass gene lists downstream to GO enrichment skill

## Sample Sheet Format

```csv
SampleID,Condition,Replicate,Batch,bamReads,Peaks,PeakCaller
ctrl_rep1,control,1,A,ctrl_rep1.bam,ctrl_rep1_peaks.narrowPeak,macs
ctrl_rep2,control,2,B,ctrl_rep2.bam,ctrl_rep2_peaks.narrowPeak,macs
treat_rep1,treated,1,A,treat_rep1.bam,treat_rep1_peaks.narrowPeak,macs
treat_rep2,treated,2,B,treat_rep2.bam,treat_rep2_peaks.narrowPeak,macs
```

Add `Donor`, `Sex`, `Time`, etc. columns as needed and include in the design formula.

## Method Selection Quick Reference

| Replicate count | Method |
|-----------------|--------|
| 1 per group | Cannot do statistical DA; report log2FC only with caveat |
| 2 per group | edgeR QL (DiffBind edgeR backend); avoid DESeq2 + apeglm |
| 3 per group | DiffBind DESeq2 (default) is fine |
| 5+ per group | DESeq2 with apeglm shrinkage; can also use limma-voom |
| 10+ heterogeneous | DESeq2 LRT for primary effect + interactions |

## Normalization Quick Reference

| Treatment biology | Normalization |
|-------------------|---------------|
| Localized differential (most experiments) | DiffBind `DBA_NORM_NATIVE` (RLE for DESeq2 / TMM for edgeR; DiffBind's own default is `DBA_NORM_LIB`) |
| Global accessibility shift (HDACi, DNMTi) | `DBA_NORM_LIB` OR exogenous spike-in |
| Outlier-prone counts | `DBA_NORM_TMM` (edgeR-style) |
| Cross-cell-type comparison | `DBA_NORM_LIB` (cell-type biology often global) |

## Effect Size Quick Reference

| Goal | Threshold |
|------|-----------|
| Standard reporting | FDR < 0.05, abs(log2FC) >= 1 |
| Conservative high-confidence | FDR < 0.01, abs(log2FC) >= 1 |
| Primary cell types (smaller effects) | FDR < 0.05, abs(log2FC) >= 0.585 (1.5x) |
| Discovery for follow-up validation | FDR < 0.1, shrunken log2FC magnitude reported |

## Tips

- DiffBind's reads-in-peaks normalization (`library=DBA_LIBSIZE_PEAKREADS`) erases global accessibility shifts; keep the full-library default (`DBA_NORM_LIB`/`DBA_LIBSIZE_FULL`) whenever the biology is whole-genome.
- Always set `summits=250` to re-center peaks on summit; otherwise peak-width differences inflate counts.
- For low replicates (n=2), edgeR QL is the safe choice. DESeq2 + apeglm shrinkage at n=2 over-shrinks.
- csaw is the only valid choice when peak structure differs dramatically between conditions; don't force a stable consensus when one doesn't exist.
- Filter low-count peaks aggressively (`filterByExpr`) before fitting; tail peaks inflate dispersion and FDR.
- Compute log2FC on shrunken estimates only when n >= 3; below that, report unshrunken with caveat.
- For peak-to-gene assignment, use `tssRegion=c(-2000, 500)` not the ChIPseeker default `(-3000, 3000)` -- the default over-counts promoter assignments by lumping in distal elements.
- Spike-in normalization is the strongest approach for global-shift biology (Reske 2020 raises spike-in as a candidate for such settings but does not test it); budget for it in experimental design.
- Reproducibility: report the consensus peakset strategy, normalization, and exact tool versions; results vary across DiffBind versions even with identical inputs.

## Related Skills

- atac-seq/atac-peak-calling - Generate input peaks per replicate
- atac-seq/consensus-peakset - Build differential-ready consensus
- atac-seq/atac-qc - Drop failing replicates before differential
- atac-seq/single-cell-atac - Pseudobulk-level differential
- atac-seq/co-accessibility - Cis-regulatory follow-up
- differential-expression/deseq2-basics - Underlying DESeq2 patterns
- differential-expression/de-results - Effect-size shrinkage details
- chip-seq/differential-binding - Same DiffBind workflow for ChIP
- pathway-analysis/go-enrichment - Downstream gene-level enrichment
