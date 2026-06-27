# Hashtag Demultiplexing and Cross-Sample Doublet Calling - Usage Guide

## Overview

Hashtag demultiplexing assigns pooled, hashed cells back to their sample of origin and calls cross-sample doublets directly from the HTO count matrix. Each sample is labeled before pooling with a unique oligo-tagged reagent (CITE-seq antibody HTO, MULTI-seq lipid/cholesterol tag, or CellPlex CMO), so a singlet is dominated by one tag above background and a cross-sample doublet shows two tags high. This skill covers Seurat HTODemux and MULTIseqDemux, scanpy hashsolo, pegasus/demuxEM, GMM-Demux, and demuxmix, when to pick each, and how hashtag demux relates to genetic demux and expression-based doublet detection.

## Prerequisites

```r
# Seurat and demuxmix (R)
install.packages('Seurat')
BiocManager::install('demuxmix')
```

```bash
# hashsolo via scanpy/solo, and demuxEM via pegasus (Python)
pip install scanpy solo-sc pegasuspy demuxEM
```

## Quick Start

Tell your AI agent what you want to do:
- "Assign my hashed cells back to their sample of origin"
- "Call cross-sample doublets from my HTO counts"
- "Run HTODemux on my Seurat HTO assay"
- "Demultiplex my MULTI-seq lipid-tagged samples"
- "I have lots of Negatives and weak staining - which demux method should I use?"

## Example Prompts

### Hashtag demultiplexing
> "Normalize my HTO assay with CLR and run HTODemux to assign samples and flag doublets"
> "Demultiplex these MULTI-seq tags with MULTIseqDemux using automated thresholding"
> "Use hashsolo in scanpy to assign samples when I only have two hashtags"

### Robustness to ambient and bad staining
> "My HTODemux call gives a huge Negative pile - rescue it with a method that models ambient background"
> "Run demuxEM on my nucleus-hashing data using empty droplets to estimate background"
> "Use demuxmix with the detected-gene count to handle uneven staining across tags"

### Choosing a modality
> "I did not hash but my samples are different donors - should I use genetic demultiplexing instead?"
> "Combine my hashtag doublet calls with expression-based doublet detection"

## What the Agent Will Do

1. Confirm the hashing chemistry (antibody HTO, MULTI-seq lipid, CellPlex CMO) and that the HTO matrix barcodes match the GEX cells
2. Normalize the HTO counts with CLR, choosing the margin deliberately (margin=2 corrects per-tag capture bias)
3. Pick a caller from the decision table: HTODemux/MULTIseqDemux for clean data, hashsolo for few hashes, demuxEM/demuxmix when staining is marginal or ambient is high
4. Classify cells into singlet (with sample), cross-sample doublet, and Negative
5. Reconcile the cross-sample doublet rate against the expected loading doublet rate to sanity-check thresholds
6. Recommend expression-based doublet detection in addition to catch within-sample doublets hashing cannot see
7. Subset to confident singlets and pass them downstream for integration and clustering

## Tips

- **The cross-sample doublet rate calibrates the total** - with two samples within- and cross-sample doublets are equally frequent, with k samples within-sample doublets fall to about 1/k of all doublets, and a near-zero hashing doublet rate signals loose thresholds.
- **Negatives are not empties** - empty droplets are removed in preprocessing; Negatives are real cells whose true tag never cleared background.
- **Choose the CLR margin deliberately** - margin=2 normalizes each tag across cells and corrects capture-efficiency differences; do not blindly accept the default.
- **Switch methods when staining is weak** - demuxEM models background from empty droplets and demuxmix regresses on detected genes; both beat a fixed quantile.
- **Nucleus hashing is harder** - lower tag capture means more Negatives; demuxEM was built for it.
- **Hashing and genetics cannot see within-sample doublets** - always add expression-based doublet detection.
- **Genetic demux cannot split same-donor samples** - use hashtag demux when several samples share a genotype.
- **Low Negatives can still mean failed staining** - if one tag captures nearly all cells the assignment is meaningless; check the per-tag singlet distribution against the expected pooling.
- **Check every tag has positives** - a single near-zero tag means a failed antibody silently dropped or misassigned that sample even when global QC looks fine.
- **Unequal pooling destabilizes minority tags** - inspect per-tag ridge plots and consider demuxmix for a rare sample; 3+ tags high signals over-loading or ambient, not ordinary doublets.

## Related Skills

- single-cell/doublet-detection - Expression-based within-sample doublet calling that complements cross-sample hashing doublets
- single-cell/preprocessing - Filter empty droplets and QC the cells before and after demultiplexing
- single-cell/batch-integration - Integrate the demultiplexed per-sample data; covers genetic demultiplexing as an alternative
- single-cell/multimodal-integration - HTOs are an ADT-like modality; the CLR normalization here parallels CITE-seq ADT handling
- single-cell/clustering - Cluster the recovered singlets after sample assignment
