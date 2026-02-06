# GRN Pipeline - Usage Guide

## Overview
Complete gene regulatory network inference workflow from processed single-cell data to regulon discovery and perturbation simulation. Supports RNA-only analysis with pySCENIC (GRNBoost2 + RcisTarget + AUCell) and multiome analysis with SCENIC+ for enhancer-driven GRNs. Includes CellOracle for in silico perturbation simulation.

## Prerequisites
```bash
# pySCENIC (RNA-only path)
pip install pyscenic arboreto ctxcore

# SCENIC+ (multiome path)
pip install scenicplus pycisTopic pycistarget

# CellOracle (perturbation simulation)
pip install celloracle

# Common dependencies
pip install scanpy pandas numpy scipy matplotlib

# Download cisTarget databases (once, ~10 GB each)
# Human: hg38_500bp_up_100bp_down and hg38_10kbp_up_10kbp_down
# Mouse: mm10 equivalents
# From: https://resources.aertslab.org/cistarget/databases/
```

**Input data:**
- Processed AnnData (QC'd, normalized, clustered) with `.raw` containing counts
- TF list for species (from cisTarget resources)
- For SCENIC+: paired scRNA-seq and scATAC-seq (multiome)

## Quick Start
Tell your AI agent what you want to do:
- "Build gene regulatory networks from my single-cell RNA-seq data"
- "Run pySCENIC to find regulons in my scRNA-seq dataset"
- "Infer enhancer-driven GRNs from my multiome data with SCENIC+"
- "Simulate a TF knockout using CellOracle"

## Example Prompts

### RNA-Only GRN
> "I have a processed scRNA-seq AnnData with cell type annotations. Run pySCENIC to discover regulons and score their activity."

> "Find the key transcription factors driving cell identity in my single-cell data."

### Multiome GRN
> "I have paired RNA and ATAC from a 10x multiome experiment. Build enhancer-driven regulatory networks with SCENIC+."

> "Use SCENIC+ to identify which enhancers are driving TF regulons in my data."

### Perturbation Simulation
> "Simulate what happens if I knock out MYC in my scRNA-seq data using CellOracle."

> "Predict cell state transitions after perturbing key transcription factors."

## What the Agent Will Do
1. Load processed AnnData and extract expression matrix
2. Run GRN inference (GRNBoost2 for RNA-only, cisTopic + pycistarget for multiome)
3. Prune TF-target links using motif enrichment (RcisTarget)
4. Score regulon activity per cell (AUCell)
5. Validate: check regulon count (50-500), known TFs present, cell type separation
6. Optionally simulate TF perturbations with CellOracle
7. Export regulon activity matrix and annotated AnnData

## Tips
- Start with pySCENIC for RNA-only data; upgrade to SCENIC+ only with multiome
- Use raw counts for GRNBoost2 (not normalized values)
- Download both 500bp and 10kbp cisTarget databases for comprehensive motif scanning
- Subsample to ~50k cells if memory is an issue during GRNBoost2
- Expect 50-500 regulons; outside this range suggests parameter tuning needed
- Verify known lineage TFs appear in your regulon list as a sanity check
- CellOracle perturbation results are predictions; validate experimentally

## Related Skills
- gene-regulatory-networks/scenic-regulons - pySCENIC implementation details
- gene-regulatory-networks/multiomics-grn - SCENIC+ enhancer-driven GRNs
- gene-regulatory-networks/perturbation-simulation - CellOracle details
- single-cell/clustering - Upstream cell type annotation
- single-cell/preprocessing - QC and normalization before GRN inference
