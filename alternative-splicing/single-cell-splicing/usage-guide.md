# Single-Cell Splicing - Usage Guide

## Overview
Analyze alternative splicing at single-cell resolution. The first decision is library chemistry - 10X 3' is fundamentally limited because RT primes from poly(A) and >70% of reads land in the 3' UTR. Plate-based full-length methods (Smart-seq3, FLASH-seq, VASA-seq, STORM-seq) and single-cell long-read (MAS-Iso-seq, scISOr-Seq2) are the chemistries that give per-cell isoform structure. Tools include MARVEL (R, Smart-seq), BRIE2 (Bayesian PSI with regulatory features), scQuint (junction-cluster, plate-based; not for 10X), SpliZ (annotation-free Z-score), Psix (graph-smoothness regulated AS), and Sierra (alternative polyadenylation, often confused with AS).

## Prerequisites
```bash
# Python
pip install brie scanpy anndata scipy

# R - both from GitHub (Sierra not on Bioconductor; MARVEL archived from CRAN 2025-10)
# devtools::install_github('VCCRI/Sierra')
remotes::install_github('wenweixiong/MARVEL')

# scQuint, SpliZ, Psix from GitHub
pip install git+https://github.com/songlab-cal/scquint
pip install git+https://github.com/salzman-lab/SpliZ
pip install git+https://github.com/lareaulab/psix
```

## Quick Start
Tell your AI agent what you want to do:
- "Will my 10X 3' data support splicing analysis?"
- "Analyze splicing in Smart-seq2 plate-based scRNA-seq with MARVEL"
- "Find cell-type-specific splicing without an event database using SpliZ"
- "Detect alternative polyadenylation in 10X data with Sierra"
- "Pseudobulk single cells by cell type and run leafcutter for differential splicing"

## Example Prompts

### Chemistry Audit
> "I have 10X Chromium 3' v3 data; can I do splicing analysis? What can I do instead?"

### Plate-Based Workflow
> "Run MARVEL on Smart-seq2 BAMs to quantify SE/A5SS/A3SS PSI, classify modality, and test differential splicing between cell types."

> "Use BRIE2 with sequence-feature prior to estimate per-cell PSI for cassette exons in low-coverage Smart-seq3 data."

### Discovery
> "Run SpliZ to find genes with cell-state-associated splicing without using an event database."

### Trajectory
> "Use Psix to detect regulated AS along a developmental pseudotime trajectory, robust to dropout."

### APA (10X 3')
> "Use Sierra to peak-call 3' ends and detect alternative polyadenylation across cell types."

### Pseudobulk
> "Aggregate cells by cluster and run leafcutter or rMATS on pseudobulk junction counts for differential splicing between cell types."

### Long-Read Single-Cell
> "I have MAS-Iso-seq + 10X 5' data; demultiplex barcodes with skera/lima, then run FLAMES for full-length isoform quantification per cell."

## What the Agent Will Do
1. Audit the library chemistry and recommend tool subset based on what's recoverable
2. For full-length plate data: quantify per-cell PSI; classify modality; test cell-state association
3. For 10X 3' data: redirect to APA (Sierra). Note: scQuint authors recommend against use on 10X data (3' bias confounds junction-cluster modeling)
4. For long-read single-cell: demultiplex, align with minimap2 splice:hq, quantify with FLAMES/Bambu/IsoQuant
5. Recommend pseudobulk for between-cell-type comparisons; per-cell only for within-cluster heterogeneity
6. Avoid PSI imputation - use graph-smoothness (Psix) or hierarchical Bayes (BRIE2) instead

## Tips
- 10X 3' chemistry CANNOT support transcriptome-wide splicing analysis; <0.1 junction read per cell per AS event vs the 5-10 needed
- 10X 5' shifts capture region but does not solve the problem
- Smart-seq2/3, FLASH-seq, VASA-seq, STORM-seq are the full-length plate options
- MAS-Iso-seq (PacBio Kinnex + 10X 5') is the practical SOTA for high-throughput single-cell isoforms in 2024-2026
- Pseudobulk first for differential testing; descend to per-cell only for full-length data
- Don't impute PSI - destroys heterogeneity; use graph-based methods
- Microexons (3-27nt) require either long-read or short-read aligners with low overhang (uLTRA, deSALT, or STAR with `--alignSJoverhangMin 6 --alignSJDBoverhangMin 1` plus a strict mismatch filter; the typical AS-pipeline value of 8/3 is too strict for microexons)
- Sierra is APA, not splicing; clearly distinguish in interpretation
- snRNA-seq IR signal is mostly nuclear-retained transcripts, not splicing dysregulation

## Related Skills

- single-cell/preprocessing - QC and normalization
- single-cell/clustering - Cell type annotation prerequisite
- single-cell/data-io - h5ad / Seurat I/O
- splicing-quantification - Bulk RNA-seq comparison
- long-read-splicing - Full-isoform analysis from MAS-Iso-seq
