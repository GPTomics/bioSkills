# STAMP / Antibody-Free RBP Profiling - Usage Guide

## Overview

Profile RBP-RNA targets without UV crosslinking or immunoprecipitation. Express an APOBEC1-RBP fusion (STAMP / scSTAMP) or ADAR-RBP fusion (TRIBE / HyperTRIBE); the deaminase edits RNA near the binding site producing detectable C->U or A->I signatures in standard RNA-seq. Compatible with single-cell readout. Trade-off: editing is offset 0-50 nt from binding; resolution is approximate. Use STAMP/TRIBE for hypothesis generation, single-cell profiling, or when antibody is unavailable. Cross-validate with CLIP for high-resolution.

## Prerequisites

```bash
conda install -c bioconda samtools bedtools
pip install pysam scanpy anndata
# Bullseye: github.com/mekoulnik/Bullseye
# SAILOR: github.com/YeoLab/sailor
# JACUSA2: github.com/dieterich-lab/JACUSA2
# REDItools2: github.com/BioinfoUNIBA/REDItools2
```

## Quick Start

Tell your AI agent:
- "STAMP with APOBEC1-RBP fusion; analyze with Bullseye, subtract APOBEC1-only control"
- "scSTAMP single-cell from 10x library; per-cell editing rates"
- "TRIBE with ADAR; filter out ALU repeats (baseline editing)"
- "DART-seq for YTHDF1 m6A reader profiling"
- "Why is my edit count enormous? Forgot the deaminase-only control"
- "STAMP vs eCLIP concordance for top RBP targets"
- "Titrate APOBEC1-RBP expression; saturated editing degrades specificity"

## Example Prompts

### STAMP

> "Bullseye on STAMP BAM with APOBEC1-only control; output C->U edits at edit rate >= 0.1"

> "Verify fusion edits / APOBEC1-only edits > 3 for clean signal"

### scSTAMP

> "10x cellranger then Bullseye per-cell; pseudobulk by cluster"

> "Smart-seq2 alternative if 10x coverage too sparse"

### TRIBE

> "REDItools2 on TRIBE BAM; subtract ADAR-only baseline; filter ALU"

> "HyperTRIBE E488Q for higher edit rate"

### DART-seq

> "APOBEC1-YTH for m6A reader profiling; expect edits 0-50nt offset from m6A"

### Diagnostics

> "Saturated editing on every gene - reduce fusion expression with inducible promoter"

> "TRIBE all edits at ALUs - subtract ADAR-only; filter ALU repeats"

> "Edit count enormous - need deaminase-only control subtraction"

> "scSTAMP per-cell sparse - aggregate by cluster"

### Validation

> "Run STAMP + eCLIP on same RBP; concordance for top 100 targets"

## What the Agent Will Do

1. Identify experimental setup: STAMP/TRIBE/DART; bulk/single-cell
2. Standard RNA-seq pipeline (no CLIP-specific preprocessing); strand-aware
3. Edit-site detection with Bullseye (STAMP/DART), SAILOR (Yeo), or JACUSA2 (general)
4. MANDATORY: subtract deaminase-only control (APOBEC1-only for STAMP; ADAR-only for TRIBE)
5. Filter at edit rate >= 0.1 and coverage >= 10
6. For scSTAMP: per-cell editing matrix; pseudobulk by cluster if sparse
7. Cross-validate with CLIP/eCLIP for high-resolution localization
8. Flag failures: saturation, off-target, ALU baseline (TRIBE), DNA editing (STAMP)

## Tips

- **Deaminase-only control is mandatory.** Without it, edits look like signal.
- **Titrate fusion expression.** Saturated APOBEC1 edits everything.
- **STAMP for ssRNA; APOBEC1 acts on ssRNA.** For structured RBPs may miss sites.
- **TRIBE ALU baseline dominates.** Filter ALU repeats unless RBP is repeat-binder.
- **DART editing offset 0-50 nt from m6A.** Not single-base resolution.
- **Resolution is ~50 nt approximate.** CLIP is finer; use both.
- **scSTAMP needs cluster pseudobulk.** Per-cell sparse; pool similar cells.
- **Strand-specific filtering.** C->U on +sense; G->A on -sense (mate-pair convention).
- **Fusion vs control edit ratio > 3.** Below this = saturation or non-specific.
- **STAMP + eCLIP cross-validation = gold standard.** Best figure.

## Related Skills

- clip-seq/m6a-clip - DART-seq part of m6A toolkit
- clip-seq/clip-deep-learning - Computational target prediction
- clip-seq/ago-clip-mirna-targets - AGO-CLIP comparison
- single-cell/preprocessing - scSTAMP processing
- single-cell/clustering - Per-cluster pseudobulk
- single-cell/markers-annotation - Cell-type targets
- methylation-analysis/methylation-calling - Editing analogue
- epitranscriptomics/m6anet-analysis - Nanopore m6A
