# Crosslink-Site Detection - Usage Guide

## Overview

Detect single-nucleotide crosslink (CL) sites in CLIP-seq data. The detection chemistry differs by CLIP variant: iCLIP/eCLIP read 5' end is CL site -1 (truncation); HITS-CLIP shows deletions at CL; PAR-CLIP shows T->C transitions at CL. PureCLIP is the most general HMM (works on all variants); CTK CITS for truncation; CTK CIMS for deletion or substitution; PARalyzer for PAR-CLIP T->C. Single-nt sites feed mCross motif registration, BEAPR allele-specific binding, and RBPNet deep-learning.

## Prerequisites

```bash
conda install -c bioconda pureclip samtools bedtools
# CTK
git clone https://github.com/chaolinzhanglab/ctk
# PARalyzer
git clone https://github.com/ohlerlab/PARalyzer
# wavClusteR
BiocManager::install('wavClusteR')
```

## Quick Start

Tell your AI agent:
- "Run PureCLIP with SMInput on my eCLIP for single-nt CL maps"
- "CTK CITS truncation sites for iCLIP, then mCross"
- "CTK CIMS deletion for HITS-CLIP - need BWA-aln upstream"
- "PARalyzer T->C clusters with the published parameters for PUM2"
- "CTK CIMS substitution T->C for PAR-CLIP at p=0.001"
- "Reconcile PureCLIP and CTK CITS for the same iCLIP dataset"
- "Restrict PureCLIP to expressed transcripts via -iv flag"

## Example Prompts

### iCLIP / eCLIP

> "PureCLIP with -dm 8 (single-nt distance to merge), restrict to expressed tx with -iv"

> "CTK CITS - tag2cluster.pl -cs5 5 -m 1 - the truncation-based alternative"

### HITS-CLIP

> "Re-align with BWA-aln for deletion tolerance; CTK CIMS -type del at p=0.01"

### PAR-CLIP

> "PARalyzer with Hafner 2010 parameters for HuR"

> "CTK CIMS substitution T->C for single-nt; alternative to PARalyzer clusters"

### Downstream

> "Feed PureCLIP sites to mCross for CL-position-registered motif"

> "Intersect CL sites with het SNPs for BEAPR allele-specific binding"

### Diagnostics

> "Why is my CTK CIMS empty for HITS-CLIP? Aligner not deletion-tolerant"

> "Why is PureCLIP convergence failing? Sparse coverage; add -iv expressed.bed"

> "Off-by-one motif position - strand handling issue"

## What the Agent Will Do

1. Identify the CLIP variant and pick the matching detection method
2. For iCLIP/eCLIP: PureCLIP (HMM) or CTK CITS (truncation)
3. For HITS-CLIP: CTK CIMS deletion mode; ensure aligner is deletion-tolerant (BWA-aln or relaxed STAR)
4. For PAR-CLIP: PARalyzer cluster-level or CTK CIMS substitution single-nt
5. Cross-validate: 2+ tools agreeing within 5 nt = high confidence
6. Downstream: mCross PWM, BEAPR ASB, RBPNet training data
7. Flag failure modes: convergence, deletion-tolerance, strand off-by-one

## Tips

- **PureCLIP works on all CLIP variants.** Most general HMM.
- **Match method to chemistry.** Truncation for iCLIP/eCLIP; deletion for HITS-CLIP; T->C for PAR-CLIP.
- **HITS-CLIP needs deletion-tolerant aligner.** BWA-aln is the Yeo convention.
- **PAR-CLIP needs raised STAR mismatch ceiling.** 0.07 not 0.04.
- **PureCLIP is focal.** F1 ~0.2 on bulk RBPs; use for single-nt sites, not broad zones.
- **Restrict PureCLIP scope.** -iv expressed.bed to avoid genome-wide convergence issues.
- **CL site = RT stop position - 1 for truncation.** Tools should handle this; verify.
- **Strand handling matters.** Plus-strand RNA: cDNA 5' is downstream of CL; minus-strand: upstream.
- **Cross-validate PureCLIP and CTK.** Same mCross PWM = high confidence.

## Related Skills

- clip-seq/clip-preprocessing - 5' preservation critical
- clip-seq/clip-alignment - End-to-end + deletion-tolerance for HITS-CLIP
- clip-seq/clip-peak-calling - Peaks vs CL sites complementary
- clip-seq/clip-motif-analysis - mCross consumes CL sites
- clip-seq/clip-deep-learning - RBPNet trained on CL distributions
- clip-seq/m6a-clip - miCLIP2 own detection
- clip-seq/stamp-antibody-free - Editing-based, not CL-based
