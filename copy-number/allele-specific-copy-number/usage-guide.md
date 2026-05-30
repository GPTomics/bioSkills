# Allele-Specific Copy Number Usage Guide

## Overview

Allele-specific copy number callers jointly model read depth (logR) and B-allele frequency (BAF) to fit tumor purity, ploidy, and integer major/minor copy number per segment. Depth alone gives only relative copy ratio; adding BAF gives absolute, allele-specific copy number and resolves loss of heterozygosity, copy-neutral LOH, and whole-genome doubling. The core difficulty is the purity-ploidy identifiability problem: many (purity, ploidy) pairs explain the same depth profile, so the likelihood surface is multimodal and every fit must be checked against its diagnostic plot. This skill covers ASCAT, Sequenza, FACETS, PURPLE, and tumor-only PureCN.

## Prerequisites

```bash
# FACETS
conda install -c bioconda snp-pileup
R -e "install.packages('facets')"   # or remotes::install_github('mskcc/facets')
# Sequenza
conda install -c bioconda sequenza-utils
R -e "remotes::install_github('ShixiangWang/copynumber'); install.packages('sequenza')"
# ASCAT and PureCN
R -e "remotes::install_github('VanLoo-lab/ascat/ASCAT'); BiocManager::install('PureCN')"
```

Inputs: tumor (and, except for PureCN, matched-normal) BAMs; the reference FASTA; a common-SNP VCF (FACETS) or GC wiggle (Sequenza); GC/replication-timing reference files (ASCAT). PureCN needs a normal database built from process-matched normals.

## Quick Start

Tell the AI agent what to do:
- "Run FACETS on this tumor-normal panel pair and report purity, ploidy, and LOH segments"
- "Estimate tumor cellularity and ploidy from this WES pair with Sequenza"
- "Run ASCAT on this tumor-normal WGS pair and check the sunrise plot"
- "Call allele-specific copy number for this tumor-only panel sample with PureCN"
- "Reconcile a ploidy disagreement between FACETS and ASCAT on the same tumor"

## Example Prompts

### Calling

> "Run the two-pass FACETS workflow on this panel tumor-normal pair: a coarse purity run then a dipLogR-seeded sensitivity run, and report per-segment total and minor CN."

> "Estimate purity and ploidy with Sequenza for this exome pair and review the alternative solutions before accepting the fit."

### Tumor-only and WGS

> "I have no matched normal for this panel tumor. Run PureCN against a normal database and flag the fit for manual curation if purity is low."

> "Run ASCAT on this WGS pair, apply GC and replication-timing correction, and tell me whether the sunrise plot indicates an ambiguous ploidy."

### Interpretation and reconciliation

> "FACETS and ASCAT disagree on ploidy by about a factor of two for this tumor. Diagnose whether this is an integer-multiple ploidy flip and which solution to trust."

> "This near-diploid tumor returns 100% purity from ASCAT. Explain why and what to do."

## What the Agent Will Do

1. Select the caller from the data type (WGS, WES, panel, tumor-only, SNP array)
2. Prepare inputs (snp-pileup, seqz, or logR/BAF tracks)
3. Run the joint logR+BAF fit, using the two-pass workflow for FACETS
4. Inspect the fit diagnostic (sunrise plot, cellularity/ploidy contour, dipLogR)
5. Report purity, ploidy, and integer allele-specific CN; flag LOH (minor CN = 0)
6. Cross-check ploidy against odd/even CN fraction and clonal-SNV multiplicity
7. Reconcile against other callers and flag indeterminate fits rather than over-reporting

## Tips

- Never trust a bare purity/ploidy number; inspect the fit diagnostic every time.
- The likelihood surface is multimodal; a factor-of-two ploidy disagreement between callers is usually an integer-multiple flip, resolved by SNV multiplicity.
- Allele-specific calling fails below ~40% tumor purity and is impossible below ~20%.
- A near-diploid genome cannot anchor purity; ASCAT drifting to ~100% means indeterminate.
- FACETS cval ~150-300 for panels/WES, ~25-100 for WGS; too low causes hyperfragmentation.
- Sequenza needs a maintained `copynumber` fork (removed from Bioconductor 3.18+).
- LOH is minor copy number = 0; copy-neutral LOH has total CN >= 2 with minor CN = 0.
- PureCN tumor-only fits are weaker than matched-normal; use `createCurationFile` review.

## Related Skills

- copy-number/copy-ratio-segmentation - logR normalization and segmentation theory
- copy-number/subclonal-copy-number - Battenberg/TITAN subclonal CN and WGD
- copy-number/hrd-scoring - LOH/LST/TAI scars from allele-specific output
- copy-number/cnvkit-analysis - Relative depth-only calling
- copy-number/gatk-cnv - GATK somatic CNV (relative; no purity/ploidy)
- variant-calling/vcf-basics - SNV VCFs for BAF and clonal-mutation cross-checks
