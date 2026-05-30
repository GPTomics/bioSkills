# GATK CNV Workflows Usage Guide

## Overview

GATK provides two distinct, best-practices copy-number workflows: a **somatic CNV** pipeline for tumor samples (read counts, tangent-normalized denoising, joint copy-ratio/allele-fraction segmentation) and a **germline GATK-gCNV** pipeline for constitutional CNV discovery across a cohort. They share almost no tools. The somatic workflow produces *relative* copy-ratio segments and +/-/0 calls, not integer allele-specific copy number, purity, or ploidy. The germline workflow produces integer CN genotype VCFs but requires a large, technically matched cohort and QS-based filtering.

## Prerequisites

```bash
conda install -c bioconda gatk4
# gCNV additionally needs the GATK Python (gatkcondaenv) environment installed
```

Inputs: reference FASTA + `.fai` + `.dict`; an interval list (exome targets or WGS bins); a common biallelic-SNP list for allelic counts (somatic); a contig-ploidy priors table (germline). The somatic workflow needs a panel of normals; the germline cohort workflow needs >= ~100 process-matched samples.

## Quick Start

Tell the AI agent what to do:
- "Run the GATK somatic CNV workflow on this tumor-normal WES pair"
- "Build a GATK panel of normals from my tumor-free control BAMs"
- "Call rare germline CNVs from my 300-sample exome cohort with GATK-gCNV"
- "Filter my raw gCNV calls to high precision for a de novo CNV analysis"
- "Explain why a known amplification disappeared after DenoiseReadCounts"

## Example Prompts

### Somatic CNV

> "Run the full GATK somatic CNV pipeline on this tumor-normal WGS pair, including AnnotateIntervals and matched-normal allelic counts, and explain what the CallCopyRatioSegments output can and cannot tell me."

> "My GATK panel of normals has only 8 samples. Assess whether tangent normalization is reliable and recommend how to improve it."

### Germline gCNV

> "Set up GATK-gCNV cohort mode for my 250-exome neurodevelopmental cohort, including DetermineGermlineContigPloidy and FilterIntervals."

> "Apply QS and sample-level filters to my gCNV output so it is suitable for a de novo CNV burden analysis."

### Workflow selection and reconciliation

> "I need integer allele-specific copy number and tumor purity. Explain why the GATK somatic CNV workflow cannot give me that and what to use instead."

> "GATK and CNVkit disagree on a focal deletion. Walk through reconciling them."

## What the Agent Will Do

1. Determine whether the task needs the somatic or the germline gCNV workflow
2. Preprocess and annotate intervals (GC content for bias correction)
3. Collect read counts; build a panel of normals (somatic) or run contig-ploidy (germline)
4. Denoise via tangent normalization (somatic) or the Bayesian gCNV model (germline)
5. Segment and call: joint copy-ratio/allele-fraction (somatic) or CN genotypes (germline)
6. Apply QS and sample-level filtering for germline calls
7. Flag when escalation to an allele-specific caller is required and reconcile with other callers

## Tips

- GATK somatic CNV is a relative caller; it never outputs purity, ploidy, or integer allele-specific CN. Use ASCAT/Sequenza/FACETS/PureCN for those.
- Never put tumors in a somatic panel of normals; tangent normalization will subtract any copy-number pattern shared across the PoN.
- Run `AnnotateIntervals` and pass `--annotated-intervals` for explicit GC-bias correction.
- gCNV raw output is ~95% recall but only ~22% precision; QS > 1000 reaches ~96% precision. Always filter before association or de novo analysis.
- gCNV cohort mode needs >= ~100 technically matched samples; split mixed capture kits.
- Case-mode gCNV must reuse the exact cohort interval list and scatter count.
- `FilterIntervals` can drop segdup-flanked loci where genomic-disorder CNVs live; check genes of interest survive filtering.

## Related Skills

- copy-number/allele-specific-copy-number - Integer ASCN, purity, ploidy
- copy-number/copy-ratio-segmentation - Segmentation and depth-normalization theory
- copy-number/cnvkit-analysis - Read-depth CNV calling for panels and exomes
- copy-number/germline-cnv-interpretation - ACMG/ClinGen classification of germline CNVs
- copy-number/cnv-visualization - Plotting denoised ratios and modeled segments
- variant-calling/gatk-variant-calling - GATK SNV/indel pipeline
