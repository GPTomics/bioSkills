# Subclonal Copy Number and Tumor Evolution Usage Guide

## Overview

A tumor is a mixture of cell populations. When a copy-number change is present in only some cancer cells, bulk sequencing averages it into a non-integer state, and a long non-integer segment is a subclonal copy-number alteration, not noise. Resolving subclonal CN reveals the tumor's clonal architecture; calling whole-genome doubling explicitly is a prerequisite because it rescales every copy number. This skill covers Battenberg (phased clonal + subclonal CN), TITAN (HMM mixture of cell populations), whole-genome-doubling detection and timing, and MEDICC2 copy-number phylogenies.

## Prerequisites

```bash
R -e "remotes::install_github('Wedge-lab/battenberg')"   # plus 1000G impute reference
R -e "BiocManager::install('TitanCNA')"
conda install -c bioconda -c conda-forge medicc2         # pip alone misses the OpenFST dependency
```

Inputs: tumour and matched-normal WGS (or WES for TITAN); allele-specific data (logR and BAF at heterozygous SNPs); a 1000 Genomes phasing/impute reference for Battenberg. For evolution analysis, multi-region or multi-sample data.

## Quick Start

Tell the AI agent what to do:
- "Call subclonal copy number for this tumour-normal WGS pair with Battenberg"
- "Run TITAN and estimate the cellular prevalence of clonal clusters"
- "Determine whether this tumour has undergone whole-genome doubling"
- "Time whole-genome doubling relative to point mutations"
- "Build a copy-number phylogeny across my multi-region samples with MEDICC2"

## Example Prompts

### Subclonal calling

> "Run Battenberg on this tumour-normal WGS pair and report which segments are subclonal, with their major/minor states and cell fractions."

> "Run TITAN, sweep the number of clonal clusters from one to five, and select the cluster count by model fit."

### Whole-genome doubling

> "Decide whether this tumour is whole-genome-doubled using absolute allele-specific copy number, and explain why a depth-only profile cannot answer this."

> "Time whole-genome doubling for this tumour using mutation copy number."

### Evolution and reconciliation

> "Build a whole-genome-doubling-aware copy-number phylogeny across my five regions with MEDICC2."

> "Battenberg and TITAN disagree on whether a segment is subclonal. Walk through reconciling them."

## What the Agent Will Do

1. Confirm depth and purity are adequate for subclonal resolution
2. Establish whole-genome-doubling status from absolute allele-specific copy number
3. Run Battenberg or TITAN to fit clonal and subclonal copy-number states
4. Report cancer cell fractions; flag mirrored subclonal allelic imbalance
5. For evolution questions, build a copy-number phylogeny with MEDICC2
6. Require concordant evidence before declaring a subclone; flag single-segment claims

## Tips

- A long non-integer segment is a subclonal CNA, not noise; use Battenberg or TITAN, not a clonal-only caller.
- Call whole-genome doubling explicitly before interpreting any copy number; missing it halves every copy number and mis-times every mutation.
- WGD calling needs absolute allele-specific copy number; depth alone cannot distinguish a doubled genome from a non-doubled one.
- Battenberg resolves subclones to ~3% of cells only at adequate WGS depth and purity; treat low-depth subclonal calls as exploratory.
- Mirrored subclonal allelic imbalance can look balanced in bulk BAF; phasing or single-cell data is needed to detect it.
- A single subclonal segment is a hypothesis; require multiple concordant segments at a consistent cell fraction before declaring a subclone.
- Single-region sampling misses spatial subclones; use multi-region data for evolution.

## Related Skills

- copy-number/allele-specific-copy-number - Clonal allele-specific CN, purity, ploidy
- copy-number/copy-ratio-segmentation - Segmentation feeding subclonal callers
- copy-number/hrd-scoring - Whole-genome-doubling correction for the LST scar
- copy-number/recurrent-cnv - Copy-number signatures including WGD and chromothripsis
- copy-number/cnv-visualization - Visualizing subclonal segments and BAF
- variant-calling/vcf-basics - SNV calls for cancer cell fraction and WGD timing
