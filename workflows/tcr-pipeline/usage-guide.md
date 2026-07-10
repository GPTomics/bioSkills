# TCR/BCR Repertoire Pipeline - Usage Guide

## Overview

This workflow takes immune-repertoire sequencing from FASTQ to clonotypes and downstream metrics, but it is a router rather than a fixed line. A repertoire measurement is a depth- and chemistry-confounded sample of an unevenly-expanded clonal population, so the correct pipeline depends on two forks: bulk versus single-cell, and TCR versus BCR. Bulk libraries give deep but unpaired chains and route to MiXCR then depth-normalized `mixcr postanalysis` (or immunarch) diversity and overlap; single-cell 10x data gives native chain pairing and routes to scirpy for gene-expression integration. TCR uses exact CDR3-nucleotide + V + J clonotypes, but BCR undergoes somatic hypermutation, so exact clonotypes are wrong and BCR must route to Immcantation for data-derived clonal clustering, germline reconstruction, SHM, and lineage trees. Two rules gate every valid result: MiXCR 4.x needs an activated license and a chemistry-matched preset, and diversity or overlap must be compared only after downsampling samples to a common depth.

## Prerequisites

```bash
# MiXCR 4.x (requires an activated license; academic is free at platforma.bio/getlicense)
conda install -c milaboratory mixcr
mixcr activate-license                         # or export MI_LICENSE_FILE=/path/mi.license

# VDJtools (Java): optional, for pre-4.x MiXCR or non-MiXCR inputs; run RInstall once for plotting deps.
# It cannot parse raw MiXCR 4.x exportClones output -- prefer `mixcr postanalysis` there.
conda install -c bioconda vdjtools

# Immcantation (BCR clonal clustering, SHM, lineages) -- prefer the Docker suite image
#   docker pull immcantation/suite:4.x
# or R: BiocManager / install.packages(c('shazam','scoper','alakazam','dowser','tigger'))

# scirpy (single-cell VDJ + gene expression)
pip install scirpy mudata scanpy
```

## Quick Start

Tell your AI agent what you want to do:
- "Run the bulk TCR pipeline from FASTQ to diversity"
- "Analyze my BCR repertoire including somatic hypermutation and lineage trees"
- "Process my 10x single-cell VDJ data and link clonotypes to cell state"
- "Compare diversity across my samples at equal depth"

## Example Prompts

### Bulk TCR

> "Assemble clonotypes from my paired-end TCR amplicon FASTQ, then compare diversity across the cohort at equal depth"

> "Compute repertoire overlap between my pre- and post-treatment TCR samples with a depth-robust metric"

### BCR

> "Analyze my IGH repertoire: cluster clones, quantify somatic hypermutation, and build lineage trees"

> "My BCR diversity looks absurdly high with no lineages -- what am I doing wrong?"

### Single-cell

> "Run MiXCR on my 10x VDJ FASTQ, then integrate clonotypes with my scRNA-seq clusters in scirpy"

> "Filter doublets and measure clonal expansion by cell state in my single-cell TCR data"

## What the Agent Will Do

1. Confirm the two forks (bulk vs single-cell, TCR vs BCR) and pick the MiXCR 4.x preset by chemistry.
2. Activate the MiXCR license and run `mixcr analyze <preset>` to assemble clonotypes.
3. Run MiXCR QC (alignment rate, chain usage) and report the correct abundance denominator (reads, UMIs, or cells).
4. Export AIRR TSV (BCR / single-cell); for bulk metrics stay in MiXCR `postanalysis` (VDJtools cannot parse raw MiXCR 4.x exportClones output).
5. Bulk TCR: downsample all samples to a common depth, then compute diversity and overlap.
6. BCR: derive the clonal threshold from distToNearest, cluster clones, reconstruct germlines, quantify SHM, build lineages.
7. Single-cell: run chain QC, define clonotypes, and integrate with the gene-expression AnnData.
8. Generate figures and, optionally, annotate antigen specificity as a hypothesis.

## Tips

- MiXCR 4.x will not run without an activated license; set `MI_LICENSE_FILE` on HPC/Docker and whitelist the phone-home IPs on firewalled clusters.
- The preset is the analysis. The wrong preset (RNA vs DNA, floating vs rigid boundary, missing tag pattern) does not error -- it silently truncates CDR3 and mis-calls V. Audit with `mixcr exportPreset`.
- Never compare raw diversity, clonality, or Jaccard overlap across samples of unequal depth. DownSample first, or read rarefaction curves at a common depth.
- BCR is the common trap: somatic hypermutation makes exact CDR3 clonotypes wrong. Route BCR to Immcantation clonal clustering before any diversity, SHM, or lineage step.
- Report the right unit: `uniqueMoleculeCount` for UMI libraries, cells for single-cell, reads only for non-UMI bulk.
- Hold the clonotype match key (nt vs aa, +/-V, +/-J) constant across a whole study; aa-level matching inflates apparent sharing.
- A specificity-database hit is a hypothesis, not a label; public clonotypes are mostly high-Pgen convergent sequences, not antigen selection.

## Related Skills

- tcr-bcr-analysis/mixcr-analysis - V(D)J alignment and clonotype assembly
- tcr-bcr-analysis/vdjtools-analysis - Depth-normalized diversity and overlap
- tcr-bcr-analysis/immcantation-analysis - BCR clonal clustering, SHM and lineages
- tcr-bcr-analysis/scirpy-analysis - Single-cell VDJ + gene-expression integration
- tcr-bcr-analysis/repertoire-visualization - Figures for the pipeline outputs
- tcr-bcr-analysis/specificity-annotation - Optional antigen-specificity annotation
