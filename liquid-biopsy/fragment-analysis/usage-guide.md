# Fragment Analysis - Usage Guide

## Overview
Extracts cfDNA fragmentomics features (DELFI genome-wide short/long ratios, WPS nucleosome positioning, Griffin GC-corrected accessibility, end-motifs/MDS, OCF) from plasma WGS for cancer detection and tissue-of-origin. The features are four views of one nucleosome-footprint object, so GC correction and protocol-matching are non-negotiable.

## Prerequisites
```bash
pip install finaletoolkit pysam numpy pandas matplotlib
# Griffin: clone github.com/adoebley/Griffin and run its Snakemake pipeline (not pip-installable)
# DELFI is a methodology + company, not a package -- compute its features via finaletoolkit
```

## Quick Start
Tell your AI agent what you want to do:
- "Compute a genome-wide DELFI short/long fragment ratio profile for cancer detection"
- "Profile nucleosome accessibility around TF binding sites for tissue of origin"
- "Calculate end-motif frequencies and the Motif Diversity Score for my sample"
- "GC-correct my fragmentomic features before comparing samples across batches"
- "Apply in-silico size selection to enrich tumor fraction at low ctDNA"

## Example Prompts

### Detection (genome-wide ratios)
> "Run a GC-corrected DELFI score across 5 Mb bins for my plasma BAM and flag bins deviating from a healthy reference."

> "Compute a custom short(100-150 bp)/long(151-220 bp) ratio profile and explain why it is not comparable across sequencing batches."

### Tissue of origin (nucleosome profiling)
> "Set up the Griffin Snakemake pipeline to profile nucleosome accessibility around a TF site list for subtype calling."

> "Compute WPS over a promoter region and identify nucleosome positions and TF footprints."

### Nuclease signal (end motifs)
> "Extract 4-mer end-motif frequencies and compute the Motif Diversity Score, and interpret a raised MDS in the context of DNASE1L3 biology."

### Low tumor fraction
> "My sample has tumor fraction below 0.03 -- recommend a feature family and whether to apply 90-150 bp size selection."

## What the Agent Will Do
1. Choose a feature family from the question (detection vs tissue-of-origin vs nuclease signal)
2. Extract fragments from BAM/CRAM or a tabix-indexed `.frag.gz` file
3. Apply GC correction (FinaleToolkit `delfi`, Griffin, or a standalone corrector)
4. Compute the chosen features (DELFI ratios, WPS, Griffin profiles, end-motifs/MDS)
5. Compare to a co-processed healthy reference and flag protocol/batch confounders

## Tips
- GC correction is the point - an uncorrected DELFI plot is a GC plot, not a tumor signal
- Never mix ssDNA and dsDNA libraries; the prep sets the size floor and moves every feature
- DELFI ratios are entangled with copy-number alterations - treat them as a hybrid, not pure fragmentation
- Griffin is the robust choice at low tumor fraction because its GC correction removes the dominant technical signal
- The 10.4 bp periodicity below 167 bp is the cleanest check that the data are genuine nucleosome footprints
- Size selection (90-150 bp) trades tumor fraction for depth; skip it when already depth-limited
- The FinaleToolkit filter subcommand is `filter-file`, not `filter-bam`

## Related Skills
- cfdna-preprocessing - library prep determines which fragments (and features) are recoverable
- tumor-fraction-estimation - fragmentomics enables signal below the CNA-based TF floor
- methylation-based-detection - orthogonal genome-wide cfDNA signal
- atac-seq/nucleosome-positioning - shared nucleosome-footprint biology
