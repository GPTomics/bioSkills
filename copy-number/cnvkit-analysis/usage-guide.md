# CNVkit Analysis Usage Guide

## Overview

CNVkit detects copy number variants from targeted (panel, exome) and whole-genome sequencing using read depth. It combines on-target and off-target ("antitarget") reads to recover genome-wide resolution from hybrid-capture data. CNVkit is a depth-only caller: it produces **relative** log2 copy ratios, not absolute integer copy number, tumor purity, or allele-specific state. Use it as a fast screening caller and escalate to an allele-specific caller (ASCAT/Sequenza/FACETS/PureCN) when absolute CN, loss of heterozygosity, purity, or whole-genome-doubling status is needed.

## Prerequisites

```bash
conda install -c bioconda cnvkit          # pulls DNAcopy (CBS) and dependencies
pip install pomegranate                   # only if HMM segmentation is needed
```

Inputs: aligned/sorted/indexed BAM files; a target BED (capture regions); the reference FASTA used for alignment; optionally `refFlat.txt` for gene annotation. For best results, collect 5-20 process-matched normal samples (same capture kit, same lab) for a panel of normals.

## Quick Start

Tell the AI agent what to do:
- "Call CNVs from my tumor-normal exome pair with CNVkit"
- "Build a panel of normals from my 12 control exomes, then call CNVs on the tumors"
- "My panel is amplicon-based - run CNVkit in the correct mode"
- "Diagnose why my tumor-only CNVkit calls look noisy and recurrent across samples"
- "Convert CNVkit log2 segments to purity-corrected integer copy number"

## Example Prompts

### Calling

> "Run CNVkit on this tumor-normal exome pair, drop low-coverage bins, and call purity-corrected integer copy number given an estimated tumor purity of 0.6."

> "Build a pooled CNVkit reference from these 15 matched normals, then process all tumor BAMs against it and report MAD for each."

### Tool and mode selection

> "This is a 50-gene multiplex-PCR amplicon panel - configure CNVkit appropriately and explain why antitarget bins must be dropped."

> "I have no matched or pooled normal for this FFPE tumor. Tell me whether CNVkit can still be used and what the limitations of a flat reference are."

### Troubleshooting and reconciliation

> "My CNVkit profile shows the whole genome as gained. Diagnose whether this is a centering artifact and what to do about whole-genome doubling."

> "CNVkit and GATK disagree on a focal amplification. Walk through how to reconcile depth-only callers and decide which to trust."

## What the Agent Will Do

1. Choose the CNVkit mode (hybrid / amplicon / WGS) from the assay type
2. Generate target and antitarget regions, restricted to accessible genome
3. Build or apply a reference; a pooled panel of normals is strongly preferred over flat
4. Compute coverage, fix bias, and segment the log2 profile (CBS or HMM)
5. Call copy-number states, applying purity correction when purity is known
6. Run QC metrics (MAD, spread, inferred sex) and flag noisy samples
7. Reconcile against orthogonal callers and recommend escalation when allelic data is needed

## Tips

- A pooled panel of normals is the single biggest quality lever; a flat reference only corrects GC content and leaves capture bias as false focal calls.
- Pass `--drop-low-coverage` for any tumor, FFPE, or low-input sample; zero-coverage bins are otherwise miscalled as homozygous deletions.
- Use `--method amplicon` for multiplex-PCR panels; antitarget bins are pure noise there.
- Check MAD in `cnvkit.py metrics` output: < 0.5 acceptable, < 0.3 good, > 0.5 unreliable.
- Depth-only calling fails below ~40% tumor purity and on whole-genome-doubled genomes; escalate to an allele-specific caller rather than trusting CNVkit centering.
- CNVkit can overlay BAF from a VCF but does not jointly fit purity from it; that needs ASCAT, Sequenza, FACETS, or PureCN.

## Related Skills

- copy-number/copy-ratio-segmentation - CBS vs HMM choice, depth normalization, bias correction
- copy-number/allele-specific-copy-number - Purity, ploidy, integer allele-specific CN
- copy-number/gatk-cnv - GATK depth-based CNV alternative
- copy-number/cnv-annotation - Gene and clinical annotation of CNV calls
- copy-number/cnv-visualization - CNV profile plots and cohort heatmaps
- copy-number/recurrent-cnv - Cohort-level recurrent and driver CNV with GISTIC2
- alignment-files/bam-statistics - QC of input BAMs
