# Small RNA Preprocessing - Usage Guide

## Overview

Preprocess small RNA-seq data by removing the kit-specific 3' adapter, handling any UMI or 4N degenerate bases, size-selecting to the target class window, and collapsing identical reads before quantification or discovery. The central facts are that the 3' adapter is sequenced through on every real insert (so `--discard-untrimmed` is correct), that ligation bias makes absolute cross-miRNA abundance untrustworthy within a sample, and that the post-trim read-length histogram is the primary library-quality readout. Correct handling depends on the kit: invariant-end kits get adapter trimming only, NEXTflex 4N spacers must be stripped, and QIAseq true UMIs must be extracted and used for deduplication.

## Prerequisites

```bash
pip install cutadapt umi_tools
conda install -c bioconda cutadapt fastp seqkit
```

## Quick Start

Tell your AI agent:
- "Trim the TruSeq adapter and size-select my small RNA-seq to 18-30 nt"
- "My library is NEXTflex with 4N adapters - trim and strip the random bases"
- "My library is QIAseq with UMIs - extract the UMI and deduplicate"
- "Plot the read-length distribution and tell me if the library is good"
- "Collapse identical reads to a counted FASTA for miRDeep2"

## Example Prompts

### Adapter Trimming and Size Selection

> "Trim Illumina TruSeq small RNA adapters and keep reads 18-30 nt, discarding any read without an adapter"

> "Set up a piRNA-focused window (24-35 nt) instead of the miRNA window"

> "My reads are still 60 nt after trimming - what went wrong with the adapter?"

### Kit-Specific Degenerate Bases and UMIs

> "This is a NEXTflex library - trim the adapter, then strip the 4 random nucleotides from each end"

> "This is a QIAseq miRNA library - extract the 12-nt UMI and deduplicate after alignment"

> "Should I PCR-deduplicate this TruSeq small RNA library?"

### Quality Control

> "Plot the read-length distribution after trimming and interpret the peak"

> "What fraction of my reads are adapter dimers?"

> "My length histogram has no 22 nt peak and a broad 30+ nt smear - what does that mean?"

## What the Agent Will Do

1. Identify the kit and its exact 3' adapter, and whether it carries a 4N spacer or a true UMI
2. Trim the 3' adapter with cutadapt or fastp, discarding untrimmed reads and applying the size window
3. Strip 4N degenerate bases (NEXTflex) or extract and retain the UMI (QIAseq), in the correct order relative to adapter removal
4. Collapse identical reads to a counted FASTA (preserving the `_xN` count) when the downstream tool expects it
5. Plot and interpret the read-length distribution as the primary library-quality check

## Tips

- The 3' adapter is on every real read, so `--discard-untrimmed` is the correct default - a no-adapter read is not a complete small RNA
- Ligation bias means absolute, cross-miRNA abundance within a sample is unreliable; only compare the same miRNA across samples, and never merge counts across different kits
- The read-length histogram is QC: a sharp 21-23 nt peak is a good miRNA library; a 26-32 nt peak is piRNA; a broad 30+ nt smear is degradation or tRNA/rRNA contamination
- Do not PCR-deduplicate small RNA without a true UMI - distinct molecules share sequence and position, so position-based dedup deletes real signal
- Distinguish the NEXTflex 4N spacer (discard it) from the QIAseq 12-nt UMI (keep and use it); conflating them corrupts quantification
- Order matters: trim the adapter first, then strip 4N or extract the UMI
- Use DV200, not RIN, to judge small-RNA input quality - small RNAs survive degradation that destroys RIN
- Run miRTrace before quantifying: it reports RNA-class composition and fingerprints clade-specific miRNAs to catch cross-species/reagent contamination that a good mapping rate hides
- The "long/broad = degradation" rule assumes a standard ligation prep; for PANDORA-seq/phospho-RNA-seq the broad tRF/rRF distribution is the expected signal, not a failure
- For plasma/serum, check hemolysis with the miR-451a:miR-23a-3p ratio (red-cell miR-451a contamination inflates many miRNAs) and use cel-miR-39 spike-ins for low-biomass normalization

## Related Skills

- mirdeep2-analysis - Novel miRNA discovery from collapsed reads
- mirge3-analysis - Fast known-miRNA and isomiR quantification
- trf-pirna-profiling - tRF and piRNA profiling, where end chemistry and wider windows matter
- read-qc/adapter-trimming - General adapter trimming
- read-qc/umi-processing - UMI extraction and deduplication
