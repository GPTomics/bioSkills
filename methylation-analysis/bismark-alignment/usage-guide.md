# Bismark Alignment - Usage Guide

## Overview
Bismark aligns bisulfite-converted short reads (WGBS, RRBS, PBAT) and enzymatic-conversion reads (EM-seq) to an in-silico C->T/G->A-converted reference, recovering methylation by comparing the original read to the original reference. The idea underneath everything: methylation is never sequenced directly - the run reads which cytosines survived deamination, against a 3-letter genome stripped of cytosines. Every call rests on two silent bets: that conversion completed in both directions, and that the read mapped despite throwing its cytosines away. This skill covers the genome index, the directional/non-directional/PBAT strand flag, deduplication, and conversion QC; it stops at a deduplicated, M-bias-aware BAM, handing extraction to methylation-calling.

## Prerequisites
```bash
conda install -c bioconda bismark bowtie2 hisat2 samtools trim-galore fastqc

bismark --version
bowtie2 --version
trim_galore --version
```

Conceptual prerequisites:
- The bisulfite index is built once per genome FASTA with a chosen backend (`--bowtie2` or `--hisat2`); the backend must match between genome preparation and alignment.
- The library protocol determines the strand flag: directional WGBS/EM-seq/RRBS use the default, PBAT/scBS-seq need `--pbat`, non-directional libraries need `--non_directional`. Guessing wrong silently drops reads.
- Conversion QC needs spike-in controls in BOTH directions: unmethylated lambda for under-conversion, CpG-methylated pUC19 for over-conversion. A single control bounds only one error direction.
- EM-seq uses the identical aligners and flags as bisulfite; only the upstream chemistry and the (flatter, higher) coverage/efficiency expectations differ.
- Standard bisulfite and standard EM-seq report 5mC and 5hmC summed; separating them requires an oxBS/TAB pairing.

## Quick Start
Tell your AI agent what you want to do:
- "Prepare a bisulfite genome index for my hg38 reference"
- "Align my paired-end WGBS reads and deduplicate them"
- "Align my RRBS data without deduplication"
- "Align my PBAT single-cell library"
- "Check my bisulfite conversion efficiency with lambda and pUC19 spike-ins"
- "My Bismark mapping rate is very low, help me diagnose it"

## Example Prompts

### Genome Preparation
> "Build a Bismark bisulfite index from my hg38.fa using HISAT2 so it uses less memory on this large genome."

### WGBS / EM-seq Alignment
> "Trim adapters, align my paired-end EM-seq reads to the prepared genome, deduplicate, then sort and index the BAM."

### RRBS
> "I have RRBS FASTQ files. Trim them with the MspI fill-in handling, align with Bismark, and do not deduplicate. Explain why dedup is skipped."

### PBAT / Single-Cell
> "Align my scBS-seq PBAT library to mm10 with the correct strand flag and aggressive 5' clipping."

### Conversion QC
> "Align my lambda and pUC19 spike-ins separately and report the under-conversion and over-conversion rates."

### Troubleshooting
> "My Bismark mapping efficiency is near zero. Walk through the diagnosis in order before relaxing any alignment parameters."

> "Help me figure out whether my library is directional, non-directional, or PBAT, and which flag to pass."

## What the Agent Will Do
1. Prepares the bisulfite genome index once with the chosen backend (if not already built), keeping the backend consistent for alignment.
2. Trims adapters and library-specific end artifacts with Trim Galore, adding `--rrbs`, `--non_directional`, or PBAT 5' clipping as the protocol requires.
3. Aligns the trimmed reads with Bismark, passing the strand flag that matches the library type.
4. Deduplicates WGBS/EM-seq output on the by-name BAM, and skips deduplication for RRBS.
5. Sorts and indexes the BAM for downstream tools and visualization.
6. Quantifies conversion efficiency from lambda and pUC19 spike-ins (both error directions) and reads mapping efficiency and per-context methylation from the Bismark report.
7. Flags M-bias at read ends so the extraction step can clip appropriately.

## Tips
- Always prepare the genome index once before aligning, with the same backend used at alignment time.
- Trim before Bismark; Bismark does not trim. Use `--rrbs` for RRBS and the PBAT 5' clip for random-priming bias.
- Skip deduplication for RRBS: MspI fixed fragment ends are not PCR duplicates. Deduplicate WGBS/EM-seq.
- A 50-70% WGBS mapping efficiency is normal (3-letter reduced complexity); EM-seq is usually higher.
- Diagnose a low mapping rate in order: library-type flag first, then trimming, then reference correctness, then real biological divergence. `-N 1` is a last resort that increases mis-mapping.
- Report conversion in both directions: lambda (under-conversion, target <=1% residual) and pUC19 (over-conversion, ~96-98% methylated). The spike-in is an optimistic floor; GC-rich regions under-convert more.
- A FastQC per-base sequence content or GC FAIL is expected for converted libraries (C is depleted) and is not a defect.
- For precious low-input DNA (cfDNA, FFPE, single-cell), prefer EM-seq or TAPS upstream; bisulfite degrades 84-96% of input.
- Treat methylation at C/T-polymorphic sites as a hypothesis until a SNP-aware caller is used.

## Related Skills

- methylation-calling - Extract per-CpG methylation from the aligned BAM
- methylkit-analysis - Downstream import, filtering, normalization
- read-qc/adapter-trimming - Trim Galore before Bismark (RRBS/PBAT handling)
- read-qc/quality-reports - FastQC (expect per-base C-depletion FAIL on converted libraries)
- alignment-files/sam-bam-basics - BAM manipulation after alignment
- sequence-io/read-sequences - FASTQ handling before alignment
- long-read-sequencing/nanopore-methylation - Native long-read MM/ML modification calling (out of scope here)
- workflows/methylation-pipeline - End-to-end bisulfite pipeline
