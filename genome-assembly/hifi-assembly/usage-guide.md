# HiFi Assembly - Usage Guide

## Overview

This skill builds haplotype-resolved diploid and telomere-to-telomere (T2T) genome assemblies from PacBio HiFi reads. HiFi reads are simultaneously long (15-25 kb) and accurate (~Q30), which made *phased diploid* assembly routine and shifted the field's deliverable from one collapsed reference to two phased haplotypes (and onward to the pangenome). The core tool is hifiasm (HiFi-only, Hi-C, or trio phasing); verkko handles true T2T when HiFi is paired with ultralong ONT and trio/Hi-C. The recurring judgment is that a "primary" assembly is a haplotype mosaic that exists in no cell, and that switch errors - the defining error of phasing - are invisible to N50, BUSCO, and base QV, so k-mer/trio QC (Merqury hap-mer blob plots, switch/hamming error) is the only window onto them.

## Prerequisites

```bash
# Core assembler + trio k-mer tool
conda install -c bioconda hifiasm yak

# GFA-to-FASTA extraction
conda install -c bioconda gfatools

# T2T (HiFi + ultralong ONT); heavy, Snakemake-orchestrated
conda install -c bioconda verkko

# Phasing QC (hap-mers, switch/hamming) - see assembly-qc
conda install -c bioconda merqury meryl
```

Profile the sample first (k-mer spectrum via GenomeScope) to learn genome size and heterozygosity - this drives the purge setting and whether phasing is even meaningful.

## Quick Start

Tell your AI agent what you want to do:
- "Assemble my diploid HiFi reads with hifiasm"
- "Phase my assembly with Hi-C since I do not have parental data"
- "Run a trio-binned assembly with both parents' Illumina reads"
- "My sample is an inbred line - assemble it without purging real duplications"
- "I have HiFi plus ultralong ONT and want a T2T assembly"
- "Convert the hifiasm GFA output to FASTA"

## Example Prompts

### HiFi-only and Hi-C phasing
> "Assemble my diploid plant HiFi reads with hifiasm, then explain why the primary assembly is not one of the two haplotypes and which output I should hand to a phased variant caller."
> "I have HiFi plus a Hi-C library and no parents. Run hifiasm Hi-C phasing and tell me whether the resulting hap1/hap2 are globally or partially phased."

### Trio phasing
> "Build yak k-mer databases from the parental Illumina reads at k=31 and run a trio-binned hifiasm assembly of the child. Note that the output uses the .dip. prefix."

### Inbred / coverage diagnostics
> "This is a doubled-haploid line. Assemble it with purging disabled so hifiasm does not delete real segmental duplications, and explain why default purging would be wrong here."
> "My hap1 came out 1.4x the size of hap2. Diagnose this as a coverage-estimate problem and re-run with --hom-cov set from the k-mer histogram."

### T2T
> "I have 40x HiFi, ultralong ONT with a 90 kb N50, and parental reads. Decide whether verkko or hifiasm --ul is the right call for a T2T assembly and justify it."

### Phasing QC handoff
> "I have a phased assembly but only N50 and BUSCO were reported. Explain why those cannot see switch errors and what k-mer QC I need to actually validate the phasing."

## What the Agent Will Do

1. Establish the sample type (outbred diploid, inbred/doubled-haploid, mole) and what data exist (parents, Hi-C, ultralong ONT) - this determines everything downstream.
2. QC the HiFi reads (length, QV, contamination) before assembling.
3. Select the phasing mechanism: HiFi-only (partial), Hi-C (global, no parents), or trio (gold standard); or verkko/hifiasm `--ul` if T2T is the genuine goal and ultralong reads exist.
4. Set the purge level appropriately (`-l0` for inbred/homozygous; default otherwise) and pin `--hom-cov` if haplotype sizes come out unbalanced.
5. Run hifiasm (or verkko) and extract contigs from GFA `S` lines to FASTA.
6. Identify the output by its prefix (`.bp.` = HiFi-only/Hi-C, `.dip.` = trio) and report which file is the mosaic primary vs the two haplotypes.
7. Route phasing QC to k-mer/trio metrics (Merqury hap-mer blob plot, switch/hamming) and contiguity/completeness/QC to assembly-qc.

## Tips

- The primary assembly is a haplotype mosaic - never use it for allele-aware analysis (phased variants, allele-specific expression, HLA/KIR typing). Use phased hap1/hap2.
- "hap1" in a filename does not certify global phasing. HiFi-only hap1/hap2 are partially phased (switch errors between blocks); only trio or Hi-C makes them globally consistent.
- Read the filename prefix: `.bp.` (balanced phasing) = HiFi-only or Hi-C; `.dip.` (diploid) = trio.
- hifiasm purges duplicates by default (`-l3` HiFi-only). That is correct for outbred diploids but DELETES real segmental duplications in inbred/homozygous/mole samples - use `-l0` there.
- An unbalanced hap1/hap2 size ratio is a mis-estimated coverage alarm, not biology. Set `--hom-cov` from the k-mer histogram.
- HiFi coverage sweet spot is ~30-40x diploid (>=13x per haplotype). Past ~40x is not free quality and can worsen false duplications.
- T2T is a project, not a flag. Ultralong ONT N50 - not HiFi coverage - is the binding constraint; verkko/`--ul` without long UL reads will not close centromeres.
- Do not reflexively polish a HiFi assembly; it is already ~Q30+ and over-polishing can lower QV.
- The k-mer size convention splits: k=31 for yak/hifiasm trio, k=30 for verkko/Merqury hap-mers. Do not mix them.
- GFA is not FASTA. Extract `S` lines (`gfatools gfa2fa` or awk) before any downstream tool, but keep the GFA for its bubble/alternate-path information.

## Related Skills

- assembly-qc - Merqury hap-mer blob plots and switch/hamming are the only QC that sees phasing
- assembly-polishing - HiFi is already accurate; deciding whether to polish at all
- scaffolding - Hi-C used here for phasing; chromosome-scale scaffolding is the adjacent step
- contamination-detection - Screen contigs for foreign sequence after assembly
- long-read-sequencing/long-read-qc - HiFi read length/QV/contamination before assembly
- workflows/genome-assembly-pipeline - End-to-end profile -> assemble -> phase -> QC -> scaffold
