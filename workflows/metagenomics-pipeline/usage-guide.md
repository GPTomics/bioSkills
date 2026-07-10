# Metagenomics Pipeline - Usage Guide

## Overview

This workflow takes shotgun metagenomic FASTQ to taxonomic and functional profiles, orchestrating controls and host depletion, Kraken2+Bracken classification, MetaPhlAn marker profiling, and HUMAnN functional profiling. Every result is a position in a choice-chain (extraction -> depletion -> depth -> classifier -> database -> normalization), not a direct observation of the community. The workflow holds that chain constant across a study and reports each link, so the central discipline is consistency, not any single "best" tool. Kraken2 read counts are not abundances, and MetaPhlAn cell fractions are not comparable to Bracken read fractions.

## Prerequisites

```bash
# CLI tools
conda install -c bioconda kraken2 bracken metaphlan humann bowtie2 samtools fastp hostile nonpareil

# Taxonomic database: the database DEFINES what can be detected, so it is the dominant batch variable.
# Use a standard/PlusPF-style DB and PIN its version for the whole study; a custom DB only for a
# well-characterized environment. Include the host in the DB (or deplete host first) so host reads are
# not misassigned. Build the bowtie2/T2T-CHM13 host index for depletion.
kraken2-build --standard --db kraken2_db        # or download a pre-built, version-pinned DB
bracken-build -d kraken2_db -k 35 -l 150        # -k MUST equal the Kraken2 build k; -l the read length
```

Carry extraction blanks and a mock community through the entire workflow - they are the only way to tell a real low-biomass signal from the kitome (see metagenomics/contamination-controls).

## Quick Start

Tell your AI agent what you want to do:
- "Profile my shotgun metagenomes from FASTQ to taxonomic and functional tables"
- "Classify my reads at a sensible confidence, then estimate species abundance"
- "Profile metabolic pathway potential across my gut microbiome samples"
- "Hold the pipeline constant across my study and report each sample the same way"

## Example Prompts

### Taxonomic profiling
> "Classify my shotgun reads at a confidence threshold and estimate species-level abundance"

> "Profile my samples with a marker-gene method and tell me why its percentages are not comparable to read-count abundances"

### Functional profiling
> "Profile the metabolic pathway potential in my metagenomes and keep the unmapped fraction in the denominator"

### Controls and interpretation
> "Run blanks and a mock through the whole pipeline and flag any taxon that could be kitome"

> "Check whether my sequencing depth is adequate before I interpret a non-detection"

## Input Requirements

| Input | Format | Description |
|-------|--------|-------------|
| FASTQ files | .fastq.gz | Paired-end shotgun metagenomic reads |
| Extraction blanks + mock | .fastq.gz | Negative and positive controls carried through every step |
| Kraken2 database | Directory | Version-pinned standard/custom DB; defines what is detectable |
| Host index | Bowtie2/T2T-CHM13 | For host depletion before classification |

## What the Workflow Does

0. **Controls, QC and Host Removal** - Carry blanks/mock through; trim with fastp; deplete host against T2T-CHM13; confirm depth with Nonpareil before trusting any non-detection.
1. **Taxonomic Classification** - Kraken2 with a raised `--confidence` and `--minimum-hit-groups 2` (the default 0 over-classifies), then Bracken for read-fraction estimates; or MetaPhlAn 4 marker profiling with a pinned `--index` for cell fractions.
2. **Functional Profiling** - HUMAnN for pathway and gene-family potential; keep UNMAPPED/UNINTEGRATED as the denominator (dropping them inflates everything else).
3. **Reporting** - Taxonomic and functional tables, interpreted relative to the fixed pipeline, never as a direct community observation. Resistome and strain questions run from the same reads via their own skills.

## Choosing a Classifier: Kraken2/Bracken vs MetaPhlAn

| Feature | Kraken2/Bracken | MetaPhlAn 4 |
|---------|-----------------|-------------|
| Method | Minimizer + spaced-seed LCA over a whole-genome DB | Clade-specific marker genes |
| Speed | Very fast | Moderate |
| Error profile | Higher recall; precision depends on DB + confidence | Higher precision; FP-conservative |
| Output quantity | Read fraction (reads redistributed by Bracken) | Cell fraction (estimated organism fraction) |
| Comparability | Counts are not abundances; not comparable to MetaPhlAn % | Cell fractions; NOT comparable to Bracken read fractions |

Frame the choice by precision/recall and by what quantity you need, never as "accuracy". The two outputs answer different questions and must not be merged into one table.

## Tips

- **Controls are not optional**: blanks and a mock are the denominator for every presence call; low-biomass samples can be entirely kitome.
- **Confidence**: raise Kraken2 `--confidence` above the default 0 (0.1-0.4 on a standard/nt DB) and require `--minimum-hit-groups 2` to cut single-region false positives.
- **Pin the database**: the DB version is a batch variable; pin Kraken2's DB and MetaPhlAn's `--index` for the whole study.
- **Host removal**: critical for human-associated samples; report the reads removed - depletion can halve usable depth.
- **Depth**: confirm coverage adequacy (Nonpareil) before interpreting any absence; absence is "not detectable by this chain", not "not present".
- **Compositionality**: classifier tables are compositional - transform (CLR) before distances or differential abundance; do not compare raw percentages across samples.
- **Replicates**: needed for differential abundance; biological replication, not deeper sequencing of one sample.

## Related Skills

- read-qc/fastp-workflow - Adapter trimming and quality filtering up front
- metagenomics/contamination-controls - Host depletion, blanks/decontam, depth checks
- metagenomics/kraken-classification - Kraken2 confidence, hit-groups, minimizer evidence
- metagenomics/metaphlan-profiling - MetaPhlAn 4 parameters and cell-fraction interpretation
- metagenomics/abundance-estimation - Bracken options and compositional handling
- metagenomics/functional-profiling - HUMAnN workflow and UNMAPPED handling
- metagenomics/metagenome-visualization - Community statistics and plotting
- metagenomics/amr-detection - Community resistome from the same reads
- metagenomics/strain-tracking - Strain resolution from the same reads
- genome-assembly/metagenome-assembly - Assembly-based route when reads cannot resolve the question
