# AGO-CLIP miRNA Targets - Usage Guide

## Overview

Identify direct miRNA-target interactions experimentally. Three approaches: chimeric methods (CLEAR-CLIP, chimeric eCLIP / miR-eCLIP, CLASH) ligate miRNA-target chimeras for direct pairing; standard AGO HITS-CLIP / eCLIP + computational seed-matching (TargetScan, miRDB) for indirect pairing; HEAP (Halo-Ago2 mouse) for in vivo. Chimeric methods are gold-standard for direct identification; standard AGO-CLIP cannot assign miRNAs without computational prediction. ~50% of miRNA-target interactions are non-canonical (3' compensatory) and missed by seed-only methods.

## Prerequisites

```bash
conda install -c bioconda samtools bedtools
# Hyb pipeline for chimeras
git clone https://github.com/gkudla/hyb
# miR-eCLIP analysis from Yeo lab
```

## Quick Start

Tell your AI agent:
- "Identify direct miRNA-target pairs with chimeric eCLIP / miR-eCLIP"
- "AGO HITS-CLIP + TargetScan 7mer-m8 + miRNA expression filter for computational pairing"
- "miR-eCLIP for hsa-miR-21 deep targets with probe enrichment"
- "Run Hyb pipeline on CLEAR-CLIP data to extract chimeras"
- "Why does my AGO-CLIP have 50 candidate miRNAs per peak? No miRNA assignment without chimera"
- "Cross-reference TargetScan predictions with AGO-CLIP peaks for high-confidence"
- "HEAP mouse strain for in vivo - won't work in human"

## Example Prompts

### Chimeric Methods

> "Hyb pipeline with bowtie2 alignment mode for short miRNA chimeras"

> "miR-eCLIP for hsa-miR-21 with probe enrichment 30-175x"

### Standard AGO-CLIP + Computational

> "Call peaks from AGO eCLIP, then scan with TargetScanHuman 8.0 conserved predictions"

> "Filter for 7mer-m8 / 8mer seeds; drop 6mer (too noisy)"

### Expression Filter

> "Filter miRNA-target predictions by miRNAs > 100 TPM in matched small-RNA-seq"

### Non-Canonical Pairing

> "Chimeric methods capture 3' compensatory pairing that TargetScan misses"

> "RNAhybrid for full miRNA-target duplex prediction"

### Diagnostics

> "Chimera rate < 1% - enrich with miR-eCLIP probes for specific miRNAs"

> "TargetScan predicts thousands per miRNA - require CLIP peak overlap"

> "HEAP results don't translate to human - mouse-only method"

## What the Agent Will Do

1. Identify experimental setup: chimeric eCLIP / standard AGO-CLIP / HEAP / CLEAR-CLIP
2. For chimeric: Hyb pipeline (bowtie2 mode) to extract miRNA-mRNA chimeras
3. For standard AGO-CLIP: CLIPper / Skipper peaks + TargetScan seed-matching + miRNA expression filter
4. Apply canonical seed rules: 7mer-m8 / 7mer-A1 / 8mer; report 6mer separately
5. For non-canonical / 3' compensatory: chimeric method or RNAhybrid full duplex
6. Cross-validate: chimeric direct + TargetScan computational + matched small-RNA-seq
7. Report per-miRNA target counts; flag rare miRNAs (< 100 TPM)
8. For specific miRNA deep profiling: miR-eCLIP probe enrichment

## Tips

- **Chimeric methods are gold standard for direct pairing.** Other methods need computational inference.
- **miR-eCLIP enriches specific miRNAs 30-175x.** For deep targets of one or a few miRNAs.
- **Hyb pipeline needs bowtie2 mode for short miRNAs.** BLAST is too stringent.
- **Standard chimera rate is 1-5%.** Sequence deep (200M+) for global chimera profiling.
- **TargetScan + CLIP overlap is high-confidence.** TargetScan alone over-predicts.
- **Filter by miRNA expression > 100 TPM.** Many miRNAs in databases not expressed in your cell type.
- **6mer matches dominate but are weak.** Report separately; require 7mer-m8 / 8mer for high confidence.
- **3' compensatory misses in seed-only methods.** ~50% of miRNA-target interactions are non-canonical.
- **HEAP is mouse-only.** Use eCLIP / chimeric eCLIP for human.

## Related Skills

- clip-seq/clip-peak-calling - AGO CLIP peaks
- clip-seq/binding-site-annotation - 3' UTR annotation
- clip-seq/clip-motif-analysis - Seed motif scan
- clip-seq/differential-clip - miRNA perturbation
- clip-seq/m6a-clip - DART-seq APOBEC1 fusion
- small-rna-seq/target-prediction - TargetScan / miRDB
- small-rna-seq/differential-mirna - miRNA expression
- small-rna-seq/mirdeep2-analysis - miRNA discovery
