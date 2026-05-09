# Splicing QC - Usage Guide

## Overview
Assesses RNA-seq data quality specifically for alternative splicing analysis. Splicing analysis is more demanding than DGE on read length, depth, library prep, alignment strategy, and annotation choice. Failures silently bias PSI estimates and inflate novel-junction false positives. QC layers include experimental design audit (library prep, read length, depth, replicates), STAR 2-pass alignment, junction saturation curves, novel-vs-known junction ratio, splice-site strength scoring (MaxEntScan, SpliceAI), strandedness, and rRNA contamination.

## Prerequisites
```bash
# Python
pip install rseqc maxentpy pysam pandas matplotlib spliceai

# CLI
conda install -c bioconda regtools fastq_screen

# Reference files
# - GENCODE basic.gtf
# - BED of canonical transcripts (for RSeQC)
# - Genome FASTA matching alignment reference
```

## Quick Start
Tell your AI agent what you want to do:
- "Audit my RNA-seq design for splicing analysis suitability"
- "Compute junction saturation curves and check whether sequencing depth is sufficient"
- "Classify junctions as known vs novel and flag samples with low known-junction fraction"
- "Score 5' and 3' splice sites with MaxEntScan and SpliceAI"
- "Verify library strandedness and recommend STAR 2-pass strategy"

## Example Prompts

### Pre-Sequencing Design Review
> "Will PE 100nt poly(A)-selected libraries at 30M reads/sample support splicing analysis? If not, what should change?"

### Post-Alignment QC
> "Run junction_saturation.py and junction_annotation.py across all my BAMs and report samples below acceptable thresholds."

### STAR 2-Pass
> "Configure cohort-style STAR 2-pass alignment for 12 samples - collect SJ.out.tab from pass 1 and re-align everything in pass 2 with the augmented junction set."

### Splice-Site Scoring
> "For my candidate cryptic splice sites, compute MaxEntScan donor and acceptor scores plus SpliceAI in-vivo usage probability."

### Diagnostic
> "I'm seeing high novel-junction rate (>40%); diagnose whether it's annotation gaps, mapping artifacts, contamination, or biology (TDP-43, SF3B1)."

## What the Agent Will Do
1. Audit experimental design vs splicing analysis requirements
2. Run RSeQC junction_saturation, junction_annotation, infer_experiment per sample
3. Compute per-junction read counts and overhang distributions with pysam
4. Score splice sites with MaxEntScan + SpliceAI
5. Check rRNA contamination
6. Flag samples failing thresholds; recommend specific fixes (re-sequencing, library re-prep, annotation update)

## Tips
- For IR analysis, rRNA-depleted libraries are mandatory; poly(A) selection biases against pre-mRNA
- 50nt single-end is poor for AS - junction-spanning reads need >=8nt overhang on each side
- Cohort-style STAR 2-pass beats per-sample basic 2-pass for differential splicing (consistent junction sets)
- High novel-junction rate may be biology (TDP-43 cryptic exons in ALS, SF3B1 cryptic 3'ss in MDS) not artifact
- MaxEntScan = sequence intrinsic strength; SpliceAI = context-aware in-vivo probability - report both
- GENCODE basic = canonical (rMATS-friendly); comprehensive = discovery (DTU-friendly, higher noise)
- Wrong strandedness halves usable junction reads; always verify with infer_experiment.py
- Junction overhang <8nt flags weakly-supported, often false-positive junctions

## Related Skills

- splicing-quantification - PSI estimation after QC passes
- read-alignment/star-alignment - STAR 2-pass detail
- read-qc/quality-reports - General sequencing QC
- read-qc/contamination-screening - rRNA / contamination
- splice-variant-prediction - SpliceAI / Pangolin for variant impact
- long-read-splicing - When short-read QC is fundamentally limiting
