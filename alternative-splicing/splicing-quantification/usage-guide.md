# Splicing Quantification - Usage Guide

## Overview
Quantify alternative splicing events from RNA-seq data as PSI (percent spliced in) values. Supports the five canonical event classes (SE, A5SS, A3SS, MXE, RI), special classes (microexons, exitrons, AFE/ALE), and intron retention subtypes (canonical RI vs detained introns). Tools span event-based (rMATS-turbo, SUPPA2, VAST-TOOLS), LSV-based (MAJIQ V3), junction-cluster (leafcutter), and the 2025 SOTA Shiba.

## Prerequisites
```bash
# Python tools
pip install suppa pandas

# rMATS-turbo + regtools (CLI via bioconda; rMATS-turbo is not on PyPI)
conda install -c bioconda rmats regtools

# leafcutter R package (GitHub, not Bioconductor)
# devtools::install_github('davidaknowles/leafcutter/leafcutter')

# MAJIQ V3 (academic license at majiq.biociphers.org)
# VAST-TOOLS (perl-based; conda install -c bioconda vast-tools)
```

## Quick Start
Tell your AI agent what you want to do:
- "Quantify exon skipping events from my RNA-seq BAMs with rMATS-turbo"
- "Calculate PSI values for all splicing events using my Salmon TPM"
- "Run MAJIQ to quantify local splice variations including complex events"
- "Detect intron retention separately from cassette exon events"
- "Quantify microexons that standard tools miss"

## Example Prompts

### Tool Choice
> "I have STAR-aligned BAMs and a GENCODE GTF; quantify SE/A5SS/A3SS/MXE/RI events with rMATS-turbo and filter for >=10 junction reads per replicate."

> "Use SUPPA2 to compute event PSI from my existing Salmon transcript TPMs without re-aligning."

> "For complex multi-junction events, run MAJIQ V3 with the splice-graph local-splice-variation framework."

### Special Event Classes
> "Detect microexons (3-27nt) in my brain RNA-seq using VAST-TOOLS or rMATS with reduced anchor length."

> "Distinguish detained introns (nuclear-retained, regulated) from canonical NMD-targeted intron retention."

> "Separate alternative first/last exon events from spliceosomal AS - they're often promoter or APA driven."

### QC and Validation
> "Compute per-event coverage statistics from rMATS output and filter for reliable PSI estimates."

> "Cross-validate rMATS event PSI against MAJIQ LSV PSI for the same comparison."

## What the Agent Will Do
1. Choose tools based on input (BAMs vs TPM) and question (event-level vs LSV vs DTU)
2. Generate event annotations from GTF (SUPPA2 generateEvents, MAJIQ build, or use VastDB)
3. Compute per-event PSI with effective-length normalization (IncFormLen / SkipFormLen for rMATS)
4. Filter by junction read support, replicate consistency, and dynamic-range thresholds
5. Output PSI matrices, per-replicate counts, and event coordinates for downstream analysis

## Tips
- Run two complementary tools (rMATS + leafcutter) and reconcile rather than relying on one
- Use rRNA-depleted libraries for IR analysis; poly(A) selection biases against pre-mRNA
- STAR 2-pass alignment is essential - 1-pass loses ~14% of novel junctions
- For SF3B1-mutant cancer samples, expect cryptic 3'ss ~10-30nt upstream of canonical (Darman 2015)
- For TDP-43-loss ALS/FTD samples, use annotation-free leafcutter or MAJIQ denovo for cryptic exons
- JC vs JCEC: prefer JC for clean cassette analysis; JCEC adds power for short alternative exons
- AFE/ALE are usually promoter or APA driven, not spliceosomal - flag separately

## Related Skills

- differential-splicing - Test PSI differences between conditions
- isoform-switching - DTU framework with functional consequences
- splicing-qc - Sequencing depth, library, alignment QC
- splice-variant-prediction - SpliceAI/Pangolin for variant impact
- long-read-splicing - Full-isoform PSI without anchor-length limits
- read-alignment/star-alignment - STAR 2-pass setup
- rna-quantification/alignment-free-quant - Salmon/kallisto for SUPPA2
