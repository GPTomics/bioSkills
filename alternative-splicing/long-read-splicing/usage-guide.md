# Long-Read Splicing - Usage Guide

## Overview
Analyze alternative splicing from PacBio Iso-Seq (HiFi, Kinnex/MAS-Iso-seq) and Oxford Nanopore (direct cDNA, direct RNA, R10.4.1+) long-read RNA-seq with full-isoform resolution. Tools include FLAIR (correct/collapse/quantify/diffSplice), IsoQuant (de-novo or annotation-guided isoform discovery, 2024 SOTA), Bambu (annotation-aware Bayesian discovery + quantification with NDR), SQANTI3 (isoform classification), rMATS-long (rMATS event calling on long-read isoforms), and minimap2 (-ax splice:hq for HiFi; -ax splice -k14 for ONT cDNA; add -uf only for direct RNA / stranded cDNA). Solves microexon detection, recursive splicing, complex multi-exon isoforms, and DTU without transcript-quantification uncertainty.

## Prerequisites
```bash
# Python tools
pip install flair-brookslab isoquant
conda install -c bioconda minimap2 samtools

# R tools
BiocManager::install(c('bambu', 'IsoformSwitchAnalyzeR'))

# SQANTI3
git clone https://github.com/ConesaLab/SQANTI3
# Plus dependencies (cDNA_Cupcake, kallisto, etc.)

# rMATS-long (separate from rMATS-turbo)
git clone https://github.com/Xinglab/rMATS-long
```

## Quick Start
Tell your AI agent what you want to do:
- "Discover and quantify isoforms from my PacBio Iso-Seq HiFi data"
- "Run FLAIR pipeline (correct -> collapse -> quantify -> diffSplice) on ONT direct cDNA"
- "Detect microexons and recursive splicing that short-read tools miss"
- "Use Bambu for annotation-aware discovery + quantification with calibrated NDR"
- "Process MAS-Iso-seq + 10X 5' single-cell long-read data for full isoform per cell"

## Example Prompts

### PacBio Iso-Seq Discovery
> "Align HiFi reads with minimap2 -ax splice:hq, run IsoQuant for de-novo discovery, classify with SQANTI3."

### ONT Direct cDNA (unstranded by default; PCS-114, PCB-114)
> "Align ONT direct cDNA with minimap2 -ax splice -k14 (omit -uf for unstranded cDNA), run FLAIR correct->collapse->quantify->diffSplice between conditions. Add -uf only for direct RNA (RNA004) or stranded cDNA preps."

### Discovery + Quant in One
> "Use Bambu in R for joint isoform discovery and quantification with NDR=0.1 for balanced novel discovery."

### Event-Level on Long Reads
> "Use rMATS-long to identify SE/A5SS/A3SS/MXE/RI events from long-read isoform GTFs between two groups."

### Single-Cell Long-Read
> "Demultiplex MAS-Iso-seq Kinnex array with skera, run lima -> isoseq3 refine -> FLAMES for barcode rescue and per-cell isoform counts."

### DTU on Long-Read Counts
> "Run DRIMSeq+DEXSeq+stageR DTU on FLAIR transcript counts (no Salmon Gibbs samples needed since long-read counts are read-level identities)."

### Microexon-Specific
> "Confirm microexon inclusion in neural samples using long-read coverage; cross-validate with VAST-TOOLS short-read calls."

## What the Agent Will Do
1. Choose alignment preset based on platform (`splice:hq` for HiFi; `splice -uf -k14` for ONT)
2. Run isoform discovery (IsoQuant for novel-heavy; Bambu for annotation-aware; FLAIR for end-to-end)
3. Classify isoforms with SQANTI3 (FSM/ISM/NIC/NNC + artifact flags)
4. Filter SQANTI3 RT-switching and intra-priming flags
5. Run differential analysis (FLAIR diffSplice; or DTU pipeline; or rMATS-long for events)
6. Cross-validate against short-read calls when available

## Tips
- minimap2 `-ax splice:hq` is the HiFi preset (NOT plain `-ax splice` - that's for ONT)
- ONT R10.4.1+ is required for reliable splicing; pre-R10 had ~85-90% junction concordance and false-positive novel junctions
- Bambu's NDR (Novel Discovery Rate) replaces per-sample heuristic thresholds with one calibrated parameter
- IsoQuant memory requirement: >=64 GB for atlas-scale runs
- SQANTI3 is the standard long-read isoform-curation/QC tool; filter intra-priming and RT-switching flags before reporting
- Direct RNA sequencing preserves modifications (m6A) and true 5'/3' termini but is lower throughput
- For single-cell long-read, MAS-Iso-seq + 10X 5' is the practical SOTA in 2024-2026
- Long-read DTU bypasses Salmon EM uncertainty - counts are read-level isoform identities
- rMATS-long is a separate tool from rMATS-turbo; the older lr2rmats only augmented short-read GTFs
- IsoformSwitchAnalyzeR v2 has explicit long-read input mode

## Related Skills

- splicing-quantification - Short-read PSI for cross-validation
- isoform-switching - DTU framework on long-read counts
- single-cell-splicing - MAS-Iso-seq + 10X integration
- long-read-sequencing/isoseq-analysis - PacBio Iso-Seq general pipeline
- long-read-sequencing/long-read-alignment - minimap2 splice:hq details
- long-read-sequencing/long-read-qc - QC for long-read data
