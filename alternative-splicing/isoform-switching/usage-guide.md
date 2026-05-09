# Isoform Switching - Usage Guide

## Overview
Identify shifts in which transcript a gene predominantly uses between conditions, with functional consequence prediction (NMD sensitivity via the 50nt rule, ORF disruption, protein domain loss/gain, signal peptide changes, IDR alterations, coding-potential shifts). Operationally distinct from differential gene expression (DGE) and differential transcript expression (DTE). Tools include IsoformSwitchAnalyzeR v2 (auto-selects satuRn for >5 reps else DEXSeq), the manual DRIMSeq+DEXSeq+stageR pipeline, and fishpond/swish for inferential-uncertainty-aware testing.

## Prerequisites
```bash
# R Bioconductor
BiocManager::install(c(
 'IsoformSwitchAnalyzeR', 'DRIMSeq', 'DEXSeq', 'satuRn',
 'stageR', 'fishpond', 'tximeta', 'tximport'
))

# Salmon for transcript quantification (with --numGibbsSamples 20 for swish)
conda install -c bioconda salmon

# External annotators (run outside R)
# CPC2, Pfam, SignalP, IUPred2A or DeepTMHMM
```

## Quick Start
Tell your AI agent what you want to do:
- "Identify isoform switches between conditions and predict NMD/domain consequences"
- "Run the DTU pipeline with proper stageR multi-stage FDR control"
- "Use swish on Salmon Gibbs samples to incorporate quantification uncertainty"
- "Find genes where a poison exon switch reduces functional protein"

## Example Prompts

### IsoformSwitchAnalyzeR Workflow
> "Import Salmon quantification, pre-filter for expressed isoforms, run isoformSwitchTestSatuRn with dIF cutoff 0.1, and annotate consequences with CPC2/Pfam/SignalP/IUPred2A."

> "Generate switch plots for the top 25 hits and compute global consequence enrichment."

### Manual DTU Pipeline
> "Run DRIMSeq filter -> DEXSeq exon-bin test -> stageR two-stage FDR for proper gene+transcript-level multiple testing."

### Inferential-Uncertainty-Aware
> "Re-run Salmon with --numGibbsSamples 20, then use fishpond/swish for DTE that propagates quantification uncertainty."

### NMD-Focused
> "Find isoform switches where the alternative form is NMD-sensitive (PTC >50nt upstream of last exon-exon junction)."

> "Identify poison-exon switches in SR/hnRNP genes that autoregulate via AS-NMD."

### Long-Read Integration
> "Use IsoformSwitchAnalyzeR v2 long-read input mode with PacBio Iso-Seq counts to bypass quantification uncertainty."

## What the Agent Will Do
1. Import transcript-level quantification (Salmon/kallisto/RSEM or long-read counts)
2. Pre-filter low-expression isoforms and single-isoform genes
3. Apply appropriate statistical test (satuRn for >5 reps, DEXSeq otherwise)
4. Extract sequences for external annotation (CPC2, Pfam, SignalP, IUPred2A)
5. Re-import annotations and run analyzeSwitchConsequences
6. Summarize consequence types globally and per-gene
7. Generate switchPlot for top hits

## Tips
- Skipping stageR inflates transcript-level FDR - gene-level multiple testing burden ignored
- Increased PSI of a poison exon decreases protein - always check NMD direction
- IsoformSwitchAnalyzeR v2 default: satuRn if any condition has >5 replicates, else DEXSeq
- Long-read input bypasses Salmon EM uncertainty entirely; preferred for genes with many similar isoforms (TTN, MAPT, NEFM)
- Switch annotation depends on GENCODE basic vs comprehensive - document version
- Run swish only with Salmon's `--numGibbsSamples 20+`; Gibbs samples are required
- For single-cell, use satuRn (DEXSeq does not scale)
- Check disease-relevant tissue: TDP-43 cryptic exons appear in post-mortem ALS brain, not blood

## Related Skills

- differential-splicing - Event-level (rMATS, leafcutter, MAJIQ) complement to DTU
- splicing-quantification - PSI is a 1D projection of DTU shifts
- splice-variant-prediction - Connects SpliceAI predictions to specific isoforms
- long-read-splicing - Full-isoform DTU bypasses transcript-quant uncertainty
- pathway-analysis/go-enrichment - Pathway enrichment of switching genes
- rna-quantification/alignment-free-quant - Salmon with `--numGibbsSamples` for swish
