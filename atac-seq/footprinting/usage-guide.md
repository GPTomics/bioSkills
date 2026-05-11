# TF Footprinting - Usage Guide

## Overview

Detect transcription factor binding footprints in ATAC-seq data: short DNA stretches of reduced Tn5 cleavage within accessible regions where bound TFs physically protect DNA. Covers Tn5 sequence-bias correction (the single most important step), per-base footprint scoring, and motif-anchored bound/unbound classification. Compares TOBIAS, HINT-ATAC, Wellington, PIQ, and scprinter; documents per-TF-family failure modes (CTCF clean, nuclear receptors transient, ZBTBs unfootprintable).

## Prerequisites

```bash
conda install -c bioconda tobias rgt-hint pydnase samtools bedtools
pip install scprinter pyBigWig    # scprinter for multi-scale / single-cell
```

Inputs required:
- Deduplicated, MAPQ-filtered, chrM-stripped BAM (>= 50M nuclear reads recommended)
- Consensus peakset (BED or narrowPeak)
- Reference genome FASTA (matched build to BAM)
- ENCODE blacklist
- Motif database: JASPAR 2024 CORE vertebrates PFM (default) or HOCOMOCO v12

```bash
wget https://jaspar.genereg.net/download/data/2024/CORE/JASPAR2024_CORE_vertebrates_non-redundant_pfms.txt
mv JASPAR2024_CORE_vertebrates_non-redundant_pfms.txt JASPAR2024.pfm
```

## Quick Start

Tell your AI agent what you want to do:
- "Run the TOBIAS three-step pipeline (ATACorrect, ScoreBigwig, BINDetect) on two conditions"
- "Verify Tn5 bias correction worked by plotting CTCF aggregate footprint"
- "Use scprinter for multi-scale footprinting that handles both CTCF and nuclear receptors"
- "Filter to NFR fragments (< 100 bp) before footprinting"
- "Cross-validate predicted FOXA1 footprints against published ChIP-seq peaks"

## Example Prompts

### Differential TOBIAS
> "Run TOBIAS ATACorrect on each condition's BAM, compute footprint scores per condition, then BINDetect with `--cond_names treated control` to get per-TF differential bound counts and p-values."

### Single-Tool Sanity Check
> "Plot aggregate footprint at JASPAR CTCF MA0139.1 motif sites using TOBIAS PlotAggregate. Confirm clean V-shape -- if shallow, depth or bias correction failed."

### Multi-Scale TF Activity
> "Run scprinter on this 80M-read library to get multi-scale footprints; compare CTCF (long footprint) and FOXA1 (short pioneer footprint) activity."

### Single-Cell Footprinting
> "Run scprinter on the scATAC fragments file with cell-type clusters; report per-cluster TF activity for hematopoietic lineage TFs (PU.1, GATA1, KLF1)."

### Tool Reconciliation
> "Run both TOBIAS and HINT-ATAC; report only TFs called bound by both tools at >= 50% site overlap as high-confidence."

### NFR-Only Pipeline
> "Filter BAM to fragments < 100 bp, then run TOBIAS three-step on the NFR BAM; compare to full-fragment BAM to see footprint sharpness gain."

## What the Agent Will Do

1. Verify input depth (>= 50M nuclear reads) and library QC
2. Optionally filter BAM to NFR fragments for sharper footprints
3. Run Tn5 bias correction (TOBIAS ATACorrect or HINT's built-in dinucleotide model)
4. Validate bias correction by checking aggregate at CTCF (should show clean V)
5. Score per-base footprints across consensus peakset
6. Classify motif sites as bound/unbound with motif database (JASPAR/HOCOMOCO)
7. For differential: run BINDetect with `--cond_names`; report per-TF differential scores and p-values
8. Cross-validate against ChIP-seq peaks where available (especially for nuclear receptors and pioneer TFs)
9. Report per-TF-family caveats (transient binding, asymmetric pioneer footprints, dynamic ZBTBs)

## Tool Decision Quick Reference

| Goal | Tool |
|------|------|
| Standard differential TF activity, two conditions | TOBIAS three-step |
| Single condition, find bound TFs | TOBIAS or HINT-ATAC |
| Multi-scale (CTCF + nuclear receptors simultaneously) | scprinter |
| Single-cell ATAC | scprinter |
| Plant / non-model | TOBIAS with custom motif database (CIS-BP) |
| Lower depth (< 50M) | PIQ (less depth-sensitive); or pool replicates first |

## TF Family Quick Reference

| Family | Footprint quality | Notes |
|--------|------------------|-------|
| CTCF | Excellent (gold standard) | Use for QC validation of bias correction |
| AP-1 (FOS/JUN) | Good but composite | Use HOCOMOCO for specific dimers |
| Pioneer (FOXA1, GATA, OCT4) | Asymmetric | Stranded mode in HINT-ATAC; expect one-shoulder V |
| Nuclear receptors (ER/AR/GR) | Shallow / transient | Validate with ChIP-seq; pool replicates |
| Forkhead/Homeobox short motifs | At resolution limit | scprinter multi-scale; aggregate over many sites |
| ZBTB / dynamic TFs | Often unfootprintable | Use ChIP-seq instead |

## Tn5 Bias Correction Quick Reference

| Method | Tool | Pros | Cons |
|--------|------|------|------|
| 9-bp k-mer PWM | TOBIAS ATACorrect | Default; well-validated | Genome-wide; single-process |
| Dinucleotide HMM | HINT-ATAC | Built into HINT pipeline | Less control over output |
| seqOutBias | seqOutBias | Independent of footprinting tool | Manual integration |
| Empirical from naked DNA | (research-grade) | Most accurate | Requires input library |

## Tips

- Tn5 bias correction is mandatory. Skipping it produces "footprints" reflecting Tn5 sequence preference, not TF binding.
- Validate bias correction by inspecting the aggregate footprint at CTCF motifs first. Clean V-shape means correction worked; inverted-V or flat means it failed.
- 50M nuclear reads is the practical floor; below that, weak / transient binding TFs cannot be reliably called.
- For nuclear receptors, footprint absence does not mean absence of binding. Validate with ChIP-seq.
- ATAC footprints DO NOT replicate ChIP-seq for all TFs. Concordance is good for stable binders (CTCF ~70%), poor for transient (nuclear receptors ~30%).
- Pioneer factor footprints are asymmetric (one face of DNA on histone). HINT-ATAC has stranded mode for this.
- For differential, identical peakset and identical blacklist must be used per condition; otherwise differential scores are not comparable.
- BINDetect differential score magnitudes 0.1-0.5 are typical for biologically relevant changes; > 1.0 is rare and warrants investigation.
- JASPAR is conservatively curated; HOCOMOCO has more motifs with quality grades; CIS-BP for non-model organisms.
- NFR-only filtering (< 100 bp fragments) sharpens footprints but discards mononucleosome-borne information; do it for footprinting only, not for general accessibility.
- TOBIAS PlotAggregate is the most useful diagnostic; always run it for top-candidate TFs.

## Related Skills

- atac-seq/atac-peak-calling - Generate consensus peakset; NFR peaks are sharper for footprinting
- atac-seq/atac-qc - Confirm depth >= 50M and TSS enrichment before footprinting
- atac-seq/motif-deviation - Complementary chromVAR for accessibility variability (different question)
- atac-seq/single-cell-atac - scprinter for single-cell footprinting
- chip-seq/peak-annotation - ChIP-seq cross-validation
- sequence-manipulation/motif-search - Underlying motif-scanning details
- gene-regulatory-networks/scenic-regulons - Downstream TF -> gene regulatory networks
