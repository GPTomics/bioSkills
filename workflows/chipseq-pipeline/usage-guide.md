# ChIP-seq Pipeline - Usage Guide

## Overview

End-to-end ChIP-seq workflow from FASTQ to annotated peaks. Covers QC, alignment with Bowtie2, deduplication, peak calling with MACS3, ENCODE QC metrics (FRiP, NSC, RSC, library complexity), IDR or naive overlap for replicate concordance, and peak annotation with ChIPseeker. Handles both narrow peaks (transcription factors and sharp histone marks like H3K4me3) and broad peaks (H3K27me3, H3K9me3, H3K36me3). For CUT&RUN/CUT&Tag protocols use chip-seq/cut-and-run-tag instead. For experiments with expected global signal shifts (HDACi, BETi, EZH2i) add spike-in normalization per chip-seq/spike-in-normalization.

## Prerequisites

```bash
# CLI tools
conda install -c bioconda fastp bowtie2 samtools macs3 deeptools bedtools \
    phantompeakqualtools idr

# R packages
BiocManager::install(c('ChIPseeker', 'TxDb.Hsapiens.UCSC.hg38.knownGene',
                       'org.Hs.eg.db', 'DiffBind', 'ChIPQC'))
```

## Quick Start

Tell your AI agent what you want to do:
- "Run the ChIP-seq pipeline on my IP and Input samples"
- "Call narrow peaks for H3K4me3 from FASTQ with ENCODE QC"
- "Call broad peaks for H3K27me3 with naive overlap across 3 replicates"
- "Process my transcription factor ChIP-seq with IDR on biological replicates"
- "Run the pipeline and apply Drosophila spike-in normalization for cross-condition comparison"
- "Annotate peaks against ENCODE cCREs (PLS, pELS, dELS, CTCF-only)"

## Example Prompts

### Starting from FASTQ
> "I have ChIP-seq FASTQ files for 2 IP replicates and 2 input controls. Run fastp, Bowtie2, MarkDuplicates, then MACS3 narrow peaks with `-q 0.01`, then ChIPseeker annotation."

> "Call broad peaks for my H3K27me3 ChIP-seq with `--broad --broad-cutoff 0.1` and use naive overlap (>=40% reciprocal) across 3 replicates."

> "Run ChIP-seq analysis with mouse genome (mm10) and the deepTools effective genome size for 100bp reads."

### QC and validation
> "Run the full ENCODE QC battery (FRiP, NSC, RSC, NRF, PBC1, PBC2, plotFingerprint) and grade against thresholds before peak calling."

> "Apply the ENCODE Nself/Nt consistency rule across true replicates and pseudoreplicates."

> "Check for hyper-ChIPable artifacts by intersecting peaks against top-1% input signal regions."

### Annotation and downstream
> "Annotate my peaks with nearby genes using ChIPseeker host-gene convention (overlap='all')."

> "Cross-reference my peaks against ENCODE cCRE atlas via SCREEN BED and report cCRE class distribution."

> "Compare peaks between treatment and control with DiffBind background-bin normalization."

## Input Requirements

| Input | Format | Description |
|-------|--------|-------------|
| IP FASTQ | .fastq.gz | ChIP samples (immunoprecipitated chromatin) |
| Input FASTQ | .fastq.gz | Input control (sonicated chromatin, same library prep batch) |
| Reference | FASTA + Bowtie2 index | Reference genome |
| Annotation | GTF | Custom annotation OR pre-built TxDb |
| Blacklist | BED | ENCODE blacklist v2 (Amemiya 2019) |

## What the Workflow Does

1. **Pre-flight checks**: antibody validation, sequencing depth adequacy (20M+ unique mapped for TF/sharp histone; 40M+ for broad histone)
2. **Quality control with fastp**: adapter trimming, length filter, Phred quality filter
3. **Alignment with Bowtie2**: ENCODE-style parameters; `samtools view -F 1804 -q 30` filter for uniquely-mapped non-duplicate non-secondary reads
4. **BAM processing**: MarkDuplicates (do not remove yet; needed for NRF computation), chrM removal, sort and index
5. **Library complexity QC**: compute NRF, PBC1, PBC2 on PRE-deduplication BAM (post-dedup metrics are uninformative)
6. **Fragment-size diagnostic**: extract distribution; classify as sub-nucleosomal (TF), mono-nucleosomal (histone), or over-sonicated (rescue impossible)
7. **Cross-correlation QC**: phantompeakqualtools for NSC, RSC, fragment-length estimate
8. **Peak calling with MACS3**: `-f BAMPE` for paired-end, `-g` from deepTools effective-genome-size table, `--keep-dup all` since dedup applied upstream
9. **FRiP and enrichment QC**: deepTools plotFingerprint (JS distance, AUC, synthetic JS); FRiP >1% (TF) or >5% (histone)
10. **Replicate concordance**: IDR with Nself/Nt rule for TFs; naive overlap with >=40% reciprocal for histones
11. **Hyper-ChIPable filter**: build top-1% input signal blacklist; flag peaks in rRNA / mtDNA / housekeeping artifact regions
12. **Peak annotation**: ChIPseeker with custom GTF or pre-built TxDb; ENCODE cCRE cross-reference
13. **Visualization**: bigWig with appropriate normalization (CPM within-sample; RPGC cross-sample; spike-in-scaled for global shifts), heatmap, profile plot
14. **Output**: annotated peaks TSV + BED, IDR results, QC report, bigWigs, heatmaps

## Narrow vs Broad Peaks

| Target | Mode | Why |
|--------|------|-----|
| Transcription factors (CTCF, p53, GATA1) | Narrow (default) | Sharp motif binding |
| H3K4me3, H3K27ac (at promoters/enhancers) | Narrow | Localized at regulatory elements |
| H3K4me1, H3K9ac | Narrow or broad depending on cell type | Variable; check published data |
| H3K27me3, H3K9me3, H3K36me3, H4K20me3 | `--broad --broad-cutoff 0.1` | Broad domains spanning 10-100 kb |
| RNA Pol II elongation | Broad option | Gene body coverage |

## Replicate Handling

ENCODE-style replicate concordance differs by mark type. Use the matching method:

| Mark type | Method | Threshold |
|-----------|--------|-----------|
| Transcription factors | IDR (Li 2011) | True reps IDR <= 0.05; pseudoreps IDR <= 0.10 |
| Sharp histone (H3K4me3, H3K27ac) | IDR or naive overlap | Same as TF or >=40% reciprocal |
| Broad histone (H3K27me3, H3K36me3) | Naive overlap | >=40% reciprocal, present in >=2 of N replicates |

Apply the ENCODE Nself/Nt consistency rule for IDR: both `max(N_t, max(N_self)) / min(...)` and per-replicate self-consistency ratios must be <= 2.

## When to Use Spike-In Normalization

Standard reads-in-peaks normalization (DiffBind default) assumes most peaks are unchanged between conditions. For experiments where this assumption is violated, spike-in is required:

| Experimental design | Spike-in required? |
|---------------------|---------------------|
| HDAC inhibitor (global H3K27ac increase) | Yes; Drosophila ChIP-Rx |
| BET inhibitor (global BRD4 / H3K27ac decrease) | Yes |
| EZH2 inhibitor (global H3K27me3 loss) | Yes |
| Target factor knockdown / degron | Yes (or matched-input subtraction) |
| Standard TF perturbation, local rebinding | No |

For spike-in workflow, see chip-seq/spike-in-normalization.

## Tips

- **Always include matched input controls**: sonicated input from the same library prep batch, same fragmentation. IgG is not equivalent for histone marks.
- **20-50M reads for TF / sharp histone, 40-60M for broad histone**: ENCODE 2012 standard; under-sequencing produces sparse peaks and unreliable FRiP.
- **Apply the ENCODE BAM filter consistently**: `samtools view -F 1804 -q 30` for uniquely-mapped non-duplicate non-secondary; do this AFTER MarkDuplicates and AFTER NRF/PBC computation.
- **Filter against ENCODE blacklist v2**: Amemiya 2019; catches repeat-driven artifacts but NOT hyper-ChIPable transcribed genes (rRNA, tRNA, HIST cluster, mtDNA).
- **Build a custom hyper-ChIPable blacklist**: top-1% input signal regions for the cell type; intersect-out before downstream interpretation.
- **FRiP < 1% is a hard fail**: increase sequencing depth does not fix antibody / IP / fragmentation problems.
- **Use deepTools effective genome size, not `-g hs/mm` shorthand**: read-length-matched values are more accurate.
- **For broad histones use naive overlap, not IDR**: IDR is too conservative for histone signal dynamic range.
- **Document antibody catalog + lot + validation method**: peer reviewers ask; lot-to-lot variation is common.
- **For CUT&RUN/CUT&Tag protocols, use the dedicated cut-and-run-tag skill**: different peak callers (SEACR), spike-in (E. coli), and QC thresholds.

## Troubleshooting

### Few peaks
- Check FRiP, NSC, RSC; verify antibody validation
- Inspect fragment-size distribution for over-sonication
- Confirm input control matches sample library prep
- Tighten `-q` only if other metrics are good

### Many peaks (>500k)
- Did not deduplicate or remove chrM
- `-q` too loose
- Hyper-ChIPable artifacts dominating

### Peaks shifted from motif by ~75 bp
- `--shift` not set with `-f BAM`; or aligner pre-applied shift (chromap)

### IDR returns 0 reproducible peaks
- Sorted by wrong column (use column 8 p.value, not column 7 signalValue)
- Library size imbalance > 2x
- One replicate failed; do not average

### Cross-condition results have wrong sign
- Composition bias from reads-in-peaks normalization; switch to background-bin or spike-in

## Related Skills

- chip-seq/chipseq-qc - FRiP, NSC/RSC, library complexity, hyper-ChIPable detection
- chip-seq/peak-calling - MACS3/MACS2/HOMER/SPP, IDR, per-tool failure modes
- chip-seq/peak-annotation - ChIPseeker, ENCODE cCRE classification, GREAT
- chip-seq/differential-binding - DiffBind, DESeq2, csaw with composition / trended / global shift handling
- chip-seq/chipseq-visualization - bigWig normalization, heatmaps, browser tracks
- chip-seq/motif-analysis - HOMER, MEME-ChIP, monaLisa
- chip-seq/super-enhancers - ROSE/ROSE2/LILY for SE calling
- chip-seq/cut-and-run-tag - CUT&RUN/CUT&Tag protocols (different peak caller and spike-in)
- chip-seq/spike-in-normalization - Drosophila ChIP-Rx for global-shift experiments
- chip-seq/chromatin-state-segmentation - ChromHMM multi-mark chromatin states
- chip-seq/chip-deep-learning - BPNet/chromBPNet for variant effect prediction
- chip-seq/allele-specific-binding - WASP/BaalChIP for allele-specific TF binding
- read-qc/fastp-workflow - Upstream adapter trimming
- read-alignment/bowtie2-alignment - Standard ChIP-seq aligner
- alignment-files/duplicate-handling - MarkDuplicates upstream
- atac-seq/atac-peak-calling - Parallel ATAC peak calling (no input control)
