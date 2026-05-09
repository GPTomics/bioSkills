# ATAC-seq Pipeline - Usage Guide

## Overview

End-to-end ATAC-seq workflow from raw FASTQ files to accessible chromatin regions, ENCODE 4 quality control, fixed-width consensus peaksets, differential accessibility, transcription factor footprinting, motif accessibility variability, and nucleosome positioning. Optional extensions for single-cell ATAC, deep-learning variant effect prediction, enhancer-gene linking, and allele-specific accessibility.

## Prerequisites

```bash
# Core CLI tools
conda install -c bioconda fastp bowtie2 samtools macs3 deeptools bedtools \
    tobias picard idr subread

# Optional: deep learning, ENCODE-rE2G, WASP
pip install chrombpnet tangermeme tfmodisco-lite
git clone https://github.com/broadinstitute/ABC-Enhancer-Gene-Prediction
git clone https://github.com/bmvdgeijn/WASP
```

```r
# R packages
BiocManager::install(c('DiffBind', 'csaw', 'DESeq2', 'ChIPseeker', 'chromVAR',
                       'motifmatchr', 'JASPAR2024', 'ATACseqQC'))
```

## Quick Start

Tell your AI agent what you want to do:
- "Run the ATAC-seq pipeline on my samples following ENCODE 4 standards"
- "Call accessibility peaks from my ATAC-seq data with proper Tn5 shift correction"
- "Find differential accessibility between treatment and control with spike-in normalization"
- "Build a Corces 2018 fixed-width consensus peakset before differential testing"
- "Run TF footprinting with TOBIAS and verify CTCF aggregate shows clean V-shape"

## Example Prompts

### Standard Bulk Pipeline
> "Process my paired-end ATAC-seq FASTQs through fastp, Bowtie2 with `-X 2000 --very-sensitive`, MAPQ>=30 + chrM strip, Picard MarkDuplicates, alignmentSieve `--ATACshift`, MACS3 with `-f BAM --shift -75 --extsize 150 -p 0.01` for ENCODE-style peaks."

### ENCODE-Compliant Quality Control
> "Compute ENCODE 4 QC: TSS enrichment using pyTSSe, FRiP, NRF/PBC1/PBC2 from raw mapped BAM, fragment-size periodicity, mitochondrial fraction. Grade each metric PASS/WARN/FAIL."

### Consensus and Differential
> "Build a Corces 2018 iterative-overlap consensus peakset (501 bp fixed-width). Run DiffBind on the consensus with `summits=250` and DESeq2 backend. If the treatment globally compacts chromatin, switch to `DBA_NORM_LIB` instead of `DBA_NORM_NATIVE`."

### Footprinting
> "Run TOBIAS three-step (ATACorrect, ScoreBigwig, BINDetect) for differential TF footprints between two conditions. Validate with CTCF aggregate plot."

### Single-Cell Extension
> "I have 10X scATAC instead of bulk. Switch to Signac (or ArchR for >100K cells), follow per-cell QC thresholds, do TF-IDF + LSI dims=2:30, leiden clustering, then per-cluster pseudobulk peak calling."

### Variant Interpretation
> "I have GWAS lead SNPs in my ATAC peaks. Run chromBPNet variant scoring (atac-seq/deep-learning-atac); cross-validate with allele-specific accessibility from WASP-filtered BAMs (atac-seq/allele-specific-accessibility)."

## Input Requirements

| Input | Format | Description |
|-------|--------|-------------|
| FASTQ files | .fastq.gz | Paired-end reads |
| Reference | FASTA + Bowtie2 index | Genome with index |
| Chrom sizes | .chrom.sizes | UCSC chrom-sizes file |
| Blacklist | BED | ENCODE blacklist v2 (Amemiya 2019) |
| Motifs (optional) | JASPAR PFM / HOCOMOCO | For footprinting / motif analysis |
| Hi-C / Micro-C (optional) | cooler / hic | For ABC enhancer-gene linking |

## What the Workflow Does

1. **Quality Control** - Trim Nextera adapters; per-base Q30 fraction
2. **Alignment** - Map reads with Bowtie2 (or bwa-mem2 / chromap)
3. **BAM Processing** - Remove chrM, MAPQ filter, deduplicate, Tn5 shift
4. **ENCODE QC** - TSS enrichment, FRiP, NRF/PBC1/PBC2, mt fraction
5. **Peak Calling** - MACS3 with ENCODE-style shift-extend; IDR + pseudoreplicates
6. **Consensus Peakset** - Corces 2018 iterative overlap (501 bp fixed-width)
7. **Differential** - DiffBind / DESeq2 / csaw; normalization based on biology
8. **Footprinting** - TOBIAS three-step; per-TF differential
9. **Optional**: motif variability (chromVAR), nucleosome positioning, scATAC, variant scoring, enhancer-gene linking

## ATAC-seq vs ChIP-seq Processing

| Aspect | ATAC-seq | ChIP-seq |
|--------|----------|----------|
| Adapters | Nextera | TruSeq |
| Control | None | Input required |
| Tn5 shift | Yes (+4/-5 bp via alignmentSieve) | No |
| chrM | High, must remove | Low |
| Peak type | Narrow (regulatory elements) | Narrow or broad |
| Normalization | Reads-in-peaks default; library-size for global shifts | Library-size standard |

## Tips

- Mitochondrial reads: expect 20-50% chrM in standard ATAC, <5% in Omni-ATAC; always filter before peak calling.
- Tn5 +4/-5 shift is essential for footprinting. Use `alignmentSieve --ATACshift` (deepTools) for the canonical shift, or apply manually in BED conversion.
- ENCODE 4 TSS enrichment threshold is >= 7 (ideal); 5-7 is acceptable; < 5 fails.
- MACS3 `-f BAMPE` silently ignores `--shift/--extsize`; use `-f BAM` for the ENCODE pattern.
- Footprinting requires >= 50M nuclear reads; weaker libraries cannot reliably call transient TFs.
- For consensus peaksets used in differential or ML, always re-center on summits with fixed width (501 bp Corces 2018 standard).
- Spike-in normalization (Reske 2020) is required when treatment causes global accessibility shift.
- For scATAC, use Signac/ArchR/SnapATAC2 instead; pipeline above is bulk-specific.

## Related Skills

- atac-seq/atac-peak-calling - MACS3 / Genrich / HMMRATAC details, ENCODE 4 IDR
- atac-seq/atac-qc - ENCODE 4 thresholds with citations
- atac-seq/consensus-peakset - Corces 2018 iterative overlap
- atac-seq/differential-accessibility - DiffBind / csaw / DESeq2; spike-in
- atac-seq/footprinting - TOBIAS three-step; per-TF failure modes
- atac-seq/motif-deviation - chromVAR motif variability
- atac-seq/nucleosome-positioning - V-plot, NucleoATAC, +1 nucleosome
- atac-seq/single-cell-atac - scATAC alternative
- atac-seq/co-accessibility - Cicero cis-regulatory
- atac-seq/enhancer-gene-linking - ABC, ENCODE-rE2G
- atac-seq/deep-learning-atac - chromBPNet variant effects
- atac-seq/allele-specific-accessibility - WASP + caQTL
- read-alignment/bowtie2-alignment - Upstream alignment details
- alignment-files/duplicate-handling - Pre-call dedup
- chip-seq/peak-annotation - Annotate ATAC peaks to genes
