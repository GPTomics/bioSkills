# CLIP-seq Pipeline - Usage Guide

## Overview

Complete CLIP-seq workflow from raw FASTQ to ENCODE-compliant binding sites, single-nucleotide crosslink maps, annotated peaks, motifs, and (optionally) differential binding. Orchestrates the 12 decision-grade clip-seq skills with ENCODE-standard parameters: 10 nt UMIs for eCLIP / NNNXXXXNN for iCLIP, 3'-only adapter trimming (preserving R2 5' truncation = crosslink site -1), STAR end-to-end alignment with 0.04 mismatch ceiling (0.07 for PAR-CLIP T->C), UMI-unique deduplication, CLIPper or Skipper peak calling against SMInput, stringent log2 FC >= 3 AND -log10 p >= 3 thresholds, and IDR rescue + self-consistency both < 2.

## Prerequisites

```bash
# Preprocessing + alignment + QC
conda install -c bioconda umi_tools cutadapt fastp star samtools bowtie2 preseq picard multiqc rseqc bedtools

# Peak calling and downstream (CLIPper is NOT on bioconda -- the bioconda `clipper` is an
# unrelated R package. Install YeoLab CLIPper from source:
#   git clone https://github.com/YeoLab/clipper && cd clipper && conda env create -f environment3.yml && pip install .)
conda install -c bioconda pureclip homer idr meme

# R / Bioconductor for annotation and differential
BiocManager::install(c('ChIPseeker','DEWSeq','TxDb.Hsapiens.UCSC.hg38.knownGene'))
```

**Input data:**
- Paired IP and size-matched input (SMInput) FASTQ, plus a STAR index and the genome FASTA
- A chromosome-sizes file and a BED of expressed 3'UTRs, for the GC-matched motif background. Motif discovery is SKIPPED without them: HOMER run without an explicit `-fasta` background first-order-scrambles the input and reports a spuriously U/AU-rich logo on any UV-CLIP library.

## Quick Start

Tell your AI agent what you want to do:
- "Run the full eCLIP pipeline from FASTQ to ENCODE-stringent peaks for my RBP and SMInput pair"
- "Adapt for PAR-CLIP - raise STAR mismatch ceiling to 0.07 and use PARalyzer downstream"
- "Run the iCLIP2 single-end variant; demultiplex NNNXXXXNN first"
- "Repeat-binding RBP (MATR3) - use multi-mapper alignment and CLAM peak rescue"
- "Add the differential binding step between treatment and control replicates with DEWSeq"
- "End the pipeline at single-nucleotide crosslink sites for mCross motif registration"
- "Run the QC five-gate diagnostic first; tell me which gate fails"

## Example Prompts

### Full pipeline

> "Complete eCLIP run for RBFOX2: preprocess with 10 nt UMI + two-pass trim, align with STAR ENCODE block, dedupe with --method=unique, CLIPper + SMInput normalization, stringent log2 FC >= 3 / -log10 p >= 3, IDR across two replicates"

> "Same but use Skipper for FASTKD2 - chrM peaks; CLIPper misses them"

### Protocol-specific variants

> "PAR-CLIP for HuR: raise STAR --outFilterMismatchNoverReadLmax to 0.07; use PARalyzer with Hafner 2010 parameters"

> "iCLIP2 multiplexed library: demultiplex NNNXXXXNN library barcode first, then UMI extract the 5 random Ns"

> "Repeat-binding MATR3: STAR --outFilterMultimapNmax 100 --outSAMmultNmax -1; CLAM EM rescue downstream"

### Step gates

> "Just run the five QC gates on my existing dedup BAM - preprocessing retention, alignment rate, library complexity, FRiP, IDR"

> "End at single-nt crosslink sites with PureCLIP; I'll feed them to mCross myself"

### Across conditions

> "After per-condition pipeline, run DEWSeq differential with type + condition + type:condition interaction"

### Cross-skill switches

> "Switch from RBP target identification to m6A profiling - use clip-seq/m6a-clip (miCLIP2 + m6Aboost or GLORI)"

> "Antibody unavailable - use clip-seq/stamp-antibody-free with STAMP and APOBEC1-only control"

## What the Agent Will Do

1. Identify the CLIP variant from library prep documentation
2. Choose UMI pattern (10 nt eCLIP, NNNXXXXNN iCLIP, 4 nt PAR-CLIP)
3. Apply protocol-specific preprocessing with 3'-only adapter trim and `-q 6 -m 18`
4. STAR align with ENCODE block (end-to-end, mismatch 0.04 or 0.07 for PAR-CLIP)
5. MAPQ >= 10 filter, UMI dedupe with `--method=unique`
6. Run five-gate QC (preprocessing retention, alignment rate, library complexity, FRiP, IDR)
7. Call peaks with CLIPper + SMInput log2 norm (stringent thresholds) OR Skipper
8. Detect single-nt CL sites with PureCLIP (optional)
9. Annotate with ChIPseeker (tight TSS region, transcript level)
10. Motif analysis with HOMER (GC-matched background) and mCross (CL-registered)
11. Optional differential binding with DEWSeq interaction-term design

## Tips

- **SMInput is mandatory.** Not IgG, not RNA-seq. The size-matched input captures CLIP-specific biases.
- **Preserve the R2 5' end.** It is the crosslink site -1; 5'-end trimming destroys nucleotide resolution.
- **CLIP 40-70% duplication is normal.** The IP enriches a small molecule pool. Unique fragments matter, not duplication rate.
- **PAR-CLIP exception:** Raise STAR mismatch ceiling to 0.07; T->C is signal, not error.
- **Skipper for rare-transcript binders.** CLIPper misses chrM and other rare-binding RBPs.
- **Multi-mapper rescue for repeat binders.** STAR + CLAM for MATR3, HNRNPK, FUS at LINE-1.
- **ChIPseeker tssRegion=c(-100,100) for CLIP.** The default 6 kb window is wrong for RNA biology.
- **mCross requires single-nt CL sites.** Passing a peak BED breaks the model.
- **DEWSeq interaction term is the critical decision.** `~ type + condition + type:condition` tests differential binding; `~ condition` confounds with expression.
- **FRiP >= 0.005 for narrow-binding RBPs.** Atypical (rare-transcript) binders exempt.
- **IDR rescue + self-consistency both < 2 is the ENCODE pass criterion.**

## Related Skills

- clip-seq/clip-preprocessing - UMI extraction and adapter trimming
- clip-seq/clip-alignment - STAR ENCODE block
- clip-seq/clip-qc - Five-gate QC
- clip-seq/clip-peak-calling - CLIPper / Skipper / PureCLIP
- clip-seq/crosslink-site-detection - Single-nt CL maps
- clip-seq/binding-site-annotation - ChIPseeker + RBP-Maps
- clip-seq/clip-motif-analysis - HOMER + mCross + RBNS
- clip-seq/differential-clip - DEWSeq cross-condition
- clip-seq/m6a-clip - For m6A profiling
- clip-seq/stamp-antibody-free - For antibody-free
- clip-seq/ago-clip-mirna-targets - For miRNA targets
- clip-seq/clip-deep-learning - For variant-effect prediction
- read-qc/quality-reports - Upstream FastQC / MultiQC
- reporting/automated-qc-reports - MultiQC aggregation (a triage snapshot, not a gate)
