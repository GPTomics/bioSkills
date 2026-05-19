# CLIP-seq Preprocessing - Usage Guide

## Overview

Turn raw CLIP-seq FASTQ (eCLIP, iCLIP, iCLIP2, iCLIP3, irCLIP, PAR-CLIP, FLASH) into UMI-deduplicated, alignment-ready FASTQ. The defining constraint of CLIP preprocessing is preservation of the 5' read end where the reverse transcriptase truncated at the protein-RNA adduct: that base is one nucleotide downstream of the crosslink and is the foundation of all single-nucleotide-resolution analyses. Quality trimming and adapter trimming must operate on the 3' end only.

## Prerequisites

```bash
conda install -c bioconda umi_tools cutadapt fastp samtools bowtie2 preseq picard
pip install pysam
```

## Quick Start

Tell your AI agent what you want to do:
- "Preprocess my eCLIP data with 10 nt UMIs"
- "Demultiplex multiplexed iCLIP2 with the NNNXXXXNN pattern"
- "Run two-pass adapter trimming for read-through events"
- "Dedupe by UMI using ENCODE convention"
- "Estimate library complexity with preseq"
- "Preprocess PAR-CLIP without losing T-to-C mutations"

## Example Prompts

### Protocol-Specific UMI Extraction

> "Extract 10-nt UMIs from R1 of my paired-end eCLIP library"

> "My iCLIP2 has the NNNXXXXNN barcode pattern - demultiplex by the 4-base library barcode first, then extract the surrounding 5 nt as UMI"

> "PAR-CLIP from the Tuschl lab - 4 nt random barcode; verify pattern by inspecting first 12 bases of reads"

### Adapter Trimming

> "Run cutadapt with -q 6 -m 18 and only the 3' adapter (-a/-A), never 5' (-g/-G on R1)"

> "Two-pass eCLIP trimming: pass 1 is standard 3' adapter, pass 2 strips the read-through 5' adapter from R2 only"

> "Raise STAR mismatch ceiling to 0.07 for PAR-CLIP so T-to-C reads are not filtered as errors"

### Deduplication and Complexity

> "After alignment, dedupe with `umi_tools dedup --method=unique` per ENCODE convention"

> "Estimate library complexity with preseq lc_extrap; flag if predicted unique < 1M at 100M reads"

> "Fall back to picard MarkDuplicates if my older library has no UMIs"

### Reconciliation

> "Compare retained read counts between fastp and cutadapt and pick the one matching ENCODE eCLIP"

> "My library has 95% PCR duplication - is this normal CLIP, or did the IP fail?"

> "Pre-map to rRNA before genome alignment to avoid STAR multi-mapper tangles"

## What the Agent Will Do

1. Inspect the library prep documentation OR the first 12 bases of reads to identify protocol and UMI pattern
2. Run `umi_tools extract` with the correct `--bc-pattern` (10 nt eCLIP R1; NNNXXXXNN iCLIP/iCLIP2; demultiplexed UMI for iCLIP3; 4 nt PAR-CLIP)
3. For iCLIP2 / multiplexed libraries: demultiplex by the inline library barcode BEFORE UMI extraction
4. Run cutadapt with the matching 3' adapter (and 5' adapter on R2 only for two-pass eCLIP) at `-q 6 -m 18`
5. (Optional) Pre-map to rRNA with bowtie2 to remove the abundant ribosomal background
6. Align with STAR (ENCODE parameters - see clip-alignment skill)
7. UMI-dedupe with `umi_tools dedup --method=unique`
8. Compute library complexity with preseq and report ENCODE compliance

## Tips

- **The 5' end is sacred.** R2 5' end in paired-end eCLIP is the crosslink site -1. Never quality-trim or adapter-trim there.
- **eCLIP uses 10 nt UMIs on R1.** Older protocols used 5 or 7 nt; verify before assuming 10.
- **iCLIP/iCLIP2 UMI is NNNXXXXNN (3+4+2).** The 4 middle bases are sample barcodes; demultiplex first, then dedupe the 5 surviving Ns.
- **PAR-CLIP needs higher mismatch tolerance.** T->C mutations are SIGNAL, not error. Raise `--outFilterMismatchNoverLmax` from 0.04 to 0.07.
- **CLIP duplication rates are 40-70% normally.** The IP enriches a small molecule pool. Low duplication means failed IP, not a good library.
- **`umi_tools dedup --method=unique` is the ENCODE default.** Directional clustering is more conservative but adds runtime.
- **18 nt minimum is mandatory.** Reads < 18 nt multi-map at > 50% rates.
- **Library complexity is the true QC.** preseq predicts plateau; ENCODE eCLIP needs >= 1M unique at sequenced depth.
- **Pre-map rRNA is performance, not biology.** Skipping it costs ~30 minutes of STAR time but does not change downstream peak calling if done correctly.

## Related Skills

- clip-seq/clip-alignment - Downstream alignment with ENCODE parameters
- clip-seq/clip-qc - QC of the preprocessed library
- clip-seq/crosslink-site-detection - Why 5' preservation matters
- clip-seq/stamp-antibody-free - STAMP uses RNA-seq preprocessing, not CLIP preprocessing
- read-qc/umi-processing - General UMI handling
- read-qc/adapter-trimming - General adapter trimming patterns
