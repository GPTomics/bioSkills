# CLIP-seq Peak Calling - Usage Guide

## Overview

Call protein-RNA binding sites from a deduplicated CLIP-seq BAM. Three statistical frameworks compete: coverage-based (CLIPper, Piranha, MACS3), single-nucleotide-crosslink-based (PureCLIP, CTK CIMS/CITS, PARalyzer), and window-based (Skipper). The Yeo lab CLIPper + SMInput log2 normalization is the ENCODE canon; Skipper 2023 is the modern high-sensitivity replacement; PureCLIP delivers single-nt CL maps for motif registration and allele-specific binding. ENCODE stringent thresholds are log2(IP/SMInput) >= 3 AND -log10 p >= 3.

## Prerequisites

```bash
conda install -c bioconda clipper pureclip piranha bedtools samtools idr
# Skipper - Snakemake workflow
git clone https://github.com/algaebrown/skipper && cd skipper && conda env create -f environment.yaml
# CTK
git clone https://github.com/chaolinzhanglab/ctk
# CLAM
pip install CLAM
```

## Quick Start

Tell your AI agent what you want to do:
- "Run CLIPper on my eCLIP and normalize against SMInput per ENCODE convention"
- "Apply ENCODE stringent thresholds: log2 FC >= 3, -log10 p >= 3"
- "Call single-nucleotide crosslink sites with PureCLIP"
- "Use Skipper for maximum sensitivity with proper input normalization"
- "Rescue repeat peaks with CLAM"
- "Call PAR-CLIP T-to-C sites with PARalyzer or CTK CIMS"
- "Run IDR for reproducibility between my two replicates"

## Example Prompts

### ENCODE-comparable CLIPper workflow

> "CLIPper on rep1 and rep2 of my RBFOX2 eCLIP, normalize against SMInput, filter at stringent log2 >= 3 / -log10 p >= 3, IDR threshold 0.05"

### Skipper (maximum sensitivity)

> "Use Skipper for FASTKD2 - mitochondrial RBP; CLIPper misses chrM peaks"

> "Compare Skipper site count to CLIPper for PUM2 - expect 2-3x more sites in 3' UTRs"

### Single-nt crosslink sites

> "Run PureCLIP with the IP BAM and SMInput BAM; output single-nt sites for mCross motif registration"

> "Get CTK CITS truncation sites for my iCLIP - this is what PureCLIP would call but using empirical FDR"

### Repeat-binding RBPs

> "MATR3 binds LINE-1 - run CLAM on the multi-mapper BAM from STAR with --outFilterMultimapNmax 100"

### PAR-CLIP

> "PARalyzer on PAR-CLIP data with T->C transitions; output high-confidence kernel-density clusters"

> "CTK CIMS with -substitution T C for PAR-CLIP - empirical FDR alternative to PARalyzer"

### No SMInput

> "I only have IP BAM, no SMInput - call peaks with Piranha ZTNB or PureCLIP without -ibam"

### Differential / replicate pooling

> "Pool rep1 and rep2 CLIPper calls; IDR threshold 0.05 for true reps and 0.10 for pseudoreps"

> "Skipper output is windowed - downstream Flipper for differential binding across conditions"

## What the Agent Will Do

1. Choose the caller based on (a) CLIP variant, (b) SMInput availability, (c) RBP binding mode, (d) goal (ENCODE comparison vs sensitivity vs single-nt)
2. Run the caller with appropriate parameters (CLIPper `--FDR 0.05 --superlocal`; PureCLIP with -ibam; Skipper Snakemake; Piranha `-d ZTNB -l` with covariates)
3. For CLIPper: normalize against SMInput with the Yeo lab scripts to produce log2 FC + -log10 p
4. Apply ENCODE stringent or lenient thresholds
5. Run IDR for replicate reproducibility; ENCODE rule: max/min of Nt and Nself <= 2
6. Output QC: peak count, mean width, FRiP, chrM peaks, top-peak transcript distribution
7. Flag failure modes (housekeeping bias from Piranha, focal PureCLIP miss of broad zones, CLIPper chrM blind spot)

## Tips

- **CLIPper alone is not "ENCODE peaks".** Always normalize against SMInput.
- **Stringent: log2 FC >= 3 AND -log10 p >= 3.** Lenient: log2 >= 1 AND -log10 >= 2. Both apply to IDR-passing peaks.
- **Skipper is the 2023 sensitivity king.** 210-320% more sites than CLIPper for mRNA-binding RBPs; requires SMInput.
- **PureCLIP is focal.** Only ~4 sites per CLIP on benchmark test sets. Pair with CLIPper for broad sites.
- **Piranha biases to high-expression.** Use SMInput as covariate via `-c -l` or switch to Skipper.
- **CLIPper misses chrM.** Switch to Skipper for mitochondrial RBPs (FASTKD2, LRPPRC).
- **PureCLIP misses broad zones.** Use it for single-nt CL maps; not for broad sites.
- **CLAM rescues repeat peaks.** Required for MATR3, ZFP36, HNRNPK at SINE/LINE, FUS at LINE-1.
- **CTK CIMS for HITS-CLIP/PAR-CLIP; CTK CITS for iCLIP/eCLIP.** Mutation vs truncation methods.
- **IDR rank by signalValue or log2 FC, not raw score.** CLIPper score is sparse and tied.
- **FRiP >= 0.005 for ENCODE narrow.** Below this = failed library or atypical-binding RBP.

## Related Skills

- clip-seq/clip-preprocessing - Upstream preprocessing
- clip-seq/clip-alignment - Upstream alignment
- clip-seq/clip-qc - FRiP, IDR, library complexity
- clip-seq/crosslink-site-detection - Single-nt CL sites
- clip-seq/differential-clip - Cross-condition differential
- clip-seq/binding-site-annotation - Annotate peaks
- clip-seq/clip-motif-analysis - Motifs on peaks
- chip-seq/peak-calling - DNA-protein comparison
