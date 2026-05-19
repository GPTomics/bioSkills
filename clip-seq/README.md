# clip-seq

## Overview

Map and analyze protein-RNA interactions from CLIP-seq family methods (eCLIP, iCLIP, iCLIP2, iCLIP3, irCLIP, PAR-CLIP, FLASH, miCLIP2, STAMP, chimeric eCLIP). Covers preprocessing, alignment, peak calling, single-nucleotide crosslink-site detection, motif analysis, annotation, QC, differential binding, m6A profiling, antibody-free profiling, miRNA target identification, and deep learning. Decision-grade: each skill encodes failure modes, reconciliation rules, ENCODE thresholds, and per-tool nuance needed by an agent or a postdoc to act independently.

**Tool type:** mixed | **Primary tools:** CLIPper, PureCLIP, Skipper, STAR, umi_tools, HOMER, mCross, DEWSeq, m6Aboost, GLORI, STAMP, Hyb, RBPNet, RNAProt

## Skills

| Skill | Description |
|-------|-------------|
| clip-preprocessing | Protocol-specific UMI extraction (eCLIP 10nt, iCLIP/iCLIP2 NNNXXXXNN, PAR-CLIP), 3'-only adapter trimming (5' end is the truncation = CL -1), two-pass for eCLIP read-through, ENCODE-style retention QC |
| clip-alignment | STAR ENCODE parameter block (`--alignEndsType EndToEnd`, `--outFilterMultimapNmax 1`, `--outFilterMismatchNoverReadLmax 0.04`), PAR-CLIP override 0.07 for T->C, HISAT2 low-memory fallback, CLAM multi-mapper rescue for repeat-binding RBPs, WASP for allele-specific |
| clip-peak-calling | CLIPper + SMInput log2 normalization (ENCODE stringent: log2 FC >= 3, -log10 p >= 3), Skipper beta-binomial 100bp windows (210-320% more sites), PureCLIP HMM single-nt (focal), Piranha ZTNB, omniCLIP, CTK CIMS/CITS, CLAM repeat rescue, IDR rescue + self-consistency < 2 |
| crosslink-site-detection | Truncation-based CTK CITS for iCLIP/eCLIP; CIMS deletion for HITS-CLIP; CIMS T->C / PARalyzer for PAR-CLIP; PureCLIP HMM general; chemistry of RT stop at CL-1; off-by-one strand handling |
| clip-motif-analysis | HOMER with GC-matched background; STREME (fast successor to MEME); mCross jointly registering motif + CL position (Feng 2019); PEKA positional k-mer (no input needed); RBPamp affinity; RBNS Kd cross-validation (Dominguez 2018 78 RBPs); uracil crosslink bias correction |
| binding-site-annotation | ChIPseeker `tssRegion=c(-100,100) level='transcript'` (default is for ChIP); RBP-Maps splicing regulatory metagene (Yeo 1400nt cassette); RCAS for ncRNA-aware; RepeatMasker overlap as separate axis; metagene 5'UTR / CDS / 3'UTR; per-RBP-class expectations |
| clip-qc | preseq lc_extrap library complexity (>= 1M unique); FRiP per RBP class (>= 0.005 narrow); IDR rescue + self-consistency < 2 (ENCODE); SMInput vs IgG rationale; RSeQC read distribution metagene; antibody validation on KD lysate WB |
| differential-clip | DEWSeq sliding-window NB GLM (Schwarzl + Hentze) with interaction design `~ type + condition + type:condition`; Flipper Skipper-companion; edgeR / limma-voom on consensus peakset; KD validation by WB; aggregation of adjacent significant windows |
| m6a-clip | miCLIP2 + m6Aboost ML (eCLIP-pipeline compatible, Kortel 2021); GLORI antibody-free stoichiometric (Liu 2023, new gold standard); DART-seq APOBEC1-YTH; m6Anet nanopore direct RNA; MeRIP-seq; DRACH motif constraint; cross-method discordance |
| stamp-antibody-free | STAMP / scSTAMP (APOBEC1-RBP, C->U editing) and TRIBE / HyperTRIBE (ADAR, A->I) for antibody-free profiling; Bullseye / SAILOR / JACUSA2 edit-site detection; deaminase-only control mandatory; per-cell pseudobulk for scSTAMP |
| ago-clip-mirna-targets | Chimeric eCLIP / miR-eCLIP (30-175x miRNA enrichment, Manakov 2022); CLEAR-CLIP chimeras (Moore 2015); HEAP Halo-Ago2 mouse; Hyb pipeline with bowtie2 mode; canonical 7mer-m8 / 8mer seeds + non-canonical 3' compensatory; miRNA expression filter |
| clip-deep-learning | RBPNet sequence-to-CL distribution at single-nt (Jens 2024); RNAProt RNN classifier (AUC 87-89%); GraphProt2 GCN with structure; DeepCLIP; DeepRiPe; chromosome-split prevents leakage; variant-effect prediction; saliency + TF-MoDISco interpretation |

## Example Prompts

- "Preprocess my eCLIP from FASTQ; UMI extract 10nt from R1; cutadapt with -q 6 only; never trim 5' of R2"
- "Align my iCLIP2 with STAR ENCODE parameters, end-to-end, unique mappers only"
- "MATR3 binds LINE-1 - emit multi-mappers from STAR and rescue with CLAM"
- "Call peaks with CLIPper, normalize against SMInput, apply ENCODE stringent log2 FC >= 3 / -log10 p >= 3"
- "Use Skipper for FASTKD2 - CLIPper misses chrM peaks"
- "Detect single-nucleotide crosslink sites with PureCLIP for mCross motif registration"
- "Run mCross for CL-registered motif; compare to RBNS Kd from Dominguez 2018"
- "Annotate peaks with ChIPseeker; tssRegion=c(-100,100), level=transcript"
- "Build the RBP-Maps splicing regulatory map for PTBP1"
- "QC my library: preseq, FRiP, IDR rescue + self-consistency"
- "DEWSeq window-level differential between treatment and control with interaction term"
- "Detect m6A with miCLIP2 + m6Aboost, then triangulate with GLORI"
- "STAMP with APOBEC1-RBP and APOBEC1-only control - subtract baseline editing"
- "Chimeric eCLIP for direct miRNA-target pairs with Hyb pipeline (bowtie2 mode)"
- "RBPNet variant-effect for heterozygous SNPs in eCLIP regions"
- "Why is my CLIP duplication rate 60%? Normal for CLIP - check unique fragments instead"

## Requirements

```bash
# Preprocessing + alignment
conda install -c bioconda umi_tools cutadapt fastp star hisat2 bowtie2 samtools

# Peak calling
conda install -c bioconda clipper pureclip piranha bedtools idr
git clone https://github.com/algaebrown/skipper  # Skipper Snakemake workflow
git clone https://github.com/chaolinzhanglab/ctk  # CTK CIMS/CITS

# Motif analysis
conda install -c bioconda homer meme
git clone https://github.com/chaolinzhanglab/mCross
pip install peka

# Annotation
BiocManager::install(c('ChIPseeker','RCAS','TxDb.Hsapiens.UCSC.hg38.knownGene'))
git clone https://github.com/YeoLab/rbp-maps

# QC
conda install -c bioconda preseq picard multiqc rseqc

# Differential
BiocManager::install(c('DEWSeq','DESeq2','edgeR','limma'))
pip install htseq-clip

# m6A
pip install m6anet  # m6Aboost: github.com/ZarnackGroup/m6Aboost

# STAMP / DART
# Bullseye: github.com/mekoulnik/Bullseye
# SAILOR: github.com/YeoLab/sailor

# AGO miRNA
git clone https://github.com/gkudla/hyb

# Deep learning
pip install rbpnet rnaprot graphprot2 deepclip biopython
```

## CLIP Variant Decision Tree

| Need | Method |
|------|--------|
| Robust ENCODE comparability | eCLIP / seCLIP |
| Highest motif specificity | iCLIP2 or iiCLIP_d |
| Photoactivatable nucleoside available; HEK293/K562 | PAR-CLIP |
| Non-radioactive iCLIP | irCLIP or iCLIP3 |
| Fast (1.5 day) | FLASH |
| No antibody / no UV / in vivo | STAMP or TRIBE |
| Single-cell RBP profiling | scSTAMP or scTRIBE |
| m6A modification map (stoichiometric) | GLORI |
| m6A modification map (eCLIP-compatible) | miCLIP2 + m6Aboost |
| m6A reader (YTHDF) targets | DART-seq |
| Direct miRNA-target pairs | chimeric eCLIP / miR-eCLIP |
| In vivo mouse AGO | HEAP |
| Isoform-resolved binding | dirCLIP (long-read, 2026) or m6Anet |
| Subcellular localization-resolved | coCLIP or RBProximity-CLIP (2025) |
| Multi-RBP co-binding | Re-CLIP / irCLIP-RNP (2024) |
| Highly multiplexed RBPs | SPIDR (2023) |

## Related Skills

- **chip-seq** - DNA-protein binding analogue; peak calling and motif methods translate
- **read-qc** - General UMI handling and adapter trimming
- **read-alignment** - General STAR / bowtie2 / HISAT2
- **alignment-files** - SAM/BAM manipulation
- **genome-intervals** - BED/GTF operations for annotation
- **alternative-splicing** - Cassette-exon tables for RBP-Maps splicing regulatory analysis
- **differential-expression** - DESeq2 / edgeR engines for differential CLIP
- **epitranscriptomics** - MeRIP-seq + m6A-specific tools (cross-referenced from m6a-clip)
- **small-rna-seq** - TargetScan / miRDB integration for AGO-CLIP miRNA targets
- **single-cell** - scSTAMP / scTRIBE / coCLIP downstream
- **long-read-sequencing** - dirCLIP and m6Anet nanopore CLIP
- **causal-genomics** - Variant-effect prediction from clip-deep-learning feeds MR / fine-mapping
- **machine-learning** - Train/test methodology and saliency for clip-deep-learning
