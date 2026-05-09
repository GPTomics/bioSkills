# atac-seq

## Overview

Decision-grade ATAC-seq workflows: peak calling (MACS, Genrich, HMMRATAC), ENCODE 4 quality control with citation-anchored thresholds, fixed-width consensus peakset construction (Corces 2018 iterative overlap), differential accessibility across replicate-count regimes (DiffBind, csaw, DESeq2, edgeR), TF footprinting with Tn5 bias correction (TOBIAS, HINT-ATAC, scprinter), nucleosome positioning, motif accessibility variability (chromVAR), single-cell ATAC (Signac/ArchR/SnapATAC2), and cis-regulatory co-accessibility (Cicero, SCENIC+).

**Tool type:** mixed | **Primary tools:** MACS3, DiffBind, csaw, chromVAR, TOBIAS, scprinter, NucleoATAC, ArchR, Signac, SnapATAC2, Cicero, ABC, ENCODE-rE2G, chromBPNet, BPNet, scBasset, EnFormer, WASP, GATK ASEReadCounter, RASQUAL

## Skills

| Skill | Description |
|-------|-------------|
| atac-peak-calling | MACS / Genrich / HMMRATAC peak calling; ENCODE 4 IDR; pseudoreplicate self-consistency; chromap shift; super-enhancer ROSE/LILY |
| atac-qc | ENCODE 4 thresholds (TSS enrichment, FRiP, NRF/PBC1/PBC2, mt fraction); fragment-size periodicity; preseq lc_extrap; sex-chr / cell-cycle / spike-in QC |
| consensus-peakset | Corces 2018 iterative overlap; DiffBind summits=250; per-condition union; cross-build liftOver; ENCODE-rE2G cCRE |
| differential-accessibility | DiffBind / csaw / DESeq2 / edgeR by replicate count; normalization choice; spike-in (Reske 2020); permutation; Hi-C-loop-anchored |
| footprinting | TOBIAS three-step + Tn5 +4/-5 dual-cut; bias correction alternatives (chromBPNet, seqOutBias, naked-DNA); per-TF failure modes |
| motif-deviation | chromVAR matched-background z-scores; scBasset / EnFormer alternatives; DecoupleR multi-method TF activity |
| nucleosome-positioning | NucleoATAC / DANPOS3 / ATACseqQC; V-plot interpretation; +1 calling; H2A.Z detection; Fiber-seq long-read |
| single-cell-atac | 10X scATAC + Multiome (cellranger-arc); Signac / ArchR / SnapATAC2; AMULET; cell-cycle + sex-chr QC; scArches |
| co-accessibility | Cicero / ArchR getCoAccessibility; alpha math; HiChIP integration; Hi-C concordance; (ABC: see enhancer-gene-linking) |
| deep-learning-atac | chromBPNet, BPNet, scBasset, EnFormer; in silico variant effect at GWAS SNPs; DeepLIFT + TF-MoDISco motif discovery |
| enhancer-gene-linking | ABC (Fulco 2019, Nasser 2021); ENCODE-rE2G; HiChIP H3K27ac; CRISPRi-FlowFISH validation framework |
| allele-specific-accessibility | WASP reference-bias correction; GATK ASEReadCounter; RASQUAL joint caQTL; per-peak ASE aggregation |

## Example Prompts

- "Call ATAC-seq peaks following the ENCODE 4 pipeline with IDR across replicates"
- "Compute TSS enrichment using the ENCODE pyTSSe formula and grade against ENCODE 4 thresholds"
- "Build a Corces 2018 iterative-overlap consensus peakset (501 bp fixed-width)"
- "Run DiffBind with edgeR backend because I only have 2 reps per condition"
- "Switch DiffBind normalization to library-size because the treatment globally compacts chromatin"
- "Run TOBIAS three-step footprinting; verify CTCF aggregate shows clean V-shape"
- "Compare TOBIAS and HINT-ATAC predictions; report two-tool concordance"
- "Run chromVAR with matched-background correction; report top variable motifs and limma differential"
- "Diagnose flat fragment-size distribution -- over-transposition vs degraded chromatin?"
- "Generate V-plot at TSSs to verify nucleosome positioning is recoverable"
- "Process 10X scATAC with Signac (TF-IDF + LSI dims 2:30 to skip depth)"
- "Use ArchR for a 200K-cell dataset because Signac will OOM"
- "Run AMULET doublet detection plus ArchR doublet score; intersection is high-confidence"
- "Run Cicero on Signac scATAC to get peak-pair co-accessibility; filter to coaccess > 0.25"
- "Compare Cicero connections against published Hi-C loops to estimate concordance"
- "Run SCENIC+ on Multiome data for TF -> enhancer -> gene networks"
- "Score 100 GWAS SNPs for chromatin effects with pre-trained chromBPNet K562 model"
- "Run TF-MoDISco on chromBPNet DeepLIFT contributions to discover de novo motifs"
- "Run ABC pipeline on K562 ATAC + H3K27ac + Micro-C; threshold ABC.Score >= 0.02"
- "Use WASP to correct reference-allele mapping bias before any allele-specific analysis"
- "Run GATK ASEReadCounter on WASP-filtered BAM and aggregate ASE within peaks"
- "Run RASQUAL joint cis-caQTL on a 50-individual cohort for max statistical power"
- "Cross-validate chromBPNet variant effect predictions against observed allelic imbalance"
- "Detect H2A.Z-containing nucleosomes from ATAC fragment-size shifts at H2A.Z ChIP peaks"
- "Combine ABC + ENCODE-rE2G + HiChIP -- intersection is high-confidence enhancer-gene calls"

## Requirements

```bash
# Core CLI
conda install -c bioconda macs3 macs2 genrich samtools bedtools idr \
    deeptools picard tobias rgt-hint multiqc subread bedops chromap preseq fithichip

# Single-cell (Python)
pip install snapatac2 amulet-py scprinter scenicplus pycistopic scvi-tools

# QC and visualization
pip install pysam pyBigWig numpy pandas matplotlib pybedtools

# Deep learning (chromBPNet, BPNet, scBasset, tangermeme, TF-MoDISco)
pip install tensorflow torch chrombpnet bpnet-lite tangermeme tfmodisco-lite captum kipoi

# Allele-specific (WASP, GATK, RASQUAL)
git clone https://github.com/bmvdgeijn/WASP
conda install -c bioconda gatk4 bcftools shapeit5
git clone https://github.com/natsuhiko/rasqual

# Enhancer-gene linking (ABC, ENCODE-rE2G)
git clone https://github.com/broadinstitute/ABC-Enhancer-Gene-Prediction
git clone https://github.com/EngreitzLab/ENCODE_rE2G
```

```r
BiocManager::install(c(
    'DiffBind', 'DESeq2', 'edgeR', 'csaw', 'limma',
    'ChIPseeker', 'sva', 'RUVSeq',
    'chromVAR', 'motifmatchr', 'JASPAR2024', 'TFBSTools',
    'ATACseqQC', 'GenomicRanges', 'GenomicAlignments', 'GenomicInteractions',
    'BSgenome.Hsapiens.UCSC.hg38', 'TxDb.Hsapiens.UCSC.hg38.knownGene',
    'EnsDb.Hsapiens.v86', 'Signac', 'Seurat', 'cicero', 'monocle3', 'scDblFinder'
))
remotes::install_github('GreenleafLab/ArchR', ref='master', repos=BiocManager::repositories())
```

External resources:
- ENCODE blacklist v2: github.com/Boyle-Lab/Blacklist
- JASPAR 2024 CORE: jaspar.genereg.net
- HOCOMOCO v12: hocomoco12.autosome.ru

## Related Skills

- **read-alignment** - Upstream ATAC alignment (bowtie2, bwa-mem2, chromap)
- **alignment-files** - BAM preprocessing (dedup, MAPQ filter, chrM strip)
- **chip-seq** - Same DiffBind/MACS frameworks; ChIP uses input control; super-enhancers (ROSE)
- **hi-c-analysis** - 3D contact validation for co-accessibility; ABC contact input
- **single-cell** - General sc patterns; scatac-analysis cross-reference
- **gene-regulatory-networks** - SCENIC / SCENIC+ for TF -> target inference
- **genome-intervals** - BED/narrowPeak file operations
- **pathway-analysis** - GO enrichment of differential accessibility gene lists
- **data-visualization** - Genome tracks, heatmaps for peak signals
- **causal-genomics** - Fine-mapping with caQTL and chromBPNet variant effects
- **variant-calling** - VCF inputs for allele-specific accessibility
- **phasing-imputation** - Phasing required before WASP allele-specific analysis
- **machine-learning** - General ML patterns for deep-learning-atac models
