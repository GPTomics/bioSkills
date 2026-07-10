# workflows

## Overview

End-to-end bioinformatics pipelines that orchestrate multiple skills into complete analysis workflows. Each workflow provides a primary recommended path plus alternatives, with QC checkpoints between major steps.

**Tool type:** mixed | **Primary tools:** Various (workflow-specific)

## Skills

| Skill | Description |
|-------|-------------|
| rnaseq-to-de | FASTQ to an annotated DE table: release-matched Salmon/STAR quantification, tximport gene-level collapse, raw-counts DESeq2/edgeR/limma-voom with batch-in-design, apeglm shrinkage, and VST visualization |
| fastq-to-variants | Germline FASTQ to a filtered, normalized, benchmarked VCF: reference-genome commitment, canonical step order (normalize before annotate), engine + joint-calling selection, cohort-size filtering, stratified hap.py/vcfeval benchmarking |
| chipseq-pipeline | ChIP-seq reads to blacklist-filtered, annotated, reproducible peaks: build+blacklist+effective-genome-size commitment, pre-dedup NRF/PBC, control-matched MACS3 (narrow/broad), IDR on per-replicate peaks, spike-in-aware tracks |
| scrnaseq-pipeline | 10x data to annotated cells: made-once counting commitments, ambient-before-doublet-before-normalize, per-sample QC before merge, integrate-then-cluster (never test on integrated values), pseudobulk DE + differential abundance |
| atacseq-pipeline | ATAC-seq reads to differential accessibility and footprints: no-input shift-extend background, chrM-before-calling, single Tn5 +4/-5 shift, MACS3/Genrich, Corces fixed-width consensus, DiffBind/csaw, TOBIAS |
| methylation-pipeline | Bisulfite/EM-seq to DMRs: conversion-rate gate (lambda + pUC19) before any beta, Bismark align/dedup, --no_overlap calling, M-bias trim, coverage-filtered count-model testing, selection-aware dmrseq/DSS |
| metagenomics-pipeline | Metagenomic reads to taxonomic and functional profiles, controls-first |
| expression-to-pathways | DE results to functional enrichment (GO, KEGG, Reactome, GSEA) with prokaryotic support and multi-condition comparison |
| genome-assembly-pipeline | Profile, QC, assemble (short/long/HiFi/meta), polish, decontaminate, scaffold, and three-axis QC from reads |
| longread-sv-pipeline | Long-read SV workflow: basecalling to platform-matched alignment, Sniffles2/cuteSV/pbsv calling, two-step .snf cohort merge, and Truvari benchmarking (Tier1 + CMRG) |
| gwas-pipeline | Genotypes to genome-wide associations: ancestry-matched panel + strand/allele harmonization gate, controls-only HWE, joint imputation to dosages, long-range-LD-excluded PCA, structure-driven engine (GLM/LMM/SAIGE), LDSC-intercept diagnostics |
| cnv-pipeline | BAM to segmented, integer-called CNVs: germline-vs-somatic fork (CNVkit / GATK gCNV), build+BED+assay-matched-PoN commitment, reference-before-segmentation, purity/ploidy-before-integer-calls, center-before-GISTIC2 |
| spatial-pipeline | Spatial transcriptomics end-to-end: platform-class fork first (imaging-segment vs sequencing-deconvolve), physical-space neighbors, FDR-gated SVGs, real domain methods, handing off to deconvolution and communication |
| hic-pipeline | Hi-C FASTQ to compartments, TADs, and loops: depth-dictates-feature, independent -SP5M pairing, long-range-cis QC, mask-before-ICE-balancing, GC-phased eigenvector, deep-map-only loop calling, HiChIP/PLAC routing |
| multiome-pipeline | Joint scRNA + scATAC (Cell Ranger ARC) to a WNN-integrated annotated object: shared-barcode intersection, per-modality QC + AMULET doublets before WNN, drop depth-correlated LSI, RNA-based annotation, LinkPeaks |
| somatic-variant-pipeline | End-to-end tumor-normal somatic SNV/indel + SV/CNV pipeline (Mutect2/Strelka2) with PoN, contamination, orientation-bias filtering, and AMP/ASCO/CAP tier interpretation |
| proteomics-pipeline | MaxQuant/FragPipe/DIA-NN to differential protein abundance: search-DB + DDA/DIA commitment, FDR at PSM+peptide+protein-group (not just PSM), contaminant-removal-before-normalize, TMT IRS bridge, MNAR modeled not downshift-imputed, limma/DEqMS/MSstats |
| microbiome-pipeline | 16S amplicon to differential taxa with DADA2 and ALDEx2 |
| crispr-screen-pipeline | FASTQ to hit genes: library+control-class + baseline (plasmid/Day-0/vehicle) commitment, six-stage QC with CEGv2-PR-AUC positive control, copy-number-correction-before-hit-calling (cancer lines), design-matched MAGeCK/BAGEL2/drugZ/Chronos, tier consensus |
| metabolomics-pipeline | Raw MS to differential metabolites via XCMS and pathway mapping |
| imc-pipeline | Imaging mass cytometry MCD to patient-level spatial analysis: panel+segmentation+pixel-size commitment, channel-spillover-on-pixels-before-segmentation vs REDSEA-lateral-after, arcsinh-cofactor-1-not-5, aggregate-to-patient (cells/ROIs are not replicates) |
| cytometry-pipeline | FCS files to differential populations via CATALYST/diffcyt |
| multi-omics-pipeline | Vertical bulk multi-omics integration (MOFA2/DIABLO/SNF): vertical-not-horizontal sample key, per-block-normalize + per-view-variance-equalize before stacking, correct-batch-once, held-out validation (in-cohort CV at n<<p is biased) |
| tcr-pipeline | TCR/BCR repertoire from FASTQ to clonotype diversity |
| smrna-pipeline | Small RNA-seq from FASTQ to differential miRNAs: kit-aware trimming, miRge3 quantification or miRDeep2 discovery, compositionally-aware DESeq2 (normalizer as the dominant choice, size-factor inspection), and expression-filtered targets |
| riboseq-pipeline | Ribo-seq FASTQ to ORF detection and translation efficiency: harvest-method commitment (CHX invalidates dwell-time), UMI-only dedup, 3-nt periodicity hard gate, per-length P-site offsets, count-GLM TE |
| merip-pipeline | MeRIP-seq FASTQ to differential m6A peaks: NO dedup for non-UMI (duplicates are signal), IP-vs-Input enrichment, exomePeak2 four-BAM differential, DRACH set-level sanity, Guitar stop-codon go/no-go |
| clip-pipeline | CLIP-seq (eCLIP/iCLIP/iCLIP2/iCLIP3/irCLIP/PAR-CLIP) from FASTQ to ENCODE-stringent binding sites (log2 FC >= 3 AND -log10 p >= 3), single-nt crosslink maps, ChIPseeker annotation, motif registration (HOMER + mCross), with optional DEWSeq differential binding |
| neoantigen-pipeline | Somatic variants to ranked vaccine candidates: binding-is-single-digit-PPV, full-resolution HLA + LOHHLA gating, proximal-variant phasing (--phased-proximal-variants-vcf), CCF-not-raw-VAF clonality, rank-within-patient |
| outbreak-pipeline | Pathogen isolates to transmission networks: bacterial (snippy/Gubbins/IQ-TREE/TreeTime/TransPhylo) vs viral (Nextstrain) fork, one-reference commitment, mandatory recombination-masking, TempEst temporal-signal gate, pathogen-specific cluster thresholds |
| crispr-editing-pipeline | Target to CRISPR constructs with branching strategies |
| metabolic-modeling-pipeline | Protein FASTA to flux predictions: reconstruction-tool locks the identifier namespace, medium-IS-the-model (set before FBA), biomass-objective drives every flux, MEMOTE-is-well-formedness-not-correctness, iterative curation with energy-generating-cycle removal |
| biomarker-pipeline | End-to-end biomarker discovery from expression to validated panels: group-aware (subject-level) splitting, in-fold feature selection (Boruta/LASSO), leakage-safe/nested CV, SHAP audit, and discrimination + calibration validation |
| splicing-pipeline | Bulk short-read alternative splicing from FASTQ to differential splicing: cohort-consistent STAR 2-pass (shared junction DB), event-level rMATS/leafcutter and isoform-level DTU (dtuScaledTPM + stageR), reconciled, with sashimi plots |
| liquid-biopsy-pipeline | cfDNA plasma-to-monitoring: pre-analytics-is-the-sensitivity-ceiling, UMI/duplex error-suppression before calling, tumor-naive-vs-informed fork, ichorCNA TF (~3% floor) / VarDict low-VAF, CHIP subtraction vs matched WBC, VAF-needs-input-GE |
| genome-annotation-pipeline | Assembled contigs to functional annotation: pro-vs-eukaryotic fork + genetic-code-table-from-taxonomy, annotate-only-decontaminated-QC-passed-assembly (CheckM2 gate), RNA-seq+protein evidence for BRAKER3, tool+DB-version pinning |
| grn-pipeline | Single-cell data to regulons + in-silico perturbation (pySCENIC/SCENIC+/CellOracle): GRNs are undirected-by-default (report the evidence tier), species/DB-vintage matching, raw-counts-of-cleaned-cells, ctx-buys-directionality |
| causal-genomics-pipeline | GWAS summary statistics to triangulated causal inference via MR (with CHP-aware sensitivity), colocalization, fine-mapping, mediation, TWAS, cis-pQTL drug-target MR, effector-gene prioritization, heritability partitioning, genetic correlation, and GenomicSEM common-factor GWAS |
| timecourse-pipeline | Expression matrix to temporal gene modules: temporal DE (limma-splines/DESeq2-LRT), z-scored soft clustering (Mfuzz/tslearn), GAM trajectories, temporal-gene-background enrichment, and a design-gated circadian rhythm branch |
| edna-pipeline | eDNA amplicons to community ecology via OBITools3/DADA2, iNEXT, and vegan |
| clinical-trial-pipeline | CDISC data to ICH E9(R1) estimand-driven analysis: FDA 2023 marginal-vs-conditional logistic, MMRM/reference-based MI, modern HTE subgroup methods, graphical multiplicity, survival, CONSORT 2025 reporting |

## Example Prompts

- "I have FASTQ files, how do I find differentially expressed genes?"
- "Run the complete RNA-seq pipeline from raw reads to DE results"
- "Process my ChIP-seq data from FASTQ to annotated peaks"
- "Analyze my 10X single-cell data end to end"
- "Call variants from my whole genome sequencing data"
- "Find structural variants from my Nanopore reads"
- "Run GWAS on my case-control study"
- "Detect CNVs from my exome sequencing"
- "Analyze my Visium spatial transcriptomics data"
- "Process my Hi-C data to find TADs and loops"
- "Analyze my CRISPR screen from FASTQ to hit genes"
- "Run metabolomics analysis from raw MS data to pathways"
- "Process my imaging mass cytometry data with spatial analysis"
- "Analyze my flow cytometry data end to end"
- "Integrate my RNA-seq, proteomics, and metabolomics data"
- "Run the TCR repertoire pipeline from FASTQ to diversity"
- "Analyze my small RNA-seq for differential miRNAs"
- "Process my Ribo-seq to translation efficiency"
- "Run m6A analysis from MeRIP-seq data"
- "Find RBP binding sites from my CLIP-seq data"
- "Find neoantigens from my somatic VCF for vaccine design"
- "Investigate this outbreak with genomic data"
- "Design CRISPR guides to knock out my target gene"
- "Build a metabolic model from my genome annotation"
- "Analyze differential splicing between my conditions"
- "Estimate tumor fraction from my plasma cfDNA"
- "Run a complete liquid biopsy pipeline for my samples"
- "Annotate my newly assembled genome from scratch"
- "Build gene regulatory networks from my single-cell data"
- "Run post-GWAS causal inference on my summary statistics"
- "Analyze my time-course expression experiment end to end"
- "Process my eDNA water samples through the full biodiversity pipeline"
- "Analyze my clinical trial data from CDISC files to odds ratios and forest plots"

## Requirements

Requirements vary by workflow. See individual skill files for specific dependencies.

```bash
# Common tools
conda install -c bioconda samtools bcftools bwa-mem2 star salmon fastp

# R/Bioconductor
BiocManager::install(c('DESeq2', 'Seurat', 'clusterProfiler'))
```

## Related Skills

- **database-access** - Download public data (SRA / GEO / NCBI Datasets / Ensembl / UniProt / interaction DBs) before running any pipeline that starts from public records
- **read-qc** - Quality control and preprocessing (first step in most workflows)
- **read-alignment** - Alignment tools used by many workflows
- **differential-expression** - DE analysis details
- **single-cell** - Single-cell analysis details
- **variant-calling** - Variant calling details
- **alternative-splicing** - Splicing analysis skills
- **liquid-biopsy** - cfDNA analysis skills
- **genome-annotation** - Genome annotation skills
- **gene-regulatory-networks** - GRN inference skills
- **causal-genomics** - Causal inference from GWAS
- **temporal-genomics** - Circadian and differential rhythmicity, periodicity, temporal clustering, trajectory modeling, dynamic GRN
- **ecological-genomics** - eDNA metabarcoding, biodiversity metrics, community ecology
- **clinical-biostatistics** - Clinical trial statistical methods
