# workflows

## Overview

End-to-end bioinformatics pipelines that orchestrate multiple skills into complete analysis workflows. Each workflow provides a primary recommended path plus alternatives, with QC checkpoints between major steps.

**Tool type:** mixed | **Primary tools:** Various (workflow-specific)

## Skills

| Skill | Description |
|-------|-------------|
| rnaseq-to-de | FASTQ to differential expression via Salmon/STAR and DESeq2 |
| fastq-to-variants | DNA sequencing to variant calls via BWA and bcftools/GATK |
| chipseq-pipeline | ChIP-seq reads to annotated peaks |
| scrnaseq-pipeline | 10X data to clustered and annotated cells |
| atacseq-pipeline | ATAC-seq reads to accessibility analysis |
| methylation-pipeline | Bisulfite-seq to differentially methylated regions |
| metagenomics-pipeline | Metagenomic reads to taxonomic profiles |
| expression-to-pathways | DE results to functional enrichment |
| genome-assembly-pipeline | Reads to polished assembly with QC |
| longread-sv-pipeline | Long reads to structural variants |
| gwas-pipeline | VCF to genome-wide associations |
| cnv-pipeline | BAM to copy number variants |
| spatial-pipeline | Spatial transcriptomics end-to-end |
| hic-pipeline | Hi-C data to compartments, TADs, and loops |
| multiome-pipeline | Joint scRNA + scATAC analysis |
| somatic-variant-pipeline | Tumor-normal somatic calling with Mutect2/Strelka2 |
| proteomics-pipeline | MaxQuant to differential protein abundance with MSstats/limma |
| microbiome-pipeline | 16S amplicon to differential taxa with DADA2 and ALDEx2 |
| crispr-screen-pipeline | FASTQ to hit genes via MAGeCK counting and analysis |
| metabolomics-pipeline | Raw MS to differential metabolites via XCMS and pathway mapping |
| imc-pipeline | Imaging mass cytometry to spatial cell analysis |
| cytometry-pipeline | FCS files to differential populations via CATALYST/diffcyt |
| multi-omics-pipeline | Multi-omics integration via MOFA2/mixOmics |
| tcr-pipeline | TCR/BCR repertoire from FASTQ to clonotype diversity |
| smrna-pipeline | Small RNA-seq from FASTQ to differential miRNAs |
| riboseq-pipeline | Ribo-seq from FASTQ to translation efficiency |
| merip-pipeline | MeRIP-seq from FASTQ to m6A peaks |
| clip-pipeline | CLIP-seq from FASTQ to binding sites and motifs |
| neoantigen-pipeline | Somatic variants to ranked vaccine candidates |
| outbreak-pipeline | Pathogen isolates to transmission networks |
| crispr-editing-pipeline | Target to CRISPR constructs with branching strategies |
| metabolic-modeling-pipeline | Genome to flux predictions with iterative curation |
| biomarker-pipeline | End-to-end biomarker discovery from expression to validated panels |
| splicing-pipeline | Alternative splicing from FASTQ to differential splicing with sashimi plots |
| liquid-biopsy-pipeline | cfDNA analysis for tumor fraction estimation and mutation detection |

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

## Requirements

Requirements vary by workflow. See individual skill files for specific dependencies.

```bash
# Common tools
conda install -c bioconda samtools bcftools bwa-mem2 star salmon fastp

# R/Bioconductor
BiocManager::install(c('DESeq2', 'Seurat', 'clusterProfiler'))
```

## Related Skills

- **read-qc** - Quality control and preprocessing (first step in most workflows)
- **read-alignment** - Alignment tools used by many workflows
- **differential-expression** - DE analysis details
- **single-cell** - Single-cell analysis details
- **variant-calling** - Variant calling details
- **alternative-splicing** - Splicing analysis skills
- **liquid-biopsy** - cfDNA analysis skills
