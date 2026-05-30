# copy-number

## Overview

Detect, interpret, and characterize copy number variants from whole-genome, exome, and targeted sequencing. Covers read-depth calling, copy-ratio segmentation, allele-specific copy number with tumor purity and ploidy, recurrent and driver CNV, clinical interpretation, homologous recombination deficiency, focal amplification architecture, and subclonal tumor evolution.

**Tool type:** mixed | **Primary tools:** CNVkit, GATK4, ASCAT, Sequenza, FACETS, GISTIC2, scarHRD, AmpliconArchitect

## Skills

| Skill | Description |
|-------|-------------|
| cnvkit-analysis | Read-depth CNV calling from panels, exomes, and WGS with CNVkit |
| gatk-cnv | GATK somatic CNV and germline gCNV best-practices workflows |
| copy-ratio-segmentation | CBS/HMM/HaarSeg segmentation and depth-bias correction |
| allele-specific-copy-number | Integer allele-specific CN, purity, ploidy (ASCAT/Sequenza/FACETS/PURPLE/PureCN) |
| recurrent-cnv | Cohort recurrent and driver CNV with GISTIC2; copy-number signatures |
| cnv-annotation | Gene, dosage-sensitivity, and clinical-database annotation of CNVs |
| germline-cnv-interpretation | ACMG/ClinGen points-based constitutional CNV classification |
| hrd-scoring | Homologous recombination deficiency from LOH/LST/TAI genomic scars |
| focal-amplification-ecdna | ecDNA vs BFB vs HSR amplicon architecture with AmpliconArchitect |
| subclonal-copy-number | Subclonal CN, whole-genome doubling, and copy-number evolution |
| cnv-visualization | Profile, allele-specific, and cohort CNV figures |

## Example Prompts

- "Call CNVs from my tumor-normal exome pair and decide whether I need an allele-specific caller"
- "Estimate tumor purity, ploidy, and integer allele-specific copy number with FACETS"
- "Choose a segmentation algorithm for my shallow whole-genome data"
- "Find recurrently amplified regions across my cohort with GISTIC2"
- "Classify this constitutional CNV with the ACMG/ClinGen points framework"
- "Compute an HRD score from copy-number scars for PARP-inhibitor eligibility"
- "Determine whether this focal amplification is ecDNA or a chromosomal HSR"
- "Detect whole-genome doubling and call subclonal copy number for this tumor"
- "Annotate my CNV calls with affected genes and dosage-sensitivity scores"
- "Diagnose why my tumor-only CNV profile looks noisy and recurrent across samples"

## Requirements

```bash
# Read-depth calling
conda install -c bioconda cnvkit gatk4 bedtools samtools
# Allele-specific copy number
conda install -c bioconda snp-pileup sequenza-utils
R -e "install.packages(c('facets')); BiocManager::install(c('DNAcopy', 'PureCN'))"
R -e "remotes::install_github(c('VanLoo-lab/ascat/ASCAT', 'ShixiangWang/copynumber'))"
R -e "install.packages('sequenza')"
# Recurrent CNV, HRD, subclonal, annotation
# GISTIC 2.0: download the compiled binary + MATLAB Compiler Runtime from the Broad
R -e "remotes::install_github('sztup/scarHRD'); BiocManager::install(c('TitanCNA', 'CINSignatureQuantification'))"
R -e "remotes::install_github('Wedge-lab/battenberg')"
conda install -c bioconda annotsv ampliconsuite
pip install pandas pybedtools matplotlib seaborn SigProfilerAssignment
# ClassifyCNV is GitHub-only: git clone https://github.com/Genotek/ClassifyCNV
```

## Related Skills

- **alignment-files** - Input BAM processing and QC before CNV calling
- **variant-calling** - SNV/indel calling; allelic counts and clonal-mutation cross-checks
- **clinical-databases** - ClinVar, gnomAD, somatic signatures, variant prioritization
- **long-read-sequencing** - Structural-variant and complex-amplicon resolution from long reads
- **liquid-biopsy** - Low-pass and cfDNA copy number, tumor fraction estimation
- **single-cell** - Single-cell CNV inference from scRNA-seq and scDNA-seq
