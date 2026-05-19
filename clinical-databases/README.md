# clinical-databases

## Overview

Decision-grade clinical and population genetics database skills covering variant interpretation, ACMG/AMP classification, pharmacogenomics, polygenic risk, somatic signatures, TMB / MSI ICI biomarkers, and HLA typing. Skills encode the ClinGen SVI 2020+ specifications, Tavtigian Bayesian point system, Pejaver 2022 PP3/BP4 calibration, Vega 2021 TMB harmonization, Caudle 2020 CYP2D6 activity scores, COSMIC v3.4 signatures, and current FDA / ESMO / NCCN guidance.

**Tool type:** mixed | **Primary tools:** myvariant, requests, cyvcf2, PharmCAT, Cyrius, OptiType, T1K, HLA-LA, SigProfilerAssignment, MSIsensor-pro, LDpred2, PRS-CSx, pgsc_calc, InterVar, GeneBe

## Skills

| Skill | Description |
|-------|-------------|
| myvariant-queries | BioThings aggregator for ClinVar / gnomAD / dbSNP / dbNSFP / COSMIC / CIViC with _meta versioning |
| clinvar-lookup | VCV / SCV / RCV hierarchy, 2024 XML schema, ClinGen VCEPs, star ratings, ClinGen Allele Registry CA IDs |
| dbsnp-queries | Build 156 JSON, RsMergeArch multi-hop merge chains, SPDI canonical, ALFA vs gnomAD frequencies |
| gnomad-frequencies | v4 / v3 / v2.1.1 decision matrix, grpmax FAF95, LOEUF constraint, MID ancestry, Whiffin max-credible-AF |
| acmg-classification | Tavtigian point system, PVS1 decision tree, Pejaver REVEL/AlphaMissense calibration, Brnich OddsPath, Walker splicing, cancer AMP Tier I-IV |
| variant-prioritization | Trio rare-disease pipeline: DeNovoGear, WhatsHap compound het, Exomiser hiPHIVE, ClinGen gene-disease validity, ACMG SF v3.2 (81 genes) |
| pharmacogenomics | CPIC vs DPWG, PharmVar, Caudle 2020 CYP2D6 activity (*10 = 0.25), Cyrius for CYP2D6 SVs, DPYD 2024 AS, PharmCAT, HLA-drug |
| polygenic-risk | LDpred2-auto, SBayesRC, MegaPRS, PRS-CSx, PROSPER, MUSSEL, BridgePRS, JointPRS, PRSmix, Ding 2023 continuous ancestry, Hingorani 2023 critique |
| somatic-signatures | COSMIC v3.4 (SBS40a/b/c split; SBS17b 5-FU; SBS10a-d POLE/POLD1), MuSiCal mvNMF, HRDetect, FFPE-as-SBS30 correction, Petljak A3A dominance |
| tumor-mutational-burden | Vega 2021 calibration (FoundationOne 0.8 Mb scored; TSO500 7.8; Oncomine 8.4), FDA pembrolizumab pan-tumor, McGrail 2021 tumor-type exclusions, HLA-LOH |
| msi-detection | MSIsensor-pro tumor-only, MSIsensor-ct cfDNA, MANTIS, Lynch syndrome workflow (IHC + MSI + MLH1 methylation), POLE-exo vs MMR-D disambiguation |
| hla-typing | T1K (class I + II + KIR), HLA-LA, OptiType, arcasHLA (RNA-seq), StarPhase / FuFiHLA long-read, HIBAG SNP-array imputation, CIWD v3.0.0, TCE3 core/non-core |

## Example Prompts

- "Resolve these HGVS variants to ClinGen Allele Registry CA IDs and pull ClinVar + gnomAD v4 grpmax FAF95"
- "Apply ACMG/AMP classification with Tavtigian point system + Pejaver 2022 REVEL calibration; check VCEP CSpec for this gene"
- "Calculate CYP2D6 activity score for *4xN/*10 using Caudle 2020 (note *4xN is clinically silent)"
- "Run PharmCAT on this VCF + Cyrius for CYP2D6 SVs; generate CPIC-compliant report"
- "Compute LDpred2-auto PRS for CAD with sample-overlap detection via EraSOR; apply Hingorani 2023 absolute-risk transform"
- "Extract de novo SBS signatures from this 100-sample WGS cohort with SigProfilerExtractor stability gates"
- "Calculate TMB on TSO500 panel; apply Vega 2021 calibration equivalent threshold (7.8/Mb)"
- "Run MSIsensor-pro tumor-only; apply Lynch syndrome workflow if MSI-H + retained IHC"
- "Type HLA from WGS with T1K (class I + II + KIR); screen B*57:01 for abacavir, B*15:02 for carbamazepine"
- "Filter trio exome to candidate Mendelian variants: rare grpmax FAF95 < 0.0001 + de novo + Exomiser hiPHIVE ranking"
- "Cross-check ACMG SF v3.2 (Miller 2023; 81 genes including new CALM1/2/3) for incidental findings"
- "Reconcile MSI-H + TMB-H + HLA-LOH for ICI eligibility decision"
- "For African-ancestry warfarin patient, apply IWPC algorithm explicitly INCLUDING CYP2C9 *5/*6/*8/*11 (COAG paradigmatic ancestry-algorithm failure)"
- "Apply Walker 2023 ClinGen Splicing Subgroup: SpliceAI DS_max + SpliceVault aberrant transcript prediction"
- "Distinguish POLE-exo hypermutator (SBS10a/10b; typically MSI-stable) from MMR-D (SBS6/15/26/44; MSI-H)"

## Requirements

```bash
pip install myvariant requests pandas cyvcf2 hail SigProfilerMatrixGenerator SigProfilerExtractor SigProfilerAssignment

# PharmCAT (Java)
# Download: https://pharmcat.org

# CYP2D6 SV calling
pip install cyrius

# MSI detection
conda install -c bioconda msisensor-pro msisensor

# HLA typing
conda install -c bioconda optitype t1k samtools

# Variant prioritization tools
# Exomiser, DeNovoGear, WhatsHap from bioconda
conda install -c bioconda denovogear whatshap exomiser

# Polygenic risk (R + CLI)
# install.packages('bigsnpr')  # LDpred2
# pgsc_calc nf-core: nextflow pull pgscatalog/pgsc_calc
# PRS-CS / PRS-CSx: github.com/getian107/PRS{,Csx}
# SBayesRC / GCTB: cnsgenomics.com/software/gctb/

# Variant classification
# InterVar: github.com/WGLab/InterVar
# GeneBe API: genebe.net
```

## Related Skills

- **variant-calling** - Upstream variant calling and annotation
- **variant-calling/clinical-interpretation** - Clinical reporting workflow
- **database-access** - General Entrez query patterns
- **causal-genomics** - Mendelian randomization (uses gnomAD + ClinVar)
- **chemoinformatics/admet-prediction** - Drug metabolism prediction
- **immunoinformatics/mhc-binding-prediction** - Downstream neoantigen prediction from HLA + TMB
- **machine-learning/biomarker-discovery** - PRS + TMB + MSI in biomarker pipelines
- **population-genetics/population-structure** - Ancestry inference for PRS / HLA imputation
- **workflows/neoantigen-pipeline** - End-to-end neoantigen workflow
