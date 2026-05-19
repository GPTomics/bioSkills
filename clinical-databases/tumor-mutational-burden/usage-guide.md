# Tumor Mutational Burden - Usage Guide

## Overview

Calculate tumor mutational burden from somatic VCFs with Friends of Cancer Research (Vega 2021) harmonization equations for per-assay calibration to the FDA pembrolizumab 10 mut/Mb pan-tumor threshold (TSO500 = 7.8, Oncomine = 8.4). Covers FoundationOne CDx convention of including synonymous variants, hypermutator / ultra-hypermutator tiering (POLE+MMR), tumor-type-specific cutoffs per McGrail 2021 + ESMO 2024 (NOT endorsed for breast / prostate / glioma), blood TMB caveats (BFAST Cohort C failed), and integration with MSI-H, HLA-LOH (LOHHLA / DASH), and neoantigen quality (Luksza fitness, McGranahan clonality).

## Prerequisites

```bash
pip install cyvcf2 pandas numpy

# VEP for annotation (or snpEff / Funcotator)
# Ensembl VEP: docker pull ensemblorg/ensembl-vep
# Or: conda install -c bioconda ensembl-vep

# HLA-LOH
# LOHHLA: https://github.com/mskcc/lohhla
# DASH: https://github.com/RamratnamLab/DASH

# Friends of Cancer Research TMB harmonization tool
# https://friendsofcancerresearch.org/tmb/
```

## Quick Start

Tell the agent what to do:
- "Calculate TMB from this FoundationOne-style somatic VCF using 0.8 Mb scored region (NOT 1.1 Mb panel total)"
- "Apply Vega 2021 cross-panel calibration: report TMB-H equivalent thresholds for TSO500 (7.8), Oncomine (8.4)"
- "For this MSK-IMPACT cohort, compute TMB excluding synonymous (vs FoundationOne which includes them); document convention"
- "Integrate TMB + MSI-H + HLA-LOH (LOHHLA) for ICI eligibility decision"
- "For breast / prostate / glioma cohorts, apply tumor-type-specific Samstein 2019 cutoff (NOT universal 10/Mb per ESMO 2024)"

## Example Prompts

### Standard TMB Calculation

> "Compute TMB from this VEP-annotated Mutect2 VCF. Use FoundationOne CDx conventions: scored region 0.8 Mb, include synonymous, VAF >= 5%, exclude hotspots, gnomAD grpmax FAF95 < 0.5%."

> "For my TSO500 panel cohort, compute TMB with 1.3 Mb scored region; apply Vega 2021 calibration equivalent threshold of 7.8/Mb for TMB-H classification."

### Cross-Panel Comparison

> "Reconcile TMB values across FoundationOne CDx (0.8 Mb scored, includes synonymous), TSO500 (1.3 Mb scored), MSK-IMPACT v4 (1.22 Mb scored). Document counting conventions."

> "Audit these published TMB values: are they comparable across studies using different panels?"

### Cancer-Type-Specific Cutoffs

> "For this breast cancer cohort, do NOT apply universal 10/Mb cutoff. Use tumor-type-specific cutoff per Samstein 2019 (top 20% within tumor type). Cite ESMO 2024 + McGrail 2021 for breast / prostate / glioma exclusion."

### MSI + HLA-LOH Integration

> "For each tumor: compute TMB, run MSI testing, run LOHHLA. Report ICI eligibility with: MSI-H supersedes TMB-H; TMB-H + HLA-LOH-positive has reduced neoantigen presentation."

### Blood TMB

> "Compute bTMB from cfDNA panel; flag low tumor fraction (BFAST Cohort C failure mechanism); recommend tissue TMB confirmation."

### Hypermutator Characterization

> "For samples with TMB > 100 mut/Mb, identify mechanism: MMR-D (SBS6/15/26/44 + ID1/2) vs POLE-exo (SBS10a/10b). For TMB > 500 mut/Mb, suspect POLE+MMR concurrent."

## What the Agent Will Do

1. Verify VCF has VEP / snpEff / Funcotator consequence annotations.
2. Apply filters: PASS-only, VAF >= 5% (FoundationOne) or higher, depth >= 100, tumor-only germline filter gnomAD grpmax FAF95 <= 0.5%, exclude COSMIC hotspots.
3. Count nonsynonymous coding variants (configurable to include synonymous for FoundationOne CDx compatibility).
4. Divide by SCORED region in Mb (NOT total panel size).
5. Apply Vega 2021 calibration: assay-specific equivalent threshold for FoundationOne 10/Mb sensitivity.
6. Classify: TMB-low; TMB-H (>= cutoff); hypermutator (>= 100); ultra-hypermutator (>= 500).
7. For breast / prostate / glioma, apply tumor-type-specific Samstein 2019 cutoff and document ESMO 2024 exclusion.
8. Integrate MSI status: MSI-H is the primary biomarker (Sha 2020); TMB-H is not additive.
9. Run LOHHLA / DASH for HLA-LOH; flag tumors with reduced neoantigen presentation.
10. For research-grade neoantigen quality, compute Luksza fitness + McGranahan clonality.

## Tips

- FoundationOne CDx **SCORED region** is 0.8 Mb (NOT the 1.1 Mb total panel). Using panel total inflates TMB ~37%.
- FoundationOne CDx INCLUDES synonymous variants in TMB; MSK-IMPACT and most academic pipelines exclude. Match the comparison reference.
- Cross-panel TMB is not directly comparable without Vega 2021 calibration. FoundationOne 10/Mb = TSO500 7.8 = Oncomine 8.4.
- bTMB is research-grade in tissue-naive contexts; BFAST Cohort C failed primary endpoint due to low ctDNA shed.
- TMB-H is NOT endorsed for breast, prostate, glioma per ESMO 2024 + McGrail 2021. Use tumor-type-specific cutoffs.
- MSI-H + TMB-H is statistical tautology (MSI-H is uniformly hypermutator). MSI-H is primary biomarker.
- HLA-LOH (~17% pan-cancer; >30% HNSCC / NSCLC / cervical) abolishes neoantigen presentation; run LOHHLA / DASH for high-TMB tumors.
- Hypermutator (>=100/Mb) signature stratification: MMR-D = SBS6/15/26/44 + ID1/2; POLE-exo = SBS10a/10b + SBS28.
- Ultra-hypermutator (>=500/Mb) is typically POLE+MMR concurrent; ICI excellent response.
- COSMIC hotspots (BRAF V600E, KRAS G12C) excluded by convention (not random; biases TMB upward).
- Population-stratified gnomAD AF (grpmax FAF95) for ancestry-diverse cohorts; EUR-only filter inflates AFR/EAS TMB.
- Tumor purity floor: FoundationOne >= 20%, MSK-IMPACT >= 30%; below this, VAF-based filtering unreliable.
- Panel size minimum: >= 0.8 Mb workable; >= 1.0 Mb preferred; < 0.5 Mb unreliable (Predicine ATLAS borderline).

## Related Skills

- clinical-databases/somatic-signatures - Mutational signatures (SBS3 HRD, MMR-D, POLE)
- clinical-databases/msi-detection - MSI-H related ICI biomarker
- clinical-databases/hla-typing - HLA typing for neoantigen + LOH
- variant-calling/variant-calling - Mutect2 / Strelka2 somatic upstream
- variant-calling/clinical-interpretation - ACMG/AMP cancer framework
