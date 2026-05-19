# ACMG Classification - Usage Guide

## Overview

Apply the ACMG/AMP 2015 framework with ClinGen SVI specifications using the Tavtigian 2018/2020 Bayesian point system (the engine inside every modern automated classifier). Cover the Abou Tayoun 2018 PVS1 decision tree, Pejaver 2022 PP3/BP4 calibrated thresholds for REVEL/BayesDel/AlphaMissense, Brnich 2020 PS3/BS3 OddsPath framework, Walker 2023 ClinGen Splicing Subgroup for SpliceAI/SpliceVault, VCEP-specific CSpec overrides, and AMP/ASCO/CAP 2017 cancer somatic Tier I-IV.

## Prerequisites

```bash
pip install requests pandas

# Automated ACMG classifiers
# GeneBe (open-source): https://genebe.net (API key optional)
# InterVar: git clone https://github.com/WGLab/InterVar.git
# VarSome (commercial, hosted)
# Franklin (Genoox, commercial, hosted)
# ClinGen VCI (gold standard for expert curation): https://curation.clinicalgenome.org/

# AutoPVS1 for VCEP-specific PVS1 decision trees
# https://autopvs1.bgi.com/
```

## Quick Start

Tell the agent what to do:
- "Apply ACMG/AMP framework to this variant; sum Tavtigian points; classify P/LP/VUS/LB/B"
- "Use Pejaver 2022 PP3/BP4 calibrated thresholds for REVEL: PP3_Strong >= 0.932, BP4_Strong <= 0.016"
- "Apply Abou Tayoun PVS1 decision tree for this nonsense variant; check NMD prediction + critical region"
- "Check VCEP CSpec for this gene before classification (e.g., Hearing Loss VCEP for OTOF)"
- "Apply AMP/ASCO/CAP 2017 cancer Tier I-IV to this somatic variant; cross-check OncoKB + CIViC"

## Example Prompts

### Standard Classification

> "Classify this missense variant: REVEL = 0.85, SpliceAI DS_max = 0.05, gnomAD grpmax FAF95 = 1e-6, in BRCA1 (Definitive validity). Apply Pejaver 2022 calibration."

> "Apply ACMG/AMP to chr17:43044294:G>A. Use one predictor (REVEL); check ClinGen ENIGMA-BRCA1 VCEP CSpec."

### PVS1 Decision Tree

> "Apply Abou Tayoun PVS1 decision tree to this BRCA1 c.5946delT frameshift. Check NMD prediction, critical region (BRCT domain), >10% coding removed."

> "PVS1 for nonsense in SCN5A: distinguish Brugada (LoF mechanism, PVS1 applies) vs LQT3 (GoF, PVS1 does NOT apply)."

### Splicing

> "Apply Walker 2023 ClinGen Splicing Subgroup framework: SpliceAI DS_max + SpliceVault predicted aberrant transcript + NMD prediction."

> "Variant has SpliceAI DS_max = 0.65; check SpliceVault 300K-RNA Top-4 events to determine which exon skips. Apply PVS1 if aberrant transcript triggers NMD."

### VCEP CSpec

> "Look up active VCEP CSpec for OTOF (Hearing Loss VCEP). Apply VCEP-specific rules: PM2 default = supporting, BA1 = 0.5%, PS3 thresholds for OTOF."

> "Check CSpec registry at cspec.genome.network for InSiGHT MMR VCEP. Apply gene-specific PM2 / PVS1 rules."

### Cancer Somatic

> "Classify this BRCA1 somatic variant in TNBC per AMP/ASCO/CAP 2017. Cross-check OncoKB Level + CIViC + CGI for therapy tier."

### Reconciliation

> "GeneBe classifies LP; VarSome says P. Identify the criteria difference and apply ClinGen VCI gold-standard for expert curation."

> "ClinVar 2018 says P; my evidence now suggests LP. Apply Tavtigian point system; report current classification with star rating."

## What the Agent Will Do

1. Pull aggregated annotations via myvariant.info (ClinVar / gnomAD / dbNSFP / VEP).
2. Compute gene-specific BS1 max-credible-AF (Whiffin 2017); apply against grpmax FAF95.
3. Apply Pejaver 2022 PP3/BP4 calibrated thresholds (REVEL/BayesDel/VEST4); ONE predictor only.
4. AlphaMissense as supporting evidence only (ClinGen not endorsed PP3 calibration).
5. Apply Abou Tayoun 2018 PVS1 decision tree for predicted LoF variants; check NMD + critical region + isoform.
6. Apply Walker 2023 Splicing Subgroup framework for splice variants.
7. Check VCEP-specific CSpec at `cspec.genome.network`; apply gene-specific overrides.
8. Sum Tavtigian points; classify P / LP / VUS / LB / B.
9. For cancer somatic, apply AMP/ASCO/CAP 2017 Tier I-IV; cross-check OncoKB / CIViC.
10. Report criteria + evidence trail + classification + ClinGen VCEP affiliation.

## Tips

- The Tavtigian point system is the engine inside every modern automated classifier (InterVar, GeneBe, VarSome, Franklin). The 2015 combining rules are subsumed.
- Final categories: P >= 10 points; LP 6-9; VUS 0-5; LB -1 to -6; B <= -7.
- REVEL PP3_Strong >= 0.932; BP4_Strong <= 0.016; the two thresholds to memorize.
- AlphaMissense developer threshold 0.564 is NOT the Pejaver PP3 calibration; treat as supporting only.
- Do NOT stack REVEL + BayesDel + VEST4; they share training data per Pejaver 2022.
- PM2 was downgraded to PM2_Supporting in SVI 2020 (1 point, not 2). Old papers with PM2_Moderate need re-classification.
- PS3 default Strong is over-strengthened; apply Brnich 2020 OddsPath framework.
- VCEP CSpec overrides generic ACMG. Check `cspec.genome.network/cspec/ui/svi/all` for active VCEP before classification.
- Abou Tayoun PVS1 decision tree refines PVS1 from binary to graded (VeryStrong / Strong / Moderate / Supporting). >15 VCEP-specific PVS1 trees as of 2024.
- PVS1 should NOT apply in GoF mechanism genes (e.g., SCN5A LQT3, KCNH2 LQT2).
- SpliceAI DS_max thresholds (Walker 2023): >=0.5 + corroborating evidence = PP3_Strong; >=0.20 minimum for ANY PP3; <0.1 = BP4_Moderate.
- SpliceVault (Dawes 2023) predicts WHICH aberrant transcript; critical for PVS1 applied to splice variants because NMD depends on the actual mis-spliced product.
- Cancer somatic = AMP/ASCO/CAP 2017 Tier I/II/III/IV; OncoKB Levels 1-4 map loosely; cross-check CIViC.
- ClinGen VCI is the gold standard for expert curation; automated tools approximate.
- This skill is for CLASSIFICATION. For rare-disease variant FILTERING + DNV + compound-het pipeline, use `clinical-databases/variant-prioritization`.

## Related Skills

- clinical-databases/variant-prioritization - Filtering pipeline (this skill is downstream classification)
- clinical-databases/clinvar-lookup - ClinVar evidence aggregation
- clinical-databases/gnomad-frequencies - Whiffin FAF95 for BS1/BA1
- clinical-databases/myvariant-queries - Aggregated annotation pull
- variant-calling/clinical-interpretation - Clinical reporting workflow
