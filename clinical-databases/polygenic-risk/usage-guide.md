# Polygenic Risk Scores - Usage Guide

## Overview

Construct and validate polygenic risk scores using the 2026 method landscape (LDpred2-auto, SBayesRC, MegaPRS, PRS-CS, PROSPER, MUSSEL, BridgePRS, JointPRS, PRSmix) and the PGS Catalog Calculator. Covers cross-ancestry portability with Ding 2023 continuous-ancestry framing, Hingorani 2023 BMJ Medicine performance critique, PRS-RS reporting standards, EraSOR sample-overlap detection, and the 11-item reviewer audit. Regulatory state: no general FDA PRS guidance exists as of 2026; August 2025 Federal Register Cancer Predisposition Risk Assessment System device classification is the operative text.

## Prerequisites

```bash
# Production multi-PGS pipeline
nextflow pull pgscatalog/pgsc_calc

# LDpred2 (R)
# install.packages('bigsnpr')

# PRS-CS / PRS-CSx
git clone https://github.com/getian107/PRScs.git
git clone https://github.com/getian107/PRScsx.git
pip install scipy h5py

# SBayesRC (GCTB)
# Download from https://cnsgenomics.com/software/gctb/
# Need UKB eigenvalue-decomposed LD + baselineLD_v2.2 annotations

# Multi-ancestry
# PROSPER: https://github.com/Jingning-Zhang/PROSPER
# MUSSEL: https://github.com/Jin93/MUSSEL
# BridgePRS: https://github.com/clivehoggart/BridgePRS
# JointPRS: https://github.com/LeqiXu/JointPRS

# Sample overlap
# Choi 2023 EraSOR: https://github.com/choishingwan/EraSOR
# LDSC: pip install ldsc
```

## Quick Start

Tell the agent what to do:
- "Compute LDpred2-auto PRS for CAD using Aragam 2022 GWAS sumstats. Run EraSOR for sample-overlap detection first."
- "Use pgsc_calc to compute PGS000004 + PGS001775 + PGS000019 on this cohort with auto-liftover and ancestry projection"
- "Multi-ancestry breast-cancer PRS using PRS-CSx with EUR + AFR + EAS sumstats. Tune phi for sparse vs polygenic architecture."
- "For my T2D PRS, apply Ding 2023 continuous-ancestry framing: report R^2 vs PC distance, NOT vs discrete ancestry boxes."
- "Audit this PRS manuscript against the 11-item reviewer checklist (PRS-RS items 13, 16, 19, 21 priority)."

## Example Prompts

### Single-Ancestry Construction

> "Compute LDpred2-auto PRS for breast cancer in this UKB-NFE cohort. Use UKB-LD reference (not 1KG-EUR); pin allow_jump_sign=FALSE, shrink_corr=0.95; report `s` LD-mismatch parameter."

> "Use SBayesRC with baselineLD v2.2 annotations + UKB eigen-LD for this lipid panel cohort. Run --impute-summary first."

> "Run MegaPRS on T2D sumstats with the BLD-LDAK heritability model; select the best auto-model on training set."

### Multi-Ancestry

> "Compute PRS-CSx joint posterior for CAD using EUR (n=200k) + AFR (n=30k) + EAS (n=50k) sumstats. Phi=1e-2 for highly polygenic; report cross-ancestry R^2."

> "For this admixed cohort, run MUSSEL with multivariate spike-and-slab + super-learner ensemble. Compare to BridgePRS for low-h^2 settings."

> "JointPRS-meta with small tuning set; compare to PRS-CSx auto-tuning."

### PGS Catalog Production

> "Run pgsc_calc on this cohort VCF with PGS000004 (CAD-Khera-2018), PGS001775 (CAD-Aragam-2022), PGS003725 (BC-Mavaddat-2019). Auto-liftover GRCh37 -> 38; project to HGDP+1KG references for ancestry."

### Calibration and Reporting

> "For my CAD PRS, integrate over BOADICEA-like age-conditional incidence to produce absolute lifetime risk; flag the Hingorani 2023 detection-rate-vs-FPR caveat."

> "Apply Ding 2023 continuous-ancestry framework: for each test sample, compute genetic distance to discovery cohort and report R^2 binned by quintile of distance."

### Audit / Reproducibility

> "Audit this PRS paper against PRS-RS 22-item checklist. Priority items: cohort independence (#13), confounder adjustment (#16), absolute-risk reporting (#19), ancestry composition of validation (#21)."

> "Run EraSOR / bivariate LDSC intercept on this PRS application; report whether |intercept| > 0.05 flags substantial overlap."

## What the Agent Will Do

1. Choose method by data: SBayesRC for EUR + functional annotations; LDpred2-auto for EUR generic; PRS-CSx / PROSPER / MUSSEL / BridgePRS / JointPRS for multi-ancestry; PRSice-2 for legacy C+T comparison.
2. Match variants strand-aware (LDpred2 snp_match with frequency-flip; PRSice-2 drops A/T C/G).
3. Use UKB-LD reference (not 1KG-EUR) for EUR; ancestry-matched for non-EUR.
4. Run EraSOR / bivariate LDSC intercept to detect sample overlap; require |intercept| < 0.05 with target n >= 1000.
5. Exclude HLA region (chr6 28-34 Mb) from main PRS; model classical HLA alleles separately for autoimmune traits.
6. Recalibrate ancestry-conditional Z using TEST cohort PCs (not discovery PCs).
7. Transform PRS percentile to absolute risk via external age-conditional incidence curve.
8. Apply Hingorani 2023 caveats: HR per SD 1.3-1.7 is comparable to family history alone; PRS does not meet population-screening utility threshold.
9. Cite methods correctly: PRSmix and MUSSEL are *Cell Genomics*; PROSPER is *Nat Commun*; Mavaddat 2023 is CEBP; Mullins 2021 is bipolar (not MDD).

## Tips

- 2026 ranking: SBayesRC ~ MegaPRS ~ LDpred2-auto > PRS-CS-auto > lassosum2 > C+T for EUR; PROSPER ~ MUSSEL ~ JointPRS > PRS-CSx > BridgePRS for multi-ancestry.
- Bayesian methods need ~1 hour CPU for typical traits; SBayesRC + UKB-LD requires ~50 GB RAM.
- LDpred2-auto requires pinning allow_jump_sign=FALSE and shrink_corr=0.95 (per Privé 2022 misspecification paper).
- Compute `s` LD-mismatch parameter routinely; `s > 0.05` means LD reference does not match GWAS.
- `pgsc_calc` is the nf-core Nextflow workflow; handles liftover + ancestry + normalization in one pass.
- HLA region (chr6 28-34 Mb) MUST be excluded or modeled separately; PRS-based on SNP dosages does not capture HLA-allele effects.
- Strand-ambiguous SNPs (A/T, C/G): LDpred2 snp_match flips by MAF (0.4-0.6 exclusion); PRSice-2 drops by default.
- Sex chromosomes routinely dropped; if included use sex-stratified dosage encoding.
- Cryptic relatedness inflates PRS performance; remove individuals with KING > 0.0884.
- LDpred2 LD: prefer UKB LD (n=40k+) over 1KG-EUR (n=489).
- Sample overlap test (EraSOR / bivariate LDSC intercept) is mandatory; |intercept| > 0.05 with target n>=1000 = problem.
- HR per SD typically 1.3-1.7 for most diseases; this is similar to family history; PRS alone does not justify population screening per Hingorani 2023.
- Citation traps: PRSmix and MUSSEL are *Cell Genomics*; PROSPER is *Nat Commun*; Hingorani 2023 is *BMJ Medicine* not main *BMJ*; Mavaddat 2023 update is *CEBP*; Mullins 2021 is bipolar disorder.
- Top 2.5% PRS detection rate ~7% of CAD cases captured; 5% false-positive rate (Hingorani 2023 PGS Catalog analysis).
- No general FDA PRS draft guidance exists as of May 2026; Aug 2025 Federal Register Cancer Predisposition Risk Assessment Class II classification is the relevant regulatory text.
- LDT Final Rule (May 2024) was vacated March 2025; LDT pathway is back to pre-2024 enforcement discretion.

## Related Skills

- clinical-databases/gnomad-frequencies - Population AF for QC + INFO filtering
- clinical-databases/variant-prioritization - Rare-variant filtering context
- clinical-databases/clinvar-lookup - Variant pathogenicity for monogenic + PRS contrast
- causal-genomics/mendelian-randomization - PRS as instrumental variable
- population-genetics/population-structure - Ancestry inference for PC computation
- machine-learning/biomarker-discovery - PRS in biomarker pipelines
