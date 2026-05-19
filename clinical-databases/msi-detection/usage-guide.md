# MSI Detection - Usage Guide

## Overview

Detect microsatellite instability from WES / WGS / targeted-panel data using MSIsensor (paired tumor-normal), MSIsensor-pro (tumor-only), MSIsensor-ct (ctDNA / liquid biopsy), MANTIS, or MSIngs. Covers FDA pembrolizumab MSI-H / dMMR pan-tumor approval (Le 2015 + KEYNOTE-016/164/158/177), Lynch syndrome universal screening (IHC + MSI), the MSI-H + TMB-H tautology (Sha 2020), and distinguishing POLE-exo hypermutator (typically MSI-stable) from MMR-D (MSI-H).

## Prerequisites

```bash
# MSIsensor-pro (tumor-only; the modern standard)
conda install -c bioconda msisensor-pro

# MSIsensor (paired normal; WES gold standard)
conda install -c bioconda msisensor

# MSIsensor-ct (liquid biopsy)
# git clone https://github.com/niu-lab/msisensor-ct

# MANTIS (alternative paired WES)
# git clone https://github.com/OSU-SRLab/MANTIS

# mSINGS (background-panel tumor-only)
# https://bitbucket.org/uwlabmed/msings/
```

## Quick Start

Tell the agent what to do:
- "Run MSIsensor-pro on this tumor-only BAM; report MSI-H if >= 20% unstable loci"
- "Paired tumor-normal WES: run MSIsensor; compare with MANTIS as orthogonal evidence"
- "For this ctDNA panel, run MSIsensor-ct; flag if tumor fraction < 3%"
- "Apply Lynch syndrome workflow: IHC results + MSI status + MLH1 methylation; classify sporadic vs Lynch suspect"
- "Reconcile MSI-H tumor with retained IHC: check POLE-exo via signatures + germline MMR"

## Example Prompts

### Standard MSI Calling

> "Run MSIsensor-pro on this tumor-only BAM with the pre-computed baseline; classify MSI-H if >= 20% unstable loci."

> "For paired tumor-normal WES, run MSIsensor; threshold MSI-H = >= 20% unstable; compare with MANTIS step-wise difference > 0.4."

### ctDNA / Liquid Biopsy

> "Estimate tumor fraction via ichorCNA; if >= 3%, run MSIsensor-ct on this cfDNA panel; otherwise fall back to tissue."

> "For longitudinal ctDNA MSI monitoring in MSI-H patient on ICI, track unstable locus proportions over time."

### Lynch Syndrome Workflow

> "Apply universal Lynch screening protocol: IHC + MSI + MLH1 methylation. For MLH1 loss + methylation -> sporadic. For MSH2/6/PMS2 loss or MLH1 loss without methylation -> germline testing."

> "MSI-H CRC patient < 70 yr: confirm with IHC (MLH1/MSH2/MSH6/PMS2); if MLH1 loss -> methylation test; if unmethylated -> germline panel for Lynch syndrome."

### ICI Eligibility

> "For this dMMR CRC patient on first-line pembrolizumab (KEYNOTE-177), confirm MSI-H + IHC + report TMB but note not additive (Sha 2020)."

> "For TMB-H breast cancer with MSS: do NOT recommend ICI based on TMB alone (ESMO 2024 / McGrail 2021 exclusion)."

### POLE vs MMR-D Disambiguation

> "Patient has 250 mut/Mb hypermutator. Run Sigprofiler: if SBS10a/10b dominant -> POLE-exo (typically MSI-stable); if SBS6/15/26/44 dominant -> MMR-D; if SBS14 + SBS20 present -> POLE+MMR concurrent ultra-hypermutator."

> "Confirm POLE-exo hypothesis with germline POLE / POLD1 sequencing; ICI excellent response expected."

### Reconciliation

> "PCR Bethesda MSI-H vs NGS MSS discordance: trust NGS with >=50 informative loci over 5-locus PCR panel."

> "NGS MSI-H + retained IHC: check germline POLE; check signatures; consider rare MSH6-only subtype."

## What the Agent Will Do

1. Choose tool by data: MSIsensor (paired WES), MSIsensor-pro (tumor-only), MSIsensor-ct (ctDNA), MANTIS (paired alternative).
2. For tumor-only, generate or use pre-computed baseline from cohort.
3. Apply panel-specific threshold: >= 20% unstable loci typical; calibrated per panel.
4. For Lynch screening, integrate IHC + MSI + MLH1 methylation + germline testing.
5. For ICI decisions, treat MSI-H / dMMR as PRIMARY biomarker; TMB-H is correlate not additive (Sha 2020).
6. For high mutation burden with unclear MSI, run Sigprofiler to distinguish POLE-exo (SBS10a/10b; typically MSI-stable) vs MMR-D (SBS6/15/26/44) vs POLE+MMR (SBS14, SBS20).
7. For ctDNA, estimate tumor fraction first; require >= 3% for reliable cfDNA MSI.
8. Cross-validate with IHC (MLH1, MSH2, MSH6, PMS2) when available.

## Tips

- FDA pembrolizumab MSI-H / dMMR pan-tumor approval (2017) requires MSI-H **or** dMMR; either suffices.
- MSI-H + TMB-H is statistical tautology (Sha 2020 *Cell Rep Med*); MSI-H is primary biomarker.
- TMB-H + MSS: ESMO 2024 endorses ICI but NOT for breast / prostate / glioma (McGrail 2021 exclusion).
- POLE-exo (SBS10a/10b) typically MSI-stable; pure POLE-exo 100-300 mut/Mb; POLE+MMR concurrent >500 mut/Mb.
- Lynch syndrome ~50% of MSI-H CRC; rest is sporadic (MLH1 hypermethylation).
- Universal Lynch screening per NCCN / ACG / EGAPP: all CRC <= 70 yr.
- MLH1 loss + methylation -> sporadic (not Lynch).
- MSH2 / MSH6 / PMS2 loss OR MLH1 loss without methylation -> Lynch suspect; germline testing.
- ctDNA MSI requires tumor fraction >= 3%; below that, false-negative likely.
- Bethesda 5-locus PCR panel is less sensitive than NGS >=50 loci; trust NGS for borderline.
- Tumor purity >= 20% for reliable MSI calling.
- MSI-L is clinically equivalent to MSS per FDA; apply MSI-H threshold strictly.
- MSH6-only loss is a more variable phenotype; may have MSI-H without all 4 MMR proteins lost.

## Related Skills

- clinical-databases/tumor-mutational-burden - TMB as correlate biomarker
- clinical-databases/somatic-signatures - SBS6/15/26/44 (MMR-D), SBS10a/10b (POLE-exo)
- clinical-databases/clinvar-lookup - Lynch genes MLH1 / MSH2 / MSH6 / PMS2 variant pathogenicity
- clinical-databases/variant-prioritization - Germline MMR variant prioritization
- variant-calling/clinical-interpretation - Clinical reporting
