# EWAS Design

## Overview

This skill covers the design-and-inference layer of an epigenome-wide association study (EWAS) on 450K/EPIC array or bisulfite methylation - the part that decides whether a hit is credible, which lives BEFORE the per-CpG test. It teaches the confounding hierarchy (cell composition as the dominant confounder, then batch/Sentrix chip/array position, age/sex, smoking, genetic ancestry/mQTL, and reverse causation/tissue relevance), how to randomize samples to chips at design time, how to remove unwanted variation with surrogate variable analysis (sva/SmartSVA), ComBat, and RUVm without over-correcting, how to handle genomic inflation with BACON instead of GWAS genomic control, which genome-wide significance threshold to use per array, how to power an EWAS with pwrEWAS, how to meta-analyze, how to replicate via the EWAS Catalog/Atlas, and how to build and interpret a methylation risk score (MRS).

The idea that reframes the whole analysis: an EWAS hit is a cell-composition difference until proven otherwise, and the signal (often under 2% methylation) is smaller than every confounder, so the p-value is the last and least interesting thing.

## Prerequisites

R / Bioconductor packages:

```r
BiocManager::install(c('meffil', 'sva', 'bacon', 'limma', 'missMethyl', 'minfi'))
install.packages('SmartSVA')
BiocManager::install('pwrEWAS')   # power estimation
```

Python (optional helper for layout balance checks):

```bash
pip install pandas numpy
```

Conceptual prerequisites:
- A normalized methylation matrix (betas and/or M-values) with samples in columns and CpGs in rows, plus a sample sheet carrying phenotype, age, sex, Sentrix chip/position, and plate.
- Cell-type proportions are estimated upstream (see cell-type-deconvolution) and supplied here as covariates.
- ComBat and per-CpG modeling run on M-values (logit of beta); effect sizes are reported back on the beta / delta-beta scale.
- The array version (450K vs EPIC v1 vs EPIC v2) sets the genome-wide threshold - confirm against the manifest.

## Quick Start

Tell your AI agent what you want to do:
- "Design a balanced chip layout for my case-control methylation study"
- "Pick the covariate set for a blood EWAS of my phenotype"
- "Run an EWAS and correct genomic inflation with BACON"
- "Choose the genome-wide significance threshold for an EPIC array"
- "Power an EWAS to detect a 2% methylation difference"
- "Check whether my top CpG is a known smoking/age hit before claiming novelty"
- "Use a methylation smoking score as a covariate instead of self-report"

## Example Prompts

### Confounding and covariates
> "I have whole-blood EPIC data for 200 cases and 200 controls. Build the EWAS covariate model, explain why cell composition is the dominant confounder, and tell me which covariates to include and why."

### Chip randomization at design time
> "I am about to run 384 samples on EPIC BeadChips. Randomize sample-to-chip and sample-to-position so case/control, age, and sex are balanced across chips, and verify no chip is single-group."

### Genomic inflation and BACON
> "My EWAS QQ plot looks inflated with lambda 1.3. Explain why GWAS genomic control is wrong here, and correct the bias and inflation with BACON, reporting the empirical-null bias and inflation."

### Over-correction check
> "I am adding surrogate variables with SmartSVA. Show me how to monitor the smoking positive control cg05575921 through the correction so I do not over-correct."

### Thresholds and replication
> "Set the genome-wide significance threshold for my 450K study, report both FWER and FDR hits with delta-beta effect sizes, and tell me how to triage hits against the EWAS Catalog."

### Power and meta-analysis
> "Power an EWAS to detect a 2% delta-beta at genome-wide significance in blood using pwrEWAS, and explain why we will likely need meta-analysis."

### Methylation risk scores
> "Explain the difference between a polygenic risk score and a methylation risk score, and whether I should treat my DNAm score for disease X as a cause or a consequence."

## What the Agent Will Do

1. Establish the confounding hierarchy for the tissue and phenotype, and assemble the default covariate set (age, sex, cell proportions, chip/position or SVs, and where relevant genetic PCs and a smoking score).
2. If the study is pre-bench, randomize or block sample-to-chip/position/plate assignment so technical batch is orthogonal to phenotype, and verify balance.
3. Remove unwanted variation with sva/SmartSVA, ComBat, or RUVm as appropriate, monitoring a positive control (cg05575921 AHRR) through the correction to detect over-correction.
4. Fit the per-CpG model on M-values (handing the mechanics to differential-cpg-testing), then correct genomic inflation with BACON and report the empirical-null bias and inflation.
5. Apply the array-specific genome-wide threshold for the headline claim and BH-FDR for discovery, reporting delta-beta effect sizes.
6. Estimate power with pwrEWAS using site-specific variance, and frame meta-analysis as the route to adequate power.
7. Triage hits against the EWAS Catalog/Atlas for replication and to flag generic age/smoking/cell-composition CpGs.
8. Where an MRS is involved, interpret it as a predictive biomarker (state consequence), contrast it with a PRS (germline cause), and defer weight-learning to the machine-learning category.

## Tips

- Treat any unadjusted hit as a cell-composition difference until proven otherwise; the cell-fraction covariates are the most important line in the model.
- Randomization at the bench is the single uncorrectable decision - it matters more than every analysis choice that follows.
- A clean QQ plot with a dead positive control is a broken pipeline, not a good one; never add surrogate variables until lambda hits 1.0 at the cost of AHRR.
- Lambda alone diagnoses nothing in an EWAS; pair it with the QQ-plot shape and whether the positive control survives, and correct with BACON rather than genomic control.
- Report effect size (delta-beta) plus the FWER hit plus the FDR hit plus replication status; a genome-wide-significant CpG with a 0.3% effect and no replication is fragile.
- An EWAS Catalog entry is included at P < 1e-4, so it is a lookup, not a genome-wide claim.
- A methylation risk score is reverse-causal-by-default; reserve "risk" and "cause" for prospectively- or MR-supported claims.

## Related Skills

- differential-cpg-testing - The per-site test this design layer feeds
- cell-type-deconvolution - Cell-fraction covariates (the dominant confounder)
- array-preprocessing - Normalization choice (funnorm) vs model-level batch correction
- array-qc-filtering - Probe filtering and chip/position batch diagnosis
- causal-genomics/mendelian-randomization - mQTL-based causal orientation (reverse causation)
- experimental-design/batch-design - General randomization and batch-design principles
- clinical-biostatistics/multiplicity-graphical - FWER for confirmatory trials (contrast with discovery FDR)
- workflows/methylation-pipeline - End-to-end pipeline
