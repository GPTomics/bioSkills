# Heritability Partitioning - Usage Guide

## Overview

Estimate SNP heritability (`h2_SNP`) and decompose it across functional categories, tissues, and individual loci using LD-score regression (LDSC), LDAK SumHer, HESS local h2, HDL high-definition likelihood, BOLT-REML, and GCTA-GREML. Choose between summary-statistic methods (LDSC, LDAK, HDL, HESS) and individual-level methods (BOLT-REML, GCTA), reconcile LDSC vs LDAK model-dependent enrichment estimates, prioritize trait-relevant cell types via Finucane 2018 chromatin partitioning, compute trans-ancestry rg with Popcorn, and report calibrated h2 with intercept and ratio diagnostics.

## Prerequisites

```bash
# Python 3 LDSC fork (official bulik/ldsc is Python 2.7 / unmaintained since 2019).
# belowlab/ldsc v3.0.1 broke the --h2/--rg/--h2-cts CLI per its README; for a
# working CLI use abdenlab/ldsc-python3 (v2.0.0) or run belowlab via Docker
# (`docker pull jtb114/ldsc:latest`).
git clone https://github.com/abdenlab/ldsc-python3.git
cd ldsc-python3 && conda env create -f environment.yml && conda activate ldsc

# LDAK 6+
wget http://dougspeed.com/wp-content/uploads/ldak6.linux_.zip

# Pre-computed reference resources (EUR; ~3-5 GB total one-time)
# All from https://alkesgroup.broadinstitute.org/LDSCORE/
#   eur_w_ld_chr.tar.bz2                           (univariate h2 / rg)
#   1000G_Phase3_baselineLD_v2.2_ldscores.tgz      (partitioned h2)
#   1000G_Phase3_frq.tgz                           (allele frequencies)
#   1000G_Phase3_weights_hm3_no_MHC.tgz            (regression weights)
#   Multi_tissue_chromatin_1000Gv3_ldscores.tgz    (Finucane 2018 cell-type prioritization)
```

```r
# HDL R package (genetic correlation with non-overlapping cohorts)
remotes::install_github('zhenin/HDL/HDL')
# Download HDL UKB SVD reference (~5 GB) per HDL wiki
```

Inputs: GWAS summary statistics with columns SNP, A1, A2, BETA (or Z), SE, P, N (per-SNP or column-supplied). Allele frequency column EAF strongly recommended. For case-control GWAS, also supply sample case fraction (`--samp-prev`) and population lifetime prevalence (`--pop-prev`).

## Quick Start

Tell your AI agent what you want to do:
- "Compute SNP heritability from this GWAS summary statistic file with LDSC, EUR ancestry"
- "Partition h2 across functional annotations using the baseline-LD v2.2 model"
- "Prioritize trait-relevant tissues from ENCODE chromatin marks via Finucane 2018 cell-type S-LDSC"
- "Compute genetic correlation between trait1 and trait2 from sumstats; use cross-trait LDSC because samples overlap"
- "Run HDL.rg for genetic correlation since cohorts are non-overlapping; want lower variance than LDSC"
- "Compute local heritability per LDetect locus with HESS and identify high-h2 loci for fine-mapping"
- "Reconcile LDSC and LDAK SumHer enrichment estimates for the conserved-region annotation"
- "Convert this case-control LDSC h2 from observed to liability scale; population prevalence is 0.005"
- "Run Popcorn for trans-ancestry rg between this EUR GWAS and the matched EAS GWAS"

## Example Prompts

### Total h2 from EUR Sumstats
> "Munge `gwas_T2D.tsv` with LDSC's `munge_sumstats.py` against the HapMap3 SNP list, then run `ldsc.py --h2` with `eur_w_ld_chr/` for the univariate estimate. Report intercept, mean chi-square, ratio, and h2 with SE on the liability scale (samp-prev 0.08, pop-prev 0.10)."

### Partitioned h2 with Baseline-LD
> "Partition T2D h2 across the baseline-LD v2.2 annotations. Use `--overlap-annot --print-coefficients`. Apply Bonferroni at `0.05 / N_annotations` for per-annotation enrichment claims. Report the top 5 categories with enrichment > 5x and joint p < 2e-3."

### Cell-Type / Tissue Prioritization (Finucane 2018)
> "Prioritize trait-relevant tissues using the Multi_tissue_chromatin_1000Gv3 ldcts file with `ldsc.py --h2-cts`. Apply Bonferroni at `0.05 / N_tissue` (~ 2.5e-4 for 200 tissues). Cross-reference top tissues against the disease's known biology."

### LDSC vs LDAK Reconciliation
> "Functional enrichment claim depends on the per-SNP heritability model. Run BOTH LDSC baseline-LD AND LDAK SumHer with LDAK-Thin tagging. If they disagree by > 2x, cite Hou 2019 Nat Genet 51:1244 and report enrichment as model-dependent. Prefer LDAK for conserved regions per Speed 2019."

### Cross-Trait Genetic Correlation
> "Estimate rg between trait1 and trait2 from sumstats. If sample overlap > 5%, use cross-trait LDSC (`ldsc.py --rg`) because HDL is biased under overlap. If non-overlapping, use HDL for ~60% lower variance. Report rg, SE, p, and cross-trait intercept."

### Local Heritability
> "Compute per-locus h2 across all 22 autosomes using HESS with the LDetect EUR partition (Berisa & Pickrell 2016). Require each locus to have >= 1000 SNPs. Identify loci with h2 > 0.001 for fine-mapping prioritization downstream."

### Liability-Scale Case-Control h2
> "Case-control GWAS for schizophrenia. Supply `--samp-prev <case_fraction>` and `--pop-prev 0.01` so LDSC reports h2 on the liability scale. Without these flags, observed-scale h2 is incomparable across studies."

### Non-EUR Ancestry h2
> "EAS GWAS for type 2 diabetes. Use ancestry-matched LD scores from the EAS folder on alkesgroup.broadinstitute.org/LDSCORE, NOT the default EUR scores. Mismatched LD biases h2 downward and inflates the intercept."

### Trans-Ancestry rg
> "Cross-population genetic correlation between EUR and EAS GWAS for the same trait. Use Popcorn-2 with per-population LD scores. Require effective N > 5000 per population for stable estimation."

### Single-Cell ATAC Annotation
> "Build cell-type-specific .ldcts from per-cluster ATAC peaks (cross-reference atac-seq/single-cell-atac). Compute per-cluster LD scores via `ldsc.py --l2`. Then run `--h2-cts` to identify which cluster's open-chromatin landscape is most heritability-enriched for the trait."

## What the Agent Will Do

1. Munge sumstats to LDSC format with HapMap3 SNP filter (`munge_sumstats.py`)
2. Compute total h2 with `--h2` and report intercept, mean chi-square, ratio jointly
3. Apply liability-scale conversion for case-control with `--samp-prev` and `--pop-prev`
4. Run baseline-LD v2.2 partitioned h2 with `--overlap-annot --print-coefficients`
5. Run Finucane 2018 cell-type S-LDSC with `--h2-cts` and Bonferroni-control per-tissue p
6. For trait pairs, decide cross-trait LDSC (overlap-tolerant) vs HDL (lower variance, no overlap)
7. Run LDAK SumHer as model-reconciliation alongside LDSC when enrichment claim is primary
8. For local h2, run HESS per chromosome with LDetect partition, aggregate genome-wide
9. For individual-level data, route to BOLT-REML (N > 100k) or GCTA-GREML (5k < N < 50k)
10. Report intercept, mean chi-square, ratio, h2 with SE, partitioned enrichment, and sensitivity over model

## Method Selection Quick Reference

| Have | N | Goal | Method |
|------|---|------|--------|
| Sumstats (EUR) | >= 50k | Total h2 | LDSC `--h2` |
| Sumstats (EUR) | >= 50k | Partitioned h2 | S-LDSC + baseline-LD v2.2 |
| Sumstats (EUR) | >= 50k | Tissue prioritization | S-LDSC `--h2-cts` (Finucane 2018) |
| Sumstats (EUR) | >= 50k each | rg, non-overlap | HDL |
| Sumstats (EUR) | >= 50k each | rg, overlap unknown | Cross-trait LDSC |
| Sumstats (non-EUR) | >= 50k | Total h2 | LDSC with ancestry-matched LD |
| Sumstats (multi-pop) | >= 5k per pop | Trans-ancestry rg | Popcorn / Popcorn-2 |
| Sumstats (EUR) | varies | Local h2 | HESS with LDetect partition |
| Individual genotypes | 5k - 50k | h2 | GCTA-GREML |
| Individual genotypes | > 100k | h2 multi-component | BOLT-REML |
| Two models matter | -- | Functional enrichment | LDSC + LDAK SumHer (report both) |

## Cell-Type / Tissue Prioritization Workflow

Finucane 2018 cell-type S-LDSC asks: which tissue's chromatin annotation contributes most to trait h2 conditional on the baseline-LD model?

1. Download `Multi_tissue_chromatin_1000Gv3_ldscores.tgz` from alkesgroup (covers Roadmap, ENCODE, GTEx)
2. Manifest `.ldcts` file lists one row per tissue: `<name> <ldscore_prefix>,<control_ldscore_prefix>`
3. `ldsc.py --h2-cts trait.sumstats.gz --ref-ld-chr-cts <manifest>.ldcts --ref-ld-chr baseline. --w-ld-chr weights.`
4. Output `*.cell_type_results.txt` has per-tissue Coefficient, SE, p-value
5. Apply Bonferroni at `0.05 / N_tissues` (~2.5e-4 for 200 tissues)
6. Top trait-relevant tissues = those passing Bonferroni AND with positive coefficient

For novel cell types (e.g. per-cluster scATAC), build a custom .ldcts by computing per-cluster LD scores from chromatin BED files via `ldsc.py --l2 --bfile <ref> --annot <cluster>.annot.gz --thin-annot --out <cluster>`.

## LDSC vs LDAK Decision Matrix

| Claim | Use | Confirm with |
|-------|-----|--------------|
| Total h2 (no enrichment) | LDSC | LDAK only if intercept anomalous |
| Functional enrichment (primary claim) | LDSC + LDAK SumHer | Report both; cite Hou 2019 if discordant |
| Cell-type prioritization | LDSC `--h2-cts` (Finucane 2018) | Optional: LDAK conditional on annotation |
| Conserved-region enrichment | LDAK SumHer (Speed 2019 recommends) | LDSC as comparison |
| Total rg, non-overlap | HDL primary | Cross-trait LDSC confirmatory |
| Total rg, with overlap | Cross-trait LDSC | -- |

## Computational Footprint

| Method | Per-trait runtime | Hardware |
|--------|------------------|----------|
| LDSC h2 (univariate) | minutes | laptop |
| Stratified LDSC + baseline-LD | 10-20 min | laptop |
| LDSC cell-type prioritization (~200 tissues) | hours, single-threaded | server |
| HESS genome-wide local h2 | hours per chromosome | cluster |
| BOLT-REML at N = 500k | days | cluster |
| HDL.rg (genetic correlation) | seconds to minutes | laptop |
| LDAK SumHer | tens of minutes | laptop |

graphREML on biobank-scale (N > 200k) typically beats S-LDSC per-trait runtime when an LDGM panel is available. Build runtime escalates with annotation count; cell-type prioritization with ~200 tissues is the most expensive per-trait step.

## Intercept Diagnostic Quick Reference

| Intercept | Ratio | Interpretation |
|-----------|-------|----------------|
| ~1.0 | ~0 | No inflation; h2 trustworthy |
| 1.0 - 1.1 | < 0.2 | Mostly polygenic; h2 trustworthy |
| 1.1 - 1.3 | 0.2 - 0.5 | Mild inflation; investigate population structure / sample overlap |
| 1.3 - 1.5 | 0.5 - 0.8 | Substantial inflation; report ratio jointly; consider re-genotype-QC |
| > 1.5 | -- | Re-run after PC adjustment or genomic control; do not interpret h2 |

## Tips

- LDSC intercept > 1 is NOT automatic evidence of confounding; report the ratio statistic jointly with mean chi-square.
- LDAK SumHer and LDSC use different per-SNP heritability priors; report both when enrichment is the primary claim (Hou 2019 Nat Genet 51:1244).
- HDL is biased when sample overlap > 5%; use cross-trait LDSC instead under any suspected overlap.
- Apply ancestry-matched LD scores for non-EUR GWAS; EUR scores on EAS / AFR sumstats yield biased h2 and inflated intercept.
- Always supply `--samp-prev` and `--pop-prev` for case-control GWAS; observed-scale h2 is not comparable across studies.
- For Finucane 2018 cell-type prioritization, Bonferroni at `0.05 / N_tissues` (~2.5e-4 for 200 tissues); report top tissues only if they pass.
- HESS requires >= 1000 SNPs per locus and an ancestry-matched LD reference; use the LDetect partition (Berisa & Pickrell 2016).
- For individual-level data with N > 100k, BOLT-REML is more precise than LDSC; for N 5k - 50k use GCTA-GREML.
- LDSC's official `bulik/ldsc` repository is Python 2.7 and unmaintained since 2019. Prefer `abdenlab/ldsc-python3` (v2.0.0) for the Python 3 CLI; the `belowlab/ldsc` v3.0.1 README states the `--h2 / --rg / --h2-cts` CLI is currently broken (Docker `jtb114/ldsc:latest` is the recommended belowlab-flavored fallback).
- Reference baseline-LD model: `baselineLD_v2.2` (Gazal 2017). Older `baselineLD_v1.1` and pre-baseline-LD `baseline_v1.1` exist but are superseded.
- Pre-compute LD scores once and reuse; downloading the EUR reference is ~3 GB but covers all standard LDSC analyses.

## Related Skills

causal-genomics/mendelian-randomization - h2 / rg-aware instrument selection and sample-overlap decisions
causal-genomics/colocalization-analysis - Per-locus shared-causal evidence complementary to HESS local h2
causal-genomics/fine-mapping - Credible-set construction at high-h2 HESS loci
causal-genomics/pleiotropy-detection - Cross-trait pleiotropy via LCV / LHC-MR using LDSC outputs
causal-genomics/genomic-sem - Genomic SEM extends LDSC rg to multivariate structural models
causal-genomics/transcriptome-wide-association - TWAS uses partitioned-h2 weights for gene-level testing
atac-seq/differential-accessibility - Per-cell-type chromatin annotations as S-LDSC input
atac-seq/single-cell-atac - scATAC peaks per cluster as .ldcts annotations
chip-seq/peak-calling - ENCODE / Roadmap chromatin marks for cell-type prioritization
population-genetics/association-testing - GWAS source summary statistics for LDSC munging
population-genetics/linkage-disequilibrium - LD reference panels for HESS / coloc.susie
workflows/gwas-pipeline - Upstream GWAS pipeline feeding sumstats to LDSC
