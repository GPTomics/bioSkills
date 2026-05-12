# Pleiotropy Detection - Usage Guide

## Overview

Detect and adjust for horizontal pleiotropy in two-sample Mendelian randomization. Distinguishes uncorrelated horizontal pleiotropy (UHP) from correlated horizontal pleiotropy (CHP) and selects a method battery whose assumptions span both regimes.

UHP vs CHP is the central distinction:

- **UHP (uncorrelated)** - Pleiotropic alpha is independent of the instrument-exposure effect gamma. The InSIDE assumption holds. MR-Egger, weighted median, weighted mode, MR-PRESSO, MR-RAPS, MR-Mix, and contamination mixture address this regime.
- **CHP (correlated)** - Pleiotropic alpha is correlated with gamma via a shared upstream factor (heritable confounder, network mediator). InSIDE is violated. IVW, Egger, PRESSO, and GSMR are all blind to CHP and return systematically biased "corrected" estimates. CAUSE, LHC-MR, and LCV are designed for CHP.

A complete sensitivity battery must include at least one CHP-aware method when LDSC rg(exposure, outcome) `>= 0.3` or biology suggests a shared upstream factor.

## Operational Decision Flow

1. Compute LDSC genetic correlation (see causal-genomics/genetic-correlation). If `|rg| > 0.3`, escalate to Step 3.
2. Run the standard battery: IVW (random-effects if Cochran Q p < 0.05), MR-Egger with NOME I^2_GX check (SIMEX correct if 0.6 <= I^2_GX < 0.9), weighted median, weighted mode, MR-PRESSO at NbDistribution `>= 10000` for stringent reporting.
3. If CHP plausible (rg > 0.3 OR PRESSO global p < 0.05 with > 50% outliers OR Egger/median/mode disagree by > 2 SE), run CAUSE (if `>= 100` sig SNPs) or LHC-MR; report ELPD delta + z + q + gamma.
4. Triangulate with pre-MR Steiger filter, bidirectional MR, LCV gcp, and the LDSC rg report.

## LCV gcp Interpretation

| gcp value | Interpretation |
|-----------|----------------|
| 0 | Pure genetic correlation; no partial causation |
| 0.5 | Partial causation; mixture |
| 0.6 | Partial causation; modestly causal direction |
| 1 | Fully causal in tested direction |
| Significant p_gcp != 0 | Directional evidence of (partial) causation |

LCV uses ALL genome-wide SNPs after LDSC-merging; it does NOT use the MR instrument set. It complements MR rather than replacing it.

## Prerequisites

```r
install.packages(c('remotes', 'TwoSampleMR', 'MendelianRandomization', 'simex'))
remotes::install_github('rondolab/MR-PRESSO')          # CRAN-never; GitHub-only
remotes::install_github('jean997/cause')               # CHP-aware (Morrison 2020)
remotes::install_github('cnfoley/mrclust')             # Mechanism heterogeneity
remotes::install_github('gqi/MRMix')                   # Mixture-of-distributions
remotes::install_github('qingyuanzhao/mr.raps')        # CRAN-archived 2025-03; install from GitHub
remotes::install_github('LizaDarrous/lhcMR')           # Heritable-confounder MR
# LCV: git clone https://github.com/lukejoconnor/LCV (R scripts, no package)
```

Inputs are typically a harmonized TwoSampleMR data.frame (beta.exposure, beta.outcome, se.exposure, se.outcome, SNP, effect_allele, eaf, etc.). LHC-MR and LCV take genome-wide GWAS sumstats with LDSC-style merged SNPs; CAUSE takes pruned signature SNPs plus full sumstats for nuisance estimation.

## Quick Start

Tell your AI agent what you want to do:
- "Run the standard MR sensitivity battery on my TwoSampleMR-harmonized data"
- "I suspect a shared heritable confounder; run CAUSE alongside IVW, Egger, and MR-PRESSO"
- "My instruments are mostly weak (mean F < 20); use MR-RAPS instead of IVW"
- "I have a polygenic exposure with few significant SNPs; run LHC-MR or LCV"
- "Cluster my instruments by causal estimate to identify heterogeneous mechanisms"
- "Run bidirectional MR with Steiger pre-filtering on both directions"
- "Apply SIMEX to my MR-Egger because I^2_GX is 0.75"
- "Build a STROBE-MR sensitivity reporting table from my MR results"

## Example Prompts

### Standard sensitivity battery (UHP-focused)

> "Run IVW, MR-Egger, weighted median, weighted mode, Cochran Q, Egger intercept, MR-PRESSO with 5000 distributions, Steiger directionality, and leave-one-out on my harmonized data and produce a comparison table."

> "Compute I^2_GX for my Egger analysis; if below 0.9 apply SIMEX correction."

### Suspected CHP (correlated pleiotropy)

> "LDSC genetic correlation between my exposure and outcome is 0.45. Run CAUSE in addition to MR-PRESSO and compare delta_ELPD against the sharing model."

> "Use LHC-MR to jointly estimate forward causal effect, reverse causal effect, and heritable-confounder contribution from genome-wide sumstats."

> "Compute LCV gcp from LDSC-merged sumstats to distinguish causation from pure genetic correlation."

### Weak-IV regime

> "Mean F-statistic across my instruments is 14. Run MR-RAPS with Huber robust loss and overdispersion modeling; report point estimate and CI alongside IVW."

> "Apply MR-Mix and contamination mixture to my weak-IV dataset; reconcile against MR-RAPS."

### Heterogeneous mechanisms

> "I suspect my LDL instruments operate through multiple lipoprotein pathways. Run MR-Clust and report per-cluster IVW estimates with biological annotation suggestions."

### Polygenic exposure

> "My exposure is polygenic with only 40 genome-wide-significant SNPs. CAUSE is underpowered. Run LHC-MR or report LCV gcp instead."

### Drug-target cis-MR

> "Run colocalization at the cis-locus instead of Egger because I only have 3 SNPs in tight LD; complement with Steiger directionality."

### Bidirectional MR

> "Instrument both directions and apply Steiger pre-filter before primary IVW; report forward and reverse IVW + Egger + median + mode side by side."

### SIMEX rescue for Egger

> "My I^2_GX is 0.78. Apply SIMEX correction to MR-Egger using the simex package; report SIMEX-corrected slope alongside naive slope."

### STROBE-MR reporting

> "Generate a STROBE-MR compliant table summarizing all methods, F-statistics, I^2_GX, Egger intercept, PRESSO global / distortion / corrected, Steiger, and CAUSE delta_ELPD with citations."

## What the Agent Will Do

1. Verify instrument selection (F-statistic distribution, LD pruning, MAF, harmonization, palindrome handling)
2. Run IVW + Egger + weighted median + weighted mode side-by-side
3. Test Egger intercept for directional UHP; compute I^2_GX for NOME validity; SIMEX-correct if needed
4. Run MR-PRESSO with `>= 5000` distributions for global + outlier + distortion tests; produce corrected estimate
5. Apply Steiger filter and report directionality test
6. Run leave-one-out for stability
7. If CHP is plausible (rg `>= 0.3` or biology suggests shared factor): add CAUSE (if sig SNPs `>= 100`) OR LHC-MR
8. If mechanisms heterogeneous: run MR-Clust
9. If weak instruments: switch primary to MR-RAPS with robust loss
10. Reconcile UHP-method agreement vs CHP-method estimate; flag discordances
11. Produce STROBE-MR table with all methods, thresholds, and limitations

## Tips

- **UHP vs CHP** - MR-PRESSO global non-significance does NOT rule out pleiotropy; it rules out UHP outliers. CHP is invisible to PRESSO.
- **NbDistribution** - PRESSO uses 5000 minimum for stable p-values; published default 1000 gives noisy global-test p.
- **Egger NOME** - Below I^2_GX 0.9, the Egger slope is biased toward null. SIMEX rescues the 0.6-0.9 range; below 0.6, use a different method (RAPS).
- **Few-SNP Egger** - With fewer than 10 SNPs, Egger intercept CI is too wide to falsify pleiotropy; do not interpret non-significant intercept as "no pleiotropy" in this regime.
- **CAUSE sample overlap** - When exposure and outcome GWAS share controls, CAUSE's rho parameter corrects for this; explicitly estimate it via `est_cause_params`.
- **CAUSE underpowered regime** - Below 100 significant SNPs, delta_ELPD CI spans zero; consider LHC-MR or LCV instead.
- **Steiger interpretation** - When exposure measurement error is high, Steiger flags can be artifactual; treat as one signal among many, especially when biology strongly supports the tested direction.
- **MR-PRESSO majority assumption** - Above 50% pleiotropic, PRESSO removes valid instruments; weighted-mode (plurality-valid) is more robust in that regime.
- **MR-RAPS install** - CRAN-archived 2025-03-01; install from GitHub `qingyuanzhao/mr.raps` and call via `TwoSampleMR::mr_raps()` (thin wrapper) or `mr.raps::mr.raps()` directly. The `MendelianRandomization` package does NOT export `mr_raps()`.
- **MR-Clust biology** - Cluster output is most informative when SNPs in each cluster annotate to distinct pathways; pure statistical clustering without biological story is weak evidence.
- **LHC-MR runtime** - Hours on full sumstats; restrict to LDSC-overlap SNPs and use `n_cores >= 4`.
- **LCV vs MR** - LCV uses genome-wide SNPs not just significant instruments; complements MR but does not replace it; gcp = 0 means rg only, no causation.
- **STROBE-MR** - 20 items + 30 subitems (Skrivankova 2021 JAMA 326:1614; explanation BMJ 375:n2233); reviewer-required at major journals since 2022.
- **Reproducibility** - Pin TwoSampleMR, MRPRESSO, CAUSE versions in the methods section; record harmonization choices.

## Related Skills

causal-genomics/mendelian-randomization - Primary causal estimation that this sensitivity battery validates
causal-genomics/genetic-correlation - LDSC rg required for Step 1 of the decision flow; CHP escalation trigger
causal-genomics/colocalization-analysis - Required for cis-MR drug-target signals where instruments are too few for Egger
causal-genomics/fine-mapping - Identify causal variants underlying instrument loci
causal-genomics/mediation-analysis - Multivariable MR for mediator-adjusted causal estimates
population-genetics/association-testing - GWAS summary statistics underlying MR instruments
clinical-biostatistics/effect-measures - Translate MR estimates to clinical effect measures
