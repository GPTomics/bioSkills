# Hit Calling - Usage Guide

## Overview

Decision-grade cross-method orchestration for calling significant hits in pooled CRISPR screens. Catalogs the seven primary methods (MAGeCK RRA, MAGeCK MLE, BAGEL2, drugZ, JACKS, Chronos, CERES), the experimental design each is built for, the failure modes outside their domain, and the reconciliation logic when methods disagree. Defines tier-based confidence stratification (Tier 1 = 3-method consensus; Tier 2 = 2-method; Tier 3 = single-method exploratory).

## Prerequisites

```bash
# All five primary tools
pip install mageck                                   # RRA + MLE
pip install bagel                                    # BAGEL2 Bayes factor
pip install jacks                                    # JACKS efficacy-aware
# drugZ via PyPI or GitHub
pip install drugz
# or
git clone https://github.com/hart-lab/drugz
# Chronos (DepMap)
pip install chronos-cn
# Custom analyses
pip install pandas numpy scipy statsmodels matplotlib seaborn
```

## Quick Start

Tell the AI agent what to call:
- "Decide which hit-calling method to run on my screen design (2-condition vs time-course vs cancer-line vs drug screen)"
- "Run MAGeCK + BAGEL2 + drugZ on the same data and build a tier-1/tier-2/tier-3 consensus hit list"
- "Reconcile MAGeCK and BAGEL2 disagreement on my essentiality screen"
- "Apply the second-best-sgRNA rule to my hit list to flag single-guide-driven false positives"
- "Pick FDR threshold by screen quality: high-quality (PR-AUC >0.85) -> FDR 0.05; low-quality (PR-AUC 0.5-0.7) -> FDR 0.01"

## Example Prompts

### Method Selection by Design

> "My screen is 5 cancer cell lines vs Day 0 controls across 14 days. Pick Chronos vs MAGeCK MLE vs JACKS and explain why."

> "I have one drug screen: vehicle vs drug at MOI 0.3, 14 days, 3 replicates each. Pick drugZ vs MAGeCK MLE."

> "I have a time-course screen with Day 0 / 7 / 14 / 21. Should I use MAGeCK MLE or run RRA pairwise?"

> "My screen is a combinatorial paired-guide library testing GIs. Decide on MAGeCK MLE with GI scoring vs custom analysis."

### Multi-Method Consensus

> "Run MAGeCK RRA and BAGEL2 on my essentiality screen. Build a 2-method consensus hit list. Use BF >6 as BAGEL2 threshold (≈ FDR 0.05 calibrated). Output Tier 1 (both methods) and Tier 2 (one method) lists separately."

> "Run MAGeCK + BAGEL2 + drugZ on my drug screen. Output the tier-1 consensus (3-method agreement) at FDR <0.05 / BF >6 across all. These hits go to arrayed validation."

> "Compute Spearman ρ between MAGeCK neg|score and BAGEL2 BF. If ρ <0.6, audit why the two methods disagree."

### Reconciliation

> "MAGeCK says ERBB2 is essential in HER2+ SK-BR-3 cells. Chronos says it isn't. Reconcile -- is this real essentiality or copy-number artifact?"

> "JACKS calls 200 hits that MAGeCK does not. Sample the top 20 disagreements and explain which to trust based on per-sgRNA efficacy."

> "BAGEL2 returns BF >12 for a gene that MAGeCK gives neg|fdr >0.5. Investigate the per-sgRNA fold changes; identify which method is correct."

### Threshold Selection

> "Suggest the right hit threshold for my screen quality grade. Pass-quality (CEGv2 PR-AUC >0.85) -> standard FDR; CN-confounded -> apply Chronos; low-quality -> use BAGEL2 with calibrated BF instead of FDR."

> "My screen has PR-AUC 0.72 (passing but not excellent). Should I tighten FDR from 0.05 to 0.01? Make the tradeoff explicit between sensitivity and specificity."

### Custom Hit Calling

> "Run a z-score hit calling using only non-targeting controls as the null. Output gene-level z, p-value, BH-adjusted FDR."

> "Apply the second-best-sgRNA rule: a gene is a hit only if the 2nd-most-extreme sgRNA also passes the threshold. Filter MAGeCK output by this rule."

## What the Agent Will Do

1. Inspect experimental design: condition count, replicates, cell lines, batches, chemistry
2. Reference screen-qc output to confirm screen is interpretable (PR-AUC >0.7)
3. Apply method-decision tree from the SKILL: primary method by design
4. Run primary method; check convergence and output sanity
5. Run secondary method on same data for reconciliation
6. (For high-stakes) run third method for tier-1 consensus
7. Apply tier thresholds: Tier 1 = 3/3, Tier 2 = 2/3, Tier 3 = 1/3
8. Apply second-best-sgRNA rule to flag single-guide-driven hits
9. For cancer-line screens, confirm copy-number correction applied (Chronos or CRISPRcleanR pre-hoc)
10. For multi-batch / multi-screen, confirm batch covariates in MAGeCK MLE or use Chronos
11. Output ranked hit list with tier annotation, per-method FDR, second-best LFC
12. Recommend orthogonal validation strategy (arrayed; orthogonal modality; cell-line panel)

## Tips

- Pick the method that matches the design, not the method you know best. RRA is great for 2-condition screens but fails on time-course; Chronos is great for cancer panels but overkill for single-line screens.
- For high-stakes hit lists (drug-target nomination, paper-level claims), require 2-of-3 or 3-of-3 method consensus. Single-method hits at FDR 0.05 will have ~5% false discoveries; consensus shrinks this dramatically.
- BAGEL2's BF >6 ≈ MAGeCK FDR 0.05 from Hart 2017 calibration. Use this when comparing across methods.
- Cancer-line screens ALWAYS need CN correction. Either pre-process with CRISPRcleanR before MAGeCK / BAGEL2, or use Chronos which models CN jointly.
- The single most common silent failure: not running an essentialome PR-AUC against CEGv2 before hit calling. A screen with PR-AUC <0.5 has no signal regardless of how many "hits" MAGeCK calls.
- For multi-cell-line studies, run per-line analysis first; pool across lines as meta-analysis downstream. Joint MLE across cell lines without indicator covariates dilutes per-line signal.
- The second-best-sgRNA rule is your friend for novel libraries with mixed efficacy: a hit driven by one extreme guide is more likely an outlier than a true effect.
- Heavy-selection drug screens (>40% guides change) break median normalization. Use BAGEL2 (reference-set-anchored, robust) or MAGeCK with `--norm-method control`.

## Decision Cheat Sheet

| If your design is... | Primary method | Backup |
|----------------------|----------------|--------|
| Single 2-condition essentiality | MAGeCK RRA | BAGEL2 |
| Time course | MAGeCK MLE | JACKS |
| Multi-cell-line cancer | Chronos | MAGeCK MLE per line |
| Drug screen | drugZ | MAGeCK MLE |
| Multi-screen, same library | JACKS | MAGeCK MLE |
| Variant function (BE/PE) | Custom + CRISPResso2 | See [[base-editing-analysis]] |
| Combinatorial | MAGeCK MLE + GI scoring | See [[combinatorial-screens]] |
| Single-cell | SCEPTRE | See [[perturb-seq-analysis]] |
| Cancer + multi-batch | Chronos | MAGeCK MLE with batch covariate |

## Confidence Tiers

| Tier | Definition | Validation requirement |
|------|------------|-------------------------|
| Tier 1 | Called by 3/3 orthogonal methods | Minimal: arrayed in matched line |
| Tier 2 | Called by 2/3 methods | Arrayed validation; orthogonal modality |
| Tier 3 | Called by 1 method | Hypothesis; needs further screens before publication |

## Related Skills

- crispr-screens/mageck-analysis - MAGeCK RRA + MLE details
- crispr-screens/bagel-essentiality - BAGEL2 details
- crispr-screens/drugz-chemogenomic - drugZ for drug screens
- crispr-screens/jacks-analysis - JACKS for multi-screen
- crispr-screens/copy-number-correction - Chronos / CRISPRcleanR for cancer-line CN bias
- crispr-screens/screen-qc - Quality gate before hit calling
- crispr-screens/combinatorial-screens - GI scoring
- crispr-screens/perturb-seq-analysis - SCEPTRE for single-cell
- pathway-analysis/go-enrichment - Functional analysis of hit lists
- pathway-analysis/gsea - GSEA on ranked gene lists
