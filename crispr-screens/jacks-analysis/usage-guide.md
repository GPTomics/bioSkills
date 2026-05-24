# JACKS Analysis - Usage Guide

## Overview

Run JACKS (Allen et al 2019 Genome Research) for joint Bayesian decomposition of CRISPR-screen log-fold-changes into per-sgRNA guide-efficacy and per-condition per-gene effect. Designed for multi-screen joint analysis where guide efficacy is shared (same library, same chemistry, cross-cell-line or cross-condition), enabling 2.5x smaller screens and 12-21% lower gene-effect error vs MAGeCK on the same data. Outputs per-gene posterior effect + posterior std and per-sgRNA efficacy + std.

## Prerequisites

```bash
pip install jacks
# or latest from GitHub (run_JACKS.py is at the repo root)
git clone https://github.com/felicityallen/JACKS && cd JACKS && pip install .
# Optional: pre-computed reference efficacy priors
# Download DepMap Brunello / Sanger Score efficacy posteriors for transfer learning
```

Required inputs:
- Count matrix (tab-separated, rows=sgRNA, columns=samples; first columns `sgRNA`, `Gene`)
- Replicate map (Sample, Experiment, Condition) -- tab-separated, no header
- Guide-to-gene map (sgRNA, Gene) -- tab-separated, no header
- Optional: reference efficacy prior file from a matched library

## Quick Start

Tell the AI agent what to do:
- "Run JACKS jointly on my 4 screens (same Brunello library, different cell lines) and output per-line gene effects + shared sgRNA efficacy"
- "Build a guide-efficacy prior from DepMap CRISPR screens and reuse it for my small custom-library screen to reduce sample-size needs"
- "Compare JACKS gene effects vs MAGeCK RRA on the same data; identify where they disagree and why"
- "Identify the low-efficacy sgRNAs (X1 <0.3) in my library for re-design in v2"
- "Decide whether JACKS or Chronos is right for my 10-cell-line cancer dependency screen"

## Example Prompts

### Multi-Screen Joint Analysis

> "I have 4 Brunello screens across HCT116, HEK293T, A375, and MCF7 cell lines, each with 3 replicate Day 0 and Day 14 samples. Run JACKS jointly with `--apply_w_hp` on. Output per-line gene effects, shared efficacy, and gene-level fdr_log10."

> "My screens were done across 6 weeks in two batches. Set up a per-batch JACKS run and compare to a single joint run; quantify whether batch sharing improves or degrades signal."

### Reference Efficacy Transfer

> "Extract the per-sgRNA efficacy posterior from a JACKS run on the DepMap Brunello panel (50 cell lines, ~10,000 screen days). Save as `brunello_efficacy_prior.tsv`. Then run JACKS on my single Brunello screen using this as `--ref_grna_efficacy_file`; verify the 2.5x sample-size reduction matches Allen 2019."

### Library Calibration / Re-design

> "Output the per-sgRNA efficacy from JACKS analysis of my custom screen. Flag sgRNAs with X1 <0.3 as low-efficacy. Group by gene; flag any gene where all guides are low-efficacy for library re-design."

> "I want to design a v2 of my custom library. Use the JACKS efficacy output to drop the bottom 25% of guides and replace with new candidates from CRISPOR."

### Comparison and Diagnostics

> "Compute Spearman ρ between JACKS X1 and MAGeCK neg|lfc on the same dataset. Investigate any rank disagreement >100 positions in the top-1000 hit list."

> "Diagnose why JACKS shows median efficacy 0.18 across my whole library. Is this a chemistry mismatch (CRISPRi screen run with Cas9 defaults), an over-shrinkage from `--apply_w_hp`, or a real library quality issue?"

> "My JACKS run varies between repeated runs (different gene rankings in top 100). Investigate convergence: ELBO trajectory, iteration count, seed."

## What the Agent Will Do

1. Verify input file formats: sgRNA names consistent between counts, guide_map, replicate_map
2. Decide whether joint analysis is appropriate: same library, same chemistry, ≥3 screens
3. If applicable, build the reference efficacy prior from a matched public dataset
4. Run JACKS via CLI `python -m jacks.run_JACKS` or programmatically via `infer_JACKS`
5. Set `--apply_w_hp` based on chemistry (on for Cas9 KO; off if CRISPRi/a requires custom prior)
6. Verify ELBO convergence; if iterations <5000, raise; if still noisy, set seed and increase further
7. Generate gene-level results (X1, X2, fdr_log10) and sgRNA efficacy results (X1, X2)
8. Apply hit thresholds: fdr_log10 <-1 (FDR<0.1); fdr_log10 <-2 (FDR<0.01)
9. Flag low-efficacy guides (X1 <0.3) and genes where all guides are weak (re-design candidates)
10. Cross-validate with MAGeCK / BAGEL2: identify high-confidence hits in agreement, single-tool hits flagged for orthogonal validation
11. Decide if Chronos is preferred (cancer-line multi-cell-line screens with CN bias)
12. Output gene_results.txt, sgrna_efficacy.txt, library_redesign_candidates.txt, comparison_with_mageck.txt

## Tips

- JACKS' core advantage is multi-screen joint analysis with shared efficacy. For a single screen with no public reference, JACKS does not outperform MAGeCK or BAGEL2.
- Always match chemistry: do not share Cas9-KO efficacy with CRISPRi/a guides; efficacy is chemistry-specific. Build separate priors for each.
- The 2.5x sample-size reduction is real but conditional on the reference being from the same library and similar cell context. Reduction is not free.
- For multi-cell-line cancer-line panels, prefer Chronos -- it models CN bias and per-screen quality jointly, which JACKS does not.
- The `X1/X2` ratio (gene-effect divided by its std) is a Bayesian z-equivalent. Use this for ranking rather than X1 alone -- a strong but uncertain effect should rank lower than a moderate but confident one.
- If running on a multi-cell-line panel, fit per-cell-line gene effects but shared efficacy. Pooling effects across cell lines is meta-analysis and should be done downstream after per-cell-line JACKS estimates are obtained.
- Library re-design: drop the bottom 25% efficacy guides; in the v2 library, every gene should have all guides at efficacy >0.4 (Brunello v2 / Avana v2 convention).
- Convergence is the silent failure mode. Always check ELBO trajectory and set seed; never trust a single non-converged run.

## Key Thresholds

| Threshold | Value | Rationale |
|-----------|-------|-----------|
| Hit fdr_log10 | <-1 (FDR <0.1); <-2 (FDR <0.01) | Allen 2019 convention |
| Effective gene signal | X1 <0 AND abs(X1/X2) >2 | Bayesian z-equivalent ≈ 95% credible |
| Low-efficacy guide flag | X1 <0.3 | Allen 2019; below this, guide non-functional |
| Joint analysis screen count | ≥3 | Below this, equivalent to MAGeCK/BAGEL2 |
| Iterations for publication | 5000+ | Verify ELBO plateau |
| Reference for prior reuse | DepMap or Sanger Score panel | ~50 cell lines, ~10k screen days |

## Decision Comparison

| Scenario | JACKS | MAGeCK | BAGEL2 | Chronos |
|----------|-------|--------|--------|---------|
| Single 2-condition screen | OK | Best | OK | N/A |
| Multi-screen joint (same library) | Best | Limited | Per-screen only | N/A |
| Multi-cell-line cancer panel | OK | OK | OK | Best (CN+quality) |
| Library re-design (efficacy info) | Best | No | No | Limited |
| Chemogenomic / drug screen | Suboptimal | OK | Suboptimal | N/A; use drugZ |
| Heavy selection (>40% guides change) | OK | RRA fails; use MLE | Robust | OK |

## Related Skills

- crispr-screens/mageck-analysis - MAGeCK MLE for joint multi-condition design (alternative)
- crispr-screens/bagel-essentiality - BAGEL2 for essentiality classification without efficacy
- crispr-screens/copy-number-correction - Chronos for cancer-line multi-screen + CN bias
- crispr-screens/screen-qc - Replicate Pearson + plasmid Gini gate JACKS use
- crispr-screens/library-design - Use JACKS efficacy output to refine library v2
- crispr-screens/hit-calling - Cross-method decision tree and reconciliation
