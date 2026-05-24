# MAGeCK Analysis - Usage Guide

## Overview

Run MAGeCK on pooled CRISPR screens for genome-wide essentiality, drug-modifier, time-course, multi-cell-line, and combinatorial designs. Covers `mageck count` from FASTQ, the alpha-RRA test (`mageck test`) for two-condition comparisons, the maximum-likelihood model (`mageck mle`) with explicit design matrices for multi-condition / time-course / batch-aware designs, sgRNA-efficiency injection, normalization strategy (median / total / control / spike-in), and downstream visualization with MAGeCKFlute (R) and MAGeCK-VISPR (Python).

## Prerequisites

```bash
pip install mageck                                  # primary CLI
pip install mageck-vispr                            # interactive QC + dashboard
R -e "BiocManager::install('MAGeCKFlute')"         # R-based downstream visualization
# Genome reference (for off-target during count if --map-genome used)
conda install -c bioconda bowtie  # MAGeCK uses for off-target check at design time
```

Required inputs: raw FASTQ (one per sample), sgRNA library CSV with `sgRNA,Gene,Sequence` header, and a sample-condition mapping. For MLE designs, additionally a design-matrix file with `Samples` plus one column per condition (each a 0/1 indicator), and the `baseline` column always set to 1.

## Quick Start

Tell the AI agent what to analyze:
- "Run MAGeCK count + test on my drug vs vehicle screen with median normalization"
- "Build an MLE design matrix for a time course (Day 0 / Day 7 / Day 14 / Day 21) and run MAGeCK MLE"
- "Compare MAGeCK MLE beta scores with JACKS efficacy-adjusted scores on the same data"
- "Diagnose why my MLE output has NaN beta scores for half the genes"
- "Pick MAGeCK RRA vs MLE vs BAGEL2 vs drugZ vs Chronos for my experimental design"
- "Apply batch-aware MLE on a multi-cell-line panel with cell-line and batch as covariates"

## Example Prompts

### Two-Condition (RRA)

> "Run MAGeCK count then test on my screen. Samples Plasmid, Day0, Veh_r1-3, Drug_r1-3. Use `--norm-method median` first; if essentialome PR-AUC is low, retry with `--norm-method control` using ntcs.txt. Output gene_summary at FDR <0.05 and LFC magnitude >1 in both directions."

> "Generate a volcano plot of my MAGeCK test output highlighting negative-selection hits (depleted) and positive-selection hits (enriched) at FDR <0.05."

### Multi-Condition (MLE)

> "Build a time-course design matrix with baseline + day7 + day14 + day21 columns. Run mageck mle, output per-condition beta + Wald-FDR. Identify genes with monotonic beta trends across timepoints."

> "I have a multi-cell-line screen across 5 cancer lines and 3 batches. Build the MLE design matrix with cell-line indicators + batch indicators + treatment indicator. Recommend Chronos as an alternative if quality / scale warrants."

### Drug-Modifier Screen

> "Run MAGeCK test on a chemogenomic screen comparing drug vs vehicle (not Day 0). Use `--norm-method control` with NTCs. Compare to drugZ output: which method is more sensitive for synthetic-lethal interactions?"

> "My drug screen has selection at MOI 0.3 then vehicle vs drug 14 days post-infection. Pick MAGeCK MLE with dose covariate or drugZ; explain the tradeoff."

### Diagnostics

> "Diagnose why my MAGeCK MLE shows NaN beta scores for 40% of genes. Check coverage, sgRNA dropout, and design-matrix singularity."

> "My MAGeCK RRA outputs every gene at FDR <0.01 in a heavy-selection screen. Diagnose and recommend a different normalization or hit-calling method."

> "Reconcile MAGeCK and BAGEL2 disagreement on my essentiality screen. Which set of hits should I trust for follow-up?"

### Visualization and Pathway Analysis

> "Run MAGeCKFlute FluteRRA on my output. Generate volcano + rank + square plot + KEGG enrichment for negative selection."

> "Launch MAGeCK-VISPR on the screen to generate an interactive QC dashboard."

## What the Agent Will Do

1. Verify library CSV format and sample-to-FASTQ mapping
2. Run `mageck count` with sample labels and `--trim-5` adapter, capturing the per-sample QC summary
3. Inspect Gini, mapping rate, % zero-count from countsummary.txt; if any fails, refer to screen-qc skill
4. Decide RRA vs MLE based on experimental design (decision tree in the SKILL)
5. For RRA: run `mageck test` with appropriate `--treatment-id`, `--control-id`, `--norm-method` (median default; control if heavy selection)
6. For MLE: build design matrix from condition assignments; ensure baseline column is all-1; one column per condition or covariate; run `mageck mle`
7. For drug screens: vehicle (not Day 0) as control; use NTCs for normalization
8. Optionally inject sgRNA efficiency from JACKS run on matching cell line
9. Filter results at FDR <0.05 and LFC threshold appropriate for screen type
10. Generate volcano + rank + replicate plots; run FluteRRA / FluteMLE for pathway analysis
11. Reconcile with alternative tools (BAGEL2 for essentiality, drugZ for chemogenomic, Chronos for cancer-line panels)
12. Output the gene_summary, sgrna_summary, an audit report explaining normalization choice, hit count at FDR thresholds, and method-comparison summary

## Tips

- Always sequence the plasmid pool and use it as baseline; comparing to Day 0 of the screen pool conflates cloning bottleneck with biology.
- For chemogenomic / drug screens, use vehicle (matched DMSO or carrier) as control, not Day 0; drug-induced shifts only mean against vehicle.
- Median normalization fails when >40% of guides change. Symptom: every gene significant. Switch to `--norm-method control` with NTCs.
- For multi-cell-line, multi-batch screens, MAGeCK MLE works but Chronos is the DepMap standard with built-in CN bias and screen-quality modeling.
- NaN beta scores in MLE indicate a fitting failure (typically <2 non-zero sgRNAs per gene per condition). Exclude these from interpretation, do not treat as zero.
- For low-effect-size screens (FDR <0.1 hits at LFC 0.3-0.5), MAGeCK results need orthogonal validation; the statistical model is sensitive enough to call hits that don't replicate.
- `--variance-estimation-samples` lets you fit dispersion from early-timepoint samples (when most guides are still neutral), then test against later samples; useful for low-replicate-number screens.
- For paired designs (each replicate from a matched donor), add a paired-replicate column to the MLE design matrix.

## Standard Thresholds

| Threshold | Value | Rationale |
|-----------|-------|-----------|
| Hit FDR (gene-level) | <0.05 | Standard publication threshold |
| Hit LFC magnitude | >1 (2-fold) | Biologically interpretable effect size |
| Replicate Pearson on log-counts | >0.85 (acceptable), >0.95 (excellent) | Pre-test QC |
| sgRNA mapping rate | >65% | Library quality |
| Gini coefficient (plasmid) | <0.1 | See [[screen-qc]] |
| sgRNAs needed for MLE per gene per condition | ≥3 with non-zero counts | Below this, NaN beta |
| Permutation iterations for RRA gene-level p-value | 1,000 default; 10,000 for tight FDR | Trades runtime |

## Method Comparison Cheat Sheet

| Screen design | Primary method | Backup / sensitivity check |
|---------------|----------------|----------------------------|
| Two-condition essentiality | MAGeCK RRA or BAGEL2 | The other |
| Time course | MAGeCK MLE | JACKS |
| Multi-cell-line panel | Chronos or MAGeCK MLE | BAGEL2 per-line then meta-analyze |
| Drug screen | drugZ | MAGeCK MLE with dose covariate |
| Combinatorial / paired guide | MAGeCK MLE with GI scoring | See [[combinatorial-screens]] |
| Single-cell screen (Perturb-seq) | SCEPTRE | See [[perturb-seq-analysis]] |

## Related Skills

- crispr-screens/screen-qc - Pre-MAGeCK QC including Gini, replicate, PR-AUC
- crispr-screens/library-design - Library design dictates normalization choice
- crispr-screens/jacks-analysis - sgRNA-efficacy-aware analysis
- crispr-screens/bagel-essentiality - Bayes-factor essentiality alternative
- crispr-screens/drugz-chemogenomic - drugZ for chemogenomic / drug-modifier screens
- crispr-screens/hit-calling - Cross-method decision tree and reconciliation
- crispr-screens/copy-number-correction - CRISPRcleanR / Chronos for cancer-line correction
- crispr-screens/batch-correction - Building batch-aware MLE design matrices
- pathway-analysis/gsea - Downstream GSEA on MAGeCK ranks
