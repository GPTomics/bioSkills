# CRISPR Screen Pipeline - Usage Guide

## Overview

End-to-end pipeline for pooled and single-cell CRISPR screens. Takes raw FASTQ files (or pre-counted matrices) through library QC, guide counting, six-stage quality control, optional copy-number correction for cancer lines, optional batch correction for multi-batch screens, design-matched hit calling across MAGeCK RRA/MLE, BAGEL2, drugZ, JACKS, and Chronos, and tier-based consensus. Branches into specialized workflows for single-cell Perturb-seq, combinatorial paralog screens, base-editor and prime-editor variant-function screens, and in vivo bottleneck-aware screens.

## Prerequisites

```bash
# Core hit-calling tools
pip install mageck mageck-vispr bagel-cas9 jacks pertpy scanpy anndata
# drugZ (clone from GitHub)
git clone https://github.com/hart-lab/drugz
# Editing analysis
pip install CRISPResso2
# Chronos for DepMap-style cancer-line CN+quality jointly
pip install chronos-cn
# R packages
R -e "BiocManager::install(c('CRISPRcleanR', 'MAGeCKFlute', 'sva', 'sceptre'))"
# Helpers
pip install pandas numpy scipy matplotlib seaborn scikit-learn biopython statsmodels
```

Required inputs:
- FASTQ files (one per sample) OR pre-counted matrix
- sgRNA library CSV with columns sgRNA, Gene, Sequence
- Sample-to-condition mapping
- For cancer-line screens: copy-number profile (WGS / SNP-array / ASCAT)
- For multi-batch: batch annotation per sample
- For in vivo: per-animal sample annotation

## Quick Start

Tell the AI agent what to run:
- "Run the end-to-end CRISPR screen pipeline on my Brunello dropout screen with vehicle and drug arms"
- "Apply CRISPRcleanR copy-number correction before MAGeCK on my HER2+ SK-BR-3 panel"
- "Pick MAGeCK RRA vs MLE vs BAGEL2 vs drugZ vs JACKS vs Chronos for my screen design and execute"
- "Build tier-based consensus across MAGeCK + BAGEL2 + drugZ for high-stakes drug-target nomination"
- "Branch the pipeline into perturb-seq-analysis for my single-cell CROP-seq dataset"
- "Diagnose CEGv2 PR-AUC under 0.5 in my screen and recommend remediation"

## Example Prompts

### Standard Pooled Screen

> "I have FASTQ from a Brunello dropout screen: plasmid + Day 0 + 3 replicates of Day 14 vehicle + 3 replicates of Day 14 drug. Run guide counting, six-stage QC, then MAGeCK test (drug vs vehicle) and drugZ in parallel. Output tier-2 consensus hits at FDR <0.05."

> "Pick the appropriate hit-calling method for my time-course screen across 5 cell lines. Run MAGeCK MLE with cell-line covariates."

### Cancer Cell Line + CN Bias

> "My screen is in SK-BR-3 (HER2+, ERBB2 ~24 copies). Run CRISPRcleanR pre-hoc, then MAGeCK on corrected counts. Verify ERBB2 is no longer top hit."

> "DepMap-style screen across 12 cancer cell lines, longitudinal sampling, matched WGS CN profiles. Use Chronos as the integrated CN + screen-quality + gene-effect model."

### Chemogenomic Drug-Modifier

> "PARPi sensitivity screen: vehicle vs olaparib at 14 days. Run drugZ with -i -o -c -x flags; expect sensitizers in DDR pathway (BRCA1/2, RAD51) and resistance in PARP1 paradox."

### Time Course / Multi-Condition

> "Build a MAGeCK MLE design matrix for Day 0 / 7 / 14 / 21 timepoints. Run mle; output per-condition beta scores; identify genes with monotonic depletion trends."

### Specialized Branches

> "Branch to perturb-seq-analysis: my data is CROP-seq with 10X 3' direct capture. Run sgRNA assignment, Mixscape escaper filter, then per-perturbation DE."

> "Branch to combinatorial-screens: Inzolia 4-guide-array Cas12a library targeting 400 paralog pairs. Compute GI scores and identify synthetic-lethal pairs at GI z <-2."

> "Branch to base-editing-analysis for my CBE saturation screen of BRCA1 exons 1-23 (Hanna 2021 methodology). Filter to >30% editing-efficient sgRNAs; aggregate to per-variant fitness."

> "Branch to in-vivo-screens: focused 2,500-gene library in B16-OVA syngeneic with 10 mice per condition. Per-animal MAGeCK + Stouffer Z meta-analysis."

### Diagnostics

> "Screen passes Gini and Pearson but CEGv2 PR-AUC is 0.45. Diagnose Cas9 selection failure vs early timepoint vs library positioning vs CN bias and recommend remediation."

> "My MAGeCK RRA says everything is significant at FDR <0.01. Diagnose heavy selection breaking median normalization; recommend switching to --norm-method control with NTCs or BAGEL2."

## What the Agent Will Do

1. Inspect library file and FASTQ filenames; verify naming consistency
2. Run mageck count with sample labels and --trim-5 adapter
3. Inspect countsummary.txt for mapping rate, Gini, % zero per sample
4. Run six-stage QC (per [[crispr-screens/screen-qc]]): plasmid Gini, Day-0 coverage, replicate Pearson/Spearman, sequencing depth, CEGv2 PR-AUC against Hart 2017
5. If cancer cell line: apply CRISPRcleanR or Chronos for CN bias (per [[crispr-screens/copy-number-correction]])
6. If multi-batch: add batch covariate to MAGeCK MLE design matrix (per [[crispr-screens/batch-correction]])
7. Pick hit-calling method by design (RRA for two-condition; MLE for multi-condition; BAGEL2 for essentiality; drugZ for chemogenomic; JACKS for multi-screen joint; Chronos for cancer panels)
8. Execute primary method and at least one orthogonal method
9. Build tier-based consensus: Tier 1 = 3-method agreement, Tier 2 = 2 of 3
10. Cross-reference cancer-line amplification database to flag CN-suspect hits
11. Run MAGeCKFlute or VISPR for visualization + KEGG/Reactome enrichment
12. For specialized branches, hand off to perturb-seq-analysis, combinatorial-screens, base-editing-analysis, prime-editing-screens, or in-vivo-screens
13. Output per-tier hit list, volcano + rank plots, QC report, recommended orthogonal validation strategy

## Tips

- The single most common failure is skipping plasmid pool sequencing. Always sequence the plasmid library before infection and use it as baseline; comparing only to Day 0 conflates cloning bottleneck with biology.
- For drug screens, always use vehicle (DMSO or carrier) as control, NOT Day 0. The drug-vs-Day-0 comparison conflates drug effect with normal proliferation.
- Cancer cell line screens ALWAYS require copy-number correction. The Aguirre 2016 / Munoz 2016 amplicon artifact is universal; ERBB2 in HER2+, MYC in MYC-amplified, FGFR1 in head-and-neck cases are textbook.
- Choose hit calling by experimental design, not by familiarity. MAGeCK RRA is great for two-condition essentiality but fails on time course; Chronos is the DepMap standard for cancer panels but overkill for single-line screens.
- For high-stakes hits (drug-target nomination), require 2-of-3 or 3-of-3 method consensus. Single-method hits at FDR 0.05 carry ~5% false discovery; tier-1 consensus shrinks this dramatically.
- BAGEL2 BF >6 corresponds to FDR 0.05 from Hart 2017 G3 calibration; use this when cross-comparing methods.
- For low-quality screens (CEGv2 PR-AUC 0.5-0.7), tighten FDR from 0.05 to 0.01 to maintain effective specificity.
- Switching from Cas9 to CRISPRi (Dolcetto library) is the cleanest way to bypass copy-number artifact in cancer lines. The tradeoff: CRISPRi knockdown is less complete than Cas9 KO.
- For single-cell screens, MOI 0.3 is non-negotiable; high MOI creates multi-perturbation cells that violate every analysis assumption.
- In vivo screens cannot achieve standard 500x coverage with most syngeneic models. Use focused libraries (3,000-15,000 sgRNAs) or CRISPR-StAR temporal activation.

## Decision Cheat Sheet

| Screen design | Primary method | Secondary check |
|---------------|----------------|-----------------|
| Two-condition essentiality | MAGeCK RRA | BAGEL2 |
| Time course | MAGeCK MLE | JACKS |
| Drug vs vehicle | drugZ | MAGeCK MLE |
| Multi-screen, same library | JACKS | MAGeCK MLE |
| Cancer-line panel | Chronos | MAGeCK MLE per line |
| Single-cell | SCEPTRE (see perturb-seq-analysis) | Mixscape pre-filter |
| Combinatorial / paralog | MAGeCK MLE + GI scoring | See combinatorial-screens |
| Variant function (BE/PE) | CRISPResso2 + MAGeCK | See base/prime-editing-screens |
| In vivo | MAGeCK with animal covariate | Per-animal Stouffer meta |

## QC Checkpoints

| Stage | Check | Pass | Action if fail |
|-------|-------|------|----------------|
| Counting | Mapping rate | >65-70% | Check adapter / library format |
| Plasmid | Gini | <0.1 | Re-amplify or re-clone |
| Endpoint | Replicate Pearson on log-counts | >0.85 | Drop outlier replicate |
| Biology | CEGv2 PR-AUC | >0.7 | Cas9 selection / timepoint / TSS issue |
| CN | Spearman LFC vs CN | abs(rho) <0.10 | Apply CRISPRcleanR / Chronos |
| Depth | Reads per sgRNA | >300 | Re-sequence |
| MOI | Poisson P(>=2) | <5% | Re-infect at lower MOI |

## Output Files

| File | Source | Description |
|------|--------|-------------|
| experiment.count.txt | mageck count | Raw count matrix |
| experiment.countsummary.txt | mageck count | Per-sample QC summary |
| screen_cleanr_corrected_counts.txt | CRISPRcleanR | CN-corrected counts |
| essentiality_rra.gene_summary.txt | mageck test | Gene-level RRA scores + FDR |
| bayes_factor.txt | BAGEL2 | Per-gene Bayes Factor |
| drugz_output.txt | drugZ | Per-direction Z + FDR |
| jacks_out_gene_JACKS_results.txt | JACKS | Gene effect + sgRNA efficacy |
| tier_consensus.csv | Pipeline aggregation | Tier-1/2/3 hit assignment |
| volcano.png, rank.png | Visualization | Standard hit plots |

## Related Skills

- crispr-screens/library-design - Library composition gates downstream analysis
- crispr-screens/screen-qc - Six-stage QC with stage-specific thresholds
- crispr-screens/mageck-analysis - MAGeCK RRA + MLE detail
- crispr-screens/bagel-essentiality - BAGEL2 Bayes factor classifier
- crispr-screens/drugz-chemogenomic - drugZ for chemogenomic screens
- crispr-screens/jacks-analysis - JACKS for joint multi-screen analysis
- crispr-screens/hit-calling - Cross-method decision tree + tier consensus
- crispr-screens/copy-number-correction - CRISPRcleanR / Chronos for cancer lines
- crispr-screens/batch-correction - Multi-batch covariate design
- crispr-screens/crispresso-editing - Editing quantification
- crispr-screens/base-editing-analysis - Variant-function BE screens
- crispr-screens/prime-editing-screens - PE screens with PRIDICT2
- crispr-screens/perturb-seq-analysis - Single-cell screens
- crispr-screens/combinatorial-screens - Cas12a paralog screens
- crispr-screens/in-vivo-screens - Bottleneck-aware in vivo design
- pathway-analysis/go-enrichment - Hit-list functional enrichment
- pathway-analysis/gsea - Pre-ranked GSEA
