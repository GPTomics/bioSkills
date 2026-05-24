# CRISPR Screen QC - Usage Guide

## Overview

Decision-grade quality control for pooled CRISPR screens. Covers six bottleneck stages (plasmid pool, Day-0 infection, selection, endpoint, biological signal, copy-number artifact), with stage-specific Gini, skew, replicate Pearson/Spearman, sequencing depth, MOI verification, and essentialome PR-AUC against CEGv2 (Hart 2017). Outputs a DepMap-style composite score and recommends remediation for each failing stage.

## Prerequisites

```bash
pip install mageck pandas numpy scipy matplotlib seaborn scikit-learn
# MAGeCK QC dashboard
pip install mageck-vispr
# R dashboard (optional)
R -e "BiocManager::install('MAGeCKFlute')"
```

Required inputs: MAGeCK count output (`screen.count.txt`), plasmid-pool counts (separate file or first sample), known copy-number profile per cell line (from WGS / SNP-array / ASCAT / matched cell-line database), and CEGv2 / NEGv1 reference gene sets (Hart 2017; `hart-lab/bagel` repository).

## Quick Start

Tell the AI agent what to audit:
- "Audit my Brunello-screen counts: plasmid Gini, Day-0 to endpoint dropout, replicate Pearson and Spearman, CEGv2 PR-AUC, depth in reads per sgRNA, MOI verification"
- "Diagnose why hits include ERBB2 in HER2+ SK-BR-3 -- is this copy-number bias or real essentiality"
- "Decide whether my screen passes DepMap quality thresholds, fails, or is salvageable"
- "Recommend which hit-calling method (MAGeCK RRA / MLE / BAGEL2 / Chronos / drugZ) my screen quality grade supports"
- "Generate a composite quality score across all my replicate pairs to gate which conditions enter hit calling"

## Example Prompts

### Library and Plasmid Pool QC

> "Compute Gini coefficient and skew ratio (p90/p10) on my plasmid-pool sequencing. Apply Joung 2017 thresholds: pass at Gini <0.1, skew <2, zero-count <0.5%. Decide whether to proceed."

> "My plasmid Gini is 0.18 and skew is 4.2. Diagnose: PCR over-amplification, synthesis defect, or cloning bottleneck? Recommend remediation."

### Replicate and Depth Audit

> "Compute pairwise Pearson on log10(counts+1) and Spearman on raw ranks between all replicates within each condition. Flag pairs below Pearson 0.85 or Spearman 0.7."

> "One of my treatment replicates shows Pearson 0.78 vs the other two replicates which are >0.95 with each other. Decide whether to drop or rescue the outlier replicate."

> "Verify sequencing depth: reads per sgRNA per sample. Pass at 300x (DepMap), warn at 200x (Joung 2017 minimum), fail below."

### Biological Signal (Essentialome Recovery)

> "Compute PR-AUC against CEGv2 essentials and NEGv1 non-essentials (Hart 2017). My screen passes only if PR-AUC >0.7. Tell me whether this screen has interpretable biology before I run hit calling."

> "PR-AUC is 0.45 even though Gini and Pearson pass. Diagnose: Cas9 not selected pre-screen, premature timepoint, or library targets wrong TSS?"

### Copy-Number Artifact Diagnostic

> "Run Aguirre 2016 copy-number bias diagnostic. Compute Spearman ρ between gene-level LFC and copy number from the matched WGS profile. Flag CN bias if abs(ρ) >0.1 and recommend CRISPRcleanR / CERES / Chronos correction."

> "Hits include ERBB2 in SK-BR-3 and FGFR1 in head-and-neck lines. Confirm whether these are copy-number artifacts before publication."

### Composite Quality Gate

> "Generate the DepMap-style composite quality score across all metrics: Gini-inverse, replicate Pearson minimum, PR-AUC, log-depth, detected-fraction. Use this to decide which conditions enter downstream hit calling."

### MOI Verification

> "Verify infection MOI: from the titration plate (8 wells with 1:2 serial dilution), interpolate the infection efficiency at the volume used in the screen. Compute Poisson P(≥2 sgRNAs/cell) at the resulting MOI."

## What the Agent Will Do

1. Identify each sequencing-sample stage (plasmid, Day 0, selection, endpoint) from metadata
2. Run library_representation() per sample: % zero, % low-count, skew ratio
3. Compute Gini per sample with stage-specific thresholds
4. Compute replicate Pearson on log-counts and Spearman on raw ranks per condition
5. Compute CEGv2 PR-AUC at the endpoint timepoint, vs Day-0 or plasmid
6. Verify MOI from titration data or qPCR of integrated copies
7. Detect copy-number artifact via gene-level LFC vs copy-number correlation
8. Run PCA to visualize sample clustering (condition vs batch)
9. Generate composite quality score and per-sample grade
10. Recommend downstream hit-calling method based on quality grade (high quality -> MAGeCK MLE or Chronos; low quality -> RRA or drugZ; cancer-line -> Chronos with CN correction; in-vivo -> bottleneck-adjusted thresholds)

## Tips

- Plasmid-pool sequencing (Gini <0.1, ≥99% guide detection at >25 reads/guide) is non-negotiable; everything downstream is normalized against this baseline. A screen with un-sequenced plasmid pool is uninterpretable.
- The single most diagnostic metric is CEGv2 PR-AUC; if it passes >0.7 the screen has biology even if individual sample metrics look weak. If it fails <0.5, no remediation in software will fix it -- the screen has no essentiality signal.
- For drug screens, expect Gini to drift up (endpoint Gini 0.3-0.5) as biology drives selection; this is normal. Compare vs vehicle, not Day 0, for any chemogenomic interpretation.
- Cancer-cell-line screens with focal amplification ALWAYS require copy-number correction. Apply CRISPRcleanR / Chronos / CERES preemptively; the Aguirre/Munoz copy-number artifact is universal, not conditional.
- MOI 0.3 is non-negotiable; there is no analytical correction for high-MOI confounding. Re-run rather than try to model around it.
- "Passing" QC at each stage is necessary but not sufficient; the final gate is essentialome PR-AUC.
- High Pearson with low Spearman = a few outlier guides dominate; switch hit calling to RRA or drugZ which use ranks.
- For CRISPRi/a, expect lower per-gene PR-AUC than Cas9 because not all essentials respond to knockdown the way they do to knockout; calibrate against the DepMap CRISPRi sub-essentialome rather than CEGv2.

## QC Decision Reference

| Stage | Primary metric | Pass | Action if fail |
|-------|----------------|------|----------------|
| Plasmid | Gini, skew | <0.1, <2 | Re-clone or re-amplify |
| Day 0 | Pearson with plasmid | >0.9 | Diagnose infection issue |
| Endpoint | Replicate Pearson | >0.85 | Drop outlier replicate |
| Biology | CEGv2 PR-AUC | >0.7 | Cas9 selection / timepoint / TSS |
| CN | Spearman LFC vs CN | abs(ρ) <0.1 | CRISPRcleanR / Chronos |
| Depth | Reads/sgRNA | >300 | Re-sequence |
| MOI | Poisson P(≥2) | <5% | Re-infect |

## Related Skills

- crispr-screens/library-design - Design library with adequate skew margin
- crispr-screens/mageck-analysis - Generate count files for QC
- crispr-screens/copy-number-correction - Remediate CN artifact
- crispr-screens/batch-correction - Address inter-batch drift
- crispr-screens/hit-calling - Pick analysis method by quality grade
- crispr-screens/in-vivo-screens - In-vivo-specific bottleneck QC
- read-qc/quality-reports - General NGS QC upstream
