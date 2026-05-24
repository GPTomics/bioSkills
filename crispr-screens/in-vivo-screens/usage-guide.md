# In Vivo CRISPR Screens - Usage Guide

## Overview

Decision-grade design and analysis of in vivo CRISPR screens. Covers the bottleneck math (cell number limits force focused libraries); Manguso 2017 immune-evasion methodology; Chen 2015 tumor screens; CRISPR-StAR temporal activation (Uijttewaal 2025 Nat Biotechnol 43:1848) for escaping early bottlenecks; syngeneic vs xenograft vs PDX choice; per-animal clonal variability; tumor-explant DNA extraction; per-animal meta-analysis; in vivo CEGv2 calibration limitations.

## Prerequisites

```bash
pip install mageck pandas numpy scipy
# Optional: focused library annotations from Manguso 2017, Joung 2017
# Tumor genomic DNA extraction kit (proteinase K + column purification)
```

Required inputs:
- Animal cohorts (n=10+ per condition for hit calling)
- Focused library (3,000-15,000 sgRNAs covering 500-3,000 genes)
- Cas9+ cell line (selected by FACS before infection)
- Plasmid pool sequencing as baseline
- Per-animal tumor DNA + sgRNA amplification primers

## Quick Start

Tell the AI agent what to do:
- "Design a focused in vivo CRISPR library targeting immune-evasion biology following Manguso 2017 methodology: 2,000 kinases + surface proteins + immune factors, 4 sgRNAs/gene"
- "Compute bottleneck math: my B16 syngeneic model can take 2M cells; pick library size for 100x coverage maintainable through bottleneck"
- "Apply CRISPR-StAR temporal activation to enable a 15,000-gene library in vivo: cells implanted at Day 0, Cas9 induced at Day 5"
- "Run MAGeCK MLE on my in vivo screen with animal-as-batch covariate; output per-condition beta scores after batch adjustment"
- "Diagnose: why is my in vivo CEGv2 PR-AUC only 0.45 despite passing all in vitro QC?"

## Example Prompts

### Library Design for In Vivo

> "Design a focused library targeting immune-evasion biology in syngeneic B16-OVA mouse melanoma. Include: kinases, cell surface proteins, immune-regulatory genes (Manguso 2017 selection criteria). Total ~2,000 genes x 4 sgRNAs = 8,000 sgRNAs. Verify maintainable at 100x coverage with 2M cells implanted (= 5x cells/sgRNA)."

> "For a syngeneic colorectal cancer model with maximum 3M cells implantable: design library at 5,000 sgRNAs (= 600x effective coverage). Pick genes by pathway relevance (metabolism, immune, proliferation)."

### CRISPR-StAR

> "Set up CRISPR-StAR (Uijttewaal 2025) screen with Tet-on Cas9 in syngeneic model. Cells infected with library at MOI 0.3; implanted in 5 mice per condition; doxycycline induction Day 5 post-implant; tumor harvest Day 21. Compute expected coverage maintenance vs constitutive-Cas9 screen."

### Per-Animal Analysis

> "Run MAGeCK on each animal vs plasmid pool separately; meta-analyze with Stouffer's Z; identify genes consistent across animals. Compare to single combined-animal analysis."

> "Run MAGeCK MLE with animal as batch covariate. Output per-condition beta after batch adjustment. Compare to per-animal RRA + meta-analysis."

### Diagnostics

> "In vivo CEGv2 PR-AUC is 0.45 despite in vitro on same line showing 0.85. Diagnose: context-specific essentialome, clonal dominance, or Cas9 selection failure?"

> "My screen has 100x effective coverage at endpoint but only 60% of library detected. Investigate: bottleneck at engraftment or PCR amplification?"

### Cross-Validation

> "Validate top in vivo hits in vitro (matched cell line, no animal context). Identify hits that are in-vivo-specific (require tumor microenvironment) vs cell-intrinsic."

> "Arrayed validation of top 10 in vivo hits with n=10 mice each. Confirm consistent effect across biological replicates."

## What the Agent Will Do

1. Determine model: syngeneic vs xenograft vs PDX
2. Calculate bottleneck-adjusted library size: max cells implantable × target coverage = total sgRNAs
3. Design focused library following Manguso 2017 / Joung 2017 conventions
4. Plan animal cohort: n=10+ per condition for hit calling
5. Use CRISPR-StAR if genome-scale library needed (Tet-on Cas9 + post-engraftment induction)
6. Verify Cas9+ cell selection before infection (FACS or selection marker)
7. Sequence plasmid pool as baseline
8. Implant cells; allow tumor growth 12-21 days
9. Harvest tumors; extract DNA via proteinase K + column purification
10. PCR amplify sgRNA cassette + sequence
11. Run MAGeCK count + MLE with animal as batch covariate (or RRA per animal + meta-analysis)
12. QC: in vivo-adjusted thresholds (CEGv2 PR-AUC >0.4 acceptable; lower than in vitro)
13. Per-condition hit calling; cross-validate with in vitro
14. Output: per-condition gene effects, per-animal consistency, arrayed validation plan

## Tips

- The single most common failure: trying to use a genome-wide library (70k+ sgRNAs) in vivo. Even at 10M cells implanted, this gives 140x coverage which collapses to 20-50x at endpoint after bottleneck. Use focused libraries unless implementing CRISPR-StAR.
- For immune-targeting screens, syngeneic models are mandatory; xenografts have impaired immunity.
- Per-animal clonal dynamics drive enormous inter-animal variability. Use n=10+ animals per condition; meta-analyze across animals rather than treating as single experiment.
- Cas9 selection before implantation is non-negotiable. Cas9-negative cells persist with their sgRNA but no editing, diluting all signal.
- In vivo CEGv2 calibration is poor because in vitro essentialome doesn't capture tumor-microenvironment biology. Use cell-type-and-context-specific reference essentialome (a matched in vitro screen in the same cell type).
- CRISPR-StAR is the modern approach for genome-scale in vivo: delay Cas9 expression until after engraftment, then induce. Uijttewaal 2025 (Nat Biotechnol 43:1848) reports intrinsic per-clone control and outperforms conventional screens in therapy-resistant melanoma models.
- Multiple animals per condition is more important than depth per animal. n=10 mice at 100x coverage is better than n=3 at 500x.
- Tumor heterogeneity arises during growth; sample multiple regions or pool whole-tumor DNA.
- For metastasis screens, each metastatic site is a separate selection event; analyze per-site.

## Decision Cheat Sheet

| Model | Use case |
|-------|----------|
| Syngeneic | Tumor-immune interaction, checkpoint biology |
| Xenograft | Human cancer cell-intrinsic biology |
| PDX | Patient-specific drug testing |
| Humanized mouse | Tumor-immune in human context |

## Thresholds

| Threshold | Value | Rationale |
|-----------|-------|-----------|
| Cells per animal (syngeneic) | 1-5M typical | Tumor model dependent |
| Library size for focused in vivo | 3,000-15,000 sgRNAs | Maintainable coverage |
| sgRNAs per gene in library | 4-6 | Standard |
| Coverage at endpoint | >50x; ideally 100-200x | Bottleneck-adjusted |
| Animals per condition for hit calling | 10+ | Inter-animal variability |
| In vivo CEGv2 PR-AUC | >0.4 acceptable | Lower than in vitro |
| Days to harvest | 12-21 days post-implant | Time for selection |

## Validation Checklist

- [ ] Library size matches bottleneck math
- [ ] Cas9+ cells selected before infection
- [ ] Plasmid pool sequenced as baseline
- [ ] n ≥10 animals per condition
- [ ] Per-animal sequencing depth ≥100 reads/sgRNA
- [ ] MAGeCK MLE with animal-as-batch covariate
- [ ] Cross-validation with in vitro screen
- [ ] Arrayed validation of top hits in matched cohort

## Related Skills

- crispr-screens/library-design - Focused library design
- crispr-screens/mageck-analysis - MAGeCK MLE with batch covariate
- crispr-screens/hit-calling - Per-animal meta-analysis
- crispr-screens/screen-qc - In vivo-specific QC
- crispr-screens/batch-correction - Animal cohort as batch
- crispr-screens/combinatorial-screens - In vivo combinatorial
- crispr-screens/copy-number-correction - Cancer-line in vivo
- pathway-analysis/go-enrichment - Functional analysis
