# crispr-screens

## Overview

Decision-grade analysis of pooled and single-cell CRISPR screens. Covers library design and sgRNA selection (Cas9, CRISPRi/a, Cas12a multiplex, base editor, prime editor); quality control across all six bottleneck stages (plasmid, infection, selection, endpoint, biology, copy-number artifact); algorithmic taxonomy of hit-calling methods (MAGeCK RRA / MLE, BAGEL2, drugZ, JACKS, Chronos); copy-number bias correction (CRISPRcleanR / CERES / Chronos); batch correction; combinatorial / paralog screens; single-cell Perturb-seq (SCEPTRE + Mixscape); base and prime editing variant-function screens; and in vivo bottleneck-aware design.

**Tool type:** mixed | **Primary tools:** MAGeCK, BAGEL2, drugZ, JACKS, Chronos, CRISPRcleanR, CRISPResso2, Pertpy, PRIDICT2

## Skills

| Skill | Description |
|-------|-------------|
| library-design | Cas9/CRISPRi/CRISPRa/Cas12a/BE/PE library design with Rule Set 2, CFD, Dolcetto/Calabrese TSS targeting, Inzolia paralog libraries, and oligo synthesis layout |
| screen-qc | Six-stage QC (plasmid Gini, replicate Pearson, CEGv2 PR-AUC, copy-number artifact diagnostic) with DepMap-style thresholds |
| mageck-analysis | MAGeCK count + RRA + MLE with decision tree, normalization choice (median/total/control), design-matrix construction, MAGeCKFlute + VISPR |
| jacks-analysis | JACKS Bayesian decomposition of LFC into gene-effect x guide-efficacy; multi-screen joint analysis; library efficacy calibration |
| bagel-essentiality | BAGEL2 Bayes Factor classifier with CEGv2/NEGv1 reference calibration; tumor-suppressor sensitivity; precision-recall threshold selection |
| drugz-chemogenomic | drugZ bidirectional Z-score for drug-modifier screens; vehicle-vs-Day-0 reference choice; 2-3x sensitivity gain over MAGeCK on small effects |
| hit-calling | Cross-method decision tree across MAGeCK / BAGEL2 / drugZ / JACKS / Chronos; consensus tier-based confidence; second-best-sgRNA rule |
| copy-number-correction | CRISPRcleanR / CERES / Chronos for cancer-cell-line CN bias (Aguirre 2016 / Munoz 2016); DepMap quarterly standard; CRISPRi bypass alternative |
| batch-correction | ComBat / RUV / SVA / MAGeCK MLE batch covariate; screen-specific batch sources; when correction harms biology |
| crispresso-editing | CRISPResso2 single-amplicon / Batch / Pooled / WGS / Compare modes; Cas9 / BE / PE quantification with editing window math |
| base-editing-analysis | BE variant-function screens (Sanson 2020 GRACE, Hanna 2021 BRCA1/2, Cuella-Martin 2021); CBE vs ABE editing windows; bystander attribution; substitution-vs-indel diagnostic |
| prime-editing-screens | PE screens with PRIDICT2 pegRNA design; PE2/PE3/PEmax/PEAR chemistry choice; chromatin context as primary determinant; MOSAIC saturation; BE cross-validation |
| perturb-seq-analysis | Single-cell pooled CRISPR screens (Perturb-seq, CROP-seq, CITE, Multiome); SCEPTRE calibrated DE; Mixscape escaper filtering; Replogle 2022 genome-scale CRISPRi |
| combinatorial-screens | Big Papi paired Cas9 + enCas12a + in4mer/Inzolia 4-guide arrays; paralog-buffering synthetic lethals; genetic-interaction scoring |
| in-vivo-screens | In vivo bottleneck math; focused library design (Manguso 2017); CRISPR-StAR temporal activation; per-animal meta-analysis |

## Example Prompts

- "Design a CBE saturation library tiling BRCA1 exons 1-23 for variant scanning with 10-15 sgRNAs per amino acid; flag bystander Cs in the editing window"
- "Run MAGeCK MLE on my time-course CRISPR screen with explicit batch covariates; output per-condition beta scores"
- "Apply CRISPRcleanR to my HER2+ SK-BR-3 screen to remove the Aguirre 2016 ERBB2 copy-number artifact before MAGeCK"
- "Compare drugZ vs MAGeCK MLE on my PARPi chemogenomic screen; identify novel sensitizers and resistance genes"
- "Diagnose why my BAGEL2 returns no hits despite known essentials present in the screen"
- "Build an Inzolia-style Cas12a 4-guide array library targeting 400 paralog pairs; include singleton controls for GI scoring"
- "Set up a Replogle-style genome-wide Perturb-seq with CRISPRi at 1,000 cells per perturbation; use SCEPTRE for differential expression"
- "Run a CRISPR-StAR in vivo screen with focused 5,000-gene library in B16 syngeneic; analyze with per-animal MAGeCK meta-analysis"
- "Design pegRNAs for 320 ClinVar MMR gene variants; filter to PRIDICT2 efficiency >50%; validate with parallel BE screen"
- "Compute consensus across MAGeCK + BAGEL2 + drugZ on the same screen; tier-1 hits get arrayed validation"

## Requirements

```bash
# Core CRISPR screen analysis
pip install mageck mageck-vispr bagel-cas9 jacks CRISPResso2 drugz pertpy scanpy
# Single-cell analysis (Perturb-seq)
pip install anndata muon
# Multi-omics
pip install scrublet
# Python dependencies
pip install pandas numpy scipy matplotlib seaborn scikit-learn biopython statsmodels
# Library design
pip install crispor crisprDesignData

# R-based tools
R -e "BiocManager::install(c('CRISPRcleanR', 'MAGeCKFlute', 'sva', 'RUVSeq', 'sceptre'))"
R -e "install.packages(c('Seurat'))"
```

For genome-wide cancer dependency analysis with DepMap-style methodology, also install:
```bash
pip install chronos-cn
```

## Related Skills

- **read-alignment** - Align screen reads to library (BWA)
- **read-qc** - General sequencing QC upstream of screen counting
- **differential-expression** - Statistical concepts shared (NB models, FDR correction)
- **pathway-analysis** - Functional enrichment of screen hit lists (GO, GSEA, KEGG)
- **variant-calling** - Edit annotation downstream of CRISPResso2 (VEP)
- **clinical-databases** - Variant pathogenicity annotation (ClinVar, COSMIC)
- **single-cell** - General single-cell preprocessing for Perturb-seq
- **copy-number** - Source CN profiles for cancer-line correction
- **chip-seq** - dCas9-KRAB analysis related (CRISPRi mechanism)
- **machine-learning** - Biomarker validation of CRISPR-identified targets
