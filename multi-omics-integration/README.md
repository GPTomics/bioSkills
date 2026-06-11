# multi-omics-integration

## Overview

Joint statistical integration of two or more BULK omics blocks on a shared sample axis: deciding the integration strategy, MOFA2 unsupervised factor discovery, mixOmics supervised multi-block discriminant analysis (DIABLO/sPLS/MINT), Similarity Network Fusion for patient stratification, and cross-omic harmonization. The framing is decision-grade: bulk multi-omics is small-n/large-p so an unvalidated integrated signature is the default noise result, vertical and horizontal integration are different problems, MOFA factors are unsupervised so a factor need not align with the phenotype, DIABLO signatures are cross-validated hypotheses not truths, SNF subtypes are hyperparameter-dependent claims to defend, and per-omic normalization belongs upstream. Single-cell multimodal integration is a different paradigm and lives in the single-cell category.

**Tool type:** r | **Primary tools:** MOFA2, mixOmics, SNFtool

## Skills

| Skill | Description |
|-------|-------------|
| integration-design | Method selection, sample correspondence, the n<<p discipline, and the variance-imbalance diagnostic |
| data-harmonization | Cross-omic preprocessing, scaling, batch strategy, and missing-data triage |
| mofa-integration | Unsupervised latent-factor discovery across omics blocks with MOFA2 |
| mixomics-analysis | Supervised multivariate integration (sPLS, DIABLO, MINT) with mixOmics |
| similarity-network | Similarity Network Fusion for patient stratification |

## Example Prompts

- "Which integration method fits my RNA, protein, and methylation data on the same patients?"
- "Integrate my RNA-seq and proteomics with MOFA2 and tell me which view drives each factor"
- "Find a cross-omic signature that discriminates my groups with DIABLO"
- "Stratify my patients into subtypes by fusing similarity networks"
- "Harmonize my omics blocks so no single one dominates the integration"

## Requirements

```r
# R/Bioconductor
BiocManager::install(c("MOFA2", "mixOmics", "SNFtool", "MultiAssayExperiment", "sva"))

# Additional packages used in examples
install.packages(c("igraph", "pheatmap", "imputeLCMD"))
```

```bash
# Python (MOFA2 training backend)
pip install mofapy2 muon
```

## Related Skills

- **single-cell** - Single-cell multimodal (CITE-seq/Multiome) integration (different paradigm)
- **differential-expression** - Per-omic RNA-seq normalization and batch correction
- **machine-learning** - Cross-validation and biomarker-panel theory
- **pathway-analysis** - Enrichment of integrated factors and signatures
- **clinical-biostatistics** - Survival validation of discovered subtypes
