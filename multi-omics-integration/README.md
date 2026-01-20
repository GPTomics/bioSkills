# multi-omics-integration

## Overview

Statistical integration of multiple omics data types (transcriptomics, proteomics, metabolomics, etc.) for biological discovery.

**Tool type:** r | **Primary tools:** MOFA2, mixOmics, SNF

## Skills

| Skill | Description |
|-------|-------------|
| mofa-integration | Multi-Omics Factor Analysis with MOFA2 |
| mixomics-analysis | Multivariate integration with mixOmics (sPLS, DIABLO) |
| similarity-network | Similarity Network Fusion for patient stratification |
| data-harmonization | Cross-omics preprocessing and batch effect handling |

## Example Prompts

- "Integrate my RNA-seq and proteomics data with MOFA2"
- "Find features correlated across omics layers with sPLS"
- "Build a patient similarity network from multi-omics data"
- "Harmonize my transcriptomics and methylation datasets"
- "Identify multi-omics signatures using DIABLO"

## Requirements

```r
# R/Bioconductor
BiocManager::install(c("MOFA2", "mixOmics", "SNFtool", "MultiAssayExperiment"))

# Additional packages used in examples
install.packages(c("survival", "survminer", "igraph", "pheatmap"))
BiocManager::install(c("DESeq2", "sva", "clusterProfiler", "msigdbr"))
```

```bash
# Python
pip install mofapy2
```

## Related Skills

- **differential-expression** - Single-omics DE analysis
- **proteomics** - Proteomics-specific workflows
- **methylation-analysis** - Epigenomics analysis
- **pathway-analysis** - Functional interpretation of integrated results
