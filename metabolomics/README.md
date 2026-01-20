# metabolomics

## Overview

LC-MS and GC-MS metabolomics analysis from raw data to metabolite identification and pathway interpretation.

**Tool type:** r | **Primary tools:** XCMS, MetaboAnalystR, lipidr, MS-DIAL

## Skills

| Skill | Description |
|-------|-------------|
| xcms-preprocessing | Peak detection, alignment, and grouping with XCMS |
| metabolite-annotation | Metabolite identification and database matching |
| normalization-qc | QC-based normalization and batch correction |
| statistical-analysis | Univariate and multivariate statistics for metabolomics |
| pathway-mapping | Map metabolites to KEGG/Reactome pathways |
| lipidomics | Lipid-specific analysis with lipidr |
| targeted-analysis | Absolute quantification with standard curves |
| msdial-preprocessing | MS-DIAL export processing and integration |

## Example Prompts

- "Process my LC-MS metabolomics data with XCMS"
- "Identify metabolites from m/z and retention time"
- "Normalize my data using QC samples"
- "Find differentially abundant metabolites between groups"
- "Map my significant metabolites to KEGG pathways"

## Requirements

```r
# R/Bioconductor
BiocManager::install(c("xcms", "MSnbase", "CAMERA", "lipidr", "SummarizedExperiment"))

# MetaboAnalystR (from GitHub)
devtools::install_github("xia-lab/MetaboAnalystR")

# Additional packages for statistical analysis
install.packages(c("limma", "pheatmap"))

# MS-DIAL: Download from https://systemsomicslab.github.io/compms/msdial/main.html
# (GUI application, not an R package)
```

```bash
# Python
pip install pyopenms matchms
```

## Related Skills

- **proteomics** - Similar MS-based workflows
- **multi-omics-integration** - Integrate with other omics
- **pathway-analysis** - Enrichment analysis concepts
