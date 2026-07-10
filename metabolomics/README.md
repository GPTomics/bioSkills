# metabolomics

## Overview

Decision-grade LC-MS and GC-MS metabolomics, from raw data to biological interpretation. Each skill makes the field's make-or-break traps explicit: the feature table is a model-dependent artifact, an annotation is a confidence level not an identification, normalization can create or erase signal, a score plot proves nothing without permutation-validated Q2, pathway enrichment launders annotation uncertainty into biology, a lipid name overstates its structural resolution, an absolute concentration is only as good as its internal standard, and labeling reports flux while pool size does not.

**Tool type:** mixed | **Primary tools:** XCMS, MS-DIAL, matchms, SIRIUS, ropls, MetaboAnalystR, lipidr, IsoCor

## Skills

| Skill | Description |
|-------|-------------|
| xcms-preprocessing | Programmatic R feature detection (modern XcmsExperiment API); the feature table as a parameter-dependent artifact |
| msdial-preprocessing | MS-DIAL MS2Dec deconvolution, DDA/DIA, GC-EI scope, and honest import of alignment output |
| metabolite-annotation | Confidence-level identification with MSI/Schymanski levels, spectral matching, SIRIUS, the isomer wall |
| normalization-qc | QC design, drift correction, compositional-aware normalization (PQN), D-ratio filtering, mechanism-aware imputation |
| statistical-analysis | Scaling, permutation-validated PLS-DA/OPLS-DA, univariate testing, and correlation-aware multiple testing |
| pathway-mapping | ORA/MSEA on identified compounds vs mummichog/PSEA on m/z peaks; background-set control; pool vs flux |
| lipidomics | Lipid nomenclature/structural-resolution honesty, class-internal-standard quantification, in-source-fragment traps |
| targeted-analysis | MRM/PRM absolute quantification, SIL internal standards, weighted calibration, ICH M10 validation |
| isotope-tracing | Stable-isotope-resolved metabolomics: 13C/15N tracers, MID, natural-abundance correction, flux vs pool |

## Example Prompts

- "Process my raw LC-MS data into a feature table with XCMS and explain which parameters set what is detectable"
- "Annotate these features and tell me the defensible MSI confidence level for each"
- "Normalize my data using QC samples and a compositional-aware method, with drift correction"
- "Find differential metabolites and validate the OPLS-DA model with a permutation test"
- "Run mummichog on my m/z peaks using the full feature table as the background"
- "Quantify these targets against an isotope-labeled internal standard with a weighted calibration curve"
- "Correct my 13C-glucose isotopologue data for natural abundance and compute the labeling pattern"

## Requirements

```r
# R/Bioconductor
BiocManager::install(c("xcms", "Spectra", "MsExperiment", "CAMERA", "lipidr", "ropls", "pmp", "structToolbox", "imputeLCMD", "SummarizedExperiment"))

# MetaboAnalystR (from GitHub)
devtools::install_github("xia-lab/MetaboAnalystR")

# missForest, FELLA for imputation and network enrichment
install.packages(c("missForest"))
BiocManager::install("FELLA")

# MS-DIAL: download from https://systemsomicslab.github.io/compms/msdial/main.html (GUI/console application, not an R package)
# SIRIUS: download from https://bio.informatik.uni-jena.de/software/sirius/ (standalone CLI/GUI; account/license required)
```

```bash
# Python
pip install matchms numpy pandas scipy statsmodels matplotlib isocor pygoslin
```

## Related Skills

- **proteomics** - Similar MS-based acquisition, identification, and quantification workflows
- **multi-omics-integration** - Integrate metabolite features with other omics layers
- **pathway-analysis** - Gene-set ORA/GSEA concepts that parallel metabolite enrichment
- **systems-biology** - Constraint-based flux-balance modeling, complementary to empirical isotope tracing
- **experimental-design** - Randomization, blocking, and batch design upstream of any metabolomics study
