# Metabolomics Pipeline Usage Guide

## Overview

This workflow orchestrates an untargeted LC-MS metabolomics study from raw mzML to enriched pathways, chaining five component skills: xcms feature extraction, QC/drift/normalization, confidence-stratified annotation, permutation-validated statistics, and background-aware pathway mapping. It is a sequencer with honest handoffs -- each stage defers to its component skill for parameters and traps, and the workflow's job is to keep the failure modes from one stage from silently corrupting the next. Stable-isotope flux tracing is a separate branch (see metabolomics/isotope-tracing), not part of this pipeline.

## Prerequisites

```r
BiocManager::install(c('xcms', 'CAMERA', 'MetaboAnalystR', 'ropls', 'pmp', 'imputeLCMD'))
```

Conceptual prerequisites: pooled QC samples injected ~1 every 5-10 samples (drift correction is impossible without them); a sample-metadata table with condition, batch, and injection order; biological groups randomized across batches (a confounded design cannot be rescued by any algorithm); and an understanding that a feature table is a parameterized hypothesis, not raw truth.

## Quick Start

Tell your AI agent what you want to do:
- "Run the full untargeted LC-MS metabolomics pipeline on my mzML files"
- "Process my LC-MS data, correct for drift, and find differential metabolites"
- "Take my feature table through QC, statistics, and pathway mapping"
- "I have no compound IDs -- run mummichog pathway analysis on my m/z peaks"

## Example Prompts

### End-to-End
> "I have centroided mzML files plus pooled QCs from an untargeted study; run xcms preprocessing, QC-correct, find differential metabolites, and map pathways."

> "Process my LC-MS data with the modern xcms 4.x API and report differential metabolites with honest annotation confidence."

### QC and Normalization
> "Correct injection-order drift against my QC samples and confirm it did not absorb biological signal."

> "Filter my feature table by QC RSD and D-ratio, PQN-normalize, and impute the residual missing values by mechanism."

### Statistics and Pathways
> "Run a permutation-validated OPLS-DA and a univariate FDR analysis, then reconcile them."

> "My features have no compound IDs -- run mummichog pathway analysis using the full feature table as background."

## What the Agent Will Do

1. Stage 1 -- Extract features with xcms 4.x (`readMsExperiment` -> `findChromPeaks` -> `adjustRtime` -> group -> `fillChromPeaks` -> `featureValues`), aligning to pooled QC and flagging filled values.
2. Stage 2 -- Filter junk features, correct within-batch drift (QCRSC), apply RSD/D-ratio filters, PQN-normalize, and impute the sparse residual holes by missingness mechanism.
3. Stage 3 -- Annotate features and attach an MSI/Schymanski confidence level to each name, collapsing ion families first.
4. Stage 4 -- Run a univariate test with BH FDR AND a permutation-validated OPLS-DA, then reconcile.
5. Stage 5 -- Map identified compounds with ORA (assay-coverage background) or raw m/z with mummichog (full-table background), reporting coverage and confidence ceilings.

## When to Use This Pipeline

- Untargeted LC-MS/MS metabolite profiling
- Metabolic biomarker discovery and treatment-response studies
- Lipidomics (adjust peak widths and annotation; see metabolomics/lipidomics)
- Studies where the feature table must reach pathway interpretation honestly

Not for: stable-isotope flux/fluxomics (see metabolomics/isotope-tracing) or absolute targeted quantification as the primary goal (see metabolomics/targeted-analysis).

## Required Inputs

1. Raw MS data -- centroided mzML/mzXML (convert vendor formats with ProteoWizard msConvert; centroid during conversion)
2. Sample metadata -- CSV with sample, condition, batch, injection_order, and a sample_group column marking QCs
3. Pooled QC samples -- bracketing the run and injected ~1 every 5-10 samples

## Sample Metadata Format

```csv
sample,sample_group,condition,batch,injection_order
QC1.mzML,QC,QC,1,1
Sample1.mzML,Control,Control,1,2
Sample2.mzML,Treatment,Treatment,1,3
QC2.mzML,QC,QC,1,4
```

## Tips

- Pooled QC samples are critical: drift correction, RSD/D-ratio filtering, and PQN reference all depend on them. No QCs, no honest pipeline.
- Validate drift correction on held-out QCs and dilution linearity, never on "QCs cluster tighter" -- that metric is exactly what an over-flexible spline games.
- A clean PLS-DA/OPLS-DA score plot is the generic output of p>>n data; require Q2 high AND a small permutation p before believing separation.
- Never promote a database hit to an identification: a name without an MSI level is incomplete, and pathway enrichment launders that uncertainty into confident biology.
- The pathway background IS the null hypothesis: ORA needs the assay-coverage metabolome, mummichog needs the entire feature table.
- Switch the front end to MS-DIAL (metabolomics/msdial-preprocessing) for MS2Dec deconvolution, GC-EI, or DIA/SWATH, then enter at Stage 2.

## Related Skills

- metabolomics/xcms-preprocessing - Stage 1 feature extraction parameters and the feature-table-as-artifact framing
- metabolomics/normalization-qc - Stage 2 drift correction, RSD/D-ratio filtering, PQN, mechanism-aware imputation
- metabolomics/metabolite-annotation - Stage 3 MSI/Schymanski confidence levels
- metabolomics/statistical-analysis - Stage 4 permutation-validated multivariate and dependence-aware FDR
- metabolomics/pathway-mapping - Stage 5 ORA vs mummichog and background construction
- metabolomics/msdial-preprocessing - Alternative front end entering at Stage 2
- metabolomics/lipidomics - Lipid-specific peak widths and annotation
- metabolomics/targeted-analysis - Absolute quantification branch
- metabolomics/isotope-tracing - Separate stable-isotope flux branch, not a stage of this untargeted pipeline
- multi-omics-integration/mofa-integration - Integrating the feature table with other omics layers
