# proteomics

## Overview

Mass spectrometry-based bottom-up proteomics, from raw MS files and search-engine output through identification, quantification, and differential abundance, with DIA, PTM, and quality-control coverage. Each skill is decision-grade: it teaches the mechanism and the trap behind the choice, not just the tool call.

**Tool type:** mixed | **Primary tools:** pyOpenMS, DIA-NN, MaxQuant, limma, DEqMS, MSstats

## Skills

| Skill | Description |
|-------|-------------|
| data-import | Load MS formats (mzML, MaxQuant proteinGroups.txt, DIA-NN report.parquet); strip decoys/contaminants; inherit the acquisition missingness contract |
| peptide-identification | Database search, target-decoy FDR, q-value vs PEP, rescoring (Percolator/mokapot) |
| protein-inference | Protein grouping, parsimony, picked (group) FDR, razor-vs-unique quant consequences |
| quantification | LFQ (MaxLFQ), isobaric TMT/iTRAQ (ratio compression, SPS-MS3, IRS bridge), SILAC; peptide-to-protein summarization |
| proteomics-qc | Three-level QC funnel (instrument, identification, matrix); inspect before normalizing |
| differential-abundance | Moderated testing (limma/DEqMS/proDA/msqrob2/MSstats); left-censored MNAR missingness; FDR control |
| ptm-analysis | PTM site localization (FLR), protein-level adjustment (MSstatsPTM), kinase-activity inference |
| dia-analysis | DIA-NN peak-group scoring, library-free vs predicted-library, q-value contexts (run vs global) |
| spectral-libraries | DDA/chromatogram/predicted libraries (Koina/Prosit/AlphaPeptDeep), RT/CCS calibration |

## Method Selection

| Decision | Options | Guidance |
|----------|---------|----------|
| Acquisition | DDA / DIA / PRM-SRM | DDA for discovery and TMT; DIA for large reproducible cohorts; targeted (PRM/SRM) for validation with standards |
| Labeling | LFQ / TMT(pro) / SILAC | LFQ for unlimited samples; TMT up to 18-plex in one run (needs an IRS bridge channel to compare across plexes); SILAC for labelable cell culture, lowest ratio variance |
| Quant summarization | MaxLFQ / Tukey median polish / msqrob | Report it: the summarizer changes the answer more than the test does |
| Differential test | limma-trend / DEqMS / proDA / msqrob2 / MSstats | DEqMS when quant depth varies; proDA or msqrob2 to model MNAR without imputing; MSstats for feature-level or labeled designs |
| Missing values | model vs impute | Model the dropout (proDA/msqrob2/MSstats-AFT); downshift imputation manufactures false positives near the detection limit |

## Example Prompts

- "Load MaxQuant proteinGroups.txt, strip contaminants and decoys, and build an LFQ intensity matrix"
- "Find differentially abundant proteins between treatment and control while handling missing values correctly"
- "Quantify a multi-batch TMT experiment and bridge the plexes"
- "Run a library-free DIA-NN analysis and filter on the right q-value context"
- "Localize phosphosites and adjust for protein-level abundance changes"
- "Estimate protein-level FDR with a picked target-decoy approach"

## Requirements

```bash
# Python
pip install pyopenms pandas numpy scipy pyarrow

# R
install.packages(c("limma", "iq", "ashr"))
BiocManager::install(c("DEqMS", "MSstats", "MSstatsPTM", "proDA", "QFeatures"))

# DIA-NN and MaxQuant are standalone tools (not pip/CRAN); install separately.
```

## Related Skills

- **differential-expression** - Empirical-Bayes moderated testing for RNA-seq (shared limma machinery)
- **pathway-analysis** - Functional enrichment of differential protein lists
- **data-visualization** - Volcano plots, heatmaps, and dimensionality-reduction plots for proteomics results
- **machine-learning** - Protein biomarker panels and classifier validation
- **workflows** - End-to-end proteomics-pipeline orchestration
