# Outlier Splicing Detection - Usage Guide

## Overview
Detect aberrant splicing in single rare-disease patients vs a control panel - the statistical question is fundamentally different from differential splicing between groups. Tools include FRASER 2.0 (Beta-binomial autoencoder on Intron Jaccard Index, Bioconductor), OUTRIDER (autoencoder on gene expression for outlier expression), LeafcutterMD (Dirichlet-multinomial outlier mode of LeafCutter), and DROP (Snakemake pipeline integrating FRASER2 + OUTRIDER + monoallelic expression). Standard tool in EU rare-disease (Solve-RD) and NIH UDN programs.

## Prerequisites
```bash
# R Bioconductor
BiocManager::install(c('FRASER', 'OUTRIDER', 'leafcutter'))

# DROP pipeline (bioconda; not on PyPI)
mamba create -n drop_env -c conda-forge -c bioconda drop --override-channels

# CLI dependencies
conda install -c bioconda regtools snakemake star samtools
```

## Quick Start
Tell your AI agent what you want to do:
- "Run FRASER2 to detect aberrant splicing in my rare-disease patient vs cohort"
- "Set up DROP pipeline for integrated outlier detection (splicing + expression + MAE)"
- "Use OUTRIDER to find genes with aberrant expression in patient samples"
- "Cross-reference SpliceAI variant hits with FRASER2 RNA outliers in the same patient"
- "Run LeafcutterMD for annotation-free outlier intron usage"

## Example Prompts

### Single-Patient Workflow
> "I have RNA-seq from a rare-disease patient and 50 control samples; run FRASER 2.0 with Intron Jaccard Index and report aberrant junctions with padj<0.05 and |delta|>=0.1."

### Pipeline Setup
> "Configure DROP for our diagnostic cohort with patient + 80 in-house controls + GTEx-derived auxiliary controls."

### Variant Confirmation
> "I have a SpliceAI hit at chr5:1234567C>T in patient X; check whether FRASER2 detects an aberrant junction within 1kb in the same sample."

### Disease-Specific
> "For ALS post-mortem brain RNA-seq, use FRASER2 to detect TDP-43-loss cryptic exons (UNC13A, STMN2, ATG4B)."

### Hyperparameter Tuning
> "Run estimateBestQ on my cohort to determine the optimal autoencoder dimension q for FRASER2."

## What the Agent Will Do
1. Audit cohort size and tissue match (n>=50 ideal, n>=20 minimum)
2. Run FRASER2 with Intron Jaccard Index (default) on cohort + patient
3. Run OUTRIDER for gene-level outlier expression
4. Optionally configure DROP for integrated pipeline
5. Filter outliers by padj, delta, gnomAD splice constraint
6. Cross-reference with DNA variant calls (SpliceAI predictions, ClinVar)
7. Generate clinical report with PS3 functional evidence weight

## Tips
- FRASER 2.0 (Bioconductor >=1.99.0) defaults to Intron Jaccard Index; delta cutoff is 0.1 (vs 0.3 in v1.x)
- Cohort size matters: n>=50 ideal; n=20-50 acceptable; n<20 needs auxiliary controls (GTEx tissue-matched)
- Tissue match is required - clinical gene expression varies dramatically across blood/fibroblast/muscle/iPSC
- Combine in-house and GTEx controls cautiously - batch effects; use ComBat or include batch covariate
- Hyperparameter `q` should be tuned via estimateBestQ; default q=10 for ~50-100 sample cohorts
- A SpliceAI hit + concordant FRASER2 outlier = strong PS3 evidence in ACMG framework
- LeafcutterMD complements FRASER2 when novel-junction sensitivity matters or Beta-binomial fits poorly
- Treat all outliers as candidates, not pathogenic - filter against gnomAD splice constraint and disease panels

## Related Skills

- splice-variant-prediction - SpliceAI / Pangolin for in-silico prediction first
- differential-splicing - When testing multiple patients vs controls
- splicing-qc - Library / depth / tissue prerequisites
- variant-calling/clinical-interpretation - ACMG/AMP framework integration
- workflows/clinical-trial-pipeline - Trial-grade RNA-seq diagnostics
