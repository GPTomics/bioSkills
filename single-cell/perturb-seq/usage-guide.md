# Perturb-seq Analysis - Usage Guide

## Overview

Perturb-seq and CROP-seq read out a pooled CRISPR screen with single-cell transcriptomes, linking a recoverable guide identity (genotype) to a transcriptome-wide phenotype. This skill treats the hard decisions: guide assignment as a mixture problem, escaper removal with Mixscape, calibrated testing with SCEPTRE, effect size with E-distance, separating compositional shifts from within-state expression, and a clear-eyed view that perturbation-prediction foundation models do not yet beat simple baselines.

## Prerequisites

```bash
pip install pertpy scanpy anndata
pip install pydeseq2 decoupler          # pseudobulk DE
```

```r
install.packages('sceptre')             # conditional-resampling test
install.packages('Seurat')              # Mixscape (Seurat v5)
```

## Quick Start

Tell your AI agent what you want to do:
- "Assign guides with a mixture model, not a flat threshold"
- "Use Mixscape to remove cells that received a guide but were not perturbed"
- "Test each perturbation with SCEPTRE so false positives are calibrated"
- "Rank my perturbations by E-distance and run the E-test"
- "Tell me whether this perturbation moves cells across states or changes a state"
- "Is this foundation-model prediction actually better than a mean baseline?"

## Example Prompts

### Guide Assignment
> "Fit a Poisson-Gaussian mixture to my guide counts and assign by posterior, then show the NT contamination floor"
> "Gate doublets before treating multi-guide cells as combinatorial"

### Effective Perturbation
> "Run Mixscape to classify KO vs non-perturbed cells and report the perturbed fraction per target"
> "This target is all non-perturbed; is the gene non-functional or did the guide fail to edit?"

### Testing and Effect Size
> "Run SCEPTRE with a calibration check before the discovery analysis"
> "Compute pairwise E-distances in PCA space and run the permutation E-test against NT"
> "Do pseudobulk DESeq2 per replicate, summing raw counts, for the within-state program change"

### Composition vs Expression
> "Run Milo to test whether the perturbation shifts cell-state proportions"
> "Separate the compositional shift from the within-state expression change"

### Foundation Models
> "Benchmark this perturbation predictor on held-out whole perturbations against an additive baseline, scored on DE genes"

## What the Agent Will Do

1. Assign guides with a per-guide background/foreground mixture and report the perturbed fraction.
2. Gate doublets and decide the MOI regime (low for single-gene attribution, high for combinatorial).
3. Run Mixscape to remove non-perturbed escaper cells before any DE.
4. Choose a calibrated test (SCEPTRE conditional resampling) over naive Wilcoxon/NB.
5. Quantify effect size with E-distance in a pinned PCA embedding and the permutation E-test.
6. Run pseudobulk DE per biological replicate (summing raw counts) for the within-state program.
7. Run a differential-abundance test (Milo/scCODA) and report composition separately from expression.

## Decision Guidance

### Guide assignment
- Mixture posterior (default): ambient contamination is biased to abundant guides, so a flat threshold mis-assigns; require a dominant-guide UMI fraction.
- MOI: low-MOI gives clean single-gene attribution but discards most cells; high-MOI is for combinatorial designs and needs deconvolution (scMAGeCK-LR, GSFA).

### Testing
- SCEPTRE when calibration matters: the depth confounder breaks naive parametric tests.
- Pseudobulk DE for the within-state program, but only with >=2-3 biological replicates; one transfection per guide has no valid replicate inference.
- E-distance + E-test for effect-size magnitude and perturbation similarity, pinned to a fixed embedding.

### Composition vs expression
- Milo/scCODA/Augur answer "does it move cells across states?"
- Pseudobulk DE/Mixscape/GSFA answer "does it change a state's program?"
- A perturbation that only redistributes cells produces a fake pseudobulk DE signature; always report both.

## Tips

- **Assignment is a mixture, not a threshold** - ambient guide contamination scales with guide abundance; assign by posterior.
- **Assignment is not perturbation** - escapers and incomplete KO dilute effect sizes; run Mixscape first.
- **All-NP is ambiguous** - it confounds "no phenotype" with "no editing"; never call a gene non-functional from Mixscape alone.
- **Naive DE is miscalibrated** - depth confounding and pseudoreplication inflate type-I error; use SCEPTRE and pseudobulk-per-replicate.
- **Sum raw counts for pseudobulk** - not means or normalized values; filter pseudobulk samples below ~10 cells.
- **E-distance is embedding-relative** - pin the pertpy version, obsm key, and metric; do not compare across studies.
- **Separate moves-cells from changes-cells** - composition vs within-state expression are different questions needing different tools.
- **Foundation models do not yet beat baselines** - on held-out perturbations they do not exceed additive/mean predictors; always benchmark on DE genes with whole-perturbation holdout.
- **Non-targeting controls define the null** - weak or contaminated NT inflate false positives across every test.

## Related Skills

single-cell/preprocessing - scRNA-seq QC and normalization upstream of the screen
single-cell/doublet-detection - gating doublets before multi-guide analysis
single-cell/markers-annotation - interpreting per-perturbation DE genes
single-cell/batch-integration - multi-sample/replicate integration
crispr-screens/mageck-analysis - bulk CRISPR screen analysis (MAGeCK RRA/MLE)
crispr-screens/perturb-seq-analysis - related single-cell CRISPR screen workflow
differential-expression/deseq2-basics - pseudobulk DESeq2 testing on summed counts
pathway-analysis/go-enrichment - pathway interpretation of perturbation signatures
