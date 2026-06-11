# mixOmics Analysis - Usage Guide

## Overview

mixOmics provides multivariate projection methods for multi-omics integration: sPLS (sparse pairwise correlation), DIABLO/block.splsda (supervised multi-block discriminant signatures), rCCA (regularized canonical correlation), and MINT (multi-study integration). The unifying caution is that these methods maximize an association criterion - covariance for PLS/DIABLO, correlation for CCA - and not biological truth. With far more features than samples they can fit almost anything, so the deliverable is a signature whose performance was estimated on data that played no role in selecting it.

The skill owns the mixOmics-specific mechanics (the design matrix, keepX tuning, BER, the circos/network plots) and the supervised-overfitting discipline. The unsupervised factor-model alternative is mofa-integration; the generic cross-validation and data-leakage theory lives in machine-learning/model-validation; per-block scaling and batch happen upstream in data-harmonization.

## Prerequisites

```r
BiocManager::install('mixOmics')
```

Conceptual prerequisites and notes:
- mixOmics input is samples-by-features (the opposite of MOFA2). DIABLO, sPLS, and rCCA relate samples row-by-row and require the same samples in the same order across every block; verify identical rownames before fitting.
- DIABLO is `block.splsda` - there is no function named `diablo`.
- `scale=TRUE` is the default and centers and scales every feature to unit variance, which changes the PLS objective; feed appropriately transformed matrices but do not pre-standardize them yourself.
- Combining one omic across cohorts is horizontal integration and needs MINT (study as a fixed effect), not DIABLO.

## Quick Start

Tell your AI agent what you want to do:
- "Find a cross-omic signature that discriminates my subtypes with DIABLO"
- "Choose the DIABLO design matrix for my goal and explain the trade-off"
- "Tune keepX honestly so my reported accuracy is not leaked"
- "Find correlated gene-metabolite pairs with sPLS"
- "Integrate RNA-seq from three cohorts with MINT"

## Example Prompts

### Supervised signature with honest validation
> "Use DIABLO to find a cross-omic signature separating my responder and non-responder groups from RNA-seq and proteomics. Tune the sparsity with repeated cross-validation using balanced error rate, and give me a performance estimate from data that was not used to select the features."

### Design-matrix decision
> "I want a biologically coherent, inter-correlated multi-omics network rather than the best classifier. Set the DIABLO design matrix accordingly, and tell me what that choice costs me in classification accuracy."

### Pairwise correlation
> "Find a sparse set of genes and metabolites that covary across my samples with sPLS, treating the two omics symmetrically since I have no directional hypothesis."

### Multi-study integration
> "I have the same transcriptomic assay from three cohorts. Build a signature that replicates across them with MINT, modeling study as a fixed effect."

## What the Agent Will Do

1. Confirm the correspondence: matched samples across blocks for DIABLO/sPLS/rCCA, or a study factor for MINT.
2. Select the method from the question (sPLS, rCCA, DIABLO, MINT) and set `mode` or the design matrix deliberately.
3. Tune the component count on a non-sparse model, then tune keepX inside cross-validation folds using balanced error rate.
4. Fit the final model and extract the selected features per block and component as candidates.
5. Estimate performance from an external test set or nested cross-validation, not from the tuning data.
6. Visualize cross-block structure (circos, network, sample plots) and route the candidates to pathway analysis with replication caveats.

## Method Selection Guide

| Method | Supervised | Blocks | Use case |
|--------|------------|--------|----------|
| sPCA | No | 1 | Sparse dimension reduction of one omic |
| sPLS | No | 2 | Sparse correlated feature pairs |
| rCCA | No | 2 | Regularized correlation landscape |
| sPLS-DA | Yes | 1 | Classify with feature selection in one omic |
| DIABLO | Yes | 2+ | Cross-omic discriminant signature (matched samples) |
| MINT | Yes | 1 | One omic across studies (horizontal) |

## Tips

- Cross-validation must wrap feature selection; tuning keepX on the full data and scoring it by CV on the same data is leakage. Report an external test set or nested CV.
- Choose the DIABLO design matrix from the goal: near 1 for a coherent inter-correlated network, below 0.5 for a predictive signature. The tutorials' 0.1 is convention, not a recommendation.
- Use balanced error rate for imbalanced classes; overall error is dominated by the majority class.
- Set `nrepeat` to at least 10 (50 for a headline number); a single fold split is noise.
- Use `mode='canonical'` for two omics on equal footing; the default `mode='regression'` treats the second block as a response.
- Never run un-regularized CCA on more features than samples; it reports correlation 1.0 by overfitting. Use rCCA with ridge or shrinkage.
- Treat selectVar output as cohort-specific candidates, not validated biomarkers; replication is the finding.

## Related Skills

- integration-design - The method-selection and paired-vs-horizontal decision
- mofa-integration - Unsupervised factor alternative where no outcome drives the fit
- data-harmonization - Per-block scaling and batch before matrices enter mixOmics
- machine-learning/model-validation - Nested cross-validation and data-leakage theory
- machine-learning/biomarker-discovery - Biomarker-panel selection and validation
- pathway-analysis/go-enrichment - Enrichment of the selected features
- differential-expression/de-results - Single-omic differential expression
- workflows/multi-omics-pipeline - End-to-end multi-omics integration pipeline
