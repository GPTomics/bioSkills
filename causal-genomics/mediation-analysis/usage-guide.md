# Mediation Analysis - Usage Guide

## Overview

Decompose genetic effects into direct and indirect paths through mediating variables. Tests whether molecular phenotypes such as gene expression or methylation mediate the relationship between genotype and disease.

## Prerequisites

```r
install.packages('mediation')

# For high-dimensional mediation
install.packages('HIMA')

# Visualization
install.packages('ggplot2')
```

## Quick Start

Tell your AI agent what you want to do:
- "Test whether gene expression mediates the effect of this SNP on disease"
- "Run mediation analysis with methylation as the mediator between genotype and phenotype"
- "Decompose the total effect of rs12345 into direct and indirect paths through GENE_X"

## Example Prompts

### eQTL Mediation

> "I have genotype, expression, and disease data -- test if expression of BRCA1 mediates the SNP effect on breast cancer"

> "Run mediation analysis for all eQTL genes at this locus"

### Multi-Omics

> "Test the mediation chain: SNP -> methylation -> expression -> disease"

> "Use HIMA to find which CpG sites mediate the genotype-disease association"

### Sensitivity

> "How robust are my mediation results to unmeasured confounding?"

> "Run sensitivity analysis on the mediation result and plot the sensitivity curve"

## What the Agent Will Do

1. Fit mediator model (treatment -> mediator)
2. Fit outcome model (treatment + mediator -> outcome)
3. Run mediation analysis with bootstrap confidence intervals
4. Report ACME, ADE, total effect, and proportion mediated
5. Run sensitivity analysis for unmeasured confounding
6. Adjust p-values for multiple testing if testing many mediators

## Tips

- **Sequential ignorability** - The key assumption is untestable; always run sensitivity analysis
- **Bootstrap CIs** - Use at least 1000 simulations; 5000 for publication
- **Proportion mediated** - Values above 0.2 suggest meaningful mediation; above 0.8 suggests the mediator explains most of the effect
- **Binary outcomes** - Use glm with family = binomial for the outcome model
- **Multiple mediators** - Use HIMA for high-dimensional mediation with penalized regression
- **Covariates** - Include population structure PCs, age, and sex in both models

## Related Skills

- causal-genomics/mendelian-randomization - Causal inference with genetic instruments
- causal-genomics/colocalization-analysis - Confirm shared causal variants
- population-genetics/association-testing - GWAS data for mediation inputs
- multi-omics-integration/mofa-integration - Multi-omics data integration
