# Community Ecology Usage Guide

## Overview

Analyzes how environmental gradients structure species communities using constrained ordination (CCA for unimodal responses, RDA for linear responses, db-RDA for arbitrary distances), variance partitioning to quantify unique and shared contributions of predictor sets, and indicator species analysis to identify taxa characteristic of habitats or conditions. All methods are implemented with the vegan and indicspecies R packages.

## Prerequisites

- R with vegan (`install.packages('vegan')`)
- indicspecies (`install.packages('indicspecies')`)
- ggplot2 for custom ordination plots (`install.packages('ggplot2')`)

## Quick Start

Tell your AI agent what you want to do:

- "Test which environmental variables drive species composition at my sites"
- "Run a CCA on my species abundance and environmental data"
- "Partition variance between spatial and environmental predictors"
- "Find indicator species for each habitat type in my dataset"
- "Select the best environmental predictors for community composition using forward selection"

## Example Prompts

### Constrained Ordination

> "I have a species-by-site abundance matrix and environmental measurements (temperature, pH, elevation, precipitation). Run CCA with permutation tests and produce a triplot showing sites, species, and environmental vectors."

> "My DCA axis 1 length is 1.8 SD, so species responses are linear. Run RDA with Hellinger transformation, forward-select the significant predictors, and check VIF for collinearity."

### Distance-Based Methods

> "I want to use Bray-Curtis distances for my community analysis. Run a db-RDA to test which environmental variables explain composition differences."

### Variance Partitioning

> "Partition the variation in species composition between environmental variables and spatial eigenvectors (PCNM). Show me the Venn diagram and test the unique fractions."

> "I have three sets of predictors: environment, space, and land use history. Run a three-way variance partition."

### Indicator Species

> "Identify indicator species for each of my four habitat types using the IndVal index with group equalization."

> "Which species are significantly associated with disturbed vs undisturbed sites? Use point-biserial correlation."

## What the Agent Will Do

1. Assess gradient length using DCA to choose between CCA (unimodal, >3 SD) and RDA (linear, <3 SD)
2. Apply appropriate data transformations (Hellinger for RDA, raw for CCA)
3. Fit the constrained ordination model with environmental predictors
4. Run permutation tests for global significance, per-axis, and per-term (marginal) contributions
5. Perform forward selection (ordiR2step) to identify the most parsimonious predictor set
6. Check VIF for multicollinearity among selected predictors (flag VIF > 10)
7. Partition variance among predictor groups and test unique fractions
8. Run indicator species analysis with permutation tests and report significant associations
9. Produce triplots, Venn diagrams, and summary tables

## Tips

- Check gradient length with DCA first: axis 1 > 3 SD suggests CCA; < 3 SD suggests RDA
- Always use Hellinger transformation for RDA on abundance data to handle the double-zero problem
- VIF > 10 indicates collinearity; remove or combine correlated predictors before final model
- Forward selection with ordiR2step prevents overfitting by stopping when adjusted R-squared stops improving
- db-RDA with Bray-Curtis is often the most flexible choice for ecological data
- Variance partitioning fractions can be negative, indicating that shared variation exceeds unique contributions
- For indicator species, use `func = 'IndVal.g'` with unbalanced group sizes (equalizes group weights)
- Set `duleg = TRUE` in multipatt to test only individual group associations, not all possible combinations

## Related Skills

- biodiversity-metrics - Alpha and beta diversity metrics
- edna-metabarcoding - Generate community data from eDNA
- landscape-genomics - Genotype-environment associations
- microbiome/diversity-analysis - Unconstrained ordination alternative
- data-visualization/ggplot2-fundamentals - Custom ordination plots
