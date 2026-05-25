# Community Ecology Usage Guide

## Overview

Analyzes how environmental gradients structure species communities using constrained ordination (CCA for unimodal responses, RDA on Hellinger-transformed data for linear responses, db-RDA for arbitrary distances), variance partitioning with adjusted R-squared (Peres-Neto 2006), indicator species analysis with the group-size-equalized IndVal.g form (De Caceres & Legendre 2009), PERMANOVA paired MANDATORILY with PERMDISP dispersion testing since the two confound each other (Anderson & Walsh 2013), Joint Species Distribution Models (HMSC, sjSDM, gjam) with explicit acknowledgement that residual covariance is NOT pure biotic interaction, phylogenetic community ecology with explicit null-model justification (Webb 2002), trait-environment via RLQ + fourth-corner with the corrected `modeltype=6` permutation (Dray 2014), bipartite network metrics with Strona 2014 curveball null randomization, and Mantel test replacements (db-RDA conditioned on spatial PCNM eigenvectors, GDM) per Legendre & Fortin 2010.

## Prerequisites

- R with vegan (`install.packages('vegan')`)
- indicspecies (`install.packages('indicspecies')`)
- Hmsc for theory-rich Bayesian JSDM (`install.packages('Hmsc')`)
- sjSDM for scalable JSDM (`install.packages('sjSDM')`)
- picante for SES_MPD/MNTD (`install.packages('picante')`)
- ade4 for RLQ + fourth-corner (`install.packages('ade4')`)
- adespatial for PCNM spatial eigenvectors (`install.packages('adespatial')`)
- ggplot2 for custom plots (`install.packages('ggplot2')`)

## Quick Start

Tell your AI agent what you want to do:

- "Run PERMANOVA AND betadisper on my community data and report both"
- "Test community-environment relationships with RDA or CCA based on gradient length"
- "Replace my Mantel test with dbRDA conditioned on spatial PCNM eigenvectors"
- "Fit an HMSC joint species distribution model with phylogenetic and trait effects"
- "Find indicator species using IndVal.g (group-size-corrected) not basic IndVal"
- "Test trait-environment relationships with RLQ + fourth-corner using modeltype=6"

## Example Prompts

### Constrained Ordination

> "I have a species-by-site abundance matrix and environmental variables. Check DCA axis-1 gradient length first; if > 3 SD use CCA, if < 3 SD use Hellinger-RDA. Apply forward selection, check VIF, and produce a triplot."

### PERMANOVA with Dispersion Companion

> "Test community composition differences across my four habitat types with PERMANOVA (adonis2). MANDATORY: also run betadisper to test dispersion heterogeneity. If betadisper is significant, the PERMANOVA conclusion must be qualified (Anderson & Walsh 2013)."

### Mantel Replacement

> "I have community distances and environmental distances. Don't run Mantel; replace with dbRDA using spatial PCNM eigenvectors as Condition() (Legendre & Fortin 2010). Test marginal environmental effect after partialling out space."

### Joint Species Distribution Models

> "Fit HMSC for my 200-species community with environmental predictors, traits, phylogeny, and spatial random effect. Report variance partitioning into fixed effects, random effects, traits, and phylogeny. Explicitly disclaim that residual covariance is not pure biotic interaction (Pollock 2014; Poggiato 2021)."

> "I have a 2000-ASV community matrix from metabarcoding. HMSC is too slow at this S; use sjSDM with elastic-net regularization (Pichler & Hartig 2021)."

### Indicator Species

> "Find indicator species for each of my four habitat types using `multipatt` with `func='IndVal.g'` (NOT 'IndVal' which is biased by group size). Use 999 permutations and report results sorted by p-value."

### Phylogenetic Community Structure

> "Compute SES_MPD and SES_MNTD across my communities with the `taxa.labels` null model. Then test whether traits track phylogeny with Blomberg's K and Pagel's lambda before interpreting clustering as environmental filtering; cite Mayfield & Levine 2010."

### Trait-Environment

> "Test how species traits respond to environmental gradients using RLQ ordination followed by the fourth-corner permutation test with `modeltype=6` (Dray 2014). Default `modeltype=2/4` gives inflated Type I error."

## What the Agent Will Do

1. Inspect species-environment data and assess gradient length with DCA before choosing CCA or RDA; apply Hellinger transformation before RDA on abundance data
2. Run PERMANOVA via `adonis2` (NOT the deprecated `adonis`) with `by='margin'` for unbalanced designs
3. ALWAYS run `betadisper` + `permutest` alongside PERMANOVA; explicitly report when the location conclusion is confounded by dispersion heterogeneity
4. Replace Mantel/partial Mantel for landscape data with db-RDA conditioned on spatial PCNM eigenvectors or GDM
5. Choose JSDM by problem scale: HMSC for S < 500 with rich priors; sjSDM for S > 500; gjam for cross-data-type
6. Explicitly disclaim "residual covariance equals biotic interaction" in any JSDM report
7. Use `func='IndVal.g'` (group-size-corrected) for indicator species; never the basic 'IndVal'
8. Compute SES_MPD/SES_MNTD with an EXPLICITLY documented null model (taxa.labels, independentswap, richness); test trait conservatism before inferring filtering vs competition (Mayfield & Levine 2010)
9. Run RLQ + fourth-corner with `modeltype=6` for trait-environment hypotheses
10. Apply Strona 2014 curveball null for bipartite network nestedness/modularity testing
11. Report VIF for collinearity, adjusted R-squared (not raw R-squared) for variance explained, and bootstrap CIs (>= 999 permutations)

## Tips

- DCA gradient length on axis 1: < 3 SD favors linear methods (RDA); > 3 SD favors unimodal (CCA); 2-3 SD is a gray zone where either is defensible; report which was chosen and why
- Hellinger transformation before RDA is mandatory for community data (Legendre & Gallagher 2001); without it, RDA is dominated by sample-size effects
- `adonis2()` is the modern PERMANOVA; `adonis()` was deprecated in vegan 2.6, so tutorials citing `adonis()` are outdated
- `by = 'margin'` in `adonis2` gives the marginal-SS PERMANOVA appropriate for unbalanced designs; `by = 'terms'` gives sequential SS (the default) which is order-dependent
- PERMANOVA + PERMDISP is non-negotiable: report BOTH. If both p < 0.05, the location-difference conclusion is not supported by these tests alone
- ANOSIM is worse than PERMANOVA for the same use case; biased by unbalanced N and dispersion (Anderson & Walsh 2013); do not use
- Pairwise PERMANOVA with Bonferroni correction is statistically incorrect for permutation tests; use `pairwiseAdonis::pairwise.adonis2` with FDR correction
- Mantel and partial Mantel have low power and biased Type I error under spatial autocorrelation; replace with db-RDA + spatial covariates or GDM
- HMSC MCMC needs thin >= 100, transient >= 1000, samples >= 1000 for publication-quality posteriors; check convergence with `convertToCodaObject(m)` and Gelman-Rubin diagnostics
- sjSDM trains via PyTorch; if no GPU, set `device='cpu'` (slower but works)
- JSDM residual covariance is NOT biotic interaction (Pollock 2014; Zurell 2018; Poggiato 2021); always disclaim this interpretation
- `IndVal.g` is the group-equalized indicator value; basic `IndVal` is biased toward larger groups; do not use 'IndVal' for unbalanced designs
- For SES_MPD/MNTD null-model choice: `taxa.labels` randomizes which species are in samples (holds richness constant); `independentswap` preserves both row and column sums; `richness` randomizes within samples; the choice answers different ecological questions
- Mayfield & Levine 2010 showed competition can produce phylogenetic clustering when traits track phylogeny; do not infer process from SES sign alone
- Curveball algorithm (Strona 2014) is the modern null for bipartite binary network randomization; exponentially faster than `r2dtable`/quasiswap; vegan/bipartite adopted as default circa 2020
- RLQ + fourth-corner with `modeltype=6` (Dray 2014); `modeltype=2/4` alone gives inflated Type I error
- Variance partitioning fractions use adjusted R^2 (Peres-Neto 2006); raw R^2 is biased upward

## Related Skills

- ecological-genomics/biodiversity-metrics - Alpha/beta diversity and Hill numbers prior to ordination
- ecological-genomics/edna-metabarcoding - Generate community data from environmental samples
- ecological-genomics/landscape-genomics - Genotype-environment associations (genetic analog of GEA)
- microbiome/diversity-analysis - Unconstrained ordination alternative for 16S microbiome data
- data-visualization/ggplot2-fundamentals - Customize triplots, ordination plots, indicator-species visualizations
