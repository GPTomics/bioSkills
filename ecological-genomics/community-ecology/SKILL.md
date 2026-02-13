---
name: bio-ecological-genomics-community-ecology
description: Analyzes community composition using constrained ordination (CCA, RDA, db-RDA), variance partitioning (varpart), indicator species analysis (indicspecies multipatt), and distance-based environmental gradient methods with vegan. Links species composition to environmental explanatory variables. Use when testing how environmental gradients structure species communities, identifying habitat indicator taxa, or partitioning explained variation among predictors. Not for basic unconstrained ordination and PERMANOVA (see microbiome/diversity-analysis).
tool_type: r
primary_tool: vegan
---

## Version Compatibility

Reference examples tested with: ggplot2 3.5+, vegan 2.6+

Before using code patterns, verify installed versions match. If versions differ:
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Community Ecology

**"Test how environmental variables structure my species communities"** → Link species composition to environmental gradients using constrained ordination (CCA, RDA, db-RDA), partition explained variation among predictor groups with varpart, and identify habitat indicator taxa with multipatt.
- R: `vegan::cca()`, `vegan::rda()`, `vegan::dbrda()` for constrained ordination
- R: `indicspecies::multipatt()` for indicator species analysis

Analyzes species-environment relationships using constrained ordination, variance partitioning, and indicator species analysis.

## Canonical Correspondence Analysis (CCA)

CCA models unimodal species responses along environmental gradients. Appropriate when species have optima along gradients (bell-shaped response curves):

```r
library(vegan)

# species_matrix: sites x species abundance matrix
# env_data: sites x environmental variables data frame
cca_result <- cca(species_matrix ~ temperature + precipitation + soil_pH + elevation,
                  data = env_data)

# Global permutation test (999 permutations by default)
anova(cca_result, permutations = 999)

# Per-axis significance
anova(cca_result, by = 'axis', permutations = 999)

# Per-term (marginal) significance
anova(cca_result, by = 'margin', permutations = 999)

# Proportion of variance explained
# Total inertia = constrained + unconstrained
summary(cca_result)$tot.chi
summary(cca_result)$constr.chi
prop_explained <- summary(cca_result)$constr.chi / summary(cca_result)$tot.chi

# Triplot: sites, species, and environmental vectors
plot(cca_result, display = c('sites', 'species', 'bp'), scaling = 2)
```

## Redundancy Analysis (RDA)

RDA models linear species responses. Preferred when species respond linearly to gradients (short gradients, <= 3 SD):

```r
# Hellinger transformation recommended for abundance data with RDA
# Reduces the weight of rare species and double-zero problem
species_hell <- decostand(species_matrix, method = 'hellinger')

rda_result <- rda(species_hell ~ temperature + precipitation + soil_pH + elevation,
                  data = env_data)

# Adjusted R-squared (corrected for number of predictors)
RsquareAdj(rda_result)

# Permutation test
anova(rda_result, permutations = 999)
anova(rda_result, by = 'margin', permutations = 999)
```

## Forward Selection with ordiR2step

Selects the most parsimonious set of environmental predictors:

```r
# Null model (intercept only)
rda_null <- rda(species_hell ~ 1, data = env_data)

# Full model (all predictors)
rda_full <- rda(species_hell ~ ., data = env_data)

# Forward selection using adjusted R-squared criterion
# Stops when adding variables no longer improves adjusted R^2
# or when p > 0.05 (permutation test)
rda_sel <- ordiR2step(rda_null, scope = formula(rda_full),
                      direction = 'forward', permutations = 999)

# Check VIF for collinearity in selected model
# VIF > 10: indicates problematic multicollinearity; remove or combine variables
vif.cca(rda_sel)
```

## Distance-Based RDA (db-RDA)

Extends RDA to any distance metric using principal coordinates:

```r
# Bray-Curtis distance for abundance data
bray_dist <- vegdist(species_matrix, method = 'bray')

# capscale (legacy) or dbrda (modern)
dbrda_result <- dbrda(bray_dist ~ temperature + precipitation + soil_pH,
                      data = env_data)

# Negative eigenvalues handled automatically (add correction if needed)
dbrda_result <- dbrda(bray_dist ~ temperature + precipitation + soil_pH,
                      data = env_data, add = 'lingoes')

anova(dbrda_result, permutations = 999)
anova(dbrda_result, by = 'margin', permutations = 999)
```

## Variance Partitioning

**Goal:** Quantify the unique and shared contributions of environmental vs spatial predictors to community composition.

**Approach:** Run varpart with the species matrix and predictor groups to decompose variation into unique and shared fractions, then test each unique fraction's significance with partial RDA conditioned on the other group.

Quantifies the unique and shared contributions of 2-4 sets of explanatory variables:

```r
# Partition variance between spatial and environmental predictors
# spatial_data: PCNM/MEM spatial eigenvectors or coordinates
# env_data: environmental variables
vp <- varpart(species_hell, env_data, spatial_data)

# Fractions: [a] = unique environment, [b] = shared, [c] = unique spatial, [d] = residual
vp

# Venn diagram of explained variation
plot(vp, digits = 2, bg = c('skyblue', 'tomato'))

# Test significance of unique fractions (testable fractions only)
# [a] unique environment
anova(rda(species_hell ~ . + Condition(as.matrix(spatial_data)), data = env_data),
      permutations = 999)

# [c] unique spatial
anova(rda(species_hell ~ . + Condition(as.matrix(env_data)), data = as.data.frame(spatial_data)),
      permutations = 999)
```

### Three-Way Partition

```r
# Environment, space, and time
vp3 <- varpart(species_hell, env_data, spatial_data, temporal_data)
plot(vp3, digits = 2, bg = c('skyblue', 'tomato', 'gold'))
```

## Indicator Species Analysis

Identifies species significantly associated with site groups using the IndVal index:

```r
library(indicspecies)

# site_groups: factor or vector of group assignments per site
# multipatt tests association of each species with group combinations
# func='IndVal.g': group-equalized indicator value (handles unbalanced groups)
# duleg=TRUE: only tests individual groups, not combinations
mp <- multipatt(species_matrix, site_groups, func = 'IndVal.g',
                duleg = TRUE, control = how(nperm = 999))

summary(mp)

# Extract significant indicators
# p.value < 0.05: standard significance threshold
sig_indicators <- mp$sign[mp$sign$p.value < 0.05, ]
sig_indicators <- sig_indicators[order(sig_indicators$p.value), ]
sig_indicators
```

### Point-Biserial Correlation

```r
# Alternative: r.g function for point-biserial correlation coefficients
# More interpretable effect sizes than IndVal
pb <- multipatt(species_matrix, site_groups, func = 'r.g',
                duleg = TRUE, control = how(nperm = 999))
summary(pb)
```

## CCA vs RDA Decision

| Criterion | Use CCA | Use RDA |
|-----------|---------|---------|
| Gradient length | > 3 SD (DCA axis 1) | <= 3 SD (DCA axis 1) |
| Species response | Unimodal (bell curves) | Linear |
| Data transformation | Raw abundances OK | Hellinger recommended |
| Chi-square distance | Appropriate | Not appropriate |

Check gradient length with DCA:

```r
dca <- decorana(species_matrix)
# Axis 1 length > 3 SD: unimodal methods (CCA)
# Axis 1 length <= 3 SD: linear methods (RDA)
dca
```

## Publication-Quality Triplot

**Goal:** Create a journal-ready CCA triplot showing sites, species, and environmental vectors in a single ordination figure.

**Approach:** Extract site, species, and biplot scores from the CCA result, combine into data frames, and overlay points, text labels, and arrow segments in ggplot2.

```r
library(ggplot2)

# Extract scores for ggplot
site_scores <- scores(cca_result, display = 'sites', scaling = 2)
species_scores <- scores(cca_result, display = 'species', scaling = 2)
env_scores <- scores(cca_result, display = 'bp', scaling = 2)

site_df <- as.data.frame(site_scores)
site_df$group <- env_data$habitat_type

sp_df <- as.data.frame(species_scores)
env_df <- as.data.frame(env_scores)
env_df$label <- rownames(env_df)

ggplot() +
    geom_point(data = site_df, aes(x = CCA1, y = CCA2, color = group), size = 3) +
    geom_text(data = sp_df, aes(x = CCA1, y = CCA2, label = rownames(sp_df)),
              size = 2, color = 'grey50', alpha = 0.7) +
    geom_segment(data = env_df, aes(x = 0, y = 0, xend = CCA1, yend = CCA2),
                 arrow = arrow(length = unit(0.2, 'cm')), color = 'red') +
    geom_text(data = env_df, aes(x = CCA1 * 1.1, y = CCA2 * 1.1, label = label),
              color = 'red', size = 3.5) +
    theme_bw() +
    labs(title = 'CCA Triplot')
```

## Related Skills

- biodiversity-metrics - Alpha and beta diversity metrics
- edna-metabarcoding - Generate community data from eDNA
- landscape-genomics - Genotype-environment associations
- microbiome/diversity-analysis - Unconstrained ordination alternative
- data-visualization/ggplot2-fundamentals - Custom ordination plots
