---
name: bio-ecological-genomics-biodiversity-metrics
description: Calculates species richness, diversity, and turnover using the Hill number framework with iNEXT coverage-based rarefaction/extrapolation, asymptotic diversity estimation, and beta diversity partitioning (betapart turnover vs nestedness). Compares assemblages using coverage-standardized rather than size-standardized rarefaction. Use when quantifying biodiversity from species abundance or incidence data, comparing diversity across sites, or constructing rarefaction curves. Not for clinical 16S microbiome alpha/beta diversity (see microbiome/diversity-analysis).
tool_type: r
primary_tool: iNEXT
---

## Version Compatibility

Reference examples tested with: ggplot2 3.5+, vegan 2.6+

Before using code patterns, verify installed versions match. If versions differ:
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Biodiversity Metrics

**"Calculate species diversity for my ecological samples"** → Compute Hill number diversity (richness, Shannon, Simpson) with coverage-based rarefaction/extrapolation using iNEXT, and partition beta diversity into turnover and nestedness components with betapart.
- R: `iNEXT::iNEXT()` for coverage-based rarefaction and extrapolation
- R: `betapart::beta.multi()` for beta diversity partitioning

Calculates alpha diversity using Hill numbers, coverage-based rarefaction/extrapolation, and beta diversity partitioning into turnover and nestedness components.

## Hill Numbers Framework

Hill numbers unify diversity indices into a single parametric family controlled by order q:

| Order (q) | Name | Sensitivity | Equivalent Index |
|-----------|------|-------------|------------------|
| q = 0 | Species richness | Rare species only | S (species count) |
| q = 1 | Shannon diversity | Equal weight all | exp(H') Shannon exponential |
| q = 2 | Simpson diversity | Dominant species | 1/D Simpson inverse |

Higher q values downweight rare species. q = 1 weights species by frequency and is the most balanced measure.

## iNEXT Coverage-Based Rarefaction

**Goal:** Compare diversity across assemblages standardized by sampling completeness rather than sample size.

**Approach:** Run iNEXT coverage-based rarefaction/extrapolation for Hill numbers q=0,1,2 and visualize with ggiNEXT.

Coverage-based rarefaction standardizes by completeness rather than sample size, enabling fair comparison of assemblages sampled with different effort (Chao & Jost 2012).

```r
library(iNEXT)

# abundance data: named list of species abundance vectors per site
abundance_data <- list(
    site_A = c(100, 45, 23, 12, 8, 5, 3, 2, 1, 1),
    site_B = c(80, 60, 40, 30, 20, 15, 10, 8, 5, 3, 2, 1),
    site_C = c(200, 10, 5, 2, 1, 1, 1)
)

# q=c(0,1,2): compute all three Hill numbers
# datatype='abundance': raw counts (use 'incidence_freq' for detection/non-detection)
# nboot=200: bootstrap replicates for confidence intervals
result <- iNEXT(abundance_data, q = c(0, 1, 2), datatype = 'abundance', nboot = 200)

# Rarefaction/extrapolation curves
ggiNEXT(result, type = 1) + theme_bw()

# Sample completeness profiles (coverage vs sample size)
ggiNEXT(result, type = 2) + theme_bw()

# Coverage-based rarefaction curves (diversity vs coverage)
ggiNEXT(result, type = 3) + theme_bw()
```

## Point Estimates at Standardized Coverage

**Goal:** Estimate diversity at a common coverage level for fair cross-site comparison.

**Approach:** Use estimateD with a target coverage of 0.95 to produce standardized Hill number estimates.

```r
# Estimate diversity at a common coverage level
# coverage=0.95: typical target for adequate sampling
# 0.95 means 95% of individuals in the community belong to detected species
est <- estimateD(abundance_data, datatype = 'abundance',
                 base = 'coverage', level = 0.95)
est
```

## Asymptotic Diversity Estimation

**Goal:** Estimate true community diversity if sampling were complete.

**Approach:** Extract asymptotic estimates (Chao1, Chao-Shannon, Chao-Simpson) from iNEXT output.

```r
# Chao1 (q=0), Chao-Shannon (q=1), Chao-Simpson (q=2)
# Asymptotic = estimated true diversity if sampling were complete
asymptotic <- iNEXT(abundance_data, q = c(0, 1, 2), datatype = 'abundance')
asymptotic$AsyEst
```

## iNEXT.3D: Taxonomic, Phylogenetic, and Functional Diversity

**Goal:** Compute diversity across three dimensions: taxonomic, phylogenetic (Faith's PD generalized), and functional (trait space).

**Approach:** Use iNEXT3D with appropriate diversity facet and supplementary data (phylogenetic tree or trait distance matrix).

```r
library(iNEXT.3D)

# Taxonomic diversity (TD) - standard Hill numbers
td <- iNEXT3D(abundance_data, diversity = 'TD', q = c(0, 1, 2), datatype = 'abundance')

# Phylogenetic diversity (PD) - requires ultrametric tree
# PD generalizes Faith's PD via Hill numbers on branch lengths
pd <- iNEXT3D(abundance_data, diversity = 'PD', q = c(0, 1, 2),
              datatype = 'abundance', PDtree = phylo_tree)

# Functional diversity (FD) - requires species distance matrix
# FD captures trait space coverage via Hill numbers
fd <- iNEXT3D(abundance_data, diversity = 'FD', q = c(0, 1, 2),
              datatype = 'abundance', FDdistM = trait_dist_matrix, FDtype = 'AUC')
```

## Classic Diversity with vegan

**Goal:** Calculate standard alpha diversity indices and rarefaction curves from a community matrix.

**Approach:** Compute Shannon, Simpson, and inverse Simpson indices with vegan::diversity, then rarefy to the minimum sample size.

```r
library(vegan)

# community_matrix: sites (rows) x species (columns)
richness <- specnumber(community_matrix)
shannon <- diversity(community_matrix, index = 'shannon')
simpson <- diversity(community_matrix, index = 'simpson')
invsimpson <- diversity(community_matrix, index = 'invsimpson')

# Classic rarefaction to minimum sample size
# Rarefies to the smallest sample, discarding extra reads
raremin <- min(rowSums(community_matrix))
rarecurve(community_matrix, step = 20, sample = raremin)
rare_richness <- rarefy(community_matrix, sample = raremin)
```

## Beta Diversity Partitioning with betapart

**Goal:** Decompose total beta diversity into turnover (species replacement) and nestedness (richness difference) components.

**Approach:** Apply betapart pairwise and multi-site decomposition on presence/absence and abundance matrices using Sorensen, Jaccard, and Bray-Curtis families.

Decomposes total beta diversity (Sorensen or Jaccard) into turnover (species replacement) and nestedness (richness difference) components:

```r
library(betapart)

# Presence/absence matrix required
pa_matrix <- ifelse(community_matrix > 0, 1, 0)

# --- Pairwise decomposition ---
# Sorensen family: beta.sim (turnover) + beta.sne (nestedness) = beta.sor (total)
pair_sor <- beta.pair(pa_matrix, index.family = 'sorensen')
pair_sor$beta.sim   # turnover component
pair_sor$beta.sne   # nestedness component
pair_sor$beta.sor   # total beta diversity

# Jaccard family: beta.jtu (turnover) + beta.jne (nestedness) = beta.jac (total)
pair_jac <- beta.pair(pa_matrix, index.family = 'jaccard')

# --- Multi-site decomposition ---
multi <- beta.multi(pa_matrix, index.family = 'sorensen')
multi$beta.SIM   # multi-site turnover
multi$beta.SNE   # multi-site nestedness
multi$beta.SOR   # multi-site total

# --- Abundance-based beta diversity ---
pair_abund <- beta.pair.abund(community_matrix, index.family = 'bray')
pair_abund$beta.bray.bal   # balanced variation (analogous to turnover)
pair_abund$beta.bray.gra   # abundance gradient (analogous to nestedness)
```

## Interpreting Beta Diversity Components

| Dominance | Ecological Meaning |
|-----------|-------------------|
| Turnover >> Nestedness | Species replacement along gradients; distinct communities |
| Nestedness >> Turnover | Poor sites are subsets of rich sites; nested pattern |
| Both similar | Mixed processes driving community differences |

## Visualization

**Goal:** Visualize the relative contribution of turnover vs nestedness to total beta diversity.

**Approach:** Plot turnover proportion against total beta diversity for each site pair using ggplot2.

```r
library(ggplot2)

# Beta diversity triangle plot (turnover vs nestedness proportions)
beta_df <- data.frame(
    turnover = as.vector(pair_sor$beta.sim),
    nestedness = as.vector(pair_sor$beta.sne),
    total = as.vector(pair_sor$beta.sor)
)
beta_df$turn_prop <- beta_df$turnover / beta_df$total

ggplot(beta_df, aes(x = total, y = turn_prop)) +
    geom_point(alpha = 0.5) +
    geom_hline(yintercept = 0.5, linetype = 'dashed') +
    labs(x = 'Total beta diversity (Sorensen)',
         y = 'Turnover proportion') +
    theme_bw()
```

## Related Skills

- edna-metabarcoding - Generate species tables from eDNA data
- community-ecology - Constrained ordination of assemblages
- microbiome/diversity-analysis - 16S microbiome diversity metrics
- data-visualization/ggplot2-fundamentals - Customize diversity plots
